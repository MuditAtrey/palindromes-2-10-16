/* main.c - command line driver.
 *
 * Commands (see usage() below):
 *   search     the original orbit DFS over a range of bit lengths, with checkpointing and threads
 *   hexdfs     the production search: outside-in over hex digit pairs (search.c), all certified runs used this
 *   bench      node counts and timings per bit length, pruning on and off
 *   estimate   projected run time for a range
 *   selftest   brute force below 2^k, the OEIS A007632 ground truth, and agreement between search modes
 *   verify     check given numbers in bases 2, 10 and 16 (independent of the search code)
 *   orbits     print the orbit structure of one bit length (development aid)
 *
 * Every feature switch of the search itself is an environment variable read in search.c (TP_R5, TP_DUAL, TP_CVP,
 * TP_OFFSET, TP_WEAKPAL and the self-check switches); README.md in this folder lists them. */
/* tripal - command line driver */
#define _GNU_SOURCE
#include "tripal.h"
#include <pthread.h>
#include <stdatomic.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

/* Wall-clock seconds (monotonic), for timings. */
double tp_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Format a large count such as 2^nfree for printing (plain, or as 2^x when huge). */
void tp_fmt_space(double v, char *buf, size_t n)
{
    if (v < 1e6)          snprintf(buf, n, "%.0f", v);
    else if (v < 1e9)     snprintf(buf, n, "%.3fM", v / 1e6);
    else if (v < 1e12)    snprintf(buf, n, "%.3fG", v / 1e9);
    else if (v < 1e15)    snprintf(buf, n, "%.3fT", v / 1e12);
    else if (v < 1e18)    snprintf(buf, n, "%.3fP", v / 1e15);
    else                  snprintf(buf, n, "%.3fE", v / 1e18);
}

/* ================================================================== */
/* configuration                                                      */
/* ================================================================== */

typedef struct {
    int        hex_from, hex_to;   /* -1 unset                        */
    int        bits_from, bits_to;
    unsigned   inv;
    int        threads;
    int        cut;            /* hexdfs: parallel cut depth              */
    int        units;
    int        checkpoint;
    int        verbose;
    int        trace;          /* debug: trace DFS decisions to this depth */
    uint64_t   node_budget;    /* per bit length; 0 = unlimited            */
    long       limit_units;    /* stop after N units (0 = no limit)        */
    double     time_budget;    /* stop after S seconds     (0 = no limit)  */
    const char *workdir;
} cfg_t;

typedef struct {
    pthread_mutex_t mu;
    tp_ckpt        *ckpt;
    int             nw;
    int             collect;
    uint64_t       *out;
    int             nout;
    int             maxout;
    int             verbose;
} hitctx_t;

/* Called by the search for every solution: print it, log it to the checkpoint, and optionally collect it for tests. */
static void hit_cb(void *ud, int B, const bn_t *n, const char *dec, int len)
{
    (void)len;
    hitctx_t *hc = (hitctx_t *)ud;
    char hex[TP_HEXBUF];
    bn_to_hex(n, hc->nw, hex);

    pthread_mutex_lock(&hc->mu);
    if (hc->collect && hc->nout < hc->maxout)
        hc->out[hc->nout++] = n->w[0];
    if (hc->verbose) {
        printf("\n*** HIT  B=%d  n=%s  0x%s ***\n", B, dec, hex);
        fflush(stdout);
    }
    pthread_mutex_unlock(&hc->mu);

    tp_ckpt_hit(hc->ckpt, B, dec, hex);
}

/* ================================================================== */
/* one bit length                                                     */
/* ================================================================== */

typedef struct {
    _Atomic int next;
    _Atomic int ndone;
    _Atomic int stop;
    long        limit_units;
    double      deadline;
} sched_t;

typedef struct {
    tp_search      s;
    const tp_unit *units;
    int            nunits;
    sched_t       *sch;
    tp_ckpt       *ckpt;
    hitctx_t      *hc;
    int            B;
    int            depth;
    int            verbose;
} worker_t;

/* Thread body for `search`: take the next unfinished work unit, run it, mark it done, stop on deadline or limit. */
static void *worker_main(void *arg)
{
    worker_t *w = (worker_t *)arg;
    sched_t  *sc = w->sch;
    for (;;) {
        if (atomic_load(&sc->stop)) break;
        int i = atomic_fetch_add(&sc->next, 1);
        if (i >= w->nunits) break;
        tp_search_run_unit(&w->s, &w->units[i], hit_cb, w->hc);
        tp_ckpt_mark(w->ckpt, w->B, w->depth, w->units[i].path);
        int nd = atomic_fetch_add(&sc->ndone, 1) + 1;
        if ((sc->limit_units > 0 && nd >= sc->limit_units) ||
            (sc->deadline > 0.0 && tp_now() > sc->deadline))
            atomic_store(&sc->stop, 1);
        if (w->verbose) {
            fprintf(stderr, "  [B=%d] unit %d/%d done  nodes=%llu\n", w->B,
                    i + 1, w->nunits, (unsigned long long)w->s.st.nodes);
        }
    }
    return NULL;
}

/* a += b, field by field (merging per-thread statistics). */
static void stats_add(tp_stats *a, const tp_stats *b)
{
    a->nodes    += b->nodes;
    a->prunes   += b->prunes;
    a->leaves   += b->leaves;
    a->hits     += b->hits;
    a->dec_conv += b->dec_conv;
    a->aborted  |= b->aborted;
    a->straddle += b->straddle;
    if (b->maxm > a->maxm) a->maxm = b->maxm;   /* a maximum, not a sum */
}

static int aborted_B[256], aborted_n = 0;

/* returns 0 ok, -1 skipped/inconsistent */
static int run_bit_length(const cfg_t *cf, int B, tp_ckpt *ck, hitctx_t *hc,
                          tp_stats *out, int report, int *pending)
{
    tp_orbits o;
    if (tp_orbits_build(&o, B, cf->inv) != 0) {
        if (report)
            fprintf(stderr, "B=%-4d skipped (%s)\n", B, o.reason);
        tp_orbits_free(&o);
        return -1;
    }

    const int nw    = tp_words(o.nsites);
    const int split = tp_split_depth(o.nfree, cf->units);
    int       cap   = (split >= 30) ? (1 << 30) : (1 << split);
    if (cap < 1) cap = 1;

    tp_unit  *units = (tp_unit *)malloc((size_t)cap * sizeof(tp_unit));
    if (pending) *pending = 0;
    tp_search en;
    tp_search_init(&en, &o, nw, split, units, cap);
    en.trace_depth = cf->trace;
    en.budget      = cf->node_budget;
    double t0 = tp_now();
    tp_search_enum_units(&en);

    int kept = 0;
    for (int i = 0; i < en.nunits; i++)
        if (!tp_ckpt_is_done(ck, B, split, units[i].path))
            units[kept++] = units[i];
    const int skipped = en.nunits - kept;

    const uint64_t enum_nodes = en.st.nodes;   /* split overhead, per run */
    tp_stats total = en.st;

    worker_t   *wk  = NULL;
    pthread_t  *th  = NULL;
    sched_t     sc;
    memset(&sc, 0, sizeof sc);
    atomic_init(&sc.next, 0);
    atomic_init(&sc.ndone, 0);
    atomic_init(&sc.stop, 0);
    sc.limit_units = cf->limit_units;
    sc.deadline    = (cf->time_budget > 0.0) ? tp_now() + cf->time_budget : 0.0;

    if (kept > 0) {
        int nt = cf->threads > kept ? kept : cf->threads;
        if (nt < 1) nt = 1;
        wk = (worker_t *)calloc((size_t)nt, sizeof(worker_t));
        th = (pthread_t *)calloc((size_t)nt, sizeof(pthread_t));
        for (int i = 0; i < nt; i++) {
            tp_search_init(&wk[i].s, &o, nw, INT_MAX, NULL, 0);
            wk[i].s.trace_depth = cf->trace;
            wk[i].s.budget      = cf->node_budget ? cf->node_budget / (uint64_t)nt : 0;
            wk[i].units   = units;
            wk[i].nunits  = kept;
            wk[i].sch     = &sc;
            wk[i].ckpt    = ck;
            wk[i].hc      = hc;
            wk[i].B       = B;
            wk[i].depth   = split;
            wk[i].verbose = cf->verbose;
        }
        for (int i = 0; i < nt; i++)
            pthread_create(&th[i], NULL, worker_main, &wk[i]);
        for (int i = 0; i < nt; i++)
            pthread_join(th[i], NULL);
        for (int i = 0; i < nt; i++)
            stats_add(&total, &wk[i].s.st);
    }
    const int ran = atomic_load(&sc.ndone);
    if (pending) *pending = kept - ran;
    const int incomplete = total.aborted;

    double t1 = tp_now() - t0;

    if (report) {
        char sp[32];
        tp_fmt_space(tp_orbits_space(&o), sp, sizeof sp);
        fprintf(stderr,
                "B=%-4d L=%-3d f=%-3d space=%-9s units=%-5d skip=%-4d pend=%-4d "
                "enum=%-5llu nodes=%-10llu prunes=%-10llu leaves=%-6llu "
                "hits=%llu  %.3fs%s\n",
                B, o.L, o.nfree, sp, ran, skipped, pending ? *pending : 0,
                (unsigned long long)enum_nodes,
                (unsigned long long)total.nodes, (unsigned long long)total.prunes,
                (unsigned long long)total.leaves, (unsigned long long)total.hits,
                t1, incomplete ? "  INCOMPLETE(node budget)" : "");
        fflush(stderr);
    }

    if (out) stats_add(out, &total);
    if (report && incomplete) aborted_B[aborted_n < 256 ? aborted_n++ : 255] = B;

    free(wk); free(th); free(units);
    tp_orbits_free(&o);
    return 0;
}

/* ================================================================== */
/* search command                                                     */
/* ================================================================== */

/* `tripal search`: run every bit length in the requested range with the orbit DFS, threads and checkpoints. */
static int cmd_search(const cfg_t *cf)
{
    int bfrom = cf->bits_from, bto = cf->bits_to;
    if (cf->hex_to >= 0) {
        bfrom = 4 * (cf->hex_from - 1) + 1;
        bto   = 4 * cf->hex_to;
    }
    if (bfrom < 1) bfrom = 1;
    if (bto > TP_MAX_BITS) bto = TP_MAX_BITS;

    tp_ckpt *ck = tp_ckpt_open(cf->workdir, cf->checkpoint);
    if (cf->checkpoint && tp_ckpt_loaded(ck))
        fprintf(stderr, "resuming: %llu work units already done (from %s/done.log)\n",
                (unsigned long long)tp_ckpt_loaded(ck), cf->workdir);

    hitctx_t hc;
    memset(&hc, 0, sizeof hc);
    pthread_mutex_init(&hc.mu, NULL);
    hc.ckpt = ck; hc.nw = TP_MAX_LIMBS; hc.verbose = 1;

    tp_stats total;
    memset(&total, 0, sizeof total);

    double t0 = tp_now();
    int nb = 0, nskip = 0;
    long pending = 0;
    for (int B = bfrom; B <= bto; B++) {
        int pend = 0;
        int r = run_bit_length(cf, B, ck, &hc, &total, 1, &pend);
        if (r == 0) nb++; else nskip++;
        pending += pend;
        if (cf->time_budget > 0.0 && (tp_now() - t0) > cf->time_budget) {
            fprintf(stderr,
                    "--time-budget reached after B=%d; progress is checkpointed\n", B);
            break;
        }
    }
    double dt = tp_now() - t0;

    char sp[32];
    tp_fmt_space((double)total.nodes + (double)total.prunes, sp, sizeof sp);
    printf("\n================ summary ================\n");
    printf("bit lengths searched : %d (%d skipped as inconsistent)\n", nb, nskip);
    printf("work units executed  : %llu nodes + %llu pruned branches  (%s total)\n",
           (unsigned long long)total.nodes, (unsigned long long)total.prunes, sp);
    printf("complete patterns    : %llu leaves\n", (unsigned long long)total.leaves);
    printf("solutions found      : %llu   (+ n=0, trivial, not in any orbit)\n",
           (unsigned long long)total.hits);
    printf("wall clock           : %.2f s\n", dt);
    if (aborted_n > 0) {
        printf("INCOMPLETE bit lengths (node budget hit, NOT proven empty):\n   ");
        for (int i = 0; i < aborted_n; i++) printf(" %d", aborted_B[i]);
        printf("\n");
    }
    if (pending > 0)
        printf("status               : PARTIAL - %ld work units pending, rerun the same\n"
               "                       command to resume where it stopped\n", pending);
    else if (aborted_n > 0)
        printf("status               : PARTIAL - %d bit length(s) hit the node budget\n",
               aborted_n);
    else
        printf("status               : COMPLETE - every unit of every bit length was run\n");
    if (total.hits == 0 && pending == 0 && aborted_n == 0)
        printf("result               : NO triple palindrome in this range\n");
    if (cf->checkpoint)
        printf("logs                 : %s/{done.log,hits.log}\n", cf->workdir);
    printf("=========================================\n");

    pthread_mutex_destroy(&hc.mu);
    tp_ckpt_close(ck);
    return total.hits ? 0 : 0;
}

/* ================================================================== */
/* bench command                                                      */
/* ================================================================== */

/* `tripal bench`: per bit length, compare brute-force size with the pruned tree (nodes, prunes, leaves, time). */
static int cmd_bench(const cfg_t *cf, int *list, int n)
{
    printf("| %5s | %3s | %3s | %16s | %14s | %13s | %8s | %9s |\n",
           "B", "L", "f", "brute force 2^f", "DFS nodes", "prunes",
           "leaves", "time(s)");
    printf("|%s|%s|%s|%s|%s|%s|%s|%s|\n", "-------", "-----", "-----",
           "------------------", "----------------", "---------------",
           "--------", "-----------");

    cfg_t b = *cf;
    b.units      = 1;          /* pure DFS, no split overhead */
    b.checkpoint = 0;
    b.threads    = 1;

    hitctx_t hc;
    memset(&hc, 0, sizeof hc);
    pthread_mutex_init(&hc.mu, NULL);
    hc.nw = TP_MAX_LIMBS; hc.verbose = 1;

    for (int i = 0; i < n; i++) {
        int B = list[i];
        tp_orbits o;
        if (tp_orbits_build(&o, B, b.inv) != 0) {
            printf("| %5d |  -- |  -- | inconsistent: %s\n", B, o.reason);
            tp_orbits_free(&o);
            continue;
        }
        char sp[32];
        tp_fmt_space(tp_orbits_space(&o), sp, sizeof sp);
        const int f = o.nfree, L = o.L;
        tp_orbits_free(&o);

        tp_stats st;
        memset(&st, 0, sizeof st);
        double t0 = tp_now();
        run_bit_length(&b, B, NULL, &hc, &st, 0, NULL);
        double dt = tp_now() - t0;

        printf("| %5d | %3d | %3d | %16s | %14llu | %13llu | %8llu | %9.3f |%s\n",
               B, L, f, sp, (unsigned long long)st.nodes,
               (unsigned long long)st.prunes, (unsigned long long)st.leaves, dt,
               st.aborted ? "  INCOMPLETE (node budget)" : "");
        printf("        span-straddle nodes: %llu\n", (unsigned long long)st.straddle);
        printf("        longest agreed decimal prefix: %u digit(s)\n", st.maxm);
    }
    pthread_mutex_destroy(&hc.mu);
    return 0;
}


/* ================================================================== */
/* estimate: exact node counts and projected runtimes for a bit range  */
/* ================================================================== */

/* `tripal estimate`: exact node counts for small lengths and a projected run time for the range. */
static int cmd_estimate(int bfrom, int bto, int threads, double rate)
{
    if (bfrom < 1) bfrom = 1;
    if (bto > TP_MAX_BITS) bto = TP_MAX_BITS;
    if (threads < 1) threads = 1;

    printf("node rate assumed: %.3g nodes/s per core\n", rate);
    printf("'full tree' is 2^(nfree+1)-2, the size with no pruning at all.\n");
    printf("'kept' is the fraction of the split-depth units that survive the\n");
    printf("real pruner; for B = 3 (mod 4) all pruning happens above the split\n");
    printf("depth, so full*kept is the exact node count.  For other B pruning\n");
    printf("continues deeper, so full*kept is only an upper bound.\n\n");
    printf("|%6s|%6s|%5s|%5s|%14s|%7s|%14s|%10s|%s\n",
           " hex ", " B ", " mod4", " f ", " full tree ", " kept ", " est. nodes ", " 1 core ", " note ");
    printf("|%s|%s|%s|%s|%s|%s|%s|%s|%s\n", "------","------","-----","-----",
           "--------------","-------","--------------","----------","------------------");

    long double tnodes = 0.0L, tsec = 0.0L;
    for (int B = bfrom; B <= bto; B++) {
        tp_orbits o;
        if (tp_orbits_build(&o, B, TP_INV_BIN | TP_INV_HEX) != 0) continue;

        const int nw    = tp_words(o.nsites);
        const int split = tp_split_depth(o.nfree, 256);
        int cap = (split >= 30) ? (1 << 30) : (1 << split);
        if (cap < 1) cap = 1;
        tp_unit *units = (tp_unit *)malloc((size_t)cap * sizeof(tp_unit));
        tp_search en;
        tp_search_init(&en, &o, nw, split, units, cap);
        tp_search_enum_units(&en);
        const long double kept = (split < 30)
            ? (long double)en.nunits / (long double)(1 << split) : 1.0L;
        free(units);
        const int f = o.nfree;
        tp_orbits_free(&o);

        const long double full  = ldexpl(1.0L, f + 1) - 2.0L;
        const long double nodes = full * kept;
        const long double secs  = nodes / (long double)rate;
        const long double par   = secs / ((long double)threads * 0.6L);

        const char *note = "";
        if ((B & 3) == 3)  note = "exact when kept < 1";
        else if (kept < 0.001L) note = "upper bound (prunes deeper)";
        else               note = "upper bound";

        printf("|%6d|%6d|%5d|%5d|%14.4Le|%6.2Lf%%|%14.4Le|%9.2Lfs| %s\n",
               (B + 3) / 4, B, B & 3, f, full, kept * 100.0L, nodes, secs, note);
        (void)par;
        if ((B & 3) == 3) { tnodes += nodes; tsec += secs; }
    }
    printf("\nB = 3 (mod 4) rows only -- these are the ones that cost anything:\n");
    printf("    %.4Le nodes total\n", tnodes);
    printf("    %.2Lf h single-threaded, %.2Lf h at %d threads (~60%% scaling)\n",
           tsec / 3600.0L, tsec / ((long double)threads * 0.6L) / 3600.0L, threads);
    printf("\ncalibrate --rate on your own machine with: tripal bench --bits 131\n");
    return 0;
}


/* ================================================================== */
/* hexdfs: outside-in hex-digit traversal (alternative tree shape)     */
/* ================================================================== */

/* Hit callback for `hexdfs`: print the number in decimal as a HIT line (the sweep scripts parse these). */
static void hexdf_hit(void *ud, int L, const bn_t *n, const char *dec, int declen)
{
    (void)ud; (void)n; (void)L; (void)declen;
    printf("HIT %s\n", dec);
}

/* `tripal hexdfs`: the production search over hex lengths, one table row per bit length. */
static int cmd_hexdfs(int bfrom, int bto, int trace, int threads, int cut)
{
    if (bfrom < 1) bfrom = 1;
    if (bto > TP_MAX_BITS) bto = TP_MAX_BITS;
    printf("|%6s|%6s|%14s|%13s|%12s|%8s|%9s|%s\n",
           " hex ", " B ", " nodes ", " prunes ", " leaves ", " hits ", " time(s) ", "");
    unsigned long long tn = 0, tp = 0, tl = 0, th = 0;
    double ttot = 0.0;
    for (int B = bfrom; B <= bto; B++) {
        tp_stats st;
        memset(&st, 0, sizeof st);
        const double t0 = tp_now();
        /* a shard request must go through the splitting path even on one
           thread, or the machine silently searches the whole tree */
        const int rc = (threads > 1 || tp_hexdfs_sharded())
            ? tp_hexdfs_run_mt(B, TP_INV_BIN | TP_INV_HEX, threads, cut,
                               &st, hexdf_hit, NULL)
            : tp_hexdfs_run(B, TP_INV_BIN | TP_INV_HEX, &st, hexdf_hit, NULL, trace);
        if (st.overflow) {
            printf("|%6d|%6d|  TASK CAP EXCEEDED -- INCOMPLETE, raise --cut lower or cap\n",
                   (B + 3) / 4, B);
            return 2;
        }
        if (rc != 0) {
            printf("|%6d|%6d|  inconsistent\n", (B + 3) / 4, B);
            continue;
        }
        const double dt = tp_now() - t0;
        printf("|%6d|%6d|%14llu|%13llu|%12llu|%8llu|%9.3f|\n",
               (B + 3) / 4, B, (unsigned long long)st.nodes,
               (unsigned long long)st.prunes, (unsigned long long)st.leaves,
               (unsigned long long)st.hits, dt);
        tn += st.nodes; tp += st.prunes; tl += st.leaves; th += st.hits; ttot += dt;
        tp_prof_dump();
    }
    printf("\ntotal: %llu nodes, %llu prunes, %llu leaves, %llu hits, %.3f s\n",
           tn, tp, tl, th, ttot);
    return 0;
}

/* ================================================================== */
/* brute force reference                                              */
/* ================================================================== */

/* Brute-force referee for 64-bit n: palindrome in bases 2, 10 and 16, checked digit by digit. */
static inline int is_tri_pal_u64(uint64_t n)
{
    if (n == 0) return 1;
    int hi = 63 - __builtin_clzll(n);
    for (int i = 0, j = hi; i < j; i++, j--)
        if (((n >> i) & 1u) != ((n >> j) & 1u)) return 0;
    char d[24];
    int L = 0;
    uint64_t t = n;
    while (t) { d[L++] = (char)('0' + (int)(t % 10)); t /= 10; }
    for (int i = 0, j = L - 1; i < j; i++, j--)
        if (d[i] != d[j]) return 0;
    int hx = hi >> 2;
    for (int i = 0, j = hx; i < j; i++, j--)
        if (((n >> (4 * i)) & 0xFu) != ((n >> (4 * j)) & 0xFu)) return 0;
    return 1;
}


/* Independent reference enumerator: walks every binary palindrome of B bits
 * (there are 2^(ceil(B/2)-1) of them) and tests the other two bases. It shares
 * no code with the orbit engine or the pruner, so it is a genuine cross-check.
 * Returns the number of hits; stores up to maxout values in out.            */
static int ref_enum(int bfrom, int bto, uint64_t *out, int maxout)
{
    int n = 0;
    for (int B = bfrom; B <= bto; B++) {
        const int half = (B + 1) / 2;
        if (half > 62) break;
        const uint64_t end = (uint64_t)1 << (half - 1);
        for (uint64_t hi = 1; hi < (end << 1); hi += 2) {   /* bit0 = 1 */
            uint64_t v = hi;
            for (int i = half; i < B; i++)
                v |= ((hi >> (B - 1 - i)) & 1u) << i;
            if (is_tri_pal_u64(v)) {
                if (n < maxout) out[n] = v;
                n++;
            }
        }
    }
    return n;
}

/* ================================================================== */
/* selftest                                                           */
/* ================================================================== */

/* Load OEIS A007632 (palindromes in bases 2 and 10) from data/oeis/b007632.txt, keeping terms below 2^maxbits.
 * The self-test drops the base-16 condition and checks that the engine reproduces exactly this list. */
static int load_b007632(uint64_t *out, int maxn, int maxbits)
{
    const char *paths[] = { "data/oeis/b007632.txt",
                            "../data/oeis/b007632.txt",
                            "tripal/data/oeis/b007632.txt" };
    FILE *f = NULL;
    for (unsigned i = 0; i < sizeof paths / sizeof paths[0] && !f; i++)
        f = fopen(paths[i], "r");
    if (!f) return -1;
    int n = 0;
    char line[512];
    while (fgets(line, sizeof line, f) && n < maxn) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char *sp = strchr(line, ' ');
        if (!sp) continue;
        uint64_t v = strtoull(sp + 1, NULL, 10);
        if (maxbits > 0 && (v >> maxbits) != 0) continue;
        out[n++] = v;
    }
    fclose(f);
    return n;
}

/* qsort comparator for uint64_t. */
static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : x > y;
}

/* `tripal selftest`: three independent checks.
 *   1. brute force over every n < 2^brute_bits against the engine;
 *   2. the engine in base-2/10 mode (hex mirror off) against OEIS A007632;
 *   3. the reference enumerator ref_enum (shares no code with the engine) against the engine up to ref_bits.
 * The Python suite in tests/ adds the pruning-specific positive controls. */
static int cmd_selftest(const cfg_t *cf, int brute_bits, int dual_bits,
                        int ref_bits)
{
    int failures = 0;

    /* ---------------- test 1: brute force, B <= brute_bits -------------- */
    fprintf(stderr, "\n[1/3] brute-force cross-check, n < 2^%d\n", brute_bits);
    uint64_t limit = (brute_bits >= 63) ? ~(uint64_t)0
                                        : ((uint64_t)1 << brute_bits);
    uint64_t *bf = (uint64_t *)malloc(sizeof(uint64_t) * 1024);
    int nbf = 0;
    bf[nbf++] = 0;                       /* 0 is a palindrome in every base */
    double t0 = tp_now();
    for (uint64_t n = 1; n < limit; n++)
        if (is_tri_pal_u64(n)) bf[nbf++] = n;
    fprintf(stderr, "      brute force: %d triple palindromes in %.2f s\n",
            nbf, tp_now() - t0);

    cfg_t c = *cf;
    c.checkpoint = 0;
    c.units      = 256;
    c.inv        = TP_INV_BIN | TP_INV_HEX;

    hitctx_t hc;
    memset(&hc, 0, sizeof hc);
    pthread_mutex_init(&hc.mu, NULL);
    hc.nw        = TP_MAX_LIMBS;
    hc.collect   = 1;
    hc.out       = (uint64_t *)malloc(sizeof(uint64_t) * 1024);
    hc.maxout    = 1024;
    hc.verbose   = 0;
    hc.out[0]    = 0;                    /* n = 0: covered outside the orbit
                                          * search, which starts at B = 1     */
    hc.nout      = 1;

    tp_stats st; memset(&st, 0, sizeof st);
    for (int B = 1; B <= brute_bits; B++)
        run_bit_length(&c, B, NULL, &hc, &st, 0, NULL);

    qsort(hc.out, (size_t)hc.nout, sizeof(uint64_t), cmp_u64);
    qsort(bf, (size_t)nbf, sizeof(uint64_t), cmp_u64);

    int same = (hc.nout == nbf);
    if (same)
        for (int i = 0; i < nbf; i++)
            if (hc.out[i] != bf[i]) { same = 0; break; }

    printf("[1] triple-palindrome search vs brute force, B <= %d\n", brute_bits);
    printf("    brute force found : %d\n", nbf);
    printf("    search    found   : %d\n", hc.nout);
    printf("    values            :");
    for (int i = 0; i < nbf && i < 64; i++) printf(" %llu", (unsigned long long)bf[i]);
    printf("\n    verdict           : %s\n\n", same ? "PASS (exact match)" : "FAIL");
    if (!same) failures++;

    /* ---------------- test 2: A007632 (base 2 + base 10) ---------------- */
    fprintf(stderr, "[2/3] A007632 reproduction, bases {2,10}, B <= %d\n", dual_bits);
    uint64_t expect[512];
    int nexp = load_b007632(expect, 512, dual_bits);
    if (nexp < 0) {
        printf("[2] A007632 reproduction\n    verdict           : SKIP (data file not found)\n\n");
    } else {
        cfg_t d = *cf;
        d.checkpoint = 0;
        d.units      = 256;
        d.inv        = TP_INV_BIN;

        hitctx_t h2;
        memset(&h2, 0, sizeof h2);
        pthread_mutex_init(&h2.mu, NULL);
        h2.nw      = TP_MAX_LIMBS;
        h2.collect = 1;
        h2.out     = (uint64_t *)malloc(sizeof(uint64_t) * 4096);
        h2.maxout  = 4096;
        h2.verbose = 0;
        h2.out[0]  = 0;                  /* trivial solution, added explicitly */
        h2.nout    = 1;

        tp_stats st2; memset(&st2, 0, sizeof st2);
        for (int B = 1; B <= dual_bits; B++)
            run_bit_length(&d, B, NULL, &h2, &st2, 0, NULL);

        qsort(h2.out, (size_t)h2.nout, sizeof(uint64_t), cmp_u64);
        qsort(expect, (size_t)nexp, sizeof(uint64_t), cmp_u64);

        int ok = (h2.nout == nexp);
        int firstbad = -1;
        if (ok)
            for (int i = 0; i < nexp; i++)
                if (h2.out[i] != expect[i]) { ok = 0; firstbad = i; break; }

        printf("[2] A007632 (palindromic in bases 2 and 10), B <= %d\n", dual_bits);
        printf("    OEIS terms in range : %d\n", nexp);
        printf("    search      found   : %d\n", h2.nout);
        if (firstbad >= 0)
            printf("    first mismatch      : index %d  expected %llu got %llu\n",
                   firstbad, (unsigned long long)expect[firstbad],
                   (unsigned long long)h2.out[firstbad]);
        printf("    workloads           : %llu nodes, %llu prunes\n",
               (unsigned long long)st2.nodes, (unsigned long long)st2.prunes);
        printf("    verdict             : %s\n\n",
               ok ? "PASS (reproduced the full OEIS list)" : "FAIL");
        if (!ok) failures++;
        free(h2.out);
        pthread_mutex_destroy(&h2.mu);
    }

    /* ---------------- test 3: verify every hit in three bases ----------- */
    fprintf(stderr, "[3/3] independent re-verification of every hit\n");
    int vb = 0;
    for (int i = 0; i < hc.nout; i++) {
        uint64_t n = hc.out[i];
        bn_t x; int nw = tp_words(64);
        bn_zero(&x, nw);
        x.w[0] = n;
        if (!(tp_is_pal_base(&x, nw, 2) && tp_is_pal_base(&x, nw, 10) &&
              tp_is_pal_base(&x, nw, 16))) vb++;
    }
    printf("[3] every reported solution re-checked in bases 2, 10, 16\n");
    printf("    solutions examined  : %d\n", hc.nout);
    printf("    bad                 : %d\n", vb);
    printf("    verdict             : %s\n\n", vb ? "FAIL" : "PASS");
    if (vb) failures++;


    /* ---------------- test 4: independent palindrome enumeration -------- */
    fprintf(stderr, "[4/4] binary-palindrome reference enumeration, B <= %d\n",
            ref_bits);
    {
        uint64_t *ref = (uint64_t *)malloc(sizeof(uint64_t) * 1024);
        int nref = ref_enum(1, ref_bits, ref, 1024);
        qsort(ref, (size_t)nref, sizeof(uint64_t), cmp_u64);

        cfg_t e = *cf;
        e.checkpoint = 0;
        e.units      = 256;
        e.inv        = TP_INV_BIN | TP_INV_HEX;

        hitctx_t h4;
        memset(&h4, 0, sizeof h4);
        pthread_mutex_init(&h4.mu, NULL);
        h4.nw      = TP_MAX_LIMBS;
        h4.collect = 1;
        h4.out     = (uint64_t *)malloc(sizeof(uint64_t) * 1024);
        h4.maxout  = 1024;
        h4.verbose = 0;
        h4.out[0]  = 0;
        h4.nout    = 1;

        tp_stats st4; memset(&st4, 0, sizeof st4);
        for (int B = 1; B <= ref_bits; B++)
            run_bit_length(&e, B, NULL, &h4, &st4, 0, NULL);
        qsort(h4.out, (size_t)h4.nout, sizeof(uint64_t), cmp_u64);

        int ok4 = (h4.nout == nref + 1);
        if (ok4)
            for (int i = 0; i < nref; i++)
                if (h4.out[i + 1] != ref[i]) { ok4 = 0; break; }

        printf("[4] search vs. exhaustive binary-palindrome enumeration, B <= %d\n",
               ref_bits);
        printf("    enumeration found   : %d (+ n=0)\n", nref);
        printf("    search      found   : %d\n", h4.nout);
        printf("    orbit workloads     : %llu nodes, %llu prunes\n",
               (unsigned long long)st4.nodes, (unsigned long long)st4.prunes);
        printf("    verdict             : %s\n\n", ok4 ? "PASS (exact match)" : "FAIL");
        if (!ok4) failures++;
        free(ref); free(h4.out);
        pthread_mutex_destroy(&h4.mu);
    }

    printf("================ selftest %s ================\n",
           failures ? "FAILED" : "PASSED");

    free(bf); free(hc.out);
    pthread_mutex_destroy(&hc.mu);
    return failures ? 1 : 0;
}

/* ================================================================== */
/* verify command                                                     */
/* ================================================================== */

/* Parse a decimal number, or a hex number with a 0x prefix. */
static int parse_num(const char *s, bn_t *out, int nw)
{
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
        bn_zero(out, nw);
        for (s += 2; *s; s++) {
            int d;
            if (*s >= '0' && *s <= '9') d = *s - '0';
            else if (*s >= 'a' && *s <= 'f') d = *s - 'a' + 10;
            else if (*s >= 'A' && *s <= 'F') d = *s - 'A' + 10;
            else return 0;
            bn_muladd(out, nw, 16, (uint64_t)d);
        }
        return 1;
    }
    return bn_from_dec(s, nw, out);
}

/* `tripal verify n ...`: print each n in decimal and hex and say whether it is a palindrome in bases 2, 10 and 16. */
static int cmd_verify(int argc, char **argv)
{
    if (argc < 1) { fprintf(stderr, "usage: tripal verify <n>...\n"); return 2; }
    int bad = 0;
    for (int i = 0; i < argc; i++) {
        bn_t x; int nw = TP_MAX_LIMBS;
        if (!parse_num(argv[i], &x, nw)) {
            printf("%s: not a valid number\n", argv[i]);
            bad++;
            continue;
        }
        char dec[TP_DECBUF], hex[TP_HEXBUF];
        bn_to_dec(&x, nw, dec);
        bn_to_hex(&x, nw, hex);
        int p2 = tp_is_pal_base(&x, nw, 2);
        int p10 = tp_is_pal_base(&x, nw, 10);
        int p16 = tp_is_pal_base(&x, nw, 16);
        printf("n     = %s\n", dec);
        printf("hex   = 0x%s\n", hex);
        printf("base2 %s   base10 %s   base16 %s\n",
               p2 ? "YES" : "no", p10 ? "YES" : "no", p16 ? "YES" : "no");
        printf("=> %s\n\n",
               (p2 && p10 && p16) ? "TRIPLE PALINDROME" : "not a triple palindrome");
        if (!(p2 && p10 && p16)) bad++;
    }
    return bad ? 1 : 0;
}

/* ================================================================== */
/* argument parsing                                                   */
/* ================================================================== */

/* Parse a comma separated list of integers such as "48,80,128". */
static int parse_int_list(const char *s, int *out, int maxn)
{
    int n = 0;
    while (*s && n < maxn) {
        out[n++] = (int)strtol(s, (char **)&s, 10);
        while (*s == ',' || *s == ' ') s++;
    }
    return n;
}

static int cmd_orbits(int B);
static double g_rate = 12e6;

/* Print the command summary. */
static void usage(void)
{
    printf(
"tripal - triple palindrome (base 2 / 10 / 16) search\n"
"\n"
"  tripal search  [--hex-from A] [--hex-to B] [--bits-from A] [--bits-to B]\n"
"                 [--threads T] [--units N] [--workdir DIR] [--no-checkpoint]\n"
"                 [--bases 2,10,16] [--limit-units N] [--time-budget SEC]\n"
"                 [--node-budget N] [--verbose]\n"
"  tripal bench   --bits 48,80,128,184,216   [--node-budget N]\n"
"  tripal selftest [--brute-bits 26] [--dual-bits 40] [--ref-bits 46]\n"
"  tripal verify  <n> [<n> ...]         (decimal or 0x-hex)\n"
"  tripal ref     --bits-from A --bits-to B     (independent enumeration)\n"
"  tripal brute   --bits-from A --bits-to B     (dumb full scan)\n"
"  tripal hexdfs  --hex-from A --hex-to B        (outside-in, alternate tree)\n"
"                 [--threads T] [--cut C] [--shard K/N]   (N machines, this is K)\n"
"                 [--bloom U:b[,U:b]]  deep reach filter at U undecided orbits, b bits/element\n"
"                 [--r5 J]  5-adic reach filter, modulus 5^J (default 12, 0 = off;\n"
"                           memory ~ 1.1 GB at 12, 4.7 GB at 13)\n"
"  tripal estimate --hex-from A --hex-to B [--threads T] [--rate R]\n"
"                                               (exact tree sizes + runtimes)\n"
"  tripal orbits  --bits-to B                   (dump orbit structure)\n"
"\n"
"  --hex N          convenience: bit length 4N (largest B with N hex digits)\n"
"  --units N        work-unit target; sets the split depth (default 256)\n"
"  --limit-units N  stop gracefully after N units, checkpointed\n"
"  --time-budget S  stop gracefully after S seconds, checkpointed\n"
"  --node-budget N  abort a bit length past N nodes -> reported INCOMPLETE\n"
"  --trace D        print DFS decisions down to depth D (debug)\n"
);
}

/* Parse the command line and dispatch to one of the cmd_* functions above. */
int main(int argc, char **argv)
{
    if (argc < 2) { usage(); return 2; }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) { usage(); return 0; }

    cfg_t cf;
    memset(&cf, 0, sizeof cf);
    cf.hex_from   = -1;
    cf.hex_to     = -1;
    cf.bits_from  = 1;
    cf.bits_to    = 64;
    cf.inv        = TP_INV_BIN | TP_INV_HEX;
    cf.threads    = 1;
    cf.cut        = 12;
    cf.units      = 256;
    cf.checkpoint = 1;
    cf.workdir    = "runs/default";

    { /* default threads = logical cpus */
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        cf.threads = (n > 0) ? (int)n : 1;
    }

    const char *cmd = argv[1];
    int list[256], nlist = 0;
    int brute_bits = 26, dual_bits = 40, ref_bits = 50;

    for (int i = 2; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] != '-') continue;           /* positional arg, e.g. verify */
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
#define NEEDVAL() do { if (!v) { fprintf(stderr, "%s needs a value\n", a); return 2; } } while (0)
        if (!strcmp(a, "--rate"))          { NEEDVAL(); g_rate = atof(v); i++; }
        else if (!strcmp(a, "--hex-from")) { NEEDVAL(); cf.hex_from = atoi(v); i++; }
        else if (!strcmp(a, "--hex-to"))   { NEEDVAL(); cf.hex_to = atoi(v); i++; }
        else if (!strcmp(a, "--bits-from")){ NEEDVAL(); cf.bits_from = atoi(v); i++; }
        else if (!strcmp(a, "--bits-to"))  { NEEDVAL(); cf.bits_to = atoi(v); i++; }
        else if (!strcmp(a, "--threads"))  { NEEDVAL(); cf.threads = atoi(v); i++; }
        else if (!strcmp(a, "--cut"))      { NEEDVAL(); cf.cut = atoi(v); i++; }
        else if (!strcmp(a, "--sieve"))    { setenv("TP_SIEVE", "1", 1); }
        else if (!strcmp(a, "--r5"))       { NEEDVAL(); setenv("TP_R5", v, 1); i++; }
        else if (!strcmp(a, "--bloom"))    { NEEDVAL(); setenv("TP_BLOOM", v, 1); i++; }
        else if (!strcmp(a, "--dual"))     { setenv("TP_DUAL", "1", 1); }
        else if (!strcmp(a, "--cvp"))      { setenv("TP_CVP", "1", 1); }
        else if (!strcmp(a, "--shard"))    { NEEDVAL();
                                             { int k2 = 0, n2 = 1;
                                               if (sscanf(v, "%d/%d", &k2, &n2) != 2 || n2 < 1) {
                                                   fprintf(stderr, "--shard wants K/N with N >= 1\n");
                                                   return 2;
                                               }
                                               tp_hexdfs_set_shard(k2, n2); }
                                             i++; }
        else if (!strcmp(a, "--units"))    { NEEDVAL(); cf.units = atoi(v); i++; }
        else if (!strcmp(a, "--workdir"))  { NEEDVAL(); cf.workdir = v; i++; }
        else if (!strcmp(a, "--brute-bits")){ NEEDVAL(); brute_bits = atoi(v); i++; }
        else if (!strcmp(a, "--dual-bits")){ NEEDVAL(); dual_bits = atoi(v); i++; }
        else if (!strcmp(a, "--ref-bits")) { NEEDVAL(); ref_bits = atoi(v); i++; }
        else if (!strcmp(a, "--bits"))     { NEEDVAL(); nlist = parse_int_list(v, list, 256); i++; }
        else if (!strcmp(a, "--hex"))      {
            NEEDVAL();
            int tmp[256], tn = parse_int_list(v, tmp, 256);
            for (int k = 0; k < tn && nlist < 256; k++) list[nlist++] = 4 * tmp[k];
            i++;
        }
        else if (!strcmp(a, "--bases")) {
            NEEDVAL();
            unsigned inv = 0;
            const char *p = v;
            while (*p) {
                int b = (int)strtol(p, (char **)&p, 10);
                if      (b == 2)  inv |= TP_INV_BIN;
                else if (b == 16) inv |= TP_INV_HEX;
                else if (b != 10) {
                    fprintf(stderr, "only bases 2, 10 and 16 are supported\n");
                    return 2;
                }
                while (*p == ',' || *p == ' ') p++;
            }
            cf.inv = inv ? inv : (TP_INV_BIN | TP_INV_HEX);
            i++;
        }
        else if (!strcmp(a, "--node-budget")) { NEEDVAL(); cf.node_budget = strtoull(v, NULL, 10); i++; }
        else if (!strcmp(a, "--trace"))       { NEEDVAL(); cf.trace = atoi(v); i++; }
        else if (!strcmp(a, "--limit-units")) { NEEDVAL(); cf.limit_units = atol(v); i++; }
        else if (!strcmp(a, "--time-budget")) { NEEDVAL(); cf.time_budget = atof(v); i++; }
        else if (!strcmp(a, "--no-checkpoint")) cf.checkpoint = 0;
        else if (!strcmp(a, "--verbose"))       cf.verbose = 1;
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(); return 0; }
        else { fprintf(stderr, "unknown option: %s\n", a); return 2; }
#undef NEEDVAL
    }

    if (!strcmp(cmd, "search")) {
        if (cf.hex_to < 0 && cf.hex_from < 0 && cf.bits_to == 64 && cf.bits_from == 1)
            { cf.hex_from = 1; cf.hex_to = 16; }
        else if (cf.hex_to >= 0 && cf.hex_from < 0) cf.hex_from = 1;
        return cmd_search(&cf);
    }
    if (!strcmp(cmd, "bench")) {
        if (nlist == 0) { fprintf(stderr, "bench needs --bits or --hex\n"); return 2; }
        return cmd_bench(&cf, list, nlist);
    }
    if (!strcmp(cmd, "selftest")) return cmd_selftest(&cf, brute_bits, dual_bits, ref_bits);
    if (!strcmp(cmd, "verify"))   return cmd_verify(argc - 2, argv + 2);
    if (!strcmp(cmd, "hexdfs")) {
        int bfrom = cf.bits_from, bto = cf.bits_to;
        if (cf.hex_to >= 0) { bfrom = 4 * (cf.hex_from - 1) + 1; bto = 4 * cf.hex_to; }
        return cmd_hexdfs(bfrom, bto, cf.trace,
                          cf.threads > 0 ? cf.threads : 1,
                          cf.cut > 0 ? cf.cut : 12);
    }
    if (!strcmp(cmd, "estimate")) {
        int bfrom = cf.bits_from, bto = cf.bits_to;
        if (cf.hex_to >= 0) { bfrom = 4 * (cf.hex_from - 1) + 1; bto = 4 * cf.hex_to; }
        return cmd_estimate(bfrom, bto, cf.threads, g_rate);
    }
    if (!strcmp(cmd, "orbits")) {
        int b = cf.bits_to;
        for (int i = 2; i < argc; i++) if (!strncmp(argv[i], "--bits-to", 9) && i+1 < argc) b = atoi(argv[i+1]);
        return cmd_orbits(b);
    }
    if (!strcmp(cmd, "ref")) {
        uint64_t *out = (uint64_t *)malloc(sizeof(uint64_t) * 65536);
        int n = ref_enum(cf.bits_from, cf.bits_to, out, 65536);
        for (int i = 0; i < n && i < 65536; i++) {
            char hex[TP_HEXBUF];
            bn_t x; int nw = tp_words(64);
            bn_zero(&x, nw); x.w[0] = out[i];
            bn_to_hex(&x, nw, hex);
            printf("%d %llu 0x%s\n", 64 - __builtin_clzll(out[i] | 1),
                   (unsigned long long)out[i], hex);
        }
        fprintf(stderr, "reference enumeration B=%d..%d: %d hits\n",
                cf.bits_from, cf.bits_to, n);
        free(out);
        return 0;
    }
    if (!strcmp(cmd, "brute")) {
        if (cf.bits_to <= 0) cf.bits_to = 30;
        uint64_t lo = (cf.bits_from <= 1) ? 0
                    : ((uint64_t)1 << (cf.bits_from - 1));
        uint64_t hi = (cf.bits_to >= 63) ? ~(uint64_t)0
                    : ((uint64_t)1 << cf.bits_to);
        int cnt = 0;
        for (uint64_t n = lo; n < hi; n++)
            if (is_tri_pal_u64(n)) { printf("%llu\n", (unsigned long long)n); cnt++; }
        fprintf(stderr, "brute force [%llu, %llu): %d hits\n",
                (unsigned long long)lo, (unsigned long long)hi, cnt);
        return 0;
    }
    fprintf(stderr, "unknown command: %s\n", cmd);
    usage();
    return 2;
}

/* ---- development helper: dump the orbit structure ------------------- */
static int cmd_orbits(int B)
{
    tp_orbits o;
    if (tp_orbits_build(&o, B, TP_INV_BIN | TP_INV_HEX) != 0) {
        printf("B=%d invalid: %s\n", B, o.reason);
        return 1;
    }
    printf("B=%d L=%d nsites=%d ncomp=%d nfree=%d\n", B, o.L, o.nsites, o.ncomp, o.nfree);
    printf("forced: ");
    for (int c = 0; c < o.ncomp; c++)
        if (o.cval[c] >= 0) printf("{%d..%d}=%d ", o.cmin[c], o.cmax[c], o.cval[c]);
    printf("\n%-4s %-6s %-6s %-9s %s\n", "#", "cmin", "cmax", "popcount", "suffixmin");
    for (int d = 0; d < o.nfree; d++) {
        int c = o.order[d], pc = 0;
        for (int i = 0; i < o.nsites; i++) if (o.comp[i] == c) pc++;
        printf("%-4d %-6d %-6d %-9d %d\n", d, o.cmin[c], o.cmax[c], pc, o.suffixmin[d]);
    }
    tp_orbits_free(&o);
    return 0;
}
