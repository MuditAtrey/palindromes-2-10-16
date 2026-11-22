/* Two-ended DFS with decimal-interval pruning.
 *
 * Invariant maintained while descending the orbit tree:
 *
 *   lo = every free site cleared      hi = every free site set
 *
 * so every completion n of the current partial assignment satisfies
 *
 *   lo <= n <= hi .
 *
 * If lo and hi have the same decimal length D and share a prefix P of length
 * m, then *every* n in [lo,hi] starts with P.  A D-digit decimal palindrome
 * starting with P must end with reverse(P), hence
 *
 *   n = reverse(P)  (mod 10^m)   =>   n = reverse(P)  (mod 2^m).
 *
 * Independently, all sites below r = min site of any still-free orbit are
 * already fixed, so n mod 2^r is known.  With k = min(m, r, 64) the two
 * statements must agree on the low k bits, otherwise the entire subtree of
 * 2^(remaining) completions is discarded at once.
 */
#include <pthread.h>
#ifdef TP_FLINT
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <flint/fmpz_vec.h>
#endif
#include <math.h>
#include <unistd.h>
#include <stdatomic.h>
#include <stdlib.h>
#include "tripal.h"
#include "viable_outer.h"

/* ------------------------------------------------------------------ */

static void setup_state(tp_search *s)
{
    const tp_orbits *o = s->o;
    bn_zero(&s->lo, s->nw);
    bn_zero(&s->hi, s->nw);
    for (int i = 0; i < o->nsites; i++) {
        int     c = o->comp[i];
        if (o->cval[c] == 1) bn_setbit(&s->lo, i);
        if (i < o->B)        bn_setbit(&s->hi, i);
    }
}

void tp_search_init(tp_search *s, const tp_orbits *o, int nw, int stop_depth,
                    tp_unit *units, int unit_cap)
{
    memset(s, 0, sizeof *s);
    s->o          = o;
    s->nw         = nw;
    s->stop_depth = stop_depth;
    s->units      = units;
    s->unit_cap   = unit_cap;
    setup_state(s);
}

void tp_search_roots(const tp_search *s, char *dlo, int *dl, char *dhi, int *dh)
{
    *dl = bn_to_dec(&s->lo, s->nw, dlo);
    *dh = bn_to_dec(&s->hi, s->nw, dhi);
}

static uint64_t g_m0 = 0, g_calls = 0;   /* diagnostics: m==0 rate */

/* 10^k as a bigint, built lazily: used to get the exact decimal length of a
 * number from its bit length plus one or two comparisons. */
static bn_t g_pow10[TP_DECBUF];
static int  g_pow10_max = -1;
static void pow10_ensure(int k)
{
    const int nw = TP_MAX_LIMBS;
    if (g_pow10_max < 0) {
        bn_zero(&g_pow10[0], nw);
        g_pow10[0].w[0] = 1;
        g_pow10_max = 0;
    }
    if (k > TP_DECBUF - 1) k = TP_DECBUF - 1;
    for (int i = g_pow10_max + 1; i <= k; i++) {
        g_pow10[i] = g_pow10[i - 1];
        bn_muladd(&g_pow10[i], nw, 10, 0);
    }
    if (k > g_pow10_max) g_pow10_max = k;
}

static int bn_bitlen(const bn_t *a, int nw)
{
    for (int i = nw - 1; i >= 0; i--)
        if (a->w[i]) {
            int b = 64 * i;
            uint64_t v = a->w[i];
            while (v) { b++; v >>= 1; }
            return b;
        }
    return 0;
}

/* exact number of decimal digits */
static int dec_digits_of(const bn_t *a, int nw)
{
    const int b = bn_bitlen(a, nw);
    if (b <= 1) return 1;
    int d = (int)((double)(b - 1) * 0.301029995663981195) + 1;   /* a lower bound */
    for (int guard = 0; guard < 8; guard++) {
        if (d < 1) d = 1;
        if (d > TP_DECBUF - 2) d = TP_DECBUF - 2;
        pow10_ensure(d + 1);
        if (d > 1 && bn_cmp(a, &g_pow10[d], nw) < 0) { d--; continue; }
        if (bn_cmp(a, &g_pow10[d], nw) >= 0)        { d++; continue; }
        return d;
    }
    return d;
}

/* Carrying lo's decimal string down the DFS (adding the level's delta in
 * place) is a measured dead end: a hex delta has a DENSE decimal expansion
 * (6*16^30 has 33 nonzero digits of 37), so the add costs a full-width pass --
 * the same order as the conversion it replaces.  Path kept for the record. */
#define TP_R5_DEFAULT 12
static int g_r5J = 0;       /* TP_R5=J: 5-adic reach filter, modulus 5^J       */
static uint64_t g_r5_tests = 0, g_r5_kills = 0;
static int g_r5_count = 0;
static int g_carry9 = 1;   /* carry the node's decimal down (TP_NOCARRY9=1 disables) */
/* Planted-solution mode (TP_OFFSET=c, decimal, c >= 0): search for base-2/16 palindromes N such that N + c is a decimal
 * palindrome. Same statistics as the real problem, but solutions can be planted (pick N0 and a decimal palindrome P0,
 * set c = P0 - N0), giving a full-scale positive control and a way to measure discovery-mode recall. Every decimal-side
 * check compares a residue R implied by the mirrored digits of N + c with N's hex-built part, so it is enough to start
 * the decimal accumulator at c and to compare R - c with N. c = 0 is the real problem. */
static int g_off = 0; static bn_t g_offc;
static uint64_t g_off5N = 0, g_off5 = 0, g_offkN = 0, g_offk = 0;
static uint64_t off_mod(uint64_t m) { bn_t t = g_offc; return bn_divmod_small(&t, TP_MAX_LIMBS, m); }
static pthread_once_t g_off_once = PTHREAD_ONCE_INIT;
static void off_parse(void)
{
    const char *e = getenv("TP_OFFSET");
    if (!e) return;
    bn_t c; bn_zero(&c, TP_MAX_LIMBS);
    if (!bn_from_dec(e, TP_MAX_LIMBS, &c)) { fprintf(stderr, "TP_OFFSET: not a decimal number\n"); exit(2); }
    g_offc = c; g_off = 1;
}
static int g_weakJ = 0;            /* TEST ONLY (TP_WEAKPAL=J): accept N when its last J digits mirror its first J */
static int g_r5chk = 0;           /* TP_R5CHK: sample real completions and assert they are in the reach set */
static uint64_t g_r5chk_n = 0, g_r5chk_bad = 0, g_r5chk_pn = 0, g_r5chk_pbad = 0;
static int g_r5lvl[TP_MAX_BITS];  /* debug copy of lvl[] */
static int g_mod11 = 1;    /* 1: leaf divisibility sieve (11 | palindrome)  */
static int g_dsieve = 0;   /* 1: digit-sum / leading-2-digit sieve (--sieve) */
static int g_suppress_hits = 0;  /* phase-1 hits belong to shard 0 only */
static int g_badprune = 0; /* TEST ONLY: deliberately unsound over-prune     */
static int g_depthhist = 0;
static uint64_t g_dnode[256], g_dprune[256];
static int g_decchk = 0;
static inline __attribute__((always_inline)) int prune_cell(const tp_search *s, int r,
                             const char *dlo, const char *dhi, int dl)
{
    int m = 0;                                 /* shared decimal prefix    */
    while (m < dl && dlo[m] == dhi[m]) m++;
    if (m > (int)s->st.maxm) ((tp_search *)s)->st.maxm = (uint32_t)m;
    g_calls++;
    if (m == 0) { g_m0++; return 0; }

    const int k = m < r ? m : r;

    if (k > 0) {
        int dead;
        if (k <= 64) {                     /* fast path: one limb          */
            uint64_t R = 0;
            for (int i = m - 1; i >= 0; i--) R = R * 10u + (uint64_t)(dlo[i] - '0');
            const uint64_t msk = (k >= 64) ? ~(uint64_t)0 : (((uint64_t)1 << k) - 1);
            dead = ((R & msk) != (s->lo.w[0] & msk));
        } else {                           /* wide path (r > 64)           */
            const int nwk = tp_words(k);
            bn_t R;
            bn_zero(&R, nwk);
            for (int i = m - 1; i >= 0; i--)
                bn_muladd(&R, nwk, 10, (uint64_t)(dlo[i] - '0'));
            dead = 0;
            for (int i = 0; i < nwk - 1; i++)
                if (R.w[i] != s->lo.w[i]) { dead = 1; break; }
            if (!dead) {
                const int rem = k & 63;
                const uint64_t msk = rem ? (((uint64_t)1 << rem) - 1) : ~(uint64_t)0;
                dead = ((R.w[nwk - 1] & msk) != (s->lo.w[nwk - 1] & msk));
            }
        }
        if (dead) return 1;
    }

    return 0;
}

static inline __attribute__((always_inline)) int prune_cell_lo(const tp_search *s,
                                                               int r, const char *dlo, int m)
{
    if (m > (int)s->st.maxm) ((tp_search *)s)->st.maxm = (uint32_t)m;
    const int k = m < r ? m : r;
    if (k > 0) {
        if (k <= 64) {
            uint64_t R = 0;
            for (int i = m - 1; i >= 0; i--) R = R * 10u + (uint64_t)(dlo[i] - '0');
            const uint64_t msk = (k >= 64) ? ~(uint64_t)0 : (((uint64_t)1 << k) - 1);
            if ((R & msk) != (s->lo.w[0] & msk)) return 1;
        } else {
            const int nwk = tp_words(k);
            bn_t R;
            bn_zero(&R, nwk);
            for (int i = m - 1; i >= 0; i--)
                bn_muladd(&R, nwk, 10, (uint64_t)(dlo[i] - '0'));
            int dead = 0;
            for (int i = 0; i < nwk - 1; i++)
                if (R.w[i] != s->lo.w[i]) { dead = 1; break; }
            if (!dead) {
                const int rem = k & 63;
                const uint64_t msk = rem ? (((uint64_t)1 << rem) - 1) : ~(uint64_t)0;
                dead = ((R.w[nwk - 1] & msk) != (s->lo.w[nwk - 1] & msk));
            }
            if (dead) return 1;
        }
    }
    return 0;
}


static __attribute__((noinline)) int prune_straddle(const tp_search *s, int r,
                                                    const char *dlo, int dl,
                                                    const char *dhi, int dh)
{
    char nine[TP_DECBUF], one[TP_DECBUF];
    ((tp_search *)s)->st.straddle++;

    /* lower cell: [lo, 10^dl - 1] */
    memset(nine, '9', (size_t)dl);
    nine[dl] = 0;
    if (!prune_cell(s, r, dlo, nine, dl)) return 0;

    /* upper cell: [10^(dh-1), hi] */
    one[0] = '1';
    memset(one + 1, '0', (size_t)(dh - 1));
    one[dh] = 0;
    if (!prune_cell(s, r, one, dhi, dh)) return 0;

    return 1;                                   /* every cell is empty      */
}

static inline __attribute__((always_inline)) int prune_check(const tp_search *s, int r,
                              const char *dlo, int dl,
                              const char *dhi, int dh)
{
    if (dl <= 0 || dh <= 0) return 0;
    if (dl == dh) return prune_cell(s, r, dlo, dhi, dl);
    if (dh != dl + 1) return 0;                /* impossible for a 2^B range */
    return prune_straddle(s, r, dlo, dl, dhi, dh);
}

static int g_nohi = 1;     /* 1: never convert hi (see hexdfs)           */
static int g_xm = 0;       /* cross-check: recompute m and verdict the long way */
static uint64_t g_xm_m_bad = 0, g_xm_v_bad = 0, g_xm_n = 0, g_xm_fwd = 0, g_xm_back = 0;
static int g_noprune = 0;  /* debug: read once, never in the hot loop   */

#ifdef TP_PROF2
#include <x86intrin.h>
static uint64_t gP_mask, gP_lo, gP_hi, gP_prune, gP_lastconv, gP_pal;
static uint64_t gP_nmask, gP_nprune, gP_nlast;
#define PC() __rdtsc()
#else
#define PC() 0ULL
#endif
void tp_prof_dump(void);

/* ------------------------------------------------------------------ */
/* decimal state in base 10^9 chunks                                   */
/*                                                                     */
/* bn_to_dec costs ~140 ns and the search needs two of them per node;   */
/* a decimal add on a 7-chunk number costs a few ns.  Keeping the       */
/* decimal form incrementally is what makes the difference.             */
/* ------------------------------------------------------------------ */
#define TP_DB  1000000000u
#define TP_DNC (TP_DECBUF / 9 + 2)

typedef struct { uint32_t c[TP_DNC]; int nc; } dec9;   /* c[TP_DNC-1] = least significant */

static const char g_lut2[200] =
    "00010203040506070809" "10111213141516171819"
    "20212223242526272829" "30313233343536373839"
    "40414243444546474849" "50515253545556575859"
    "60616263646566676869" "70717273747576777879"
    "80818283848586878889" "90919293949596979899";

static void dec9_from_bn(const bn_t *x, int nw, dec9 *d)
{
    uint64_t tmp[TP_MAX_BITS / 64 + 2];
    if (nw > TP_MAX_BITS / 64 + 2) nw = TP_MAX_BITS / 64 + 2;
    memcpy(tmp, x->w, sizeof(uint64_t) * (size_t)nw);
    memset(d->c, 0, sizeof d->c);
    d->nc = 0;
    for (;;) {
        uint64_t rem = 0;
        int nz = 0;
        for (int i = nw - 1; i >= 0; i--) {
            const __uint128_t cur = ((__uint128_t)rem << 64) | tmp[i];
            tmp[i] = (uint64_t)(cur / TP_DB);
            rem = (uint64_t)(cur % TP_DB);
        }
        d->c[TP_DNC - 1 - d->nc] = (uint32_t)rem;
        d->nc++;
        for (int i = 0; i < nw; i++) if (tmp[i]) { nz = 1; break; }
        if (!nz || d->nc >= TP_DNC - 1) break;
    }
}

static void dec9_add(dec9 *dst, const dec9 *a, const dec9 *b)
{
    const int nc = a->nc > b->nc ? a->nc : b->nc;
    uint32_t out[TP_DNC];
    uint64_t carry = 0;
    for (int k = 0; k < nc; k++) {                 /* k from the least significant */
        const uint32_t av = (k < a->nc) ? a->c[TP_DNC - 1 - k] : 0;   /* nc may differ */
        const uint32_t bv = (k < b->nc) ? b->c[TP_DNC - 1 - k] : 0;
        const uint64_t s = (uint64_t)av + bv + carry;
        carry = (s >= TP_DB);
        out[TP_DNC - 1 - k] = (uint32_t)(carry ? s - TP_DB : s);
    }
    int total = nc;
    if (carry && nc < TP_DNC) { out[TP_DNC - 1 - nc] = 1; total = nc + 1; }
    for (int k = 0; k < total; k++) dst->c[TP_DNC - 1 - k] = out[TP_DNC - 1 - k];
    dst->nc = total;
}

/* digits of one chunk, most significant first, exactly nd of them */
/* number of leading decimal digits on which a and b agree; -1 if their digit
   counts differ (the caller then has no shared prefix) */
/* the first m digits of a, most significant first (no NUL) */
/* dl nines, and 10^dl */
static int dec9_to_str(const dec9 *d, char *out)
{
    const int base = TP_DNC - d->nc;
    int n = 0;
    {   /* top chunk: no leading zeros */
        uint32_t t = d->c[base];
        char buf[12];
        int k = 0;
        do { buf[k++] = (char)('0' + t % 10); t /= 10; } while (t);
        while (k) out[n++] = buf[--k];
    }
    for (int j = 1; j < d->nc; j++) {
        const uint32_t v = d->c[base + j];
        const int q0 = (int)(v / 10000000u),
                  q1 = (int)((v / 100000u) % 100u),
                  q2 = (int)((v / 1000u) % 100u),
                  q3 = (int)((v / 10u) % 100u),
                  q4 = (int)(v % 10u);
        out[n++] = g_lut2[2*q0]; out[n++] = g_lut2[2*q0+1];
        out[n++] = g_lut2[2*q1]; out[n++] = g_lut2[2*q1+1];
        out[n++] = g_lut2[2*q2]; out[n++] = g_lut2[2*q2+1];
        out[n++] = g_lut2[2*q3]; out[n++] = g_lut2[2*q3+1];
        out[n++] = (char)('0' + q4);
    }
    out[n] = 0;
    return n;
}

/* ---- chunk-native node kernel (no decimal string is ever built) ----------
 * A node's decimal value is carried as base-10^9 chunks.  Everything the
 * pruner reads is a handful of digits near the top, one comparison of the low
 * k digits against a plan-time threshold, and the reversed prefix, which is
 * carried as a multi-limb accumulator (the shared prefix only extends as the
 * tree descends, so the child adds just the newly revealed digits). */
static uint64_t g_fastchk_n = 0;
static int g_fastchk = 0;                        /* TP_FASTCHK: cross-check the kernel */
static int g_fastk = 1;                          /* TP_NOFASTK=1 disables       */
static uint64_t g_p10w[TP_DECBUF][TP_MAX_LIMBS]; /* 10^i mod 2^(64*TP_MAX_LIMBS) */
static const uint32_t g_p10_32[10] = {1u,10u,100u,1000u,10000u,100000u,1000000u,
                                      10000000u,100000000u,1000000000u};
static pthread_once_t g_p10_once = PTHREAD_ONCE_INIT;
static void p10w_init(void)
{
    memset(g_p10w, 0, sizeof g_p10w);
    g_p10w[0][0] = 1;
    for (int i = 1; i < TP_DECBUF; i++) {
        uint64_t cy = 0;
        for (int w = 0; w < TP_MAX_LIMBS; w++) {
            const __uint128_t t = (__uint128_t)g_p10w[i - 1][w] * 10u + cy;
            g_p10w[i][w] = (uint64_t)t; cy = (uint64_t)(t >> 64);
        }
    }
}
static inline int ilog10_u32(uint32_t v)
{
    const int t = ((32 - __builtin_clz(v | 1u)) * 1233) >> 12;
    return t - (v < g_p10_32[t]);
}
/* i-th decimal digit from the top (i = 0 is the leading digit); td = digits of the top chunk */
static inline int dec9_digit(const dec9 *x, int td, int i)
{
    const int pos = i + (9 - td);
    const uint32_t v = x->c[(TP_DNC - x->nc) + pos / 9];
    return (int)((v / g_p10_32[8 - pos % 9]) % 10u);
}
/* Rc += (digits [from,to) of x, read as a number with digit `from` least significant) * 10^from */
static inline void rev_add_digits(uint64_t *Rc, int nw, const dec9 *x, int td, int from, int to)
{
    while (from < to) {
        int e = from + 19; if (e > to) e = to;
        uint64_t t = 0;
        for (int i = e - 1; i >= from; i--) t = t * 10u + (uint64_t)dec9_digit(x, td, i);
        const uint64_t *pw = g_p10w[from];
        uint64_t cy = 0;
        for (int w = 0; w < nw; w++) {
            const __uint128_t q = (__uint128_t)t * pw[w] + Rc[w] + cy;
            Rc[w] = (uint64_t)q; cy = (uint64_t)(q >> 64);
        }
        from = e;
    }
}


/* ------------------------------------------------------------------ */
/* the same pruner, reading the base-10^9 chunk form directly          */
/*                                                                     */
/* bn_to_dec costs ~96 ns and a node needs two of them.  These         */
/* routines consume the incremental decimal state instead, so the       */
/* search performs no radix conversion at all.  They are exact          */
/* transliterations of prune_cell / prune_straddle / prune_check;       */
/* TP_XCHECK compares the verdicts call by call.                        */
/* ------------------------------------------------------------------ */
typedef struct hexplan hexplan;
typedef struct { bn_t lo, hi; int depth; int s_sum; int s_outer; uint64_t nodes; } hxtask;
#define TP_MAXLVL (TP_MAX_BITS / 4 + 2)

struct hexplan {
    int L;                 /* number of nibbles                        */
    int depth_max;         /* (L+1)/2 -- nibbles 0..depth_max-1 decide */
    /* optional parallel cut: when d == cut, record the state and stop */
    int             cut;
    hxtask         *tasks;
    int            *ntask;
    int             cap;
    int nb[TP_MAXLVL][4];             /* sites in nibble d             */
    int nn[TP_MAXLVL];                /* how many (<= 4)               */
    int trace;
    /* Viable digits per level.  Whether a digit survives depends only on the
       level, the digit and the FORCED components, none of which change during
       the search, so the 16-way digit scan that used to run at every node
       collapses to the handful that actually survive -- and their whole
       effect on the state is two masks.                                   */
    int       nviab[TP_MAXLVL];
    int       vtab[TP_MAXLVL][16];
    uint64_t *setm, *clrm;            /* [level][16][nw] words         */
    int       ddig[64];               /* [level]: decimal digits of (hi-lo)   */
    int       add11[64][16];          /* [level][digit]: N += this, mod 11     */
    char     *dstr;                   /* [level] 10^k - (hi-lo), k=ddig, k digits */
    /* R5 filter: the 5-adic reach of the not-yet-decided orbits */
    int       r5J;                    /* 0 = off; else modulus is 5^r5J       */
    uint64_t  r5N;                    /* 5^r5J                                */
    uint64_t  add5[64][16];           /* [level][digit]: N += this, mod 5^J    */
    const uint64_t *reach5[64];       /* [level]: bitset of residues reachable by
                                         the orbits decided AFTER this level    */
    uint32_t  thr[64][TP_DNC];        /* [level]: 10^k - width as low-first base-1e9 chunks */
    uint32_t  thrpart[64];            /* [level]: the k%9 leading digits of it */
    int       bK;                     /* deep filter: 0 = off, else digits used  */
    uint64_t  bN;                     /* 5^bK                                    */
    uint64_t  addk[64][16];           /* [level][digit]: N += this, mod 5^bK     */
    int       bat[64];                /* [level]: filter index + 1, 0 = none     */
    const void *bfp;                  /* the filters (bfilter[TP_BF_MAX])        */
    const struct duallev_s *dual[64]; /* [level]: dual-lattice certificates      */
    const struct cvplev_s *cvp[64];   /* [level]: exact enumeration test         */
    dec9     *dlo9;                   /* [level][16]: decimal of the bits
                                         this level sets in lo / clears
                                         from hi, first-assignment only */
};

static void hexplan_init(hexplan *p, const tp_orbits *o)
{
    p->L = (o->nsites + 3) / 4;
    if (p->L > TP_MAX_BITS / 4) p->L = TP_MAX_BITS / 4;
    p->depth_max = (p->L + 1) / 2;
    for (int d = 0; d < p->L; d++) {
        int n = 0;
        for (int k = 0; k < 4; k++) {
            int s = 4 * d + k;
            if (s < o->nsites) p->nb[d][n++] = s;
        }
        p->nn[d] = n;
    }
}

/* Collect the distinct (orbit, value) implications of setting nibble d to v.
 * Returns -1 if v is internally inconsistent, else the count. */
static int hex_implications(const tp_orbits *o, const hexplan *p, int d, int v,
                            int *orb, int *val)
{
    int n = 0;
    for (int i = 0; i < p->nn[d]; i++) {
        const int s = p->nb[d][i];
        const int c = o->comp[s];
        const int k = (v >> (s & 3)) & 1;
        int found = 0;
        for (int j = 0; j < n; j++)
            if (orb[j] == c) {
                if (val[j] != k) return -1;
                found = 1;
                break;
            }
        if (!found) { orb[n] = c; val[n] = k; n++; }
    }
    return n;
}

/* level at which component c is first assigned: min over its sites of
   min(nib, L-1-nib).  Level d is entered on every surviving path. */
static void hexplan_levels(const hexplan *p, const tp_orbits *o, int *lvl)
{
    for (int c = 0; c < o->ncomp; c++) lvl[c] = -1;
    for (int s = 0; s < o->nsites; s++) {
        const int nb = s / 4, m = p->L - 1 - nb;
        const int l  = nb < m ? nb : m;
        const int c  = o->comp[s];
        if (lvl[c] < 0 || l < lvl[c]) lvl[c] = l;
    }
}

static void hexplan_free(hexplan *p);

/* ---- R5: 5-adic reach of the undecided orbits ---------------------------
 * N = lo + sum_{c undecided} x_c * M_c with each x_c free in {0,1}, so the set of
 * residues N mod 5^J can still reach is  lo + SubsetSums({M_c mod 5^J}).  A
 * palindrome needs N mod 10^j to equal the reversal of N's top j digits, and
 * N mod 5^j is therefore known as soon as j leading digits are shared by lo
 * and hi.  If (target - lo) is not a subset sum of the undecided orbits, no
 * completion of this node can be a palindrome.  The subset-sum sets are built
 * once per bit length, from the centre outwards, as bitsets over Z/5^J.      */
static pthread_mutex_t g_r5_mu = PTHREAD_MUTEX_INITIALIZER;
static struct { int key_n, key_J; uint64_t key_h; int D; uint64_t *tab[64]; } g_r5c = { -1, -1, 0, 0, {0} };

/* the table depends on the orbit structure, not only on nsites: several bit
 * lengths share the same padded nsites but differ in their forced orbits */
static uint64_t r5_orbit_hash(const tp_orbits *o)
{
    uint64_t h = 1469598103934665603ull;
    for (int c = 0; c < o->ncomp; c++) h = (h ^ (uint64_t)(o->cval[c] + 2)) * 1099511628211ull;
    for (size_t i = 0; i < (size_t)o->ncomp * (size_t)o->nw64; i++) h = (h ^ o->mask[i]) * 1099511628211ull;
    return h ^ (uint64_t)o->ncomp ^ ((uint64_t)o->nsites << 32);
}

/* cyclic rotate-left of an N-bit bitset by M (0 <= M < N): dst[i] = src[(i - M) mod N] */
static void bs_shl(uint64_t *dst, const uint64_t *src, size_t W, uint64_t M)
{   /* dst = src << M, truncated to W words, no wrap */
    const size_t ws = (size_t)(M >> 6); const unsigned bs = (unsigned)(M & 63);
    for (size_t i = W; i-- > 0; ) {
        uint64_t v = 0;
        if (i >= ws) {
            v = src[i - ws] << bs;
            if (bs && i > ws) v |= src[i - ws - 1] >> (64 - bs);
        }
        dst[i] = v;
    }
}
static void bs_shr_or(uint64_t *dst, const uint64_t *src, size_t W, uint64_t S)
{   /* dst |= src >> S */
    const size_t ws = (size_t)(S >> 6); const unsigned bs = (unsigned)(S & 63);
    for (size_t i = 0; i + ws < W; i++) {
        uint64_t v = src[i + ws] >> bs;
        if (bs && i + ws + 1 < W) v |= src[i + ws + 1] << (64 - bs);
        dst[i] |= v;
    }
}
static void bs_rotl_or(uint64_t *cur, uint64_t *tmp, size_t W, uint64_t N, uint64_t M)
{   /* cur |= rotl(cur, M) over N bits */
    bs_shl(tmp, cur, W, M);                    /* bits shifted up (may spill past N) */
    bs_shr_or(tmp, cur, W, N - M);             /* wrap-around part                   */
    if (N & 63) tmp[W - 1] &= (1ull << (N & 63)) - 1;
    for (size_t i = 0; i < W; i++) cur[i] |= tmp[i];
}
static uint64_t bs_popcount(const uint64_t *a, size_t W)
{ uint64_t c = 0; for (size_t i = 0; i < W; i++) c += (uint64_t)__builtin_popcountll(a[i]); return c; }

static void r5_build(hexplan *p, const tp_orbits *o, const int *lvl)
{
    pthread_mutex_lock(&g_r5_mu);
    const int D = p->depth_max;
    const uint64_t oh = r5_orbit_hash(o);
    if (g_r5c.key_n != o->nsites || g_r5c.key_J != p->r5J || g_r5c.key_h != oh) {
        for (int d = 0; d < 64; d++) { free(g_r5c.tab[d]); g_r5c.tab[d] = NULL; }
        const uint64_t N = p->r5N;
        const size_t W = (size_t)((N + 63) / 64);
        uint64_t *cur = calloc(W, 8), *tmp = malloc(W * 8);
        if (!cur || !tmp) { fprintf(stderr, "R5: out of memory\n"); exit(2); }
        cur[0] = 1;
        int full = 0;
        for (int d = D - 1; d >= 0; d--) {
            /* tab[d] = subset sums of the orbits with lvl > d; NULL once that
               set is every residue (the test is then vacuous)                */
            if (!full) {
                uint64_t *t = malloc(W * 8);
                if (!t) { fprintf(stderr, "R5: out of memory\n"); exit(2); }
                memcpy(t, cur, W * 8);
                g_r5c.tab[d] = t;
                if (getenv("TP_R5DBG")) {
                    int u = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0 && lvl[c] > d) u++;
                    fprintf(stderr, "  [r5] d=%2d undecided=%2d density=%.4f\n", d, u, (double)bs_popcount(cur, W) / (double)N);
                }
            }
            if (d == 0) break;
            if (!full) {
                for (int c = 0; c < o->ncomp; c++) {
                    if (o->cval[c] >= 0 || lvl[c] != d) continue;
                    bn_t b; bn_zero(&b, o->nw64);
                    memcpy(b.w, o->mask + (size_t)c * (size_t)o->nw64, (size_t)o->nw64 * sizeof(uint64_t));
                    const uint64_t M = bn_divmod_small(&b, o->nw64, N);
                    if (M) bs_rotl_or(cur, tmp, W, N, M);
                }
                if (bs_popcount(cur, W) == N) full = 1;
            }
        }
        free(cur); free(tmp);
        g_r5c.key_n = o->nsites; g_r5c.key_J = p->r5J; g_r5c.key_h = oh; g_r5c.D = D;
    }
    for (int d = 0; d < D && d < 64; d++) p->reach5[d] = g_r5c.tab[d];
    pthread_mutex_unlock(&g_r5_mu);
}


/* ---- deep reach filter ("--bloom U:bits[,U:bits]") --------------------------
 * The exact reach bitsets cost 5^J bits, so J stops at 13.  Near the bottom of the
 * tree far more leading digits are shared than that, and the undecided orbits are
 * few enough to list.  At a chosen level the set S = { sum of a subset of the
 * orbits decided after this level, mod 5^K } has exactly 2^U elements; storing it
 * as a split-block Bloom filter costs a few bits per element and is independent
 * of the modulus, so K can be 25.  A Bloom filter has false positives, never false
 * negatives, so a child rejected here cannot be a palindrome. */
#define TP_BF_MAX 3
typedef struct { uint32_t *w; uint64_t nb; int d, u; } bfilter;
static struct { int key_n, K; uint64_t key_h; char spec[64]; int n; bfilter f[TP_BF_MAX]; }
    g_bfc = { -1, 0, 0, "", 0, {{0,0,0,0}} };
static pthread_mutex_t g_bf_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_bchk = 0;      /* TP_BLOOMCHK: verify carried residue/target and sampled members */
static uint64_t g_bchk_n = 0, g_bchk_t = 0;
static const uint32_t g_bfs[8] = { 0x47b6137bU, 0x44974d91U, 0x8824ad5bU, 0xa2b7289dU,
                                   0x705495c7U, 0x2df1424bU, 0x9efc4947U, 0x5c6bfb31U };
static inline uint64_t bf_mix(uint64_t x)
{
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return x;
}
static inline uint32_t *bf_block(const bfilter *f, uint64_t h)
{ return f->w + 8 * (size_t)(((h >> 32) * f->nb) >> 32); }
static inline int bf_has(const bfilter *f, uint64_t x)
{
    const uint64_t h = bf_mix(x);
    const uint32_t *b = bf_block(f, h);
    for (int i = 0; i < 8; i++)
        if (!((b[i] >> (((uint32_t)h * g_bfs[i]) >> 27)) & 1u)) return 0;
    return 1;
}
static inline void bf_add_h(bfilter *f, uint64_t h)
{
    uint32_t *b = bf_block(f, h);
    for (int i = 0; i < 8; i++)
        __atomic_fetch_or(&b[i], 1u << (((uint32_t)h * g_bfs[i]) >> 27), __ATOMIC_RELAXED);
}
typedef struct { bfilter *f; const uint64_t *M; int u, p; uint64_t N; uint64_t next; } bf_job;
static void *bf_worker(void *arg)
{
    bf_job *J = arg;
    const int q = J->u - J->p;
    enum { RING = 32 };
    uint64_t ring[RING]; int fill = 0, pos = 0;
    for (;;) {
        const uint64_t t = __atomic_fetch_add(&J->next, 1, __ATOMIC_RELAXED);
        if (t >> J->p) break;
        uint64_t x = 0, gray = 0;
        for (int i = 0; i < J->p; i++)
            if ((t >> i) & 1) { x += J->M[q + i]; if (x >= J->N) x -= J->N; }
        const uint64_t n = (uint64_t)1 << q;
        for (uint64_t i = 0; i < n; i++) {
            const uint64_t h = bf_mix(x);
            __builtin_prefetch(bf_block(J->f, h), 1, 0);
            if (fill == RING) bf_add_h(J->f, ring[pos]); else fill++;
            ring[pos] = h; pos = (pos + 1) % RING;
            if (i + 1 < n) {
                const int b = __builtin_ctzll(i + 1);
                const uint64_t m = J->M[b];
                if ((gray >> b) & 1) x = (x >= m) ? x - m : x + J->N - m;
                else { x += m; if (x >= J->N) x -= J->N; }
                gray ^= (uint64_t)1 << b;
            }
        }
    }
    for (int i = 0; i < fill; i++) bf_add_h(J->f, ring[i]);
    return NULL;
}

/* ---- dual-lattice reachability certificate ("--dual") ------------------------
 * For ANY integer lam, put v_c = lam*M_c mod+- N.  A subset sum q of the M_c then
 * satisfies lam*q = sum x_c v_c (mod N) with the right side an integer in
 * [neg, pos], neg/pos the sums of the negative/positive v_c.  If pos-neg < N, a
 * target whose lam*q falls outside that window is PROVABLY unreachable: this is
 * arithmetic, there are no false negatives and no table.  Short (v_c) are short
 * vectors of a lattice, found once per level by LLL.  N = 5^K fits 256 bits. */
#define TP_DUAL_TMAX 128
#define TP_DUAL_KMAX 109
#define TP_SPEC_MAX 256
typedef struct { uint64_t w[4]; } u256;
typedef struct duallev_s {
    int K, nt, u, nl; u256 N; uint64_t n0; u256 c64, c10;
    u256 add[16];
    u256 lamM[TP_DUAL_TMAX], negabs[TP_DUAL_TMAX], span[TP_DUAL_TMAX];
    int ns; double sNd, sB; u256 slam[TP_SPEC_MAX]; double scen[TP_SPEC_MAX], swt[TP_SPEC_MAX];
    int nadd; u256 ladd[TP_DUAL_TMAX][16];
    u256 coff;                                              /* TP_OFFSET: c mod N */                  /* lam_i * add[k] mod N, for sibling sharing */   /* Fourier certificate */
} duallev;
static struct { int key_n; uint64_t key_h; int K; duallev *lev[64]; } g_duc = { -1, 0, 0, {0} };
static pthread_mutex_t g_du_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_duchk = 0;
static int g_noshare = 0;                   /* TP_NOSHARE: recompute every product per child (reference path) */
static inline int u256_ge(const u256 *a, const u256 *b)
{ for (int i = 3; i >= 0; i--) { if (a->w[i] != b->w[i]) return a->w[i] > b->w[i]; } return 1; }
static inline void u256_sub(u256 *r, const u256 *a, const u256 *b)
{ uint64_t br = 0; for (int i = 0; i < 4; i++) { const __uint128_t t = (__uint128_t)a->w[i] - b->w[i] - br; r->w[i] = (uint64_t)t; br = (uint64_t)(t >> 64) & 1; } }
static inline void u256_add(u256 *r, const u256 *a, const u256 *b)
{ uint64_t c = 0; for (int i = 0; i < 4; i++) { const __uint128_t t = (__uint128_t)a->w[i] + b->w[i] + c; r->w[i] = (uint64_t)t; c = (uint64_t)(t >> 64); } }
static inline void u256_addmod(u256 *r, const u256 *a, const u256 *b, const u256 *N)
{ u256 t; u256_add(&t, a, b); if (u256_ge(&t, N)) u256_sub(r, &t, N); else *r = t; }
static inline void u256_submod(u256 *r, const u256 *a, const u256 *b, const u256 *N)
{ if (u256_ge(a, b)) u256_sub(r, a, b); else { u256 t; u256_add(&t, a, N); u256_sub(r, &t, b); } }
static inline void u256_mont(u256 *r, const u256 *a, const u256 *b, const u256 *N, uint64_t n0, int nl)
{
    uint64_t t[6] = {0,0,0,0,0,0};
    for (int i = 0; i < nl; i++) {
        __uint128_t c = 0;
        for (int j = 0; j < nl; j++) { c += (__uint128_t)a->w[j] * b->w[i] + t[j]; t[j] = (uint64_t)c; c >>= 64; }
        c += t[nl]; t[nl] = (uint64_t)c; t[nl + 1] = (uint64_t)(c >> 64);
        const uint64_t m = t[0] * n0;
        c = (__uint128_t)m * N->w[0] + t[0]; c >>= 64;
        for (int j = 1; j < nl; j++) { c += (__uint128_t)m * N->w[j] + t[j]; t[j - 1] = (uint64_t)c; c >>= 64; }
        c += t[nl]; t[nl - 1] = (uint64_t)c; c >>= 64; t[nl] = t[nl + 1] + (uint64_t)c;
    }
    u256 z = {{0,0,0,0}}; for (int i = 0; i < nl; i++) z.w[i] = t[i];
    if (t[nl] || u256_ge(&z, N)) u256_sub(r, &z, N); else *r = z;
}
static inline void mod_reduce(u256 *out, const uint64_t *x, int nw, const u256 *N, uint64_t n0, const u256 *c64, int nl)
{
    u256 a = {{0,0,0,0}}, t;
    const int small = !(N->w[1] | N->w[2] | N->w[3]);
    for (int i = nw - 1; i >= 0; i--) {
        u256_mont(&t, &a, c64, N, n0, nl);
        u256 lv = {{ small ? x[i] % N->w[0] : x[i], 0, 0, 0 }};
        u256_addmod(&a, &t, &lv, N);
    }
    *out = a;
}
/* window test for child k of a parent with residue base, sharing lam_i * base between siblings (exact) */
static inline int dual_pass_child(const duallev *L, const u256 *base, int k, u256 *zb, uint8_t *zok)
{
    for (int i = 0; i < L->nt; i++) {
        if (!zok[i]) { u256_mont(&zb[i], &L->lamM[i], base, &L->N, L->n0, L->nl); zok[i] = 1; }
        u256 z, a;
        u256_submod(&z, &zb[i], &L->ladd[i][k], &L->N);
        u256_addmod(&a, &z, &L->negabs[i], &L->N);
        if (!u256_ge(&L->span[i], &a)) return 0;
    }
    return 1;
}
static inline void dual_reduce(u256 *out, const uint64_t *x, int nw, const duallev *L)
{ mod_reduce(out, x, nw, &L->N, L->n0, &L->c64, L->nl); }
static inline int dual_pass(const duallev *L, const u256 *q)
{
    for (int i = 0; i < L->nt; i++) {
        u256 z, a;
        u256_mont(&z, &L->lamM[i], q, &L->N, L->n0, L->nl);
        u256_addmod(&a, &z, &L->negabs[i], &L->N);
        if (!u256_ge(&L->span[i], &a)) return 0;
    }
    if (L->ns) {                        /* Fourier certificate: sum (1 - cos 2 pi y_w) <= sB for every reachable q */
        double def = 0;
        for (int i = 0; i < L->ns; i++) {
            u256 z; u256_mont(&z, &L->slam[i], q, &L->N, L->n0, L->nl);
            const double v = (ldexp((double)z.w[3], 192) + ldexp((double)z.w[2], 128) + ldexp((double)z.w[1], 64) + (double)z.w[0]) / L->sNd;
            double y = v - L->scen[i]; y -= rint(y);                 /* 2 sin^2: no cancellation, safe for huge weights */
            const double sn = sin(3.141592653589793 * y);
            def += L->swt[i] * 2.0 * sn * sn;
            if (def > L->sB) return 0;
        }
    }
    return 1;
}
#ifdef TP_FLINT
static void f2u(u256 *r, const fmpz_t x)
{
    fmpz_t t, l; fmpz_init_set(t, x); fmpz_init(l);
    for (int i = 0; i < 4; i++) { fmpz_fdiv_r_2exp(l, t, 64); r->w[i] = fmpz_get_ui(l); fmpz_fdiv_q_2exp(t, t, 64); }
    fmpz_clear(t); fmpz_clear(l);
}
static void limbs2f(fmpz_t r, const uint64_t *w, int nw)
{
    fmpz_zero(r);
    for (int i = nw - 1; i >= 0; i--) { fmpz_mul_2exp(r, r, 64); fmpz_add_ui(r, r, w[i]); }
}
/* ---- learned choice of windows (TP_DUAL_POOL) ---------------------------------
 * Every integer combination of LLL rows is a valid window, so the choice can be made
 * freely without any risk to soundness: spans and lambdas are recomputed exactly.
 * A pool is grown by L1-targeted pairwise reduction, then windows are picked greedily
 * (boosting) by how many still-alive random residues each one newly rejects. */
static int cmp_dpi(const void *a, const void *b)
{ const double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
static uint64_t dl_rng(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; return *s = x; }
static void dual_learn(duallev *L, const fmpz_mat_t A, const fmpz_t N, const fmpz_t R, int d, int dbg)
{
    const int P = atoi(getenv("TP_DUAL_POOL")), u = L->u;
    int n = 0, rowid[257];
    const int rounds = getenv("TP_DUAL_ROUNDS") ? atoi(getenv("TP_DUAL_ROUNDS")) : 10;
    const int T = getenv("TP_DUAL_REDUCERS") ? atoi(getenv("TP_DUAL_REDUCERS")) : 100;
    const int Cmax = getenv("TP_DUAL_CAND") ? atoi(getenv("TP_DUAL_CAND")) : 1500;
    const int S = getenv("TP_DUAL_S") ? atoi(getenv("TP_DUAL_S")) : 2048, H = 8192;
    const int minrej = getenv("TP_DUAL_MINREJ") ? atoi(getenv("TP_DUAL_MINREJ")) : 4;
    if (P < u + 1 || u + 1 > 257) return;
    fmpz_t t, neg, pos, lam, sp; fmpz_init(t); fmpz_init(neg); fmpz_init(pos); fmpz_init(lam); fmpz_init(sp);
    double *Vb = malloc(sizeof(double) * (size_t)(u + 1) * (size_t)u), *V = malloc(sizeof(double) * (size_t)P * (size_t)u), *l1 = malloc(sizeof(double) * (size_t)P);
    double *ord = malloc(sizeof(double) * 2 * (size_t)P), *W = malloc(sizeof(double) * (size_t)u);
    int *C = calloc((size_t)P * (size_t)(u + 1), sizeof(int));
    if (!Vb || !V || !l1 || !ord || !W || !C) { fprintf(stderr, "dual: out of memory\n"); exit(2); }
    const double Nd = fmpz_get_d(N);
    for (int i = 0; i <= u; i++) {                /* usable rows: span not absurdly above N */
        int any = 0; double l1r = 0;
        for (int j = 0; j < u; j++) {
            fmpz_divexact(t, fmpz_mat_entry(A, i, j + 1), N);
            if (!fmpz_is_zero(t)) any = 1;
            const double f = fmpz_get_d(t) / Nd;     /* selection only: spans are recomputed exactly */
            Vb[(size_t)n * u + j] = f; l1r += fabs(f);
        }
        if (any && l1r < 64.0) rowid[n++] = i;
    }
    if (n < 2) goto done;
    uint64_t rng = 0xD1B54A32D192ED03ULL ^ ((uint64_t)d << 32) ^ (uint64_t)u; for (int i = 0; i < 8; i++) dl_rng(&rng);
    for (int i = 0; i < P; i++) {
        int *c = C + (size_t)i * n; double *v = V + (size_t)i * u; memset(v, 0, sizeof(double) * (size_t)u);
        if (i < n) c[i] = 1;
        else { const int w = 2 + (int)(dl_rng(&rng) % 5); for (int k = 0; k < w; k++) c[dl_rng(&rng) % (uint64_t)n] += (dl_rng(&rng) & 1) ? 1 : -1; }
        for (int b = 0; b < n; b++) if (c[b]) for (int j = 0; j < u; j++) v[j] += c[b] * Vb[(size_t)b * u + j];
    }
    for (int r = 0; r < rounds; r++) {
        for (int i = 0; i < P; i++) { double a = 0; const double *v = V + (size_t)i * u; for (int j = 0; j < u; j++) a += fabs(v[j]); l1[i] = a; ord[2 * i] = (a > 0) ? a : 1e300; ord[2 * i + 1] = (double)i; }
        qsort(ord, (size_t)P, 2 * sizeof(double), cmp_dpi);
        int changed = 0;
        for (int ti = 0; ti < T && ti < P; ti++) {
            const int pi = (int)ord[2 * ti + 1]; const double *pv = V + (size_t)pi * u; const int *pc = C + (size_t)pi * n;
            for (int i = 0; i < P; i++) {
                if (i == pi) continue;
                double *v = V + (size_t)i * u; double ip = 0, b1 = 0, b2 = 0;
                for (int j = 0; j < u; j++) ip += v[j] * pv[j];
                const double sg = ip >= 0 ? 1 : -1;
                for (int j = 0; j < u; j++) { W[j] = v[j] - sg * pv[j]; b1 += fabs(W[j]); b2 += fabs(W[j]) > 0; }
                if (b2 > 0 && b1 < l1[i] * 0.9995) {
                    memcpy(v, W, sizeof(double) * (size_t)u); l1[i] = b1; int *c = C + (size_t)i * n;
                    for (int b = 0; b < n; b++) c[b] -= (int)sg * pc[b];
                    changed++;
                }
            }
        }
        if (!changed) break;
    }
    /* candidates: smallest L1, deduplicated up to sign, recomputed exactly */
    for (int i = 0; i < P; i++) { double a = 0; const double *v = V + (size_t)i * u; for (int j = 0; j < u; j++) a += fabs(v[j]); ord[2 * i] = (a > 0) ? a : 1e300; ord[2 * i + 1] = (double)i; }
    qsort(ord, (size_t)P, 2 * sizeof(double), cmp_dpi);
    u256 *clm = malloc(sizeof(u256) * (size_t)Cmax), *cna = malloc(sizeof(u256) * (size_t)Cmax), *csp = malloc(sizeof(u256) * (size_t)Cmax);
    uint64_t *chs = malloc(sizeof(uint64_t) * (size_t)Cmax);
    if (!clm || !cna || !csp || !chs) { fprintf(stderr, "dual: out of memory\n"); exit(2); }
    int nc = 0;
    fmpz *vz = _fmpz_vec_init(u);
    for (int oi = 0; oi < P && nc < Cmax; oi++) {
        if (ord[2 * oi] >= 1e299) break;
        const int *c = C + (size_t)((int)ord[2 * oi + 1]) * n;
        fmpz_zero(lam); _fmpz_vec_zero(vz, u);
        for (int b = 0; b < n; b++) { if (!c[b]) continue;
            const int ro = rowid[b]; if (c[b] > 0) fmpz_addmul_ui(lam, fmpz_mat_entry(A, ro, 0), (ulong)c[b]); else fmpz_submul_ui(lam, fmpz_mat_entry(A, ro, 0), (ulong)(-c[b]));
            for (int j = 0; j < u; j++) { fmpz_divexact(t, fmpz_mat_entry(A, ro, j + 1), N);
                if (c[b] > 0) fmpz_addmul_ui(vz + j, t, (ulong)c[b]); else fmpz_submul_ui(vz + j, t, (ulong)(-c[b])); } }
        fmpz_zero(neg); fmpz_zero(pos);
        uint64_t hs = 0x9E3779B97F4A7C15ULL; int sgn = 0;
        for (int j = 0; j < u; j++) { const int z = fmpz_sgn(vz + j); if (!sgn && z) sgn = z; if (z < 0) fmpz_add(neg, neg, vz + j); else if (z > 0) fmpz_add(pos, pos, vz + j);
            { const uint64_t mp = 2305843009213693951ULL; uint64_t r = fmpz_fdiv_ui(vz + j, mp); if (sgn < 0 && r) r = mp - r; hs = (hs ^ r) * 0x100000001B3ULL; } }
        if (!sgn) continue;
        fmpz_sub(sp, pos, neg);
        if (fmpz_cmp(sp, N) >= 0) continue;
        int dup = 0; for (int q = 0; q < nc; q++) if (chs[q] == hs) { dup = 1; break; }
        if (dup) continue;
        fmpz_mod(lam, lam, N); fmpz_mul(t, lam, R); fmpz_mod(t, t, N);
        chs[nc] = hs; f2u(&clm[nc], t); fmpz_neg(t, neg); fmpz_mod(t, t, N); f2u(&cna[nc], t); f2u(&csp[nc], sp); nc++;
    }
    _fmpz_vec_clear(vz, u);
    /* boosting: greedy cover of random residues */
    const int ws = (S + 63) / 64;
    uint64_t *rej = calloc((size_t)(nc ? nc : 1) * (size_t)ws, sizeof(uint64_t)), *alive = malloc(sizeof(uint64_t) * (size_t)ws);
    char *used = calloc((size_t)(nc ? nc : 1), 1);
    u256 *q = malloc(sizeof(u256) * (size_t)(S + H));
    if (!rej || !alive || !used || !q) { fprintf(stderr, "dual: out of memory\n"); exit(2); }
    for (int a = 0; a < S + H; a++) { uint64_t x[4]; for (int b = 0; b < 4; b++) x[b] = dl_rng(&rng); dual_reduce(&q[a], x, 4, L); }
    int hp0 = 0; for (int a = S; a < S + H; a++) hp0 += dual_pass(L, &q[a]);
    for (int k = 0; k < nc; k++) for (int a = 0; a < S; a++) {
        u256 z, w; u256_mont(&z, &clm[k], &q[a], &L->N, L->n0, L->nl); u256_addmod(&w, &z, &cna[k], &L->N);
        if (!u256_ge(&csp[k], &w)) rej[(size_t)k * ws + a / 64] |= 1ULL << (a % 64);
    }
    memset(alive, 0, sizeof(uint64_t) * (size_t)ws); for (int a = 0; a < S; a++) alive[a / 64] |= 1ULL << (a % 64);
    const int oldnt = L->nt; int nt = 0;
    while (nt < TP_DUAL_TMAX) {
        int bk = -1, bg = minrej - 1;
        for (int k = 0; k < nc; k++) { if (used[k]) continue; int g = 0; const uint64_t *rk = rej + (size_t)k * ws; for (int w = 0; w < ws; w++) g += __builtin_popcountll(rk[w] & alive[w]); if (g > bg) { bg = g; bk = k; } }
        if (bk < 0) break;
        used[bk] = 1; L->lamM[nt] = clm[bk]; L->negabs[nt] = cna[bk]; L->span[nt] = csp[bk]; nt++;
        const uint64_t *rk = rej + (size_t)bk * ws; for (int w = 0; w < ws; w++) alive[w] &= ~rk[w];
    }
    if (nt > 0) { L->nt = nt; }
    { int hp1 = 0; for (int a = S; a < S + H; a++) hp1 += dual_pass(L, &q[a]);
      if (dbg) fprintf(stderr, "  [learn] d=%d u=%d K=%d cands=%d windows %d->%d holdout pass %.4f->%.4f\n", d, u, L->K, nc, oldnt, L->nt, (double)hp0 / H, (double)hp1 / H); }
    free(rej); free(alive); free(used); free(q); free(clm); free(cna); free(csp); free(chs);
done:
    free(Vb); free(V); free(l1); free(ord); free(W); free(C);
    fmpz_clear(t); fmpz_clear(neg); fmpz_clear(pos); fmpz_clear(lam); fmpz_clear(sp);
}
/* ---- Fourier certificate (TP_SPEC), after Aharonov-Regev ------------------------
 * For a reachable q there is x in {0,1}^u with q = sum x_c M_c mod N.  Put z = x - 1/2,
 * w = v/N for a dual vector v = lam*M mod N.  Then lam*q/N - sum(w)/2 = <w,z> mod 1, and
 * 1 - cos(2 pi y) <= 2 pi^2 y^2 for all real y, so
 *     sum_w (1 - cos(2 pi (lam_w q/N - c_w))) <= 2 pi^2 z^T (W^T W) z <= 2 pi^2 mu u/4,
 * with mu >= lambda_max(W^T W).  Exceeding that proves q unreachable, even where every
 * single window is useless.  mu is certified by a Cholesky factorisation of mu I - A. */
static int spec_chol_ok(const double *A, int u, double mu)
{
    double *L = malloc(sizeof(double) * (size_t)u * (size_t)u); int ok = 1;
    if (!L) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    for (int i = 0; i < u && ok; i++) for (int j = 0; j <= i; j++) {
        double t = (i == j ? mu : 0) - A[(size_t)i * u + j];
        for (int k = 0; k < j; k++) t -= L[(size_t)i * u + k] * L[(size_t)j * u + k];
        if (i == j) { if (t <= 0) { ok = 0; break; } L[(size_t)i * u + i] = sqrt(t); }
        else L[(size_t)i * u + j] = t / L[(size_t)j * u + j];
    }
    free(L); return ok;
}
static double spec_mu(const double *A, int u)
{
    double *x = malloc(sizeof(double) * (size_t)u), *y = malloc(sizeof(double) * (size_t)u); double mu = 0;
    if (!x || !y) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    for (int i = 0; i < u; i++) x[i] = 1.0 + 0.01 * i;
    for (int it = 0; it < 400; it++) {
        double n = 0;
        for (int i = 0; i < u; i++) { double a = 0; for (int j = 0; j < u; j++) a += A[(size_t)i * u + j] * x[j]; y[i] = a; n += a * a; }
        n = sqrt(n); if (n <= 0) break; mu = n;
        for (int i = 0; i < u; i++) x[i] = y[i] / n;
    }
    free(x); free(y);
    mu *= 1.000001; while (!spec_chol_ok(A, u, mu)) mu *= 1.0001;
    return mu;
}
static void jacobi_eig(double *A, int n, double *ev, double *V)
{
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) V[(size_t)i * n + j] = (i == j);
    for (int sw = 0; sw < 12; sw++) {
        double off = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) off += A[(size_t)i * n + j] * A[(size_t)i * n + j];
        if (off < 1e-26) break;
        for (int p = 0; p < n; p++) for (int q = p + 1; q < n; q++) {
            const double apq = A[(size_t)p * n + q]; if (fabs(apq) < 1e-30) continue;
            const double th = (A[(size_t)q * n + q] - A[(size_t)p * n + p]) / (2 * apq);
            const double t = (th >= 0 ? 1 : -1) / (fabs(th) + sqrt(th * th + 1)), c = 1 / sqrt(t * t + 1), sn = t * c;
            for (int k = 0; k < n; k++) { const double a = A[(size_t)k * n + p], b = A[(size_t)k * n + q]; A[(size_t)k * n + p] = c * a - sn * b; A[(size_t)k * n + q] = sn * a + c * b; }
            for (int k = 0; k < n; k++) { const double a = A[(size_t)p * n + k], b = A[(size_t)q * n + k]; A[(size_t)p * n + k] = c * a - sn * b; A[(size_t)q * n + k] = sn * a + c * b; }
            for (int k = 0; k < n; k++) { const double a = V[(size_t)k * n + p], b = V[(size_t)k * n + q]; V[(size_t)k * n + p] = c * a - sn * b; V[(size_t)k * n + q] = sn * a + c * b; }
        }
    }
    for (int i = 0; i < n; i++) ev[i] = A[(size_t)i * n + i];
}
/* maximise sum a_w subject to lambda_max(sum a_w w w^T) <= 1, a >= 0, by multiplicative weights on a smoothed lambda_max */
static double spec_design(const double *W, int m, int u, double *a, int iters)
{
    double *A = malloc(sizeof(double) * (size_t)u * (size_t)u), *Aw = malloc(sizeof(double) * (size_t)u * (size_t)u), *V = malloc(sizeof(double) * (size_t)u * (size_t)u), *ev = malloc(sizeof(double) * (size_t)u), *g = malloc(sizeof(double) * (size_t)m), *best = malloc(sizeof(double) * (size_t)m);
    if (!A || !Aw || !V || !ev || !g || !best) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    double bestR = 0;
    for (int i = 0; i < m; i++) { double nr = 0; for (int j = 0; j < u; j++) nr += W[(size_t)i * u + j] * W[(size_t)i * u + j]; a[i] = 1.0 / (nr + 1e-30); }
    for (int it = 0; it < iters; it++) {
        memset(A, 0, sizeof(double) * (size_t)u * (size_t)u);
        for (int k = 0; k < m; k++) { const double *w = W + (size_t)k * u; if (a[k] <= 0) continue; for (int i = 0; i < u; i++) for (int j = 0; j < u; j++) A[(size_t)i * u + j] += a[k] * w[i] * w[j]; }
        memcpy(Aw, A, sizeof(double) * (size_t)u * (size_t)u); jacobi_eig(Aw, u, ev, V);
        double mx = 0, sa = 0; for (int i = 0; i < u; i++) if (ev[i] > mx) mx = ev[i]; for (int k = 0; k < m; k++) sa += a[k];
        if (mx <= 0) break;
        const double R = sa / mx; if (R > bestR) { bestR = R; for (int k = 0; k < m; k++) best[k] = a[k] / mx; }
        /* gradient of the smoothed max eigenvalue at this point */
        const double beta = 8.0 * log((double)u + 1) / mx; double *p = malloc(sizeof(double) * (size_t)u), tot = 0;
        if (!p) exit(2);
        for (int i = 0; i < u; i++) { p[i] = exp(beta * (ev[i] - mx)); tot += p[i]; }
        for (int i = 0; i < u; i++) p[i] /= tot;
        double gbar = 0;
        for (int k = 0; k < m; k++) { const double *w = W + (size_t)k * u; double gg = 0; for (int e = 0; e < u; e++) { if (p[e] < 1e-12) continue; double x = 0; for (int j = 0; j < u; j++) x += w[j] * V[(size_t)j * u + e]; gg += p[e] * x * x; } g[k] = gg; gbar += a[k] * gg; }
        free(p);
        gbar /= sa;
        const double eta = 0.3 / (1 + it / 30.0);
        for (int k = 0; k < m; k++) a[k] *= exp(eta * (1.0 - g[k] / (gbar + 1e-30)) > 2 ? 2 : eta * (1.0 - g[k] / (gbar + 1e-30)));
    }
    memcpy(a, best, sizeof(double) * (size_t)m);
    free(A); free(Aw); free(V); free(ev); free(g); free(best);
    return bestR;
}
/* minimise mean(D) + lambda_max(A - D): the max-cut (hypercube) bound, z^T A z <= (u/4)(mean D + lambda_max(A-D)) for z in {+-1/2}^u */
static void spec_diag(const double *A, int u, double *D, int iters)
{
    double *M = malloc(sizeof(double) * (size_t)u * (size_t)u), *V = malloc(sizeof(double) * (size_t)u * (size_t)u), *ev = malloc(sizeof(double) * (size_t)u), *bestD = calloc((size_t)u, sizeof(double)), *gv = malloc(sizeof(double) * (size_t)u);
    if (!M || !V || !ev || !bestD || !gv) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    memset(D, 0, sizeof(double) * (size_t)u);
    double best = 1e300, scale = 0; for (int i = 0; i < u; i++) scale += A[(size_t)i * u + i]; scale /= u;
    for (int it = 0; it < iters; it++) {
        memcpy(M, A, sizeof(double) * (size_t)u * (size_t)u); double md = 0;
        for (int i = 0; i < u; i++) { M[(size_t)i * u + i] -= D[i]; md += D[i]; } md /= u;
        jacobi_eig(M, u, ev, V);
        double mx = -1e300; for (int i = 0; i < u; i++) if (ev[i] > mx) mx = ev[i];
        const double f = md + mx; if (f < best) { best = f; memcpy(bestD, D, sizeof(double) * (size_t)u); }
        const double beta = 20.0 / (scale + 1e-30); double tot = 0; memset(gv, 0, sizeof(double) * (size_t)u);
        for (int e = 0; e < u; e++) { const double p = exp(beta * (ev[e] - mx)); tot += p; if (p < 1e-9) continue; for (int i = 0; i < u; i++) gv[i] += p * V[(size_t)i * u + e] * V[(size_t)i * u + e]; }
        const double step = 2.0 * scale * (double)u / (1 + it / 20.0);
        for (int i = 0; i < u; i++) D[i] -= step * (1.0 / u - gv[i] / tot);
    }
    memcpy(D, bestD, sizeof(double) * (size_t)u);
    free(M); free(V); free(ev); free(bestD); free(gv);
}
static void dual_spec(duallev *L, const fmpz_mat_t A, const fmpz_t N, const fmpz_t R, int d, int dbg)
{
    const int u = L->u, n = u + 1;
    const int Mcap = getenv("TP_SPEC_M") ? atoi(getenv("TP_SPEC_M")) : 1500, S = 4096;
    const int ntrip = getenv("TP_SPEC_TRIP") ? atoi(getenv("TP_SPEC_TRIP")) : 200000;
    const double Nd = fmpz_get_d(N);
    fmpz_t t, lam, ss; fmpz_init(t); fmpz_init(lam); fmpz_init(ss);
    double *Vb = malloc(sizeof(double) * (size_t)n * (size_t)u);
    const size_t maxp = (size_t)n + 2 * (size_t)n * (size_t)n + (size_t)ntrip;
    double *ord = malloc(sizeof(double) * 2 * maxp); int (*cf)[6] = malloc(sizeof(int[6]) * maxp);
    if (!Vb || !ord || !cf) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    for (int i = 0; i < n; i++) for (int j = 0; j < u; j++) { fmpz_divexact(t, fmpz_mat_entry(A, i, j + 1), N); Vb[(size_t)i * u + j] = fmpz_get_d(t) / Nd; }
    uint64_t rng = 0xA0761D6478BD642FULL ^ ((uint64_t)d << 40) ^ (uint64_t)u; for (int i = 0; i < 8; i++) dl_rng(&rng);
    size_t np = 0;
    #define ADDC(a,sa,b,sb,c,sc) do { cf[np][0]=a; cf[np][1]=sa; cf[np][2]=b; cf[np][3]=sb; cf[np][4]=c; cf[np][5]=sc; \
        double nr = 0; for (int j = 0; j < u; j++) { double e = (sa)*Vb[(size_t)(a)*u+j]; if (sb) e += (sb)*Vb[(size_t)(b)*u+j]; if (sc) e += (sc)*Vb[(size_t)(c)*u+j]; nr += e*e; } \
        ord[2*np] = nr; ord[2*np+1] = (double)np; np++; } while (0)
    for (int a = 0; a < n; a++) ADDC(a, 1, 0, 0, 0, 0);
    for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { ADDC(a, 1, b, 1, 0, 0); ADDC(a, 1, b, -1, 0, 0); }
    for (int k = 0; k < ntrip; k++) { int a = (int)(dl_rng(&rng) % (uint64_t)n), b = (int)(dl_rng(&rng) % (uint64_t)n), c = (int)(dl_rng(&rng) % (uint64_t)n);
        if (a == b || b == c || a == c) continue; ADDC(a, 1, b, (dl_rng(&rng) & 1) ? 1 : -1, c, (dl_rng(&rng) & 1) ? 1 : -1); }
    qsort(ord, np, 2 * sizeof(double), cmp_dpi);
    /* exact lambda and centre for the shortest distinct vectors */
    u256 *lm = malloc(sizeof(u256) * (size_t)Mcap); double *cen = malloc(sizeof(double) * (size_t)Mcap), *W = malloc(sizeof(double) * (size_t)Mcap * (size_t)u);
    uint64_t *hs = malloc(sizeof(uint64_t) * (size_t)Mcap);
    if (!lm || !cen || !W || !hs) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    int m = 0;
    for (size_t oi = 0; oi < np && m < Mcap; oi++) {
        const int *c = cf[(size_t)ord[2 * oi + 1]];
        if (ord[2 * oi] <= 1e-18) continue;
        fmpz_zero(lam); fmpz_zero(ss); uint64_t h = 0x9E3779B97F4A7C15ULL; int sgn = 0; double *w = W + (size_t)m * u;
        for (int j = 0; j <= u; j++) {
            fmpz_zero(t);
            for (int q = 0; q < 3; q++) { if (!c[2 * q + 1]) continue; if (c[2 * q + 1] > 0) fmpz_add(t, t, fmpz_mat_entry(A, c[2 * q], j)); else fmpz_sub(t, t, fmpz_mat_entry(A, c[2 * q], j)); }
            if (j == 0) { fmpz_set(lam, t); continue; }
            fmpz_divexact(t, t, N); fmpz_add(ss, ss, t); w[j - 1] = fmpz_get_d(t) / Nd;
            if (!sgn) sgn = fmpz_sgn(t);
            { const uint64_t mp = 2305843009213693951ULL; uint64_t r = fmpz_fdiv_ui(t, mp); if (sgn < 0 && r) r = mp - r; h = (h ^ r) * 0x100000001B3ULL; }
        }
        if (!sgn) continue;
        int dup = 0; for (int q = 0; q < m; q++) if (hs[q] == h) { dup = 1; break; }
        if (dup) continue;
        hs[m] = h;
        fmpz_mod(t, lam, N); fmpz_mul(t, t, R); fmpz_mod(t, t, N); f2u(&lm[m], t);
        cen[m] = fmpz_get_d(ss) / Nd / 2.0;
        m++;
    }
    /* score m-prefixes on random residues; keep the prefix with the lowest pass rate */
    double *Am = calloc((size_t)u * (size_t)u, sizeof(double)), *dfc = malloc(sizeof(double) * (size_t)TP_SPEC_MAX * S);
    u256 *q = malloc(sizeof(u256) * S);
    if (!Am || !dfc || !q) { fprintf(stderr, "spec: out of memory\n"); exit(2); }
    for (int a = 0; a < S; a++) { uint64_t x[4]; for (int b = 0; b < 4; b++) x[b] = dl_rng(&rng); dual_reduce(&q[a], x, 4, L); }
    double *a_w = malloc(sizeof(double) * (size_t)(m ? m : 1));
    if (!a_w) exit(2);
    const int iters = getenv("TP_SPEC_ITERS") ? atoi(getenv("TP_SPEC_ITERS")) : 300;
    double minn = 1e300; for (int k = 0; k < m; k++) { double nr = 0; for (int j = 0; j < u; j++) nr += W[(size_t)k * u + j] * W[(size_t)k * u + j]; if (nr < minn) minn = nr; }
    const double Rdes = m ? spec_design(W, m, u, a_w, iters) : 0;
    /* keep the heaviest vectors only; everything is re-certified on what is kept */
    int m2 = 0; { double amx = 0; for (int k = 0; k < m; k++) if (a_w[k] > amx) amx = a_w[k];
      int *oi = malloc(sizeof(int) * (size_t)m); if (!oi) exit(2);
      for (int k = 0; k < m; k++) oi[k] = k;
      for (int i = 0; i < m && i < TP_SPEC_MAX; i++) { int bi = i; for (int j = i + 1; j < m; j++) if (a_w[oi[j]] > a_w[oi[bi]]) bi = j; const int tt = oi[i]; oi[i] = oi[bi]; oi[bi] = tt; }
      for (int i = 0; i < m && i < TP_SPEC_MAX; i++) if (a_w[oi[i]] > 1e-4 * amx) m2++; else break;
      double *W2 = malloc(sizeof(double) * (size_t)m2 * (size_t)u), *a2 = malloc(sizeof(double) * (size_t)m2), *c2 = malloc(sizeof(double) * (size_t)m2); u256 *l2 = malloc(sizeof(u256) * (size_t)m2);
      if (!W2 || !a2 || !c2 || !l2) exit(2);
      for (int i = 0; i < m2; i++) { memcpy(W2 + (size_t)i * u, W + (size_t)oi[i] * u, sizeof(double) * (size_t)u); a2[i] = a_w[oi[i]]; c2[i] = cen[oi[i]]; l2[i] = lm[oi[i]]; }
      memcpy(W, W2, sizeof(double) * (size_t)m2 * (size_t)u); memcpy(a_w, a2, sizeof(double) * (size_t)m2); memcpy(cen, c2, sizeof(double) * (size_t)m2); memcpy(lm, l2, sizeof(u256) * (size_t)m2);
      free(oi); free(W2); free(a2); free(c2); free(l2); }
    /* residues and phases only for the kept vectors */
    for (int k = 0; k < m2; k++) for (int a = 0; a < S; a++) {
        u256 z; u256_mont(&z, &lm[k], &q[a], &L->N, L->n0, L->nl);
        const double v = (ldexp((double)z.w[3], 192) + ldexp((double)z.w[2], 128) + ldexp((double)z.w[1], 64) + (double)z.w[0]) / Nd;
        { double y = v - cen[k]; y -= rint(y); const double sn = sin(3.141592653589793 * y); dfc[(size_t)k * S + a] = 2.0 * sn * sn; }
    }
    memset(Am, 0, sizeof(double) * (size_t)u * (size_t)u);
    for (int k = 0; k < m2; k++) { const double *w = W + (size_t)k * u; for (int i = 0; i < u; i++) for (int j = 0; j < u; j++) Am[(size_t)i * u + j] += a_w[k] * w[i] * w[j]; }
    double mu0 = m2 ? spec_mu(Am, u) : 1, B0 = 19.739208802178716 * mu0 * (double)u / 4.0;
    double *Dg = malloc(sizeof(double) * (size_t)u), *Ad = malloc(sizeof(double) * (size_t)u * (size_t)u);
    if (!Dg || !Ad) exit(2);
    double B = B0;
    if (m2 && !getenv("TP_SPEC_NODIAG")) {
        spec_diag(Am, u, Dg, 200);
        memcpy(Ad, Am, sizeof(double) * (size_t)u * (size_t)u); double md = 0;
        for (int i = 0; i < u; i++) { Ad[(size_t)i * u + i] -= Dg[i]; md += Dg[i]; } md /= u;
        /* shift so power iteration finds the top of the spectrum */
        double shift = 0; for (int i = 0; i < u; i++) shift = fmax(shift, fabs(Dg[i]));
        for (int i = 0; i < u; i++) Ad[(size_t)i * u + i] += shift;
        const double mud = spec_mu(Ad, u) - shift;
        const double Bd = 19.739208802178716 * (md + mud) * (double)u / 4.0;
        if (Bd < B) B = Bd;
    }
    B = B * (1 + 1e-7) + 1e-7;
    {   /* rounding: |2 sin^2(pi(y+e)) - 2 sin^2(pi y)| <= 2 pi^2 (2|y||e| + e^2), with |y| <= |<w,z>| <= |w| sqrt(u)/2 for reachable q */
        double cmx = 0, Bnum = 0; for (int k = 0; k < m2; k++) if (fabs(cen[k]) > cmx) cmx = fabs(cen[k]);
        const double e = 8e-16 * (1.0 + cmx);
        for (int k = 0; k < m2; k++) { double nr = 0; const double *w = W + (size_t)k * u; for (int j = 0; j < u; j++) nr += w[j] * w[j];
            const double ym = sqrt(nr * (double)u) / 2.0 + e;
            Bnum += a_w[k] * (19.739208802178716 * (2.0 * ym * e + e * e) + 1e-15); }
        B += Bnum;
    }
    free(Dg); free(Ad);
    int ok = 0; double sa = 0; for (int k = 0; k < m2; k++) sa += a_w[k];
    for (int qi = 0; qi < S; qi++) { double df = 0; for (int k = 0; k < m2; k++) { df += a_w[k] * dfc[(size_t)k * S + qi]; if (df > B) break; } if (df <= B) ok++; }
    const double pr = (double)ok / S;
    if (dbg) fprintf(stderr, "  [spec] d=%d u=%d K=%d shortest|w|=%.4f pool=%d kept=%d margin plain %.3f hypercube %.3f pass %.4f (windows nt=%d)\n", d, u, L->K, sqrt(minn), m, m2, sa / B0, sa / B, pr, L->nt);
    if (pr < 0.98) {
        int nk = m2 < TP_SPEC_MAX ? m2 : TP_SPEC_MAX;
        for (int k = 0; k < nk; k++) { L->slam[k] = lm[k]; L->scen[k] = cen[k]; L->swt[k] = a_w[k]; }
        L->ns = nk; L->sNd = Nd; L->sB = B * (1 + 1e-7);
    }
    (void)Rdes; free(a_w);
    free(Am); free(dfc); free(q); free(lm); free(cen); free(W); free(hs); free(Vb); free(ord); free(cf);
    fmpz_clear(t); fmpz_clear(lam); fmpz_clear(ss);
}
#endif
static void dual_build(hexplan *p, const tp_orbits *o, const int *lvl)
{
    memset(p->dual, 0, sizeof p->dual);
    if (!getenv("TP_DUAL")) return;
#ifndef TP_FLINT
    fprintf(stderr, "--dual needs a FLINT build\n"); exit(2);
#else
    if (!(g_fastk && g_carry9 && !g_xm && !g_decchk && g_nohi && !g_noprune && !g_badprune && !g_r5chk)) return;
    if ((o->nsites % 4 == 1 || o->nsites % 4 == 2) && !getenv("TP_DUAL_ALL")) return;
    const int D = p->depth_max, nw = o->nw64, dbg = getenv("TP_DUALDBG") != NULL;
    const int Kov = getenv("TP_DUAL_K") ? atoi(getenv("TP_DUAL_K")) : 0;
    const int maxl = getenv("TP_DUAL_LEVELS") ? atoi(getenv("TP_DUAL_LEVELS")) : 8;
    pthread_mutex_lock(&g_du_mu);
    const uint64_t oh = r5_orbit_hash(o);
    if (g_duc.key_n != o->nsites || g_duc.key_h != oh || g_duc.K != Kov) {
        for (int d = 0; d < 64; d++) { free(g_duc.lev[d]); g_duc.lev[d] = NULL; }
        fmpz_mat_t Asv[64]; int hasA[64]; memset(hasA, 0, sizeof hasA);
        const int dmin = (int)floor((double)(o->nsites - 1) * 0.30102999566) + 1;
        int seen = 0, extra = ((getenv("TP_DUAL_POOL") || getenv("TP_SPEC")) && getenv("TP_DUAL_EXTRA")) ? atoi(getenv("TP_DUAL_EXTRA")) : 0;
        for (int d = (D < 64 ? D : 64) - 2; d >= 1; d--) {
            int u = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0 && lvl[c] > d) u++;
            if (u < 8) continue;
            int K = dmin - p->ddig[d - 1] - 2;
            if (Kov) K = Kov;
            if (K > TP_DUAL_KMAX) K = TP_DUAL_KMAX;
            if (K < 2) { if (seen) break; continue; }
            fmpz_t N, t, neg, pos, lam, tmp, R; fmpz_init(N); fmpz_init(t); fmpz_init(neg); fmpz_init(pos); fmpz_init(lam); fmpz_init(tmp); fmpz_init(R);
            fmpz_ui_pow_ui(N, 5, (ulong)K);
            fmpz_mat_t A; fmpz_mat_init(A, u + 1, u + 1);
            fmpz_one(fmpz_mat_entry(A, 0, 0));
            { int j = 0;
              for (int c = 0; c < o->ncomp; c++) {
                  if (o->cval[c] >= 0 || lvl[c] <= d) continue;
                  limbs2f(t, o->mask + (size_t)c * (size_t)o->nw64, o->nw64);
                  fmpz_mod(t, t, N);
                  fmpz_mul(fmpz_mat_entry(A, 0, j + 1), t, N);
                  fmpz_mul(fmpz_mat_entry(A, j + 1, j + 1), N, N);
                  j++;
              } }
            fmpz_lll_t fl; fmpz_lll_context_init_default(fl);
            fmpz_lll(A, NULL, fl);
            duallev *L = calloc(1, sizeof(duallev));
            if (!L) { fprintf(stderr, "dual: out of memory\n"); exit(2); }
            L->K = K; L->u = u; f2u(&L->N, N);
            { int nl = (int)(fmpz_bits(N) + 64) / 64; if (nl > 4) nl = 4; L->nl = nl; }
            { uint64_t inv = 1; for (int i = 0; i < 6; i++) inv *= 2 - L->N.w[0] * inv; L->n0 = (uint64_t)0 - inv; }
            fmpz_one(R); fmpz_mul_2exp(R, R, (ulong)(64 * L->nl)); fmpz_mod(R, R, N);
            fmpz_mul_2exp(t, R, 64); fmpz_mod(t, t, N); f2u(&L->c64, t);
            fmpz_mul_ui(t, R, 10); fmpz_mod(t, t, N); f2u(&L->c10, t);
            double best = 1e9;
            for (int i = 0; i <= u; i++) {
                fmpz_zero(neg); fmpz_zero(pos);
                int any = 0;
                for (int j = 1; j <= u; j++) {
                    fmpz_divexact(t, fmpz_mat_entry(A, i, j), N);
                    if (fmpz_sgn(t) < 0) { fmpz_add(neg, neg, t); any = 1; }
                    else if (fmpz_sgn(t) > 0) { fmpz_add(pos, pos, t); any = 1; }
                }
                if (!any) continue;
                fmpz_sub(tmp, pos, neg);
                if (fmpz_cmp(tmp, N) >= 0) continue;
                const double f = exp(fmpz_dlog(tmp) - fmpz_dlog(N));
                if (f < best) best = f;
                if (L->nt >= TP_DUAL_TMAX) continue;
                fmpz_mod(lam, fmpz_mat_entry(A, i, 0), N);
                fmpz_mul(t, lam, R); fmpz_mod(t, t, N);
                u256 lm, na, sp; f2u(&lm, t); fmpz_neg(t, neg); fmpz_mod(t, t, N); f2u(&na, t); f2u(&sp, tmp);
                int at = L->nt++;
                while (at > 0 && !u256_ge(&sp, &L->span[at - 1])) {    /* tightest window first */
                    L->lamM[at] = L->lamM[at - 1]; L->negabs[at] = L->negabs[at - 1]; L->span[at] = L->span[at - 1]; at--;
                }
                L->lamM[at] = lm; L->negabs[at] = na; L->span[at] = sp;
            }
            const int ext = (L->nt == 0 && seen && extra > 0); if (ext) extra--;
            if (L->nt == 0 && !ext) { fmpz_mat_clear(A); free(L); fmpz_clear(N); fmpz_clear(t); fmpz_clear(neg); fmpz_clear(pos); fmpz_clear(lam); fmpz_clear(tmp); fmpz_clear(R); if (seen) break; continue; }
            seen = 1;
            fmpz_mat_init_set(Asv[d], A); hasA[d] = 1; fmpz_mat_clear(A);
            for (int k = 0; k < p->nviab[d]; k++) {
                const int v = p->vtab[d][k];
                bn_t set; bn_zero(&set, nw);
                for (int i = 0; i < p->nn[d]; i++) {
                    const int site = p->nb[d][i];
                    const int c = o->comp[site];
                    const int bit = (v >> (site & 3)) & 1;
                    if (o->cval[c] >= 0 || lvl[c] != d) continue;
                    if (bit) { const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64; for (int w = 0; w < nw; w++) set.w[w] |= mk[w]; }
                }
                limbs2f(t, set.w, nw); fmpz_mod(t, t, N); f2u(&L->add[k], t);
            }
            if (g_off) { limbs2f(t, g_offc.w, TP_MAX_LIMBS); fmpz_mod(t, t, N); f2u(&L->coff, t); }
            if (dbg) fprintf(stderr, "  [dual] d=%d undecided=%d K=%d tests=%d best span/N=%.4f\n", d, u, K, L->nt, best);
            g_duc.lev[d] = L;
            fmpz_clear(N); fmpz_clear(t); fmpz_clear(neg); fmpz_clear(pos); fmpz_clear(lam); fmpz_clear(tmp); fmpz_clear(R);
        }
        { int kept = 0;                  /* keep the shallowest levels: the earliest kills */
          for (int d = 1; d < 64; d++) if (g_duc.lev[d]) { if (kept >= maxl) { free(g_duc.lev[d]); g_duc.lev[d] = NULL; } else kept++; } }
        for (int d = 1; d < 64; d++) {
            if (!hasA[d]) continue;
            duallev *L = g_duc.lev[d];
            if (L) {
                const int u = L->u, K = L->K; (void)K;
                fmpz_t N, t, R; fmpz_init(N); fmpz_init(t); fmpz_init(R);
                fmpz_ui_pow_ui(N, 5, (ulong)L->K);
                fmpz_one(R); fmpz_mul_2exp(R, R, (ulong)(64 * L->nl)); fmpz_mod(R, R, N);
                if (getenv("TP_DUAL_POOL")) dual_learn(L, Asv[d], N, R, d, dbg);
                if (getenv("TP_SPEC")) dual_spec(L, Asv[d], N, R, d, dbg);
                if (g_duchk) {       /* genuine subset sums must always pass */
                    u256 Mv[256]; int j = 0;
                    for (int c = 0; c < o->ncomp; c++) {
                        if (o->cval[c] >= 0 || lvl[c] <= d) continue;
                        limbs2f(t, o->mask + (size_t)c * (size_t)o->nw64, o->nw64); fmpz_mod(t, t, N); f2u(&Mv[j++], t);
                    }
                    uint64_t rng = 0x9E3779B97F4A7C15ULL ^ (uint64_t)d, bad = 0, rp = 0;
                    const uint64_t tot = 200000;
                    for (uint64_t a = 0; a < tot; a++) {
                        u256 q = {{0,0,0,0}};
                        for (int b = 0; b < j; b++) {
                            rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                            if (rng >> 63) u256_addmod(&q, &q, &Mv[b], &L->N);
                        }
                        if (!dual_pass(L, &q)) bad++;
                        uint64_t x[4]; for (int b = 0; b < 4; b++) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; x[b] = rng; }
                        u256 r2; dual_reduce(&r2, x, 4, L); if (dual_pass(L, &r2)) rp++;
                    }
                    fprintf(stderr, "  [dualchk] d=%d: %llu of %llu real sums rejected; random-residue pass rate %.4f\n", d, (unsigned long long)bad, (unsigned long long)tot, (double)rp / (double)tot);
                    if (bad) exit(9);
                }
                (void)u;
                fmpz_clear(N); fmpz_clear(t); fmpz_clear(R);
                if (L->nt == 0 && L->ns == 0) { free(L); g_duc.lev[d] = NULL; }
                else { L->nadd = p->nviab[d];
                       for (int i = 0; i < L->nt; i++) for (int k = 0; k < L->nadd; k++) u256_mont(&L->ladd[i][k], &L->lamM[i], &L->add[k], &L->N, L->n0, L->nl); }
            }
            fmpz_mat_clear(Asv[d]);
        }
        g_duc.key_n = o->nsites; g_duc.key_h = oh; g_duc.K = Kov;
    }
    for (int d = 0; d < 64; d++) p->dual[d] = g_duc.lev[d];
    pthread_mutex_unlock(&g_du_mu);
#endif
}

/* ---- exact reachability by lattice enumeration ("--cvp") ----------------------
 * Lambda = { y in Z^u : sum y_c M_c = 0 mod N } has an LLL basis b_j with tiny entries.
 * The residue q is reachable iff the coset contains a point of {0,1}^u; every such
 * point lies on the sphere |x - 1/2|^2 = u/4, so Fincke-Pohst enumeration with that
 * radius sees all of them.  x = sum (f_j + z_j) b_j with f_j = frac(a_j q / N), so a
 * query is u modular multiplications and a short tree search.  If the search
 * completes empty, the target is provably unreachable; if it hits the node limit
 * the child is simply kept. */
#define TP_CVP_MX 128
typedef struct cvplev_s {
    int K, u, nl; u256 N; uint64_t n0; u256 c64, c10; double Nd;
    long limit;
    u256 add[16];
    u256 aM[TP_CVP_MX];
    u256 cadd[TP_CVP_MX][16];
    u256 coff;                                         /* TP_OFFSET: c mod N */                          /* a_j * add[k] mod N: children share the parent's products */
    double h[TP_CVP_MX], n2[TP_CVP_MX], wbox[TP_CVP_MX], R2;
    int nal[TP_CVP_MX], al[TP_CVP_MX][4];   /* alphabet of each coordinate */
    double bs[TP_CVP_MX][TP_CVP_MX], cen[TP_CVP_MX];   /* Gram-Schmidt vectors, box half-widths */
    double rmax[TP_CVP_MX];                            /* radius bound per level; learned when TP_CVP_PRUNE is set */
    double mu[TP_CVP_MX][TP_CVP_MX];
    double B[TP_CVP_MX][TP_CVP_MX];
} cvplev;
static struct { int key_n; uint64_t key_h; long bud; int K; cvplev *lev[64]; } g_cvc = { -1, 0, 0, 0, {0} };
static pthread_mutex_t g_cv_mu = PTHREAD_MUTEX_INITIALIZER;
static int g_cvchk = 0;
static int cmp_double(const void *a, const void *b) { const double x = *(const double *)a, y = *(const double *)b; return (x > y) - (x < y); }
static int g_cvstat = 0; static uint64_t g_cv_q = 0, g_cv_n = 0, g_cv_r = 0;
static void cv_report(void) { if (g_cvstat) fprintf(stderr, "  [cvpstat] queries %llu, enum nodes %llu (avg %.1f), rejected %llu\n", (unsigned long long)g_cv_q, (unsigned long long)g_cv_n, g_cv_q ? (double)g_cv_n/(double)g_cv_q : 0.0, (unsigned long long)g_cv_r); }
typedef struct {
    const cvplev *L; const u256 *q; double f[TP_CVP_MX]; int fdone; double R2; double y[TP_CVP_MX]; long nodes; int found, over;
    u256 *zb; uint8_t *zok; int kidx;                  /* sibling sharing: zb[j] = a_j * base, child k subtracts cadd[j][k] */
    double sum[TP_CVP_MX][TP_CVP_MX + 1];   /* sum[i][j] = sum_{k>=j} y_k mu[k][i], valid above rc[i] */
    int rc[TP_CVP_MX];
    double p[TP_CVP_MX + 1][TP_CVP_MX];     /* p[i] = projection of (x - centre) onto b*_i .. b*_{n-1} */
} cvp_st;
static inline double u256_dbl(const u256 *a)
{ return ldexp((double)a->w[3], 192) + ldexp((double)a->w[2], 128) + ldexp((double)a->w[1], 64) + (double)a->w[0]; }
static void cvp_rec(cvp_st *S, int i, double rho)
{
    const cvplev *L = S->L;
    const int top = S->rc[i];
    for (int j = top; j > i; j--) S->sum[i][j] = S->sum[i][j + 1] + S->y[j] * L->mu[j][i];
    S->rc[i] = i;
    if (i > 0 && S->rc[i - 1] < top) S->rc[i - 1] = top;   /* the level below must see these changes even if we never descend */
    const double c = L->h[i] - S->sum[i][i + 1];
    const double room = S->R2 - rho;
    double w = (room > 0) ? sqrt(room / L->n2[i]) : 0.0;
    if (w > L->wbox[i]) w = L->wbox[i];     /* x lies in the unit cube: |<x-1/2,b*_i>| <= |b*_i|_1 / 2 */
    while (S->fdone > i) {                   /* coordinates are needed top-down: most queries die in a few levels */
        S->fdone--;
        u256 z;
        if (S->zb) {                         /* a_j (base - add_k) = a_j base - a_j add_k (mod N): same residue, one multiply per parent */
            const int f = S->fdone;
            if (!S->zok[f]) { u256_mont(&S->zb[f], &L->aM[f], S->q, &L->N, L->n0, L->nl); S->zok[f] = 1; }
            u256_submod(&z, &S->zb[f], &L->cadd[f][S->kidx], &L->N);
        } else u256_mont(&z, &L->aM[S->fdone], S->q, &L->N, L->n0, L->nl);
        const double v = u256_dbl(&z) / L->Nd;
        S->f[S->fdone] = v - floor(v);
    }
    const double fi = S->f[i];
    const long lo = (long)ceil(c - w - fi), hi = (long)floor(c + w - fi);
    for (long z = lo; z <= hi; z++) {
        if (++S->nodes > L->limit) { S->over = 1; return; }
        const double yi = fi + (double)z, d = yi - c, r2 = rho + d * d * L->n2[i];
        if (r2 > L->rmax[i]) continue;
        {   /* the projection p of a point of the box satisfies |p|^2 = <x-c,p> <= sum c_k |p_k| */
            double s1 = 0, s2 = 0;
            const double *bv = L->bs[i], *pu = S->p[i + 1], *cn = L->cen;
            double *pv = S->p[i];
            for (int k = 0; k < L->u; k++) {
                const double t = pu[k] + d * bv[k];
                pv[k] = t; s2 += t * t; s1 += cn[k] * fabs(t);
            }
            if (s2 > s1 * (1.0 + 1e-9) + 1e-6) continue;
        }
        S->y[i] = yi;
        if (i == 0) {
            int ok = 1;
            for (int k = 0; k < L->u && ok; k++) {
                double v = 0; for (int j = 0; j < L->u; j++) v += S->y[j] * L->B[j][k];
                const double r = rint(v);
                if (fabs(v - r) > 1e-4) { ok = 0; break; }
                int hit = 0; for (int a = 0; a < L->nal[k]; a++) if (r == (double)L->al[k][a]) hit = 1;
                if (!hit) ok = 0;
            }
            if (ok) { S->found = 1; return; }
        } else {
            if (S->rc[i - 1] < i) S->rc[i - 1] = i;
            cvp_rec(S, i - 1, r2);
            if (S->found || S->over) return;
        }
    }
}
/* the same test for child k of a parent whose residue is base: q = base - add[k] */
static inline int cvp_keep_child(const cvplev *L, const u256 *base, int k, u256 *zb, uint8_t *zok, long *nodes)
{
    static __thread cvp_st *Sq = NULL;
    if (!Sq) { Sq = malloc(sizeof(cvp_st)); if (!Sq) { fprintf(stderr, "cvp: out of memory\n"); exit(2); } }
    cvp_st *S = Sq;
    S->L = L; S->q = base; S->fdone = L->u; S->R2 = L->R2; S->nodes = 0; S->found = 0; S->over = 0;
    S->zb = zb; S->zok = zok; S->kidx = k;
    for (int i = 0; i < L->u; i++) { S->rc[i] = L->u - 1; S->sum[i][L->u] = 0.0; S->p[L->u][i] = 0.0; }
    cvp_rec(S, L->u - 1, 0.0);
    if (nodes) *nodes = S->nodes;
    return S->found || S->over;
}
/* 1 = keep (reachable or undecided), 0 = provably unreachable */
static inline int cvp_keep(const cvplev *L, const u256 *q, long *nodes)
{
    /* the working state is ~270 KB, far too big for a worker's stack: one heap copy per thread */
    static __thread cvp_st *Sp = NULL;
    if (!Sp) { Sp = malloc(sizeof(cvp_st)); if (!Sp) { fprintf(stderr, "cvp: out of memory\n"); exit(2); } }
    cvp_st *S = Sp;
    S->L = L; S->q = q; S->fdone = L->u; S->R2 = L->R2; S->nodes = 0; S->found = 0; S->over = 0; S->zb = NULL;
    for (int i = 0; i < L->u; i++) { S->rc[i] = L->u - 1; S->sum[i][L->u] = 0.0; S->p[L->u][i] = 0.0; }
    cvp_rec(S, L->u - 1, 0.0);
    if (nodes) *nodes = S->nodes;
    return S->found || S->over;
}
#ifdef TP_FLINT
/* BKZ: strengthen the LLL basis once per level, offline.  Entries are small, so exact int64 rows are fine. */
typedef struct { int n; const double *n2; const double (*mu)[TP_CVP_MX]; double best; long c[TP_CVP_MX], bc[TP_CVP_MX]; int k, found; long nodes; } svp_t;
static void svp_rec(svp_t *T, int i, double rho, int nz)
{
    if (T->nodes > 20000000) return;
    double c = 0; for (int j = i + 1; j < T->k + T->n; j++) c -= (double)T->c[j] * T->mu[j][i];
    long v0 = lround(c), step = 0;
    for (;;) {
        const long off = (step == 0) ? 0 : ((step & 1) ? (step + 1) / 2 : -step / 2);
        const long v = v0 + off;
        const double d = (double)v - c, r2 = rho + d * d * T->n2[i];
        T->nodes++;
        step++;
        if (r2 >= T->best) { if (step > 2 && fabs((double)off) * sqrt(T->n2[i]) > sqrt(T->best) + fabs(d) + sqrt(T->n2[i])) break; if (step > 64) break; if (r2 >= T->best && step >= 2 && fabs(d) > sqrt(T->best / T->n2[i]) + 2) break; continue; }
        T->c[i] = v;
        if (i == T->k) {
            if (nz || v != 0) { T->best = r2 * 0.9999; memcpy(T->bc, T->c, sizeof T->c); T->found = 1; }
        } else svp_rec(T, i - 1, r2, nz || v != 0);
        if (step > 400) break;
    }
    T->c[i] = 0;
}
static void gso_i64(long B[][TP_CVP_MX], int u, double bs[][TP_CVP_MX], double mu[][TP_CVP_MX], double *n2)
{
    for (int i = 0; i < u; i++) {
        for (int k = 0; k < u; k++) bs[i][k] = (double)B[i][k];
        for (int j = 0; j < i; j++) {
            double a = 0; for (int k = 0; k < u; k++) a += (double)B[i][k] * bs[j][k];
            mu[i][j] = a / n2[j];
            for (int k = 0; k < u; k++) bs[i][k] -= mu[i][j] * bs[j][k];
        }
        double a = 0; for (int k = 0; k < u; k++) a += bs[i][k] * bs[i][k];
        n2[i] = a;
    }
}
static void lll_i64(long B[][TP_CVP_MX], int u)
{
    fmpz_mat_t A; fmpz_mat_init(A, u, u);
    for (int i = 0; i < u; i++) for (int k = 0; k < u; k++) fmpz_set_si(fmpz_mat_entry(A, i, k), B[i][k]);
    fmpz_lll_t fl; fmpz_lll_context_init_default(fl); fmpz_lll(A, NULL, fl);
    for (int i = 0; i < u; i++) for (int k = 0; k < u; k++) B[i][k] = fmpz_get_si(fmpz_mat_entry(A, i, k));
    fmpz_mat_clear(A);
}
static void bkz_i64(long B[][TP_CVP_MX], int u, int beta, int tours)
{
    static __thread double bs[TP_CVP_MX][TP_CVP_MX], mu[TP_CVP_MX][TP_CVP_MX]; double n2[TP_CVP_MX];
    for (int t = 0; t < tours; t++) {
        int changed = 0;
        for (int k = 0; k + 1 < u; k++) {
            int n = (k + beta <= u) ? beta : u - k;
            if (n < 2) break;
            gso_i64(B, u, bs, mu, n2);
            svp_t T; memset(&T, 0, sizeof T);
            T.n = n; T.n2 = n2; T.mu = (const double (*)[TP_CVP_MX])mu; T.best = n2[k] * 0.98; T.k = k;
            svp_rec(&T, k + n - 1, 0.0, 0);
            if (!T.found) continue;
            int jp = -1; for (int j = k + n - 1; j >= k; j--) if (T.bc[j] == 1 || T.bc[j] == -1) { jp = j; break; }
            if (jp < 0) continue;
            long v[TP_CVP_MX]; memset(v, 0, sizeof v);
            for (int j = k; j < k + n; j++) if (T.bc[j]) for (int m = 0; m < u; m++) v[m] += T.bc[j] * B[j][m];
            for (int j = jp; j > k; j--) memcpy(B[j], B[j - 1], sizeof(long) * (size_t)u);
            memcpy(B[k], v, sizeof(long) * (size_t)u);
            lll_i64(B, u);
            changed++;
        }
        if (!changed) break;
    }
}
#endif
static void cvp_build(hexplan *p, const tp_orbits *o, const int *lvl)
{
    memset(p->cvp, 0, sizeof p->cvp);
    if (!getenv("TP_CVP")) return;
    if (getenv("TP_CVPSTAT") && !g_cvstat) { g_cvstat = 1; atexit(cv_report); }
#ifndef TP_FLINT
    fprintf(stderr, "--cvp needs a FLINT build\n"); exit(2);
#else
    if (!(g_fastk && g_carry9 && !g_xm && !g_decchk && g_nohi && !g_noprune && !g_badprune && !g_r5chk)) return;
    if ((o->nsites % 4 == 1 || o->nsites % 4 == 2) && !getenv("TP_CVP_ALL")) return;
    const int D = p->depth_max, nw = o->nw64, dbg = getenv("TP_CVPDBG") != NULL;
    const long bud = getenv("TP_CVP_BUDGET") ? atol(getenv("TP_CVP_BUDGET")) : 30;
    const int Kov = getenv("TP_CVP_K") ? atoi(getenv("TP_CVP_K")) : 0;
    const int bkzb = getenv("TP_CVP_BKZ") ? atoi(getenv("TP_CVP_BKZ")) : 0;
    pthread_mutex_lock(&g_cv_mu);
    const uint64_t oh = r5_orbit_hash(o);
    if (g_cvc.key_n != o->nsites || g_cvc.key_h != oh || g_cvc.bud != bud || g_cvc.K != Kov) {
        for (int d = 0; d < 64; d++) { free(g_cvc.lev[d]); g_cvc.lev[d] = NULL; }
        const int dmin = (int)floor((double)(o->nsites - 1) * 0.30102999566) + 1;
        cvplev *best = NULL; int bestd = -1;
        for (int d = (D < 64 ? D : 64) - 2; d >= 1; d--) {
            int u = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0 && lvl[c] > d) u++;
            if (u < 12) continue;
            if (u >= TP_CVP_MX) break;
            int K = dmin - p->ddig[d - 1] - 2;
            if (Kov) K = Kov;
            else if (getenv("TP_CVP_KX")) K += atoi(getenv("TP_CVP_KX"));
            if (K > TP_DUAL_KMAX) K = TP_DUAL_KMAX;
            {   /* a modulus much beyond 2^(3.5u) adds nothing and makes the basis entries huge */
                const int cap = (int)ceil(3.5 * u / 2.3219) + 2;
                if (!Kov && K > cap) K = cap;
            }
            if (K < 2) break;
            fmpz_t N, t, R, cinv, den; fmpz_init(N); fmpz_init(t); fmpz_init(R); fmpz_init(cinv); fmpz_init(den);
            fmpz_ui_pow_ui(N, 5, (ulong)K);
            /* variables: a level whose two orbits have proportional masks becomes one digit
             * variable over {0,a,b,a+b}; anything else stays a 0/1 variable.  This halves the
             * dimension.  Coordinates are weighted so every alphabet has the same half-width. */
            fmpz *M = _fmpz_vec_init(TP_CVP_MX);
            int nvars = 0, nalv[TP_CVP_MX], alv[TP_CVP_MX][4]; long Wv[TP_CVP_MX]; double hwv[TP_CVP_MX];
            {
                const int use_digit = getenv("TP_CVP_BIT") == NULL;
                int cand[8], nc;
                fmpz_t ma, mb, g; fmpz_init(ma); fmpz_init(mb); fmpz_init(g);
                for (int j = d + 1; j <= D && nvars < TP_CVP_MX - 2; j++) {
                    nc = 0;
                    for (int c = 0; c < o->ncomp && nc < 8; c++) if (o->cval[c] < 0 && lvl[c] == j) cand[nc++] = c;
                    int done = 0;
                    if (use_digit && nc == 2) {
                        limbs2f(ma, o->mask + (size_t)cand[0] * (size_t)o->nw64, o->nw64);
                        limbs2f(mb, o->mask + (size_t)cand[1] * (size_t)o->nw64, o->nw64);
                        for (int al = 1; al <= 15 && !done; al++) for (int be = 1; be <= 15 && !done; be++) {
                            int x = al, y = be; while (y) { int z = x % y; x = y; y = z; }
                            if (x != 1) continue;
                            fmpz_mul_ui(t, ma, (ulong)be); fmpz_mul_ui(g, mb, (ulong)al);
                            if (!fmpz_equal(t, g)) continue;
                            if (!fmpz_divisible_si(ma, al)) continue;
                            fmpz_divexact_si(g, ma, al); fmpz_mod(&M[nvars], g, N);
                            alv[nvars][0] = 0; alv[nvars][1] = al; alv[nvars][2] = be; alv[nvars][3] = al + be;
                            nalv[nvars] = 4; hwv[nvars] = (al + be) / 2.0; nvars++; done = 1;
                        }
                    }
                    if (!done) for (int q = 0; q < nc && nvars < TP_CVP_MX - 2; q++) {
                        limbs2f(&M[nvars], o->mask + (size_t)cand[q] * (size_t)o->nw64, o->nw64);
                        fmpz_mod(&M[nvars], &M[nvars], N);
                        alv[nvars][0] = 0; alv[nvars][1] = 1; nalv[nvars] = 2; hwv[nvars] = 0.5; nvars++;
                    }
                }
                fmpz_clear(ma); fmpz_clear(mb); fmpz_clear(g);
                for (int i = 0; i < nvars; i++) { Wv[i] = lround(35.0 / hwv[i]); if (Wv[i] < 1) Wv[i] = 1; }
            }
            const int nu = nvars;
            int un = -1; for (int i = 0; i < nu; i++) if (fmpz_fdiv_ui(&M[i], 5)) { un = i; break; }
            if (un < 0 || nu < 6) { _fmpz_vec_clear(M, TP_CVP_MX); fmpz_clear(N); fmpz_clear(t); fmpz_clear(R); fmpz_clear(cinv); fmpz_clear(den); continue; }
            { fmpz_swap(&M[0], &M[un]);
              int ti; long tl; double td;
              ti = nalv[0]; nalv[0] = nalv[un]; nalv[un] = ti;
              for (int a = 0; a < 4; a++) { ti = alv[0][a]; alv[0][a] = alv[un][a]; alv[un][a] = ti; }
              tl = Wv[0]; Wv[0] = Wv[un]; Wv[un] = tl; td = hwv[0]; hwv[0] = hwv[un]; hwv[un] = td; }
            fmpz_invmod(cinv, &M[0], N);
            fmpz_mat_t A, Ai; fmpz_mat_init(A, nu, nu); fmpz_mat_init(Ai, nu, nu);
            fmpz_set(fmpz_mat_entry(A, 0, 0), N); fmpz_mul_si(fmpz_mat_entry(A, 0, 0), fmpz_mat_entry(A, 0, 0), Wv[0]);
            for (int i = 1; i < nu; i++) {
                fmpz_mul(t, &M[i], cinv); fmpz_neg(t, t); fmpz_mod(t, t, N);
                fmpz_mul_si(fmpz_mat_entry(A, i, 0), t, Wv[0]);
                fmpz_set_si(fmpz_mat_entry(A, i, i), Wv[i]);
            }
            fmpz_lll_t fl; fmpz_lll_context_init_default(fl);
            /* Cross-length basis cache (TP_CVP_CACHE=dir). The lattice is fixed by the first column and the weights of A
             * (the masks scaled by M_0^-1, which removes the length-dependent unit 16^c), so a stored reduced basis with the
             * identical key spans exactly the same lattice. Expensive reduction can then be paid once for all lengths. */
            const char *cdir = getenv("TP_CVP_CACHE");
            int cached = 0; char *ckey = NULL; size_t cklen = 0; char cpath[1024] = {0};
            if (cdir && !getenv("TP_CVP_SEED")) {
                FILE *ms = open_memstream(&ckey, &cklen);
                fprintf(ms, "v1 K=%d nu=%d bkz=%d dred=%d", K, nu, bkzb, getenv("TP_CVP_DRED") ? 1 : 0);
                for (int i = 0; i < nu; i++) {
                    char *e0 = fmpz_get_str(NULL, 16, fmpz_mat_entry(A, i, 0)); fprintf(ms, " %s/%ld/%d:", e0, Wv[i], nalv[i]); flint_free(e0);
                    for (int a2 = 0; a2 < nalv[i]; a2++) fprintf(ms, "%d,", alv[i][a2]);
                }
                fclose(ms);
                uint64_t hh = 1469598103934665603ULL; for (size_t z = 0; z < cklen; z++) hh = (hh ^ (unsigned char)ckey[z]) * 1099511628211ULL;
                snprintf(cpath, sizeof cpath, "%s/cvp_%016llx.mat", cdir, (unsigned long long)hh);
                FILE *cf = fopen(cpath, "r");
                if (cf) {
                    char *line = NULL; size_t cap = 0; ssize_t got = getline(&line, &cap, cf);
                    if (got > 0 && (size_t)(got - 1) == cklen && !memcmp(line, ckey, cklen)) {
                        fmpz_mat_t Cm; fmpz_mat_init(Cm, 0, 0);
                        if (fmpz_mat_fread(cf, Cm) > 0 && fmpz_mat_nrows(Cm) == nu && fmpz_mat_ncols(Cm) == nu) { fmpz_mat_set(A, Cm); cached = 1; }
                        fmpz_mat_clear(Cm);
                    }
                    free(line); fclose(cf);
                }
            }
            if (!cached) {
            if (getenv("TP_CVP_SEED")) {
                /* Discovery mode: rerandomize the basis (random unimodular row operations) before LLL, so that
                 * pruned searches with different seeds miss different solutions. The lattice is unchanged. */
                uint64_t sr = 0x2545F4914F6CDD1DULL ^ (uint64_t)strtoull(getenv("TP_CVP_SEED"), NULL, 10) * 0x9E3779B97F4A7C15ULL ^ (uint64_t)d;
                for (int it = 0; it < 4 * nu * nu; it++) {
                    sr ^= sr << 13; sr ^= sr >> 7; sr ^= sr << 17;
                    const int i = (int)(sr % (uint64_t)nu), j = (int)((sr >> 20) % (uint64_t)nu);
                    if (i == j) continue;
                    for (int k = 0; k < nu; k++) {
                        if ((sr >> 40) & 1) fmpz_add(fmpz_mat_entry(A, i, k), fmpz_mat_entry(A, i, k), fmpz_mat_entry(A, j, k));
                        else fmpz_sub(fmpz_mat_entry(A, i, k), fmpz_mat_entry(A, i, k), fmpz_mat_entry(A, j, k));
                    }
                }
            }
            fmpz_lll(A, NULL, fl);
            if (bkzb >= 2) {
                static __thread long Bi[TP_CVP_MX][TP_CVP_MX];
                int fits = 1;
                for (int i = 0; i < nu && fits; i++) for (int k = 0; k < nu; k++) { if (!fmpz_fits_si(fmpz_mat_entry(A, i, k))) { fits = 0; break; } Bi[i][k] = fmpz_get_si(fmpz_mat_entry(A, i, k)); }
                if (fits) { bkz_i64(Bi, nu, bkzb, 4); for (int i = 0; i < nu; i++) for (int k = 0; k < nu; k++) fmpz_set_si(fmpz_mat_entry(A, i, k), Bi[i][k]); }
            }
            if (getenv("TP_CVP_DRED")) {
                /* make the LAST Gram-Schmidt vectors long: LLL the dual basis, then invert and reverse.
                 * The search starts at the last index, so this shrinks the top of its tree. */
                fmpz_mat_t Dn, X; fmpz_t d1, d2; fmpz_init(d1); fmpz_init(d2);
                fmpz_mat_init(Dn, nu, nu); fmpz_mat_init(X, nu, nu);
                fmpz_mat_inv(X, d1, A);
                fmpz_mat_transpose(Dn, X);
                fmpz_lll(Dn, NULL, fl);
                fmpz_mat_inv(X, d2, Dn);
                for (int i = 0; i < nu; i++) for (int k = 0; k < nu; k++) {
                    fmpz_mul(t, fmpz_mat_entry(X, k, nu - 1 - i), d1);
                    fmpz_divexact(fmpz_mat_entry(A, i, k), t, d2);
                }
                fmpz_mat_clear(Dn); fmpz_mat_clear(X); fmpz_clear(d1); fmpz_clear(d2);
            }
            if (ckey) {                 /* store for every later length; write to a temp file, then rename (atomic) */
                char tmp[1100]; snprintf(tmp, sizeof tmp, "%s.%d.tmp", cpath, (int)getpid());
                FILE *cf = fopen(tmp, "w");
                if (cf) { fwrite(ckey, 1, cklen, cf); fputc('\n', cf); fmpz_mat_fprint(cf, A); fputc('\n', cf); fclose(cf); rename(tmp, cpath); }
            }
            }   /* !cached */
            free(ckey);
            /* A is the weighted basis; the unweighted one divides each column by its weight */
            fmpz_mat_t U; fmpz_mat_init(U, nu, nu);
            for (int i = 0; i < nu; i++) for (int k = 0; k < nu; k++) fmpz_divexact_si(fmpz_mat_entry(U, i, k), fmpz_mat_entry(A, i, k), Wv[k]);
            fmpz_mat_inv(Ai, den, U);
            {   /* doubles carry the enumeration: refuse any basis they cannot represent exactly */
                int big = 0;
                for (int i = 0; i < nu && !big; i++) for (int k = 0; k < nu; k++)
                    if (fmpz_bits(fmpz_mat_entry(A, i, k)) > 40) { big = 1; break; }
                if (big) {
                    if (dbg) fprintf(stderr, "  [cvp] d=%d undecided=%d K=%d skipped: basis entries too large\n", d, u, K);
                    fmpz_mat_clear(A); fmpz_mat_clear(Ai); fmpz_mat_clear(U); _fmpz_vec_clear(M, TP_CVP_MX);
                    fmpz_clear(N); fmpz_clear(t); fmpz_clear(R); fmpz_clear(cinv); fmpz_clear(den);
                    continue;
                }
            }
            cvplev *L = calloc(1, sizeof(cvplev));
            if (!L) { fprintf(stderr, "cvp: out of memory\n"); exit(2); }
            L->K = K; L->u = nu; f2u(&L->N, N); L->Nd = u256_dbl(&L->N);
            { int nl = (int)(fmpz_bits(N) + 64) / 64; if (nl > 4) nl = 4; L->nl = nl; }
            { uint64_t inv = 1; for (int i = 0; i < 6; i++) inv *= 2 - L->N.w[0] * inv; L->n0 = (uint64_t)0 - inv; }
            fmpz_one(R); fmpz_mul_2exp(R, R, (ulong)(64 * L->nl)); fmpz_mod(R, R, N);
            fmpz_mul_2exp(t, R, 64); fmpz_mod(t, t, N); f2u(&L->c64, t);
            fmpz_mul_ui(t, R, 10); fmpz_mod(t, t, N); f2u(&L->c10, t);
            double bs[TP_CVP_MX][TP_CVP_MX], Bw[TP_CVP_MX][TP_CVP_MX], cen[TP_CVP_MX];
            double R2 = 0;
            for (int k = 0; k < nu; k++) {
                L->nal[k] = nalv[k]; memcpy(L->al[k], alv[k], sizeof alv[k]);
                cen[k] = (double)Wv[k] * hwv[k]; R2 += cen[k] * cen[k];
            }
            L->R2 = R2 * (1.0 + 1e-9) + 1e-5;
            memcpy(L->cen, cen, sizeof(double) * (size_t)nu);
            for (int i = 0; i < nu; i++) {
                for (int k = 0; k < nu; k++) { L->B[i][k] = fmpz_get_d(fmpz_mat_entry(U, i, k)); Bw[i][k] = fmpz_get_d(fmpz_mat_entry(A, i, k)); }
                fmpz_mul(t, fmpz_mat_entry(Ai, 0, i), N); fmpz_divexact(t, t, den);
                fmpz_mul(t, t, cinv); fmpz_mod(t, t, N); fmpz_mul(t, t, R); fmpz_mod(t, t, N);
                f2u(&L->aM[i], t);
            }
            for (int i = 0; i < nu; i++) {
                memcpy(bs[i], Bw[i], sizeof(double) * (size_t)nu);
                for (int j = 0; j < i; j++) {
                    double a = 0; for (int k = 0; k < nu; k++) a += Bw[i][k] * bs[j][k];
                    L->mu[i][j] = a / L->n2[j];
                    for (int k = 0; k < nu; k++) bs[i][k] -= L->mu[i][j] * bs[j][k];
                }
                double a = 0; for (int k = 0; k < nu; k++) a += bs[i][k] * bs[i][k];
                L->n2[i] = a;
                double h = 0, l1 = 0; for (int k = 0; k < nu; k++) { h += cen[k] * bs[i][k]; l1 += cen[k] * fabs(bs[i][k]); }
                L->h[i] = h / a;
                L->wbox[i] = l1 / a + 1e-7;
                memcpy(L->bs[i], bs[i], sizeof(double) * (size_t)nu);
            }
            fmpz_mat_clear(A); fmpz_mat_clear(Ai); fmpz_mat_clear(U);
            for (int i = 0; i < nu; i++) L->rmax[i] = L->R2;
            if (getenv("TP_CVP_PRUNE")) {
                /* Learned pruning (FINDING mode, not a proof): sample genuine digit vectors, record how large
                 * their partial projected norm gets at each level, and search only up to that quantile. */
                const double q = atof(getenv("TP_CVP_PRUNE"));
                enum { NS = 40000 };
                double *rs = malloc(sizeof(double) * NS * (size_t)nu);
                if (!rs) { fprintf(stderr, "cvp: out of memory\n"); exit(2); }
                uint64_t lr = 0xC4CEB9FE1A85EC53ULL ^ (uint64_t)d;
                for (int a = 0; a < NS; a++) {
                    double v[TP_CVP_MX];
                    for (int k = 0; k < nu; k++) { lr ^= lr << 13; lr ^= lr >> 7; lr ^= lr << 17;
                        v[k] = (double)Wv[k] * alv[k][(lr >> 58) % (uint64_t)nalv[k]] - cen[k]; }
                    double acc = 0;
                    for (int i = nu - 1; i >= 0; i--) {
                        double c1 = 0; for (int k = 0; k < nu; k++) c1 += v[k] * L->bs[i][k];
                        acc += c1 * c1 / L->n2[i];
                        rs[(size_t)a * (size_t)nu + (size_t)i] = acc;
                    }
                }
                for (int i = 0; i < nu; i++) {
                    double *col = malloc(sizeof(double) * NS);
                    for (int a = 0; a < NS; a++) col[a] = rs[(size_t)a * (size_t)nu + (size_t)i];
                    qsort(col, NS, sizeof(double), cmp_double);
                    int idx = (int)(q * (NS - 1)); if (idx >= NS) idx = NS - 1;
                    L->rmax[i] = col[idx] * (1.0 + 1e-9) + 1e-6;
                    if (L->rmax[i] > L->R2) L->rmax[i] = L->R2;
                    free(col);
                }
                free(rs);
                for (int i = nu - 2; i >= 0; i--) if (L->rmax[i] < L->rmax[i + 1]) L->rmax[i] = L->rmax[i + 1];   /* monotone */
            }
            /* calibrate on random residues: cost and rejection rate */
            L->limit = bud * 40;
            uint64_t rng = 0x9E3779B97F4A7C15ULL ^ (uint64_t)d;
            long tot = 0, rej = 0, over = 0; const int ns = 64;
            for (int a = 0; a < ns; a++) {
                uint64_t x[4]; for (int b = 0; b < 4; b++) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; x[b] = rng; }
                u256 q; mod_reduce(&q, x, 4, &L->N, L->n0, &L->c64, L->nl);
                long nn = 0;
                const int kp = cvp_keep(L, &q, &nn);
                tot += nn; if (!kp) rej++; if (nn > L->limit) over++;
            }
            const double avg = (double)tot / ns;
            if (dbg) fprintf(stderr, "  [cvp] d=%d undecided=%d vars=%d K=%d avg nodes %.0f rejected %ld/%d limit-hit %ld\n", d, u, nu, K, avg, rej, ns, over);
            const int good = (avg <= (double)bud && rej >= ns * 3 / 4);
            if (good) {
                for (int k = 0; k < p->nviab[d]; k++) {
                    const int v = p->vtab[d][k];
                    bn_t set; bn_zero(&set, nw);
                    for (int i = 0; i < p->nn[d]; i++) {
                        const int site = p->nb[d][i];
                        const int c = o->comp[site];
                        const int bit = (v >> (site & 3)) & 1;
                        if (o->cval[c] >= 0 || lvl[c] != d) continue;
                        if (bit) { const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64; for (int w = 0; w < nw; w++) set.w[w] |= mk[w]; }
                    }
                    limbs2f(t, set.w, nw); fmpz_mod(t, t, N); f2u(&L->add[k], t);
                    for (int j = 0; j < nu; j++) u256_mont(&L->cadd[j][k], &L->aM[j], &L->add[k], &L->N, L->n0, L->nl);
                }
                if (g_off) { limbs2f(t, g_offc.w, TP_MAX_LIMBS); fmpz_mod(t, t, N); f2u(&L->coff, t); }
                if (g_cvchk) {          /* a genuine subset sum must always be found */
                    long miss = 0; const long nt = getenv("TP_CVPCHK_N") ? atol(getenv("TP_CVPCHK_N")) : 2000;
                    /* planted samples use their own generator, so every seed / pruning level sees the same plants */
                    uint64_t prng = 0xA24BAED4963EE407ULL ^ ((uint64_t)d << 32) ^ (uint64_t)K;
                    char *mbits = getenv("TP_CVPCHK_DUMP") ? calloc((size_t)nt + 1, 1) : NULL;
                    long enodes = 0;
                    for (long a = 0; a < nt; a++) {
                        fmpz_zero(t);
                        for (int j = 0; j < nu; j++) { prng ^= prng << 13; prng ^= prng >> 7; prng ^= prng << 17; int v = alv[j][(prng >> 60) % (uint64_t)nalv[j]]; if (v) fmpz_addmul_ui(t, &M[j], (ulong)v); }
                        fmpz_mod(t, t, N);
                        u256 q; f2u(&q, t);
                        long nn; const int kept = cvp_keep(L, &q, &nn); enodes += nn;
                        if (!kept) miss++;
                        if (mbits) mbits[a] = kept ? '1' : '0';
                    }
                    fprintf(stderr, "  [cvpchk] d=%d: %ld of %ld genuine sums rejected\n", d, miss, nt);
                    if (mbits) {                /* discovery calibration: which plants this run finds, and its enumeration cost */
                        FILE *df = fopen(getenv("TP_CVPCHK_DUMP"), "a");
                        if (df) { fprintf(df, "d=%d K=%d u=%d prune=%s seed=%s avgnodes=%.2f found=%s\n", d, K, nu,
                                          getenv("TP_CVP_PRUNE") ? getenv("TP_CVP_PRUNE") : "none", getenv("TP_CVP_SEED") ? getenv("TP_CVP_SEED") : "none",
                                          (double)enodes / (double)nt, mbits); fclose(df); }
                        free(mbits);
                    }
                    if (miss && !getenv("TP_CVP_PRUNE")) exit(9);
                }
                free(best); best = L; bestd = d;
            } else { free(L); }
            _fmpz_vec_clear(M, TP_CVP_MX); fmpz_clear(N); fmpz_clear(t); fmpz_clear(R); fmpz_clear(cinv); fmpz_clear(den);
            if (!good && best) break;              /* costs only grow as we move up the tree */
        }
        if (best) g_cvc.lev[bestd] = best;
        g_cvc.key_n = o->nsites; g_cvc.key_h = oh; g_cvc.bud = bud; g_cvc.K = Kov;
    }
    for (int d = 0; d < 64; d++) p->cvp[d] = g_cvc.lev[d];
    pthread_mutex_unlock(&g_cv_mu);
#endif
}

static void bf_build(hexplan *p, const tp_orbits *o, const int *lvl)
{
    p->bK = 0; memset(p->bat, 0, sizeof p->bat); p->bfp = NULL;
    const char *spec = getenv("TP_BLOOM");
    if (!spec || !*spec) return;
    if (!(g_fastk && g_carry9 && !g_xm && !g_decchk && g_nohi && !g_noprune && !g_badprune && !g_r5chk)) return;
    const int D = p->depth_max;
    const int dbg = getenv("TP_BLOOMDBG") != NULL;
    if (o->nsites % 4 == 1 || o->nsites % 4 == 2) return;   /* half the free bits: these trees are tiny */
    int nfree = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0) nfree++;

    /* parse the levels; a level exists only if some orbits remain undecided above it */
    int Us[TP_BF_MAX], bpi[TP_BF_MAX], nu = 0;
    { char b[64]; snprintf(b, sizeof b, "%s", spec);
      for (char *t = strtok(b, ","); t && nu < TP_BF_MAX; t = strtok(NULL, ",")) {
          int a = 0, c = 0; if (sscanf(t, "%d:%d", &a, &c) == 2 && a > 0 && c > 0) { Us[nu] = a; bpi[nu] = c; nu++; }
      } }
    int lev[TP_BF_MAX], ulev[TP_BF_MAX], n = 0, Kmin = 25;
    for (int i = 0; i < nu; i++) {
        int dl = -1, uu = 0;
        for (int d = 1; d < D && d < 64; d++) {
            int u = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0 && lvl[c] > d) u++;
            if (u <= Us[i]) { dl = d; uu = u; break; }
        }
        if (dl < 1 || uu < 8) continue;
        /* digits the parent of this level always shares: dl_min - k - 1 (one lost to a carry) */
        const int dmin = (int)floor((double)(o->nsites - 1) * 0.30102999566) + 1;
        int K = dmin - p->ddig[dl - 1] - 2;
        if (getenv("TP_BLOOM_K")) K = atoi(getenv("TP_BLOOM_K"));
        if (K > 25) K = 25;
        if (K < 2 || 2.3219 * K < uu + 3) continue;   /* set would be too dense to bite */
        if (K < Kmin) Kmin = K;
        lev[n] = dl; ulev[n] = uu; bpi[n] = bpi[i]; n++;
    }
    if (!n) return;
    const int K = Kmin;
    uint64_t N = 1; for (int i = 0; i < K; i++) N *= 5;
    p->bK = K; p->bN = N;

    pthread_mutex_lock(&g_bf_mu);
    const uint64_t oh = r5_orbit_hash(o);
    if (g_bfc.key_n != o->nsites || g_bfc.key_h != oh || g_bfc.K != K || strcmp(g_bfc.spec, spec)) {
        for (int i = 0; i < TP_BF_MAX; i++) { free(g_bfc.f[i].w); memset(&g_bfc.f[i], 0, sizeof(bfilter)); }
        g_bfc.key_n = -1; g_bfc.n = 0;
        for (int i = 0; i < n; i++) {
            uint64_t M[256]; int u = 0;
            for (int c = 0; c < o->ncomp; c++) {
                if (o->cval[c] >= 0 || lvl[c] <= lev[i]) continue;
                bn_t b; bn_zero(&b, o->nw64);
                memcpy(b.w, o->mask + (size_t)c * (size_t)o->nw64, (size_t)o->nw64 * sizeof(uint64_t));
                if (u < 256) M[u++] = bn_divmod_small(&b, o->nw64, N);
            }
            const double bytes = ldexp(1.0, u) * bpi[i] / 8.0;
            uint64_t nb = (uint64_t)(bytes / 32.0); if (nb < 16) nb = 16;
            void *mem = NULL;
            if (posix_memalign(&mem, 64, (size_t)nb * 32) || !mem) { fprintf(stderr, "bloom: out of memory\n"); exit(2); }
            memset(mem, 0, (size_t)nb * 32);
            bfilter *f = &g_bfc.f[i]; f->w = mem; f->nb = nb; f->d = lev[i]; f->u = u;
            bf_job J; J.f = f; J.M = M; J.u = u; J.p = (u > 12) ? 10 : 0; J.N = N; J.next = 0;
            long nc = sysconf(_SC_NPROCESSORS_ONLN); if (nc < 1) nc = 1; if (nc > 32) nc = 32;
            pthread_t th[32];
            for (long t = 0; t < nc; t++) pthread_create(&th[t], NULL, bf_worker, &J);
            for (long t = 0; t < nc; t++) pthread_join(th[t], NULL);
            if (dbg) fprintf(stderr, "  [bloom] d=%d undecided=%d K=%d %.2f GB\n", f->d, u, K, (double)nb * 32 / 1e9);
            if (g_bchk) {            /* every genuine subset sum must be present */
                uint64_t rng = 0x9E3779B97F4A7C15ULL, miss = 0, tot = 2000000;
                for (uint64_t a = 0; a < tot; a++) {
                    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17;
                    uint64_t x = 0, r2 = rng;
                    for (int j = 0; j < u; j++) { if (r2 & 1) { x += M[j]; if (x >= N) x -= N; } r2 = (r2 >> 1) | (r2 << 63); }
                    if (!bf_has(f, x)) miss++;
                }
                uint64_t fp = 0;
                for (uint64_t a = 0; a < tot; a++) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; if (bf_has(f, rng % N)) fp++; }
                fprintf(stderr, "  [bloomchk] d=%d: %llu of %llu real sums missing; random-residue hit rate %.4f\n",
                        f->d, (unsigned long long)miss, (unsigned long long)tot, (double)fp / (double)tot);
                if (miss) exit(9);
            }
        }
        g_bfc.key_n = o->nsites; g_bfc.key_h = oh; g_bfc.K = K; g_bfc.n = n;
        snprintf(g_bfc.spec, sizeof g_bfc.spec, "%s", spec);
    }
    p->bfp = g_bfc.f;
    for (int i = 0; i < g_bfc.n; i++) p->bat[g_bfc.f[i].d] = i + 1;
    pthread_mutex_unlock(&g_bf_mu);
}

static int hexplan_prepare(hexplan *p, const tp_orbits *o)
{
    /* every entry point (serial, threaded, sharded) builds a plan, so the
       debug / feature flags are read here -- previously the threaded path read
       only TP_NOPRUNE, so TP_NOMOD11, TP_HIDEC, TP_XM... silently did nothing
       with --threads > 1. */
    g_noprune  = getenv("TP_NOPRUNE") != NULL;
    g_nohi     = getenv("TP_HIDEC") == NULL;
    g_badprune = getenv("TP_BADPRUNE") != NULL;
    g_xm       = getenv("TP_XM") != NULL;
    g_decchk   = getenv("TP_DECCHK") != NULL;
    g_mod11    = getenv("TP_NOMOD11") == NULL;
    g_carry9 = getenv("TP_NOCARRY9") == NULL;
    g_fastk = getenv("TP_NOFASTK") == NULL;
    g_fastchk = getenv("TP_FASTCHK") != NULL;
    g_bchk = getenv("TP_BLOOMCHK") != NULL;
    g_duchk = getenv("TP_DUALCHK") != NULL;
    g_noshare = getenv("TP_NOSHARE") != NULL;
    g_cvchk = getenv("TP_CVPCHK") != NULL;
    g_depthhist = getenv("TP_DEPTH") != NULL;
    g_dsieve = getenv("TP_SIEVE") != NULL;
    {   /* default modulus 5^13 needs ~3.5 GB of reach tables; fall back to 5^12 on small machines */
        int defJ = TP_R5_DEFAULT;
        const long pg = sysconf(_SC_PHYS_PAGES), ps = sysconf(_SC_PAGESIZE);
        if (pg > 0 && ps > 0 && (double)pg * (double)ps >= 12.0e9 && o->nsites >= 372) defJ = 13;
        g_r5J = getenv("TP_R5") ? atoi(getenv("TP_R5")) : defJ;
    }
    if (g_r5J > 0) {
        /* no point in a modulus much larger than the number of free orbits can
           fill: with u free orbits the reach set has at most 2^u residues */
        int u0 = 0;
        for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0) u0++;
        int jmax = (int)((u0 + 3) / 2.3219) + 1;
        if (jmax < 2) jmax = 2;
        if (g_r5J > jmax) g_r5J = jmax;
    }
    g_r5chk = getenv("TP_R5CHK") != NULL;
    g_weakJ = getenv("TP_WEAKPAL") ? atoi(getenv("TP_WEAKPAL")) : 0;
    /* parsed exactly once: every worker thread prepares its own plan, and re-parsing here would zero the shared
     * offset while other threads are using it */
    pthread_once(&g_off_once, off_parse);
    if (g_off) { g_mod11 = 0; g_dsieve = 0; }   /* these two reason about N's own digits: off (they only ever prune) */
    if (g_r5J > 13) g_r5J = 13;
    p->r5J = g_r5J; p->r5N = 1;
    for (int i = 0; i < p->r5J; i++) p->r5N *= 5;
    const int nw = o->nw64, D = p->depth_max;
    p->setm = calloc((size_t)D * 16 * (size_t)nw, sizeof(uint64_t));
    p->clrm = calloc((size_t)D * 16 * (size_t)nw, sizeof(uint64_t));
    if (!p->setm || !p->clrm) { hexplan_free(p); return -1; }

    for (int d = 0; d < D; d++) {
        int nv = 0;
        for (int v = 0; v < 16; v++) {
            if (d == 0 && g_dsieve && p->L <= TP_VIABLE_MAX_M) {
                if (!(g_viable_outer[p->L] & (1 << v))) continue;
            }
            int orb[8], val[8];
            const int n = hex_implications(o, p, d, v, orb, val);
            if (n < 0) continue;                    /* internally inconsistent */
            int ok = 1;
            for (int i = 0; i < n; i++)
                if (o->cval[orb[i]] >= 0 && o->cval[orb[i]] != val[i]) { ok = 0; break; }
            if (!ok) continue;                      /* clashes with a forced bit */
            uint64_t *sm = p->setm + ((size_t)d * 16 + (size_t)nv) * (size_t)nw;
            uint64_t *cm = p->clrm + ((size_t)d * 16 + (size_t)nv) * (size_t)nw;
            for (int i = 0; i < n; i++) {
                const uint64_t *mk = o->mask + (size_t)orb[i] * (size_t)o->nw64;
                uint64_t *t = val[i] ? sm : cm;
                for (int w = 0; w < nw; w++) t[w] |= mk[w];
            }
            p->vtab[d][nv++] = v;
        }
        p->nviab[d] = nv;
    }

    /* Decimal form of what each (level, digit) changes, so the child's
       decimal state is one add / one subtract away.  Restricted to components
       whose FIRST assignment is this level: a component re-touched later is
       already in the value, and its OR is a no-op that must not be counted
       twice.  At the last level this is also the offset that completes the
       number, so no second endpoint and no interval test are needed there. */
    /* digits of the interval width at each level: hi - lo = 2^(nsites-4(d+1))
       - 2^(4(d+1)), a function of the level only. */
    for (int d = 0; d < D && d < 64; d++) {
        int r = 4 * (d + 1);
        if (r > o->nsites) r = o->nsites;
        int k = o->nsites - r;
        if (k < 0) k = 0;
        p->ddig[d] = (int)(0.301029995663981195 * (double)k) + 1;
    }

    int lvl[TP_MAX_BITS];
    hexplan_levels(p, o, lvl);
    if (getenv("TP_PLANDBG")) {
        fprintf(stderr, "  plan: nsites=%d L=%d depth_max=%d nw=%d\n", o->nsites, p->L, D, nw);
        for (int d = 0; d < D; d++) {
            fprintf(stderr, "   d=%2d nn=%d nviab=%2d  v:", d, p->nn[d], p->nviab[d]);
            for (int k = 0; k < p->nviab[d]; k++) fprintf(stderr, " %d", p->vtab[d][k]);
            fprintf(stderr, "\n");
        }
    }
    /* Per level: the interval width hi-lo = 2^(nsites-r) - 2^r depends only on
     * the level, so its decimal form is a constant.  Store k = digits(width)
     * and R = 10^k - width (as exactly k digits, leading zeros kept): then
     * lo + width >= 10^k  <=>  (lo mod 10^k) >= R, a cheap digit compare. */
    p->dstr = calloc((size_t)D * (TP_DECBUF + 1), 1);
    if (!p->dstr) { hexplan_free(p); return -1; }
    for (int d = 0; d < D && d < 64; d++) {
        /* exact width: lo carries 0 and hi carries 1 in every bit that is not
         * yet pinned, so hi - lo is the union of the masks of the sites whose
         * level has not been reached (padded top nibbles included). */
        bn_t delta;
        bn_zero(&delta, nw);
        for (int c = 0; c < o->ncomp; c++) {
            if (o->cval[c] >= 0 || lvl[c] <= d) continue;
            const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64;
            for (int w = 0; w < nw; w++) delta.w[w] |= mk[w];
        }
        char dbuf[TP_DECBUF];
        int kd = 1;
        if (!bn_is_zero(&delta, nw)) kd = bn_to_dec(&delta, nw, dbuf);
        else { kd = 0; dbuf[0] = 0; }
        char *R = p->dstr + (size_t)d * (TP_DECBUF + 1);
        p->ddig[d] = kd;
        /* R = 10^kd - delta, digitwise: 9's complement then +1 */
        for (int i = 0; i < kd; i++) R[i] = (char)('0' + (9 - (dbuf[i] - '0')));
        for (int i = kd - 1; i >= 0; i--) {
            if (R[i] < '9') { R[i]++; break; }
            R[i] = '0';
        }
        R[kd] = 0;
        if (d < 64) {                    /* chunked copy for the fast kernel */
            const int nfull = kd / 9, rem = kd % 9;
            for (int j = 0; j < nfull; j++) {
                uint32_t v = 0;
                for (int i = kd - 9 * (j + 1); i < kd - 9 * j; i++) v = v * 10u + (uint32_t)(R[i] - '0');
                p->thr[d][j] = v;
            }
            uint32_t v = 0;
            for (int i = 0; i < rem; i++) v = v * 10u + (uint32_t)(R[i] - '0');
            p->thrpart[d] = v;
        }
    }
    pthread_once(&g_p10_once, p10w_init);

    p->dlo9 = calloc((size_t)D * 16, sizeof(dec9));
    if (!p->dlo9) { hexplan_free(p); return -1; }
    for (int d = 0; d < D; d++) {
        for (int k = 0; k < p->nviab[d]; k++) {
            const int v = p->vtab[d][k];
            bn_t set;
            bn_zero(&set, nw);
            for (int i = 0; i < p->nn[d]; i++) {
                const int site = p->nb[d][i];
                const int c    = o->comp[site];
                const int bit  = (v >> (site & 3)) & 1;
                if (o->cval[c] >= 0 || lvl[c] != d) continue;
                const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64;
                if (bit) for (int w = 0; w < nw; w++) set.w[w] |= mk[w];
            }
            dec9_from_bn(&set, nw, &p->dlo9[(size_t)d * 16 + k]);
            {   /* N mod 11 of the same bits, for the leaf sieve */
                bn_t t = set;
                p->add11[d][k] = (int)bn_divmod_small(&t, nw, 11);
            }
            if (p->r5J) {
                bn_t t = set;
                p->add5[d][k] = bn_divmod_small(&t, nw, p->r5N);
            }
        }
    }
    if (p->r5J) r5_build(p, o, lvl);
    bf_build(p, o, lvl);
    dual_build(p, o, lvl);
    cvp_build(p, o, lvl);
    if (p->bK) for (int d = 0; d < D && d < 64; d++)
        for (int k = 0; k < p->nviab[d]; k++) {
            const int v = p->vtab[d][k];
            bn_t set; bn_zero(&set, nw);
            for (int i = 0; i < p->nn[d]; i++) {
                const int site = p->nb[d][i];
                const int c    = o->comp[site];
                const int bit  = (v >> (site & 3)) & 1;
                if (o->cval[c] >= 0 || lvl[c] != d) continue;
                if (bit) { const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64;
                           for (int w = 0; w < nw; w++) set.w[w] |= mk[w]; }
            }
            p->addk[d][k] = bn_divmod_small(&set, nw, p->bN);
        }
    memcpy(g_r5lvl, lvl, sizeof(int) * (size_t)o->ncomp);
    return 0;
}

static void hexplan_free(hexplan *p)
{
    free(p->setm);   p->setm = NULL;
    free(p->clrm);   p->clrm = NULL;
    free(p->dlo9);   p->dlo9 = NULL;
    free(p->dstr);   p->dstr = NULL;
}

/* Digit-sum sieve (ported from v13).  For an N whose hex form is a palindrome
 * with outer digit s_outer and hex digit sum S, N mod 100 is determined by
 * (s_outer mod 4) and S mod 25 (CRT).  Once the first two decimal digits of the
 * interval are pinned (lead2), they must equal the reversed last two digits. */
static inline int is_lead2_possible(int hex_len_m, int s_outer, int S_fixed, int d, int lead2)
{
    const int rem = hex_len_m - 2 * (d + 1);
    if (rem <= 0) {
        int r4 = s_outer % 4;
        int r25 = (S_fixed * (6 - 5 * hex_len_m)) % 25;
        if (r25 < 0) r25 += 25;
        int r100 = (25 * r4 + 76 * r25) % 100;
        if ((r100 % 10) % 2 == 0) return 0;
        int p = (r100 % 10) * 10 + (r100 / 10);
        return (p == lead2);
    }

    const int max_d = (s_outer == 1) ? 1 : (s_outer == 3) ? 3 : (s_outer == 5 || s_outer == 7) ? 7 : 15;
    const int min_S = S_fixed;
    const int max_S = S_fixed + rem * max_d;
    if ((max_S - min_S) >= 50) return 1;

    const int step = (rem % 2 == 0) ? 2 : 1;
    const int C = ((6 - 5 * hex_len_m) % 25 + 25) % 25;
    const int r4_term = 25 * (s_outer % 4);

    for (int S = min_S; S <= max_S; S += step) {
        int r25 = (S * C) % 25;
        int r100 = (r4_term + 76 * r25) % 100;
        int d0 = r100 % 10;
        if (d0 & 1) {
            int p = d0 * 10 + (r100 / 10);
            if (p == lead2) return 1;
        }
    }
    return 0;
}


/* TP_R5CHK: at this node, draw random completions of the undecided orbits,
 * build the REAL number with bigint arithmetic, and assert that its residue
 * mod 5^J lies in the reach set the filter is about to rely on.  */
static uint64_t g_rng = 88172645463325252ull;
static inline uint64_t rng64(void) { g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17; return g_rng; }
static void r5_selfcheck(const tp_search *s, const hexplan *p, int d)
{
    const tp_orbits *o = s->o; const int nw = s->nw;
    for (int trial = 0; trial < 6; trial++) {
        bn_t n = s->lo;
        for (int c = 0; c < o->ncomp; c++) {
            if (o->cval[c] >= 0 || g_r5lvl[c] <= d) continue;
            if (trial == 0 || (rng64() & 1)) {
                const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64;
                for (int w = 0; w < nw; w++) n.w[w] |= mk[w];
            }
        }
        const uint64_t resN = bn_divmod_small(&n, nw, p->r5N);
        const uint64_t need = (resN >= s->r5) ? resN - s->r5 : resN + p->r5N - s->r5;
        g_r5chk_n++;
        if (p->reach5[d] && !((p->reach5[d][need >> 6] >> (need & 63)) & 1)) {
            g_r5chk_bad++;
            if (g_r5chk_bad < 4) {
                bn_t t2 = s->lo; const uint64_t lr = bn_divmod_small(&t2, nw, p->r5N);
                int und = 0; for (int c = 0; c < o->ncomp; c++) if (o->cval[c] < 0 && g_r5lvl[c] > d) und++;
                fprintf(stderr, "  [r5chk] FAIL B=%d d=%d D=%d trial=%d r5=%llu lo%%N=%llu undecided=%d\n", o->nsites, d, p->depth_max, trial,
                        (unsigned long long)s->r5, (unsigned long long)lr, und);
            }
        }
    }
}

/* TP_R5CHK: every completion of this node must share the first J decimal digits
 * the filter read off lo -- that is what lets it turn them into a target.     */
static void r5_prefixcheck(const tp_search *s, const hexplan *p, int d, const char *dlos)
{
    const tp_orbits *o = s->o; const int nw = s->nw;
    for (int trial = 0; trial < 4; trial++) {
        bn_t n = s->lo;
        for (int c = 0; c < o->ncomp; c++) {
            if (o->cval[c] >= 0 || g_r5lvl[c] <= d) continue;
            if (trial == 0 || (rng64() & 1)) {
                const uint64_t *mk = o->mask + (size_t)c * (size_t)o->nw64;
                for (int w = 0; w < nw; w++) n.w[w] |= mk[w];
            }
        }
        char ds[TP_DECBUF];
        bn_to_dec(&n, nw, ds);
        g_r5chk_pn++;
        if (memcmp(ds, dlos, (size_t)p->r5J) != 0) g_r5chk_pbad++;
    }
}

static __attribute__((noinline)) void
hexdfs(tp_search *s, const hexplan *p, int d, const dec9 *dnp,
       const uint64_t *Rp, int mp, int s_sum, int s_outer,
       tp_hit_fn onhit, void *ud)
{
    const tp_orbits *o = s->o;
    const int nw = s->nw;

    if (s->budget && s->st.nodes >= s->budget) { s->st.aborted = 1; return; }

    if (p->tasks && d >= p->cut) {
        if (*p->ntask >= p->cap) { s->st.overflow = 1; return; }
        p->tasks[*p->ntask].lo    = s->lo;
        p->tasks[*p->ntask].hi    = s->hi;
        p->tasks[*p->ntask].depth = d;
        p->tasks[*p->ntask].s_sum   = s_sum;
        p->tasks[*p->ntask].s_outer = s_outer;
        (*p->ntask)++;
        return;
    }

    /* the low 4*(d+1) bits are contiguous and fixed once this level is done */
    int r = 4 * (d + 1);
    if (r > o->nsites) r = o->nsites;

    const int is_last = (d == p->depth_max - 1);
    const int nv = p->nviab[d];
    const uint64_t *smv = p->setm + (size_t)d * 16 * (size_t)nw;
    const uint64_t *cmv = p->clrm + (size_t)d * 16 * (size_t)nw;
    const dec9 *dlop = p->dlo9 + (size_t)d * 16;

    char dlos[TP_DECBUF];

    /* This node's own decimal value, taken BEFORE any child's bits are OR'd into
     * lo.  (An earlier revision built it lazily inside the loop, i.e. after lo
     * already carried the child's bits, and then added the child's delta on top:
     * the leaf's decimal was wrong and the palindrome test ran on the wrong
     * number.  TP_DECCHK now guards this.)  One conversion per parent of leaves. */
    dec9 dn;
    if (is_last) { if (g_carry9) dn = *dnp; else dec9_from_bn(&s->lo, nw, &dn); }
    const uint32_t r11p = s->r11;         /* this node's lo mod 11 */
    const uint64_t r5p  = s->r5;          /* this node's lo mod 5^J */
    const int      t5ok_p = s->t5ok;      /* parent's shared-prefix target, if it had one */
    const uint64_t t5_p   = s->t5;
    const uint64_t rk_p   = s->rk, tk_p = s->tk;
    const int      tkok_p = s->tkok;
    if (g_bchk && p->bK) {
        bn_t tt = s->lo; const uint64_t tr = bn_divmod_small(&tt, nw, p->bN);
        __atomic_fetch_add(&g_bchk_n, 1, __ATOMIC_RELAXED);
        if (tr != rk_p) { fprintf(stderr, "BLOOMCHK carried residue wrong d=%d\n", d); exit(9); }
        if (tkok_p) {
            char v2[TP_DECBUF]; (void)bn_to_dec(&s->lo, nw, v2);
            uint64_t t = 0;
            for (int i = p->bK - 1; i >= 0; i--) t = (t * 10u + (uint64_t)(v2[i] - '0')) % p->bN;
            __atomic_fetch_add(&g_bchk_t, 1, __ATOMIC_RELAXED);
            if (t != tk_p) { fprintf(stderr, "BLOOMCHK target wrong d=%d\n", d); exit(9); }
        }
    }

    /* Batched, prefetched child filter.  When this node already has a fixed
     * 5-adic target (t5ok_p) every child is decided by one probe into the reach
     * table, which lives in DRAM; issuing all the probes first lets the misses
     * overlap, and a rejected child never touches the 8-word lo/hi masks. */
    const int fastk = g_fastk && g_carry9 && !g_xm && !g_decchk && g_nohi && !g_noprune
                      && !g_badprune && !g_r5chk;
    uint8_t  fok[16];
    uint64_t fx[16];
    const int fast5 = (p->r5J && t5ok_p && !is_last && p->reach5[d] && !g_r5chk);
    if (fast5) {
        const uint64_t *rb = p->reach5[d];
        const uint64_t N5 = p->r5N;
        uint64_t fneed[16];
        for (int k = 0; k < nv; k++) {
            uint64_t x = r5p + p->add5[d][k]; if (x >= N5) x -= N5;
            fx[k] = x;
            const uint64_t need = (t5_p >= x) ? t5_p - x : t5_p + N5 - x;
            fneed[k] = need;
            __builtin_prefetch(&rb[need >> 6], 0, 0);
        }
        for (int k = 0; k < nv; k++) fok[k] = (uint8_t)((rb[fneed[k] >> 6] >> (fneed[k] & 63)) & 1);
    }

    /* dual-lattice certificate: needs only this node's lo and reversed prefix */
    uint8_t dok[16];
    const duallev *dlv = (d < 64 && fastk && !is_last && p->dual[d] && mp >= p->dual[d]->K) ? (const duallev *)p->dual[d] : NULL;
    if (dlv) {
        u256 a, b, base;
        dual_reduce(&a, Rp, nw, dlv);
        dual_reduce(&b, s->lo.w, nw, dlv);
        u256_submod(&base, &a, &b, &dlv->N);
        if (g_off) u256_submod(&base, &base, &dlv->coff, &dlv->N);
        if (g_duchk) {          /* the target must equal the one read off the decimal string */
            char v2[TP_DECBUF]; (void)bn_to_dec(&s->lo, nw, v2);
            u256 acc = {{0,0,0,0}}, tm, dg = {{0,0,0,0}}, b2;
            const int smalln = !(dlv->N.w[1] | dlv->N.w[2] | dlv->N.w[3]);
            for (int i = dlv->K - 1; i >= 0; i--) {
                u256_mont(&tm, &acc, &dlv->c10, &dlv->N, dlv->n0, dlv->nl);
                dg.w[0] = (uint64_t)(v2[i] - '0'); if (smalln) dg.w[0] %= dlv->N.w[0];
                u256_addmod(&acc, &tm, &dg, &dlv->N);
            }
            u256_submod(&b2, &acc, &b, &dlv->N);
            if (memcmp(&b2, &base, sizeof base)) { fprintf(stderr, "DUALCHK target mismatch d=%d\n", d); exit(9); }
        }
        u256 dzb[TP_DUAL_TMAX]; uint8_t dzok[TP_DUAL_TMAX]; memset(dzok, 0, (size_t)dlv->nt);
        for (int k = 0; k < nv; k++) {
            if (g_noshare || dlv->ns || k >= dlv->nadd) { u256 nd; u256_submod(&nd, &base, &dlv->add[k], &dlv->N); dok[k] = (uint8_t)dual_pass(dlv, &nd); }
            else dok[k] = (uint8_t)dual_pass_child(dlv, &base, k, dzb, dzok);
        }
    }

    /* deep filter: same two-pass shape, one cache line per child */
    uint8_t bok[16];
    const int bfi = (d < 64 && !is_last && tkok_p && p->bK) ? p->bat[d] : 0;
    if (bfi) {
        const bfilter *f = &((const bfilter *)p->bfp)[bfi - 1];
        const uint64_t N = p->bN;
        uint64_t nd[16];
        for (int k = 0; k < nv; k++) {
            uint64_t x = rk_p + p->addk[d][k]; if (x >= N) x -= N;
            nd[k] = (tk_p >= x) ? tk_p - x : tk_p + N - x;
            __builtin_prefetch(bf_block(f, bf_mix(nd[k])), 0, 0);
        }
        for (int k = 0; k < nv; k++) bok[k] = (uint8_t)bf_has(f, nd[k]);
    }

    /* exact lattice enumeration: same inputs as the certificate; run last, only on survivors */
    uint8_t cok[16];
    const cvplev *clv = (d < 64 && fastk && !is_last && p->cvp[d] && mp >= p->cvp[d]->K) ? (const cvplev *)p->cvp[d] : NULL;
    if (clv) {
        u256 a, b, base, czb[TP_CVP_MX]; uint8_t czok[TP_CVP_MX]; memset(czok, 0, (size_t)clv->u);
        mod_reduce(&a, Rp, nw, &clv->N, clv->n0, &clv->c64, clv->nl);
        mod_reduce(&b, s->lo.w, nw, &clv->N, clv->n0, &clv->c64, clv->nl);
        u256_submod(&base, &a, &b, &clv->N);
        if (g_off) u256_submod(&base, &base, &clv->coff, &clv->N);
        for (int k = 0; k < nv; k++) {
            /* cheapest tests first: a child already rejected never pays for an enumeration */
            if ((dlv && !dok[k]) || (bfi && !bok[k]) || (fast5 && !fok[k])) { cok[k] = 1; continue; }
            long nn = 0;
            if (g_noshare) { u256 nd; u256_submod(&nd, &base, &clv->add[k], &clv->N); cok[k] = (uint8_t)cvp_keep(clv, &nd, &nn); }
            else cok[k] = (uint8_t)cvp_keep_child(clv, &base, k, czb, czok, &nn);
            {
              if (g_cvstat) { __atomic_fetch_add(&g_cv_q, 1, __ATOMIC_RELAXED); __atomic_fetch_add(&g_cv_n, (uint64_t)nn, __ATOMIC_RELAXED); if (!cok[k]) __atomic_fetch_add(&g_cv_r, 1, __ATOMIC_RELAXED); } }
        }
    }

    for (int k = 0; k < nv; k++) {
        if (clv && !cok[k]) {
            s->st.nodes++; s->st.prunes++;
            if (g_depthhist) { g_dnode[d]++; g_dprune[d]++; }
            continue;
        }
        if (dlv && !dok[k]) {
            s->st.nodes++; s->st.prunes++;
            if (g_depthhist) { g_dnode[d]++; g_dprune[d]++; }
            continue;
        }
        if (bfi && !bok[k]) {
            s->st.nodes++; s->st.prunes++;
            if (g_depthhist) { g_dnode[d]++; g_dprune[d]++; }
            continue;
        }
        if (fast5 && !fok[k]) {
            s->st.nodes++; s->st.prunes++;
            if (g_depthhist) { g_dnode[d]++; g_dprune[d]++; }
            continue;
        }
        const int v = p->vtab[d][k];
        const int cur_s_outer = (d == 0) ? v : s_outer;
        const int mult = (is_last && (p->L % 2 == 1)) ? 1 : 2;
        const int cur_s_sum = (d == 0) ? ((p->L == 1) ? v : 2 * v) : (s_sum + mult * v);
        const uint64_t *sm = smv + (size_t)k * (size_t)nw;
        const uint64_t *cm = cmv + (size_t)k * (size_t)nw;
#ifdef TP_PROF2
        uint64_t t0 = PC();
#endif
        for (int w = 0; w < nw; w++) { s->lo.w[w] |= sm[w]; s->hi.w[w] &= ~cm[w]; }
#ifdef TP_PROF2
        uint64_t t1 = PC(); gP_mask += t1 - t0;
#endif
        s->st.nodes++;
        if (g_depthhist) g_dnode[d]++;
        const uint32_t r11c = (uint32_t)((r11p + (uint32_t)p->add11[d][k]) % 11u);
        s->r11 = r11c;
        if (p->r5J) { uint64_t x = r5p + p->add5[d][k]; if (x >= p->r5N) x -= p->r5N; s->r5 = x;
                      if (g_r5chk) r5_selfcheck(s, p, d);
                      /* Intervals nest, so a descendant shares every leading digit its
                         ancestor shared: once a node has J shared digits its subtree's
                         target is fixed, and the filter needs no decimal conversion. */
                      if (t5ok_p && !is_last && p->reach5[d]) {
                          const uint64_t need = (t5_p >= x) ? t5_p - x : t5_p + p->r5N - x;
                          if (!((p->reach5[d][need >> 6] >> (need & 63)) & 1)) {
                              s->st.prunes++;
                              if (g_depthhist) g_dprune[d]++;
                              for (int w = 0; w < nw; w++) { s->lo.w[w] &= ~sm[w]; s->hi.w[w] |= cm[w]; }
                              continue;
                          }
                      } }

        if (is_last) {
#ifdef TP_PROF2
            uint64_t q0 = PC();
#endif
            s->r11 = r11c;
            /* A decimal palindrome with an even number of digits is divisible
             * by 11.  N mod 11 is known from the hex digits alone (one add per
             * level), so a leaf that fails it cannot be a palindrome and never
             * needs its decimal form built at all. */
            if (g_mod11 && !g_off) {
                const uint32_t r = s->r11;
                if (g_decchk) {          /* the sieve must never reject a multiple of 11 */
                    bn_t t2 = s->lo; const uint32_t tr = (uint32_t)bn_divmod_small(&t2, nw, 11);
                    if (tr != r) { fprintf(stderr, "MOD11 carried residue wrong at d=%d: %u vs %u\n", d, r, tr); exit(9); }
                }
                if (r && (dec_digits_of(&s->lo, nw) & 1) == 0) {
                    s->st.prunes++;
                    if (g_depthhist) g_dprune[d]++;
                    for (int w = 0; w < nw; w++) { s->lo.w[w] &= ~sm[w]; s->hi.w[w] |= cm[w]; }
                    continue;
                }
            }
            dec9 dc;
            dec9_add(&dc, &dn, &dlop[k]);
            const int dl = dec9_to_str(&dc, dlos);
#ifdef TP_PROF2
            uint64_t q1 = PC(); gP_lastconv += q1 - q0; gP_nlast++;
#endif
            if (g_decchk) {
                char vv[TP_DECBUF];
                const int vl = bn_to_dec(&s->lo, nw, vv);
                if (vl != dl || memcmp(vv, dlos, (size_t)vl)) {
                    fprintf(stderr, "DECCHK last d=%d v=%d\n  got %s\n  exp %s\n",
                            d, p->vtab[d][k], dlos, vv);
                    exit(9);
                }
            }
#ifdef TP_PROF2
            q0 = PC();
#endif
            int is_pal;
            if (g_weakJ) {
                is_pal = (dl >= g_weakJ);
                for (int i = 0; is_pal && i < g_weakJ; i++) if (dlos[i] != dlos[dl - 1 - i]) is_pal = 0;
            } else is_pal = tp_is_pal_str(dlos, dl);
            if (is_pal) {
#ifdef TP_PROF2
                q1 = PC(); gP_pal += q1 - q0;
#endif
                s->st.leaves++;
                if (!g_suppress_hits) {       /* shard k>0: not this machine's */
                    s->st.hits++;
                    if (onhit) onhit(ud, o->L, &s->lo, dlos, dl);
                }
            } else {
#ifdef TP_PROF2
                q1 = PC(); gP_pal += q1 - q0;
#endif
                s->st.prunes++;
                if (g_depthhist) g_dprune[d]++;
            }
            for (int w = 0; w < nw; w++) { s->lo.w[w] &= ~sm[w]; s->hi.w[w] |= cm[w]; }
            continue;
        }

#ifdef TP_PROF2
        t0 = PC();
#endif
        dec9 dnc;
        int pruned = 0;
        int t5ok_new = t5ok_p; uint64_t t5_new = t5_p;
        int tkok_new = tkok_p; uint64_t tk_new = tk_p;
        uint64_t Rc[TP_MAX_LIMBS];
        int mchild = 0;
        if (fastk) {
            dec9_add(&dnc, dnp, &dlop[k]);
            const int td = ilog10_u32(dnc.c[TP_DNC - dnc.nc]) + 1;
            const int dl = 9 * (dnc.nc - 1) + td;
            const int kd = p->ddig[d];
            int m;
            if (kd <= 0 || dl < kd) {
                m = (kd <= 0) ? dl : 0;
            } else {
                /* carry iff (lo mod 10^kd) >= (10^kd - width): lexicographic compare, chunk-wise */
                const int nfull = kd / 9, rem = kd % 9;
                int cmp = 0;
                if (rem) {
                    const int ci = nfull;
                    const uint32_t lv = ((ci < dnc.nc) ? dnc.c[TP_DNC - 1 - ci] : 0u) % g_p10_32[rem];
                    const uint32_t tv = p->thrpart[d];
                    if (lv != tv) cmp = (lv > tv) ? 1 : -1;
                }
                for (int j = nfull - 1; !cmp && j >= 0; j--) {
                    const uint32_t lv = (j < dnc.nc) ? dnc.c[TP_DNC - 1 - j] : 0u;
                    const uint32_t tv = p->thr[d][j];
                    if (lv != tv) cmp = (lv > tv) ? 1 : -1;
                }
                if (cmp > 0) {
                    int j = dl - kd - 1;
                    while (j >= 0 && dec9_digit(&dnc, td, j) == 9) j--;
                    m = (j < 0) ? 0 : j;
                } else {
                    m = dl - kd;
                }
            }
            if (m > dl) m = dl;
            if (g_weakJ && m > g_weakJ) m = g_weakJ;
            const int m_exact = m;
            if (m > 0) {
                const int kk = (m < r) ? m : r;
                if (m > (int)s->st.maxm) s->st.maxm = (uint32_t)m;
                if (mp > 0 && mp <= m) memcpy(Rc, Rp, sizeof(uint64_t) * (size_t)nw);
                else { memset(Rc, 0, sizeof(uint64_t) * (size_t)nw); }
                rev_add_digits(Rc, nw, &dnc, td, (mp > 0 && mp <= m) ? mp : 0, m);
                mchild = m;
                if (kk > 0) {
                    const int q = kk >> 6, rm = kk & 63;
                    const uint64_t *Rx = Rc; uint64_t Ra[TP_MAX_LIMBS];
                    if (g_off) {                     /* N + c = R (mod 10^m)  <=>  N = R - c */
                        uint64_t br = 0;
                        for (int w = 0; w <= q && w < nw; w++) { const unsigned __int128 t = (unsigned __int128)Rc[w] - g_offc.w[w] - br; Ra[w] = (uint64_t)t; br = (uint64_t)(t >> 64) & 1; }
                        Rx = Ra;
                    }
                    for (int w = 0; w < q; w++) if (Rx[w] != s->lo.w[w]) { pruned = 1; break; }
                    if (!pruned && rm) {
                        const uint64_t msk = (((uint64_t)1 << rm) - 1);
                        pruned = ((Rx[q] & msk) != (s->lo.w[q] & msk));
                    }
                }
            }
            if (g_fastchk) {   /* TP_FASTCHK: recompute m and the 2-adic verdict the long way */
                char vv[TP_DECBUF];
                const int vl = bn_to_dec(&s->lo, nw, vv);
                int mr;
                if (kd <= 0 || vl < kd) mr = (kd <= 0) ? vl : 0;
                else {
                    const char *L = vv + vl - kd, *Rr = p->dstr + (size_t)d * (TP_DECBUF + 1);
                    int i = 0;
                    while (i < kd && L[i] == Rr[i]) i++;
                    if ((i < kd) && (L[i] >= Rr[i])) {
                        int j = vl - kd - 1;
                        while (j >= 0 && vv[j] == '9') j--;
                        mr = (j < 0) ? 0 : j;
                    } else mr = vl - kd;
                }
                if (mr > vl) mr = vl;
                if (g_weakJ && mr > g_weakJ) mr = g_weakJ;
                const int pr = (mr <= 0) ? 0 : prune_cell_lo(s, r, vv, mr);
                if (vl != dl || memcmp(vv, vv, 0) || mr != m_exact || pr != pruned) {
                    fprintf(stderr, "FASTCHK MISMATCH d=%d k=%d dl=%d/%d m=%d/%d pruned=%d/%d\n",
                            d, k, dl, vl, m_exact, mr, pruned, pr);
                    exit(9);
                }
                for (int i = 0; i < dl && i < 40; i++)
                    if (dec9_digit(&dnc, td, i) != vv[i] - '0') {
                        fprintf(stderr, "FASTCHK DIGIT MISMATCH d=%d i=%d\n", d, i); exit(9);
                    }
                g_fastchk_n++;
            }
            if (!pruned && g_dsieve && m_exact >= 2 && dl >= 2) {
                const int lead2 = dec9_digit(&dnc, td, 0) * 10 + dec9_digit(&dnc, td, 1);
                if (!is_lead2_possible(p->L, cur_s_outer, cur_s_sum, d, lead2)) pruned = 1;
            }
            if (!pruned && p->bK && !tkok_p && m_exact >= p->bK) {
                uint64_t t = 0;
                for (int i = p->bK - 1; i >= 0; i--) t = (t * 10u + (uint64_t)dec9_digit(&dnc, td, i)) % p->bN;
                if (g_off) { if (g_offkN != p->bN) { g_offk = off_mod(p->bN); g_offkN = p->bN; } t = (t + p->bN - g_offk) % p->bN; }
                tk_new = t; tkok_new = 1;
            }
            if (!pruned && p->r5J && !t5ok_p && m_exact >= p->r5J) {
                uint64_t t = 0;
                for (int i = p->r5J - 1; i >= 0; i--) t = t * 10u + (uint64_t)dec9_digit(&dnc, td, i);
                t %= p->r5N;
                if (g_off) { if (g_off5N != p->r5N) { g_off5 = off_mod(p->r5N); g_off5N = p->r5N; } t = (t + p->r5N - g_off5) % p->r5N; }
                t5_new = t; t5ok_new = 1;
                uint64_t need = (t >= s->r5) ? t - s->r5 : t + p->r5N - s->r5;
                const uint64_t *rb = p->reach5[d];
                if (rb && !((rb[need >> 6] >> (need & 63)) & 1)) pruned = 1;
            }
        } else {
        int dl;
        if (g_carry9) {
            dec9_add(&dnc, dnp, &dlop[k]);
            dl = dec9_to_str(&dnc, dlos);
            if (g_decchk) {
                char vv[TP_DECBUF];
                const int vl = bn_to_dec(&s->lo, nw, vv);
                if (vl != dl || memcmp(vv, dlos, (size_t)vl)) {
                    fprintf(stderr, "DECCHK carry d=%d k=%d\n  got %s\n  exp %s\n", d, k, dlos, vv);
                    exit(9);
                }
            }
        } else {
            dl = bn_to_dec(&s->lo, nw, dlos);
        }
#ifdef TP_PROF2
        t1 = PC(); gP_lo += t1 - t0; t0 = t1;
#endif
        int m_exact = 0;                     /* exact shared leading digits of lo/hi */
        if (g_noprune) {
            pruned = 0;
        } else if (g_nohi) {
            /* hi is never converted.  hi = lo + width, and the width is a
             * per-level constant, so the shared decimal prefix follows from lo
             * alone: carrying out of the low k digits bumps the digits above,
             * and the run of 9s there is exactly how far that carry reaches. */
            const int k = p->ddig[d];
            int m;
            if (k <= 0 || dl < k) {
                m = (k <= 0) ? dl : 0;
            } else {
                const char *L = dlos + dl - k, *R = p->dstr + (size_t)d * (TP_DECBUF + 1);
                int i = 0;
                while (i < k && L[i] == R[i]) i++;
                /* carry iff (lo mod 10^k) >= (10^k - width); the equality case
                 * also carries.  It cannot actually occur -- lo is odd and the
                 * width is a multiple of 16, so the two cannot be equal -- but
                 * `>=` is what the arithmetic says, so say that. */
                const int carry = (i < k) && (L[i] >= R[i]);
                if (carry) {
                    int j = dl - k - 1;              /* rightmost non-9 above */
                    while (j >= 0 && dlos[j] == '9') j--;
                    m = (j < 0) ? 0 : j;             /* hi = 10^dl: straddle */
                } else {
                    m = dl - k;
                }
            }
            if (m > dl) m = dl;
            if (g_weakJ && m > g_weakJ) m = g_weakJ;      /* TEST ONLY */
            m_exact = m;
            if (g_badprune) m += 2;      /* TEST ONLY: claims 2 more shared digits */
            pruned = (m <= 0) ? 0 : prune_cell_lo(s, r, dlos, m);
            if (g_xm) {
                /* reference: both endpoints converted, prefix compared digitwise */
                char dref[TP_DECBUF];
                const int dh = bn_to_dec(&s->hi, nw, dref);
                int mr = 0;
                while (mr < dl && mr < dh && dlos[mr] == dref[mr]) mr++;
                const int vr = (dl == dh)
                             ? (mr ? prune_cell(s, r, dlos, dref, dl) : 0)
                             : prune_check(s, r, dlos, dl, dref, dh);
                g_xm_n++;
                if (mr != m) {
                    g_xm_m_bad++;
                    if (getenv("TP_XMDBG"))
                        fprintf(stderr, "  [xm] d=%d dl=%d dh=%d k=%d m=%d mref=%d lo=%s hi=%s\n",
                                d, dl, dh, p->ddig[d], m, mr, dlos, dref);
                }
                if (vr != pruned) {
                    g_xm_v_bad++;
                    if (pruned && !vr) g_xm_fwd++;      /* we prune, ref keeps  */
                    else g_xm_back++;
                    if (pruned && !vr)
                        fprintf(stderr, "  [xm] UNSOUND? d=%d m=%d mref=%d lo=%s hi=%s\n",
                                d, m, mr, dlos, dref);
                }
            }
        } else {
            char dhis[TP_DECBUF];            /* reference path only (TP_HIDEC) */
            const int dh = bn_to_dec(&s->hi, nw, dhis);
            pruned = prune_check(s, r, dlos, dl, dhis, dh);
            if (dl == dh && dl >= 2 && dlos[0] == dhis[0] && dlos[1] == dhis[1]) m_exact = 2;
        }
        /* m_exact >= 2 means lo and hi have the same length and agree on their
           first two digits, which is v13's condition for the lead2 test. */
        if (!pruned && g_dsieve && m_exact >= 2 && dl >= 2) {
            const int lead2 = (dlos[0] - '0') * 10 + (dlos[1] - '0');
            if (!is_lead2_possible(p->L, cur_s_outer, cur_s_sum, d, lead2)) pruned = 1;
        }
        if (!pruned && p->r5J && m_exact >= p->r5J) {
            if (g_r5chk) r5_prefixcheck(s, p, d, dlos);
            /* last J digits of N must be the first J digits reversed */
            uint64_t t = 0;
            for (int i = p->r5J - 1; i >= 0; i--) t = t * 10u + (uint64_t)(dlos[i] - '0');
            t %= p->r5N;
            t5_new = t; t5ok_new = 1;
            uint64_t need = (t >= s->r5) ? t - s->r5 : t + p->r5N - s->r5;
            const uint64_t *rb = p->reach5[d];
            if (rb && !((rb[need >> 6] >> (need & 63)) & 1)) pruned = 1;
        }
        }
#ifdef TP_PROF2
        t1 = PC(); gP_prune += t1 - t0; gP_nprune++;
#endif

        if (pruned) {
            s->st.prunes++;
            if (g_depthhist) g_dprune[d]++;
            for (int w = 0; w < nw; w++) { s->lo.w[w] &= ~sm[w]; s->hi.w[w] |= cm[w]; }
            continue;
        }

        s->t5ok = t5ok_new; s->t5 = t5_new;
        s->tkok = tkok_new; s->tk = tk_new;
        if (p->bK) { uint64_t x = rk_p + p->addk[d][k]; if (x >= p->bN) x -= p->bN; s->rk = x; }
        hexdfs(s, p, d + 1, g_carry9 ? &dnc : dnp, Rc, mchild, cur_s_sum, cur_s_outer, onhit, ud);

        for (int w = 0; w < nw; w++) { s->lo.w[w] &= ~sm[w]; s->hi.w[w] |= cm[w]; }
    }
}

static void hexdfs_from(tp_search *s, const hexplan *p, int d, int s_sum, int s_outer,
                        tp_hit_fn onhit, void *ud)
{
    if (g_mod11) {                       /* seed N mod 11 at the entry node */
        bn_t t = s->lo;
        s->r11 = (uint32_t)bn_divmod_small(&t, s->nw, 11);
    }
    if (p->r5J) { bn_t t = s->lo; s->r5 = bn_divmod_small(&t, s->nw, p->r5N); }
    s->t5ok = 0; s->t5 = 0;
    s->tkok = 0; s->tk = 0;
    if (p->bK) { bn_t t = s->lo; s->rk = bn_divmod_small(&t, s->nw, p->bN); }
    dec9 d0;
    if (g_off) {
        if (!(g_fastk && g_carry9 && g_nohi) || g_xm || g_decchk || g_fastchk || g_duchk || g_r5chk || g_bchk || g_noprune || g_badprune) {
            fprintf(stderr, "TP_OFFSET works only on the default fast path (no TP_XM/DECCHK/FASTCHK/DUALCHK/R5CHK/BLOOMCHK/NOFASTK/NOCARRY9/HIDEC)\n"); exit(2);
        }
        /* lo + c over the full width. Only lo's first nw words are defined (a work unit handed to another thread
         * copies the struct, so the words above nw can hold leftovers): treat them as zero explicitly. */
        bn_t lc; unsigned __int128 cy = 0;
        for (int w = 0; w < TP_MAX_LIMBS; w++) {
            cy += (unsigned __int128)(w < s->nw ? s->lo.w[w] : 0) + g_offc.w[w];
            lc.w[w] = (uint64_t)cy; cy >>= 64;
        }
        dec9_from_bn(&lc, TP_MAX_LIMBS, &d0);
    } else dec9_from_bn(&s->lo, s->nw, &d0);
    uint64_t R0[TP_MAX_LIMBS];
    memset(R0, 0, sizeof R0);
    hexdfs(s, p, d, &d0, R0, 0, s_sum, s_outer, onhit, ud);
}

int tp_hexdfs_run(int B, unsigned inv, tp_stats *out,
                  tp_hit_fn onhit, void *ud, int trace)
{
    static int init = 0;
    if (!init) { init = 1;
        g_noprune = getenv("TP_NOPRUNE") ? 1 : 0;
        g_nohi    = getenv("TP_HIDEC") == NULL;
        g_badprune = getenv("TP_BADPRUNE") != NULL;
        g_xm      = getenv("TP_XM") ? 1 : 0;
        g_decchk  = getenv("TP_DECCHK") ? 1 : 0;
        g_mod11   = getenv("TP_NOMOD11") == NULL;
        g_carry9 = getenv("TP_NOCARRY9") == NULL;
        g_fastk = getenv("TP_NOFASTK") == NULL;
        g_fastchk = getenv("TP_FASTCHK") != NULL;
        g_bchk = getenv("TP_BLOOMCHK") != NULL;
        g_duchk = getenv("TP_DUALCHK") != NULL;
    g_noshare = getenv("TP_NOSHARE") != NULL;
        g_cvchk = getenv("TP_CVPCHK") != NULL;
        g_depthhist = getenv("TP_DEPTH") ? 1 : 0;
    }
    tp_orbits o;
    if (tp_orbits_build(&o, B, inv) != 0) return -1;
    const int nw = tp_words(o.nsites);

    tp_search s;
    tp_search_init(&s, &o, nw, 0, NULL, 0);
    s.budget = 0;

    g_noprune  = getenv("TP_NOPRUNE") != NULL;
    g_nohi     = getenv("TP_HIDEC") == NULL;
    g_badprune = getenv("TP_BADPRUNE") != NULL;
    g_xm       = getenv("TP_XM") != NULL;
    g_decchk   = getenv("TP_DECCHK") != NULL;
    g_mod11    = getenv("TP_NOMOD11") == NULL;
    g_carry9 = getenv("TP_NOCARRY9") == NULL;
    g_fastk = getenv("TP_NOFASTK") == NULL;
    g_fastchk = getenv("TP_FASTCHK") != NULL;
    g_bchk = getenv("TP_BLOOMCHK") != NULL;
    g_duchk = getenv("TP_DUALCHK") != NULL;
    g_noshare = getenv("TP_NOSHARE") != NULL;
    g_cvchk = getenv("TP_CVPCHK") != NULL;
    g_depthhist = getenv("TP_DEPTH") != NULL;
    memset(g_dnode, 0, sizeof g_dnode);
    g_m0 = 0; g_calls = 0;
    g_xm_n = g_xm_m_bad = g_xm_v_bad = 0;
    memset(g_dprune, 0, sizeof g_dprune);

    hexplan p;
    memset(&p, 0, sizeof p);
    hexplan_init(&p, &o);
    hexplan_prepare(&p, &o);
    p.trace = trace;

    hexdfs_from(&s, &p, 0, 0, 0, onhit, ud);

    if (g_xm)
        fprintf(stderr, "  [xm] %llu calls, m bad %llu, verdict bad %llu "
                        "(we-prune/ref-keeps %llu, ref-prunes/we-keep %llu)\n",
                (unsigned long long)g_xm_n, (unsigned long long)g_xm_m_bad,
                (unsigned long long)g_xm_v_bad,
                (unsigned long long)g_xm_fwd, (unsigned long long)g_xm_back);
    if (g_r5chk)
        fprintf(stderr, "  [r5chk] %llu sampled completions, %llu NOT in reach set; %llu prefix samples, %llu prefix MISMATCH\n",
                (unsigned long long)g_r5chk_n, (unsigned long long)g_r5chk_bad,
                (unsigned long long)g_r5chk_pn, (unsigned long long)g_r5chk_pbad);
    if (g_depthhist) {
        fprintf(stderr, "  [m==0] %llu of %llu prune_cell calls = %.1f%%  (Idea-1 hit rate)\n",
                (unsigned long long)g_m0, (unsigned long long)g_calls,
                g_calls ? 100.0 * (double)g_m0 / (double)g_calls : 0.0);
        fprintf(stderr, "  depth |      nodes |    prunes | prune%%\n");
        for (int d = 0; d < p.depth_max; d++)
            if (g_dnode[d])
                fprintf(stderr, "  %5d | %10llu | %9llu | %5.1f%%\n", d,
                        (unsigned long long)g_dnode[d], (unsigned long long)g_dprune[d],
                        100.0 * (double)g_dprune[d] / (double)g_dnode[d]);
    }

    if (out) *out = s.st;
    hexplan_free(&p);
    tp_orbits_free(&o);
    return 0;
}

/* ------------------------------------------------------------------ */
/* DFS                                                                */
/* ------------------------------------------------------------------ */

static void dfs(tp_search *s, int d,
                const char *dlo, int dl, const char *dhi, int dh,
                tp_hit_fn onhit, void *ud)
{
    const tp_orbits *o = s->o;

    if (s->budget && s->st.nodes >= s->budget) { s->st.aborted = 1; return; }

    if (d == s->stop_depth) {                  /* work-unit boundary       */
        if (s->nunits < s->unit_cap) {
            s->units[s->nunits].B     = o->B;
            s->units[s->nunits].depth = d;
            s->units[s->nunits].path  = s->path;
            s->nunits++;
        }
        return;
    }

    if (d == o->nfree) {                       /* fully assigned           */
        s->st.leaves++;
        if (tp_is_pal_str(dlo, dl)) {
            s->st.hits++;
            if (onhit) onhit(ud, o->B, &s->lo, dlo, dl);
        }
        return;
    }

    const int       r  = o->suffixmin[d + 1];
    const int       c  = o->order[d];
    const uint64_t *cm = o->mask + (size_t)c * o->nw64;
    const int       nw = s->nw;
    char            buf[TP_DECBUF];

    for (int v = 1; v >= 0; v--) {
        if (v) { for (int i = 0; i < nw; i++) s->lo.w[i] |=  cm[i]; }
        else   { for (int i = 0; i < nw; i++) s->hi.w[i] &= ~cm[i]; }
        s->st.nodes++;

        int bl, pruned;
        if (v) {
            bl = bn_to_dec(&s->lo, nw, buf);
            s->st.dec_conv++;
            pruned = prune_check(s, r, buf, bl, dhi, dh);
        } else {
            bl = bn_to_dec(&s->hi, nw, buf);
            s->st.dec_conv++;
            pruned = prune_check(s, r, buf, bl, dlo, dl);
        }

        if (s->trace_depth > 0 && d <= s->trace_depth) {
            int mm = 0;
            if (v) { while (mm < bl && mm < dh && buf[mm] == dhi[mm]) mm++; }
            else    { while (mm < bl && mm < dl && buf[mm] == dlo[mm]) mm++; }
            fprintf(stderr,
                    "  d=%-3d v=%d r=%-4d m=%-3d %s\n"
                    "        lo=%s\n        hi=%s\n",
                    d, v, r, mm, pruned ? "PRUNE" : "keep",
                    v ? buf : dlo, v ? dhi : buf);
        }

        if (pruned) {
            s->st.prunes++;
        } else {
            if (d < 64) {
                if (v) s->path |=  (uint64_t)1 << d;
                else   s->path &= ~((uint64_t)1 << d);
            }
            if (v) dfs(s, d + 1, buf, bl, dhi, dh, onhit, ud);
            else   dfs(s, d + 1, dlo, dl, buf, bl, onhit, ud);
        }

        if (v) { for (int i = 0; i < nw; i++) s->lo.w[i] &= ~cm[i]; }
        else   { for (int i = 0; i < nw; i++) s->hi.w[i] |=  cm[i]; }
    }
}

/* ------------------------------------------------------------------ */
/* public entry points                                                */
/* ------------------------------------------------------------------ */

int tp_split_depth(int nfree, int target_units)
{
    int d = 0;
    if (target_units < 1) target_units = 1;
    while (d < nfree && d < 62 && ((long)1 << d) < (long)target_units) d++;
    return d;
}

void tp_search_enum_units(tp_search *s)
{
    char dlo[TP_DECBUF], dhi[TP_DECBUF];
    int  dl, dh;
    setup_state(s);
    s->path   = 0;
    s->nunits = 0;
    tp_search_roots(s, dlo, &dl, dhi, &dh);
    dfs(s, 0, dlo, dl, dhi, dh, NULL, NULL);
}

void tp_search_run_unit(tp_search *s, const tp_unit *u,
                        tp_hit_fn onhit, void *ud)
{
    const tp_orbits *o = s->o;
    char dlo[TP_DECBUF], dhi[TP_DECBUF];
    int  dl, dh;

    setup_state(s);

    for (int d = 0; d < u->depth; d++) {
        const int       c  = o->order[d];
        const uint64_t *cm = o->mask + (size_t)c * o->nw64;
        if ((u->path >> d) & 1u) { for (int i = 0; i < s->nw; i++) s->lo.w[i] |=  cm[i]; }
        else                     { for (int i = 0; i < s->nw; i++) s->hi.w[i] &= ~cm[i]; }
    }
    tp_search_roots(s, dlo, &dl, dhi, &dh);
    dfs(s, u->depth, dlo, dl, dhi, dh, onhit, ud);
}

/* ------------------------------------------------------------------ */
/* parallel front end: expand to a cut depth, then farm out subtrees   */
/* ------------------------------------------------------------------ */

typedef struct {
    const hexplan   *p;
    const tp_orbits *o;
    int              nw;
    hxtask          *tasks;
    _Atomic int      next;
    int              ntask;
    tp_hit_fn        onhit;
    void            *ud;
    tp_stats        *st;          /* one slot per thread */
} hxpool;

typedef struct { hxpool *pool; int idx; } hxarg;

static void *hx_worker(void *argp)
{
    hxarg *a = argp;
    hxpool *pool = a->pool;
    _Alignas(64) tp_search s;
    tp_search_init(&s, pool->o, pool->nw, 0, NULL, 0);
    s.budget = 0;
    memset(&s.st, 0, sizeof s.st);

    /* rebuild the plan in this frame rather than copying it from an opaque
       pointer, so the compiler can see its fields and keep the recursion
       specialised (a memcpy of the same bytes is ~2x slower per node) */
    hexplan lp;
    memset(&lp, 0, sizeof lp);
    hexplan_init(&lp, pool->o);
    hexplan_prepare(&lp, pool->o);
    lp.tasks = NULL; lp.ntask = NULL; lp.cut = 0;

    for (;;) {
        const int i = atomic_fetch_add(&pool->next, 1);
        if (i >= pool->ntask) break;
        s.lo = pool->tasks[i].lo;
        s.hi = pool->tasks[i].hi;
        const uint64_t before = s.st.nodes;
        hexdfs_from(&s, &lp, pool->tasks[i].depth, pool->tasks[i].s_sum, pool->tasks[i].s_outer,
                    pool->onhit, pool->ud);
        pool->tasks[i].nodes = s.st.nodes - before;
    }
    hexplan_free(&lp);
    pool->st[a->idx] = s.st;
    return NULL;
}

/* run one bit length in parallel; returns 0 on success */
/* ---- shard splitting --------------------------------------------------
 * Phase 1 expands the tree to the cut single-threaded, so the task list is a
 * deterministic function of (B, inv, cut).  Two machines with the same binary
 * therefore build byte-identical task lists, and task i belongs to shard
 * i mod n.  Every shard reports its own task-only totals, so the full run is
 * the serial part plus the sum of the shards -- exact, not statistical.
 */
static int g_shard_k = 0, g_shard_n = 1;
int tp_hexdfs_sharded(void) { return g_shard_n > 1; }

void tp_hexdfs_set_shard(int k, int n)
{
    if (n < 1) n = 1;
    if (k < 0) k = 0;
    g_shard_n = n;
    g_shard_k = n ? (k % n) : 0;
}

static int hexdfs_mt_once(int B, unsigned inv, int nthreads, int cut,
                          tp_stats *out, tp_hit_fn onhit, void *ud)
{
    if (nthreads < 1) nthreads = 1;

    tp_orbits o;
    if (tp_orbits_build(&o, B, inv) != 0) return -1;
    const int nw = tp_words(o.nsites);

    hexplan p;
    memset(&p, 0, sizeof p);
    hexplan_init(&p, &o);
    hexplan_prepare(&p, &o);
    g_noprune = getenv("TP_NOPRUNE") != NULL;
    if (getenv("TP_BUILD_ONLY")) {
        /* Build (and, with TP_CVP_CACHE, store) every level's tables and lattices, then stop without searching.
         * Used to warm a shared basis cache, e.g. with strong reduction (TP_CVP_BKZ), before a sharded run. */
        if (out) memset(out, 0, sizeof *out);
        hexplan_free(&p); tp_orbits_free(&o);
        return 0;
    }

    tp_search s;
    tp_search_init(&s, &o, nw, 0, NULL, 0);
    s.budget = 0;

    if (cut > p.depth_max) cut = p.depth_max;
    const int cap = getenv("TP_CAP") ? atoi(getenv("TP_CAP")) : (1 << 16);
    hxtask *tasks = malloc((size_t)cap * sizeof *tasks);
    if (!tasks) { tp_orbits_free(&o); return -1; }
    int ntask = 0;

    /* phase 1: single-threaded expansion down to the cut, using the one
       and only traversal routine */
    p.cut = cut; p.tasks = tasks; p.ntask = &ntask; p.cap = cap;
    /* Every shard walks the same prefix down to the cut, so a hit found there
       would be printed by all of them; only shard 0 owns that report.        */
    g_suppress_hits = (g_shard_n > 1 && g_shard_k != 0);
    hexdfs_from(&s, &p, 0, 0, 0, onhit, ud);
    g_suppress_hits = 0;
    tp_stats total = s.st;
    const tp_stats serial = s.st;   /* phase-1 contribution, one per machine */
    if (getenv("TP_SPLIT"))
        fprintf(stderr, "  [split] cut=%d tasks=%d serial_nodes=%llu (%.1f%% of work)\n",
                cut, ntask, (unsigned long long)s.st.nodes,
                100.0 * (double)s.st.nodes / 99511142.0);

    /* If the whole tree fitted in the serial phase there is nothing to split;
       it belongs to shard 0 alone, otherwise every machine would report the
       same hits and every machine's node count would include the whole tree. */
    if (g_shard_n > 1 && ntask == 0 && g_shard_k != 0) {
        memset(&total, 0, sizeof total);
        if (out) *out = total;
        printf("SHARD %d/%d tasks=0 task_nodes=0 | serial_nodes=0 serial_prunes=0\n",
               g_shard_k, g_shard_n);
        fflush(stdout);
        free(tasks); tp_orbits_free(&o);
        return 0;
    }

    /* keep only this machine's slice of the task list */
    if (g_shard_n > 1) {
        const int nall = ntask;
        int w = 0;
        for (int i = 0; i < nall; i++)
            if (i % g_shard_n == g_shard_k) tasks[w++] = tasks[i];
        ntask = w;
        if (getenv("TP_SPLIT"))
            fprintf(stderr, "  [shard] %d/%d keeps %d of %d tasks\n",
                    g_shard_k, g_shard_n, ntask, nall);
    }

    /* phase 2: workers resume below the cut (no further task recording) */
    hexplan pw = p;
    pw.tasks = NULL; pw.ntask = NULL; pw.cut = 0;

    if (nthreads > ntask) nthreads = ntask > 0 ? ntask : 1;

    tp_stats *acc = calloc((size_t)nthreads, sizeof *acc);
    pthread_t *th  = malloc((size_t)nthreads * sizeof *th);
    hxarg     *ar  = malloc((size_t)nthreads * sizeof *ar);
    hxpool pool;
    memset(&pool, 0, sizeof pool);
    pool.p = &pw; pool.o = &o; pool.nw = nw;
    pool.tasks = tasks; pool.ntask = ntask;
    pool.onhit = onhit; pool.ud = ud; pool.st = acc;
    atomic_init(&pool.next, 0);

    for (int t = 0; t < nthreads; t++) {
        ar[t].pool = &pool;
        ar[t].idx  = t;
        th[t] = 0;
        if (t < nthreads - 1) pthread_create(&th[t], NULL, hx_worker, &ar[t]);
    }
    hx_worker(&ar[nthreads - 1]);          /* the caller is a worker too */

    for (int t = 0; t < nthreads - 1; t++)
        if (th[t]) pthread_join(th[t], NULL);

    for (int t = 0; t < nthreads; t++) {
        total.nodes  += acc[t].nodes;
        total.prunes += acc[t].prunes;
        total.leaves += acc[t].leaves;
        total.hits   += acc[t].hits;
        if (acc[t].maxm > total.maxm) total.maxm = acc[t].maxm;
    }

    /* The serial prefix down to the cut is walked by every machine, so
       crediting it to each of them makes the per-shard totals sum to more than
       the unsplit run.  Credit it to shard 0 only: then the shard totals add up
       to exactly the unsplit node count and the shard hits are the hit list. */
    if (g_shard_n > 1 && g_shard_k != 0) {
        total.nodes  -= serial.nodes;
        total.prunes -= serial.prunes;
    }
    if (out) *out = total;

    if (g_shard_n > 1) {
        uint64_t tn = 0;
        for (int i = 0; i < ntask; i++) tn += tasks[i].nodes;
        printf("SHARD %d/%d tasks=%d task_nodes=%llu | serial_nodes=%llu "
               "serial_prunes=%llu\n",
               g_shard_k, g_shard_n, ntask, (unsigned long long)tn,
               (unsigned long long)serial.nodes,
               (unsigned long long)serial.prunes);
        fflush(stdout);
    }

    if (getenv("TP_TASKS")) {
        /* exact task sizes are deterministic, so this predicts the makespan on
           any thread count without needing that many cores */
        uint64_t tot = 0, mx = 0;
        for (int i = 0; i < ntask; i++) { tot += tasks[i].nodes; if (tasks[i].nodes > mx) mx = tasks[i].nodes; }
        static int Nth[3] = {10, 4, 8};
        fprintf(stderr, "  [tasks] cut=%d n=%d total=%llu max=%llu (%.1f%% of tree)\n",
                cut, ntask, (unsigned long long)tot, (unsigned long long)mx,
                tot ? 100.0 * (double)mx / (double)tot : 0.0);
        for (int k = 0; k < 3; k++) {
            int nw = Nth[k];
            uint64_t load[64] = {0};
            for (int i = 0; i < ntask; i++) {          /* LPT: largest first */
                int best = 0;
                for (int j = 1; j < nw && j < 64; j++) if (load[j] < load[best]) best = j;
                load[best] += tasks[i].nodes;
            }
            uint64_t mk = 0;
            for (int j = 0; j < nw && j < 64; j++) if (load[j] > mk) mk = load[j];
            fprintf(stderr, "  [tasks]   ideal %2d-way speedup %.2fx (LPT makespan-bound %.2fx)\n",
                    nw, (double)nw, mk ? (double)tot / (double)mk : 0.0);
        }
    }

    free(acc); free(th); free(ar); free(tasks);
    tp_orbits_free(&o);
    return 0;
}


/* Public entry: if the requested cut generates more subtrees than the task
 * array can hold, back the cut off rather than truncate.  Truncation would
 * silently skip subtrees and report a false "no solutions". */
int tp_hexdfs_run_mt(int B, unsigned inv, int nthreads, int cut,
                     tp_stats *out, tp_hit_fn onhit, void *ud)
{
    for (int c = cut; ; c -= 2) {
        if (c < 1) c = 1;
        tp_stats st;
        memset(&st, 0, sizeof st);
        if (hexdfs_mt_once(B, inv, nthreads, c, &st, onhit, ud) != 0) return -1;
        if (!st.overflow) {
            if (out) {
                *out = st;
                out->cut_used = c;
            }
            return 0;
        }
        if (c == 1) { if (out) *out = st; return 0; }   /* cannot go shallower */
    }
}

void tp_prof_dump(void)
{
#ifdef TP_PROF2
    const double M = 1e6;
    fprintf(stderr, "  [prof] mask   %8.1f Mcyc\n", (double)gP_mask / M);
    fprintf(stderr, "  [prof] lo dec %8.1f Mcyc\n", (double)gP_lo / M);
    fprintf(stderr, "  [prof] hi dec %8.1f Mcyc\n", (double)gP_hi / M);
    fprintf(stderr, "  [prof] prune  %8.1f Mcyc  (%llu calls)\n", (double)gP_prune / M,
            (unsigned long long)gP_nprune);
    fprintf(stderr, "  [prof] lastcv %8.1f Mcyc  (%llu nodes)\n", (double)gP_lastconv / M,
            (unsigned long long)gP_nlast);
    fprintf(stderr, "  [prof] palstr %8.1f Mcyc\n", (double)gP_pal / M);
#endif
}
