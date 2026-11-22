/* Crash-safe checkpointing.
 *
 *   <workdir>/done.log   append-only, fsync'd after every unit
 *   <workdir>/hits.log   append-only, fsync'd immediately on every hit
 *
 * A work unit is identified by (B, split depth, path).  On restart the log is
 * replayed into a hash set and finished units are skipped, so a kill -9 costs
 * at most the units that were in flight.
 */
#include "tripal.h"
#include <pthread.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

typedef struct { int B; int depth; uint64_t path; } rec_t;

struct tp_ckpt {
    int        enabled;
    char       dir[1024];
    char       done_path[1200];
    char       hits_path[1200];
    FILE      *fdone;
    FILE      *fhits;
    pthread_mutex_t mu;

    rec_t     *tab;          /* open addressing hash set */
    size_t     tabcap;       /* power of two             */
    size_t     nrec;
    uint64_t   loaded;
};

static size_t hash_rec(int B, int depth, uint64_t path)
{
    uint64_t h = (uint64_t)B * 0x9E3779B97F4A7C15ULL;
    h ^= (uint64_t)(unsigned)depth * 0xC2B2AE3D27D4EB4FULL;
    h ^= path * 0x165667B19E3779F9ULL;
    h ^= h >> 29;
    h *= 0xBF58476D1CE4E5B9ULL;
    h ^= h >> 32;
    return (size_t)h;
}

static void tab_grow(tp_ckpt *c, size_t need)
{
    size_t cap = 64;
    while (cap < need * 2) cap <<= 1;
    rec_t *nt = (rec_t *)calloc(cap, sizeof(rec_t));
    if (!nt) return;
    for (size_t i = 0; i < cap; i++) nt[i].B = -1;
    if (c->tab) {
        for (size_t i = 0; i < c->tabcap; i++) {
            if (c->tab[i].B < 0) continue;
            size_t j = hash_rec(c->tab[i].B, c->tab[i].depth, c->tab[i].path) & (cap - 1);
            while (nt[j].B >= 0) j = (j + 1) & (cap - 1);
            nt[j] = c->tab[i];
        }
        free(c->tab);
    }
    c->tab    = nt;
    c->tabcap = cap;
}

static void tab_put(tp_ckpt *c, int B, int depth, uint64_t path)
{
    if (!c->tab || c->nrec * 2 >= c->tabcap) tab_grow(c, c->nrec + 1);
    if (!c->tab) return;
    size_t j = hash_rec(B, depth, path) & (c->tabcap - 1);
    while (c->tab[j].B >= 0) {
        if (c->tab[j].B == B && c->tab[j].depth == depth && c->tab[j].path == path)
            return;
        j = (j + 1) & (c->tabcap - 1);
    }
    c->tab[j].B     = B;
    c->tab[j].depth = depth;
    c->tab[j].path  = path;
    c->nrec++;
}

static int tab_has(tp_ckpt *c, int B, int depth, uint64_t path)
{
    if (!c->tab) return 0;
    size_t j = hash_rec(B, depth, path) & (c->tabcap - 1);
    while (c->tab[j].B >= 0) {
        if (c->tab[j].B == B && c->tab[j].depth == depth && c->tab[j].path == path)
            return 1;
        j = (j + 1) & (c->tabcap - 1);
    }
    return 0;
}

/* mkdir -p */
static int mkpath(const char *path)
{
    char buf[1200];
    snprintf(buf, sizeof buf, "%s", path);
    for (char *p = buf + 1; *p; p++) {
        if (*p != '/') continue;
        *p = 0;
        if (mkdir(buf, 0777) != 0 && errno != EEXIST) { *p = '/'; return -1; }
        *p = '/';
    }
    if (mkdir(buf, 0777) != 0 && errno != EEXIST) return -1;
    return 0;
}

tp_ckpt *tp_ckpt_open(const char *dir, int enable)
{
    tp_ckpt *c = (tp_ckpt *)calloc(1, sizeof *c);
    if (!c) return NULL;
    c->enabled = enable;
    pthread_mutex_init(&c->mu, NULL);
    snprintf(c->dir, sizeof c->dir, "%s", dir ? dir : ".");

    if (!enable) return c;

    if (mkpath(c->dir) != 0) {
        fprintf(stderr, "tripal: cannot create checkpoint directory '%s': %s\n",
                c->dir, strerror(errno));
        fprintf(stderr, "tripal: refusing to start a long run without durable "
                        "logs (use --no-checkpoint to override)\n");
        exit(1);
    }
    snprintf(c->done_path, sizeof c->done_path, "%s/done.log", c->dir);
    snprintf(c->hits_path, sizeof c->hits_path, "%s/hits.log", c->dir);

    /* replay done.log */
    FILE *f = fopen(c->done_path, "r");
    if (f) {
        char line[256];
        while (fgets(line, sizeof line, f)) {
            int B, depth;
            unsigned long long path;
            if (sscanf(line, "%d %d %llx", &B, &depth, &path) == 3) {
                tab_put(c, B, depth, (uint64_t)path);
                c->loaded++;
            }
        }
        fclose(f);
    }

    c->fdone = fopen(c->done_path, "a");
    c->fhits = fopen(c->hits_path, "a");
    if (!c->fdone || !c->fhits) {
        fprintf(stderr, "tripal: cannot open checkpoint logs in '%s': %s\n",
                c->dir, strerror(errno));
        exit(1);
    }
    return c;
}

void tp_ckpt_close(tp_ckpt *c)
{
    if (!c) return;
    if (c->fdone) { fflush(c->fdone); fsync(fileno(c->fdone)); fclose(c->fdone); }
    if (c->fhits) { fflush(c->fhits); fsync(fileno(c->fhits)); fclose(c->fhits); }
    free(c->tab);
    pthread_mutex_destroy(&c->mu);
    free(c);
}

int tp_ckpt_is_done(tp_ckpt *c, int B, int depth, uint64_t path)
{
    if (!c || !c->enabled) return 0;
    return tab_has(c, B, depth, path);
}

void tp_ckpt_mark(tp_ckpt *c, int B, int depth, uint64_t path)
{
    if (!c || !c->enabled) return;
    pthread_mutex_lock(&c->mu);
    if (c->fdone) {
        fprintf(c->fdone, "%d %d %016llx\n", B, depth, (unsigned long long)path);
        fflush(c->fdone);
        fsync(fileno(c->fdone));
    }
    tab_put(c, B, depth, path);
    pthread_mutex_unlock(&c->mu);
}

void tp_ckpt_hit(tp_ckpt *c, int B, const char *dec, const char *hex)
{
    if (!c || !c->enabled) return;
    pthread_mutex_lock(&c->mu);
    if (c->fhits) {
        fprintf(c->fhits, "%d %s 0x%s\n", B, dec, hex);
        fflush(c->fhits);
        fsync(fileno(c->fhits));
    }
    pthread_mutex_unlock(&c->mu);
}

uint64_t tp_ckpt_loaded(tp_ckpt *c) { return c ? c->loaded : 0; }

const char *tp_ckpt_hits_path(tp_ckpt *c)
{
    return c ? c->hits_path : "";
}
