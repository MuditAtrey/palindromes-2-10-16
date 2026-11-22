/* Generic union-find orbit construction.
 *
 * A "site" is a bit position 0 .. nsites-1 of the pattern we are building.
 * The constraints placed on the pattern are involutions on sites:
 *
 *   binary mirror   B(i) = B - 1 - i                    (i < B)
 *   hex mirror      H(4j+k) = 4(L-1-j) + k              (0 <= j < L, 0<=k<4)
 *
 * Every site of one orbit must carry the same bit value, so an orbit is a
 * single boolean degree of freedom ("free orbit") unless something forces it.
 *
 * Forced values:
 *   site B-1          := 1   (n has exactly B bits)
 *   site B .. 4L-1    := 0   (padding of the top hex digit)
 *
 * If an orbit would receive two different forced values the whole bit length
 * is inconsistent and is rejected before any search happens.
 */
#include "tripal.h"
#include <math.h>

static int uf_find(int *p, int x)
{
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

static void uf_union(int *p, int a, int b)
{
    a = uf_find(p, a);
    b = uf_find(p, b);
    if (a != b) p[a] = b;
}

static int setval(tp_orbits *o, int c, int8_t v)
{
    if (o->cval[c] < 0) { o->cval[c] = v; return 0; }
    return o->cval[c] == v ? 0 : -1;
}

/* Sort free orbits.
 *   primary key: highest site, descending  -> fixes the top bits (and, via the
 *                mirror, the bottom bits) first, which is exactly what
 *                stabilises the leading decimal digits.
 *   tie break:   lowest site, ascending    -> makes the "known low bits"
 *                frontier r grow as fast as possible.
 * Insertion sort: nfree <= a few hundred, no qsort_r portability headaches.
 */
static void sort_free(tp_orbits *o)
{
    for (int i = 1; i < o->nfree; i++) {
        int v = o->order[i];
        int j = i - 1;
        while (j >= 0) {
            int u = o->order[j];
            int better;
            if (o->cmax[u] != o->cmax[v])      better = o->cmax[v] > o->cmax[u];
            else if (o->cmin[u] != o->cmin[v]) better = o->cmin[v] < o->cmin[u];
            else                               better = v < u;
            if (!better) break;
            o->order[j + 1] = u;
            j--;
        }
        o->order[j + 1] = v;
    }
}

int tp_orbits_build(tp_orbits *o, int B, unsigned inv)
{
    memset(o, 0, sizeof *o);
    o->B     = B;
    o->L     = (B + 3) / 4;
    o->inv   = inv;
    o->valid = 0;

    const int hex = (inv & TP_INV_HEX) != 0;
    const int bin = (inv & TP_INV_BIN) != 0;
    const int nsites = hex ? 4 * o->L : B;

    if (B < 1 || B > TP_MAX_BITS) {
        snprintf(o->reason, sizeof o->reason, "bit length %d out of range", B);
        return -1;
    }
    o->nsites = nsites;
    o->nw64   = tp_words(nsites);

    int *par = (int *)malloc((size_t)nsites * sizeof(int));
    if (!par) { snprintf(o->reason, sizeof o->reason, "out of memory"); return -1; }
    for (int i = 0; i < nsites; i++) par[i] = i;

    if (hex) {
        for (int i = 0; i < nsites; i++)
            uf_union(par, i, 4 * (o->L - 1 - i / 4) + (i % 4));
    }
    if (bin) {
        for (int i = 0; i < B; i++)
            uf_union(par, i, B - 1 - i);
    }

    o->comp = (int *)malloc((size_t)nsites * sizeof(int));
    if (!o->comp) goto oom;
    int *rootid = (int *)malloc((size_t)nsites * sizeof(int));
    if (!rootid) goto oom;
    for (int i = 0; i < nsites; i++) rootid[i] = -1;

    int ncomp = 0;
    for (int i = 0; i < nsites; i++) {
        int r = uf_find(par, i);
        if (rootid[r] < 0) rootid[r] = ncomp++;
        o->comp[i] = rootid[r];
    }
    free(rootid);
    free(par);
    o->ncomp = ncomp;

    o->cval = (int8_t *)malloc((size_t)ncomp);
    o->cmin = (int *)malloc((size_t)ncomp * sizeof(int));
    o->cmax = (int *)malloc((size_t)ncomp * sizeof(int));
    o->mask = (uint64_t *)calloc((size_t)ncomp * o->nw64, sizeof(uint64_t));
    if (!o->cval || !o->cmin || !o->cmax || !o->mask) goto oom;

    for (int c = 0; c < ncomp; c++) {
        o->cval[c] = -1;
        o->cmin[c] = INT_MAX;
        o->cmax[c] = -1;
    }
    for (int i = 0; i < nsites; i++) {
        int c = o->comp[i];
        if (i < o->cmin[c]) o->cmin[c] = i;
        if (i > o->cmax[c]) o->cmax[c] = i;
        o->mask[(size_t)c * o->nw64 + (i >> 6)] |= (uint64_t)1 << (i & 63);
    }

    /* forced values */
    if (setval(o, o->comp[B - 1], 1) < 0) {
        snprintf(o->reason, sizeof o->reason,
                 "top bit shares an orbit with a forced-0 site");
        return -1;
    }
    for (int i = B; i < nsites; i++) {
        if (setval(o, o->comp[i], 0) < 0) {
            snprintf(o->reason, sizeof o->reason,
                     "top bit shares an orbit with a forced-0 site");
            return -1;
        }
    }

    o->order = (int *)malloc((size_t)ncomp * sizeof(int));
    if (!o->order) goto oom;
    int nfree = 0;
    for (int c = 0; c < ncomp; c++)
        if (o->cval[c] < 0) o->order[nfree++] = c;
    o->nfree = nfree;

    sort_free(o);

    o->suffixmin = (int *)malloc(((size_t)nfree + 1) * sizeof(int));
    if (!o->suffixmin) goto oom;
    o->suffixmin[nfree] = INT_MAX;
    for (int d = nfree - 1; d >= 0; d--) {
        int m = o->cmin[o->order[d]];
        o->suffixmin[d] = m < o->suffixmin[d + 1] ? m : o->suffixmin[d + 1];
    }

    o->valid = 1;
    return 0;

oom:
    snprintf(o->reason, sizeof o->reason, "out of memory");
    return -1;
}

void tp_orbits_free(tp_orbits *o)
{
    free(o->comp);      o->comp = NULL;
    free(o->cval);      o->cval = NULL;
    free(o->cmin);      o->cmin = NULL;
    free(o->cmax);      o->cmax = NULL;
    free(o->order);     o->order = NULL;
    free(o->suffixmin); o->suffixmin = NULL;
    free(o->mask);      o->mask = NULL;
}

double tp_orbits_space(const tp_orbits *o)
{
    return ldexp(1.0, o->nfree);
}
