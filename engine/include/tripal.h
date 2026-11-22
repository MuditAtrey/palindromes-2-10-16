/* tripal - triple palindrome (base 2 / 10 / 16) search engine
 *
 * Two-ended DFS over union-find orbits with decimal-interval pruning.
 */
#ifndef TRIPAL_H
#define TRIPAL_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* ------------------------------------------------------------------ */
/* configuration                                                       */
/* ------------------------------------------------------------------ */

#ifndef TP_MAX_LIMBS
#define TP_MAX_LIMBS 32                    /* 32 * 64 = 2048 bits          */
#endif
#define TP_MAX_BITS (TP_MAX_LIMBS * 64)    /* 512 hex digits               */
#define TP_DECBUF   (TP_MAX_BITS / 3 + 32) /* > max decimal digits         */
#define TP_HEXBUF   (TP_MAX_BITS / 4 + 8)

/* which non-decimal reflections are imposed on the bit pattern */
#define TP_INV_BIN 1u   /* n is a base-2 palindrome  (bit reversal)      */
#define TP_INV_HEX 2u   /* n is a base-16 palindrome (nibble reversal)   */

/* ------------------------------------------------------------------ */
/* fixed-width big integer                                             */
/* ------------------------------------------------------------------ */

typedef struct { uint64_t w[TP_MAX_LIMBS]; } bn_t;

static inline int tp_words(int nbits) { return (nbits + 63) >> 6; }

void     bn_zero(bn_t *a, int nw);
void     bn_setbit(bn_t *a, int i);
void     bn_clrbit(bn_t *a, int i);
int      bn_getbit(const bn_t *a, int i);
int      bn_is_zero(const bn_t *a, int nw);
void     bn_muladd(bn_t *a, int nw, uint64_t mul, uint64_t add);
uint64_t bn_divmod_small(bn_t *a, int nw, uint64_t d);
void     bn_add_u64(bn_t *a, int nw, uint64_t v);
int      bn_cmp(const bn_t *a, const bn_t *b, int nw);

int  bn_to_dec(const bn_t *a, int nw, char *out);   /* -> length          */
int  bn_to_hex(const bn_t *a, int nw, char *out);   /* -> length, no 0x   */
int  bn_from_dec(const char *s, int nw, bn_t *out);

int  tp_is_pal_str(const char *s, int len);
int  tp_is_pal_base(const bn_t *a, int nw, int base);

/* ------------------------------------------------------------------ */
/* orbit engine                                                        */
/* ------------------------------------------------------------------ */

typedef struct {
    int      B;          /* bit length of n                             */
    int      L;          /* hex digits = ceil(B/4)                      */
    int      nsites;     /* 4L with hex mirror, else B                  */
    int      ncomp;
    int      nfree;
    int      nw64;       /* words spanning nsites                       */
    unsigned inv;        /* TP_INV_* mask actually applied              */
    int      valid;
    char     reason[96];

    int      *comp;      /* comp[site]                                  */
    int8_t   *cval;      /* per comp: -1 free, 0 / 1 forced             */
    int      *cmin;      /* per comp: lowest site                       */
    int      *cmax;      /* per comp: highest site                      */
    int      *order;     /* free comps, sorted by cmax descending       */
    int      *suffixmin; /* suffixmin[d] = min cmin over order[d..]     */
    uint64_t *mask;      /* mask[comp*nw64 .. ] = site bitmap           */
} tp_orbits;

int  tp_orbits_build(tp_orbits *o, int B, unsigned inv);
void tp_orbits_free(tp_orbits *o);
/* number of complete assignments = 2^nfree (capped for display) */
double tp_orbits_space(const tp_orbits *o);

/* ------------------------------------------------------------------ */
/* search                                                              */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t nodes;    /* DFS nodes entered                             */
    uint64_t prunes;   /* branches killed by decimal-interval pruning   */
    uint64_t leaves;   /* fully assigned patterns                       */
    uint64_t hits;     /* decimal palindromes found at leaves           */
    uint64_t dec_conv; /* decimal conversions performed                */
    int      aborted;  /* node budget hit: result is INCOMPLETE        */
    uint64_t straddle; /* nodes whose interval straddles a power of 10  */
    uint32_t maxm;     /* longest decimal prefix any node ever agreed on */
    int      overflow;    /* task cap exceeded: result INCOMPLETE */
    int      cut_used;    /* cut actually used after back-off       */
} tp_stats;

typedef void (*tp_hit_fn)(void *ud, int B, const bn_t *n,
                          const char *dec, int declen);

typedef struct {
    int             B;
    int             depth;   /* number of fixed orbit choices           */
    uint64_t        path;    /* bit i = value chosen for order[i]       */
} tp_unit;

typedef struct {
    const tp_orbits *o;
    int       nw;
    bn_t      lo, hi;
    int       stop_depth;    /* record unit & return at this depth      */
    uint64_t  path;
    tp_unit  *units;
    int       nunits;
    int       unit_cap;
    int       trace_depth;     /* debug: print decisions for d <= this     */
    uint32_t  r11;             /* lo mod 11, carried down the DFS          */
    uint64_t  r5;              /* lo mod 5^J, carried down the DFS (R5 filter) */
    uint64_t  t5;              /* last J digits N must end in = top J digits reversed, mod 5^J */
    int       t5ok;            /* t5 valid: this subtree shares its first J digits */
    uint64_t  rk, tk;          /* deep filter: lo mod 5^K and the K-digit target */
    int       tkok;            /* tk valid */
    uint64_t  budget;          /* 0 = unlimited; else abort past this      */
    tp_stats  st;
} tp_search;

void tp_search_init(tp_search *s, const tp_orbits *o, int nw, int stop_depth,
                    tp_unit *units, int unit_cap);
void tp_search_roots(const tp_search *s, char *dlo, int *dl, char *dhi, int *dh);
void tp_search_enum_units(tp_search *s);          /* fill s->units            */
void tp_search_run_unit(tp_search *s, const tp_unit *u,
                        tp_hit_fn onhit, void *ud);

/* min number of split units, deterministic in nfree (checkpoint stability) */
int tp_split_depth(int nfree, int target_units);

/* outside-in hex-digit traversal: same orbit set, different tree shape */
int  tp_hexdfs_run(int B, unsigned inv, tp_stats *out,
                   tp_hit_fn onhit, void *ud, int trace);
void tp_hexdfs_set_shard(int k, int n);
int  tp_hexdfs_sharded(void);
void tp_prof_dump(void);
int  tp_hexdfs_run_mt(int B, unsigned inv, int nthreads, int cut,
                      tp_stats *out, tp_hit_fn onhit, void *ud);

/* ------------------------------------------------------------------ */
/* checkpoint                                                          */
/* ------------------------------------------------------------------ */

typedef struct tp_ckpt tp_ckpt;

tp_ckpt *tp_ckpt_open(const char *dir, int enable);
void     tp_ckpt_close(tp_ckpt *c);
int      tp_ckpt_is_done(tp_ckpt *c, int B, int depth, uint64_t path);
void     tp_ckpt_mark(tp_ckpt *c, int B, int depth, uint64_t path);
void     tp_ckpt_hit(tp_ckpt *c, int B, const char *dec, const char *hex);
uint64_t tp_ckpt_loaded(tp_ckpt *c);
const char *tp_ckpt_hits_path(tp_ckpt *c);

/* ------------------------------------------------------------------ */
/* misc                                                                */
/* ------------------------------------------------------------------ */

double tp_now(void);
void   tp_fmt_space(double v, char *buf, size_t n);

#endif /* TRIPAL_H */
