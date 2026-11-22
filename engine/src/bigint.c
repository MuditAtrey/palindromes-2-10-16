#include "tripal.h"

/* ================================================================== */
/* big integer primitives (fixed width, TP_MAX_LIMBS limbs)           */
/* ================================================================== */

void bn_zero(bn_t *a, int nw)
{
    for (int i = 0; i < nw; i++) a->w[i] = 0;
}

void bn_setbit(bn_t *a, int i)
{
    a->w[i >> 6] |= (uint64_t)1 << (i & 63);
}

void bn_clrbit(bn_t *a, int i)
{
    a->w[i >> 6] &= ~((uint64_t)1 << (i & 63));
}

int bn_getbit(const bn_t *a, int i)
{
    return (int)((a->w[i >> 6] >> (i & 63)) & 1u);
}

int bn_is_zero(const bn_t *a, int nw)
{
    uint64_t acc = 0;
    for (int i = 0; i < nw; i++) acc |= a->w[i];
    return acc == 0;
}

void bn_muladd(bn_t *a, int nw, uint64_t mul, uint64_t add)
{
    unsigned __int128 carry = add;
    for (int i = 0; i < nw; i++) {
        unsigned __int128 cur = (unsigned __int128)a->w[i] * mul + carry;
        a->w[i] = (uint64_t)cur;
        carry   = cur >> 64;
    }
}

/* in-place a /= d, returns a % d.  d must be < 2^64. */
uint64_t bn_divmod_small(bn_t *a, int nw, uint64_t d)
{
    unsigned __int128 rem = 0;
    for (int i = nw - 1; i >= 0; i--) {
        unsigned __int128 cur = (rem << 64) | a->w[i];
        a->w[i] = (uint64_t)(cur / d);
        rem     = cur % d;
    }
    return (uint64_t)rem;
}

/* ------------------------------------------------------------------ */
/* decimal / hex rendering                                            */
/* ------------------------------------------------------------------ */

#define TP_DEC_BASE 10000000000000000000ULL   /* 10^19 */
#define TP_DEC_CHUNKS (TP_DECBUF / 19 + 4)

/* Two digits per table lookup, so formatting costs 9 multiply-shifts instead
 * of 19 integer divisions, and no printf machinery at all.  bn_to_dec is the
 * single hottest function in the engine (75% of a DFS node), so this matters. */
static const char D2[201] =
    "00010203040506070809101112131415161718192021222324252627282930313233343536"
    "37383940414243444546474849505152535455565758596061626364656667686970717273"
    "7475767778798081828384858687888990919293949596979899";

/* writes v left-trimmed, returns the new end pointer */
static inline char *put_u64(char *p, uint64_t v)
{
    char tmp[20];
    int n = 0;
    while (v >= 100) {
        const uint64_t q = v / 100;
        const unsigned r = (unsigned)(v - q * 100);
        tmp[n++] = D2[2 * r + 1];
        tmp[n++] = D2[2 * r];
        v = q;
    }
    if (v >= 10) { tmp[n++] = D2[2 * v + 1]; tmp[n++] = D2[2 * v]; }
    else         { tmp[n++] = (char)('0' + v); }
    while (n) *p++ = tmp[--n];
    return p;
}

/* writes exactly 19 digits, zero padded */
static inline char *put_pad19(char *p, uint64_t v)
{
    for (int i = 18; i >= 1; i -= 2) {
        const uint64_t q = v / 100;
        const unsigned r = (unsigned)(v - q * 100);
        p[i]     = D2[2 * r + 1];
        p[i - 1] = D2[2 * r];
        v = q;
    }
    p[0] = (char)('0' + v);
    return p + 19;
}

/* t /= 10^19 in place over limbs [0..top], returns t % 10^19.
 *
 * The generic 128-by-64 division (`unsigned __int128 / d`) compiles to a call
 * into compiler-rt (__udivmodti4, ~50-100 cycles) and it ran 6 limbs x 6 chunks
 * per conversion: more than half of all run time.  10^19 has its top bit set, so
 * it is already a normalised divisor and the Moller-Granlund 2-by-1 division with
 * a precomputed reciprocal applies: one 64x64->128 multiply and a couple of
 * adjustments per limb.  v = floor((2^128 - 1) / 10^19) - 2^64.            */
#define TP_D19 10000000000000000000ULL
static uint64_t g_v19;
static int      g_v19_ready;
static void v19_init(void)
{
    unsigned __int128 all1 = ~(unsigned __int128)0;
    g_v19 = (uint64_t)(all1 / TP_D19);           /* low 64 bits == (that) - 2^64 */
    g_v19_ready = 1;
}

static inline uint64_t divrem19(uint64_t u1, uint64_t u0, uint64_t *q)
{   /* requires u1 < TP_D19 */
    const uint64_t d = TP_D19, v = g_v19;
    unsigned __int128 p = (unsigned __int128)v * u1;
    p += ((unsigned __int128)(u1 + 1) << 64) | u0;
    uint64_t q1 = (uint64_t)(p >> 64), q0 = (uint64_t)p;
    uint64_t r = u0 - q1 * d;
    if (r > q0) { q1--; r += d; }
    if (r >= d) { q1++; r -= d; }
    *q = q1;
    return r;
}

static inline uint64_t bn_divmod_d19(bn_t *a, int top)
{
    uint64_t rem = 0;
    for (int i = top; i >= 0; i--) {
        uint64_t q;
        rem = divrem19(rem, a->w[i], &q);
        a->w[i] = q;
    }
    return rem;
}

int bn_to_dec(const bn_t *a, int nw, char *out)
{
    bn_t t = *a;
    if (!g_v19_ready) v19_init();
    int top = nw - 1;
    while (top >= 0 && t.w[top] == 0) top--;
    if (top < 0) { out[0] = '0'; out[1] = 0; return 1; }

    uint64_t chunk[TP_DEC_CHUNKS];
    int nc = 0;
    while (top >= 0) {
        chunk[nc++] = bn_divmod_d19(&t, top);
        while (top >= 0 && t.w[top] == 0) top--;
        if (nc >= TP_DEC_CHUNKS) break;              /* cannot happen */
    }

    char *p = put_u64(out, chunk[nc - 1]);
    for (int i = nc - 2; i >= 0; i--) p = put_pad19(p, chunk[i]);
    *p = 0;                                   /* keep it a C string */
    return (int)(p - out);
}

int bn_to_hex(const bn_t *a, int nw, char *out)
{
    static const char hx[] = "0123456789abcdef";
    int top = nw - 1;
    while (top > 0 && a->w[top] == 0) top--;
    char *p = out;
    if (a->w[top] == 0) { *p++ = '0'; *p = 0; return 1; }

    int started = 0;
    for (int i = top; i >= 0; i--) {
        for (int s = 60; s >= 0; s -= 4) {
            int d = (int)((a->w[i] >> s) & 0xF);
            if (!started) { if (d == 0) continue; started = 1; }
            *p++ = hx[d];
        }
    }
    *p = 0;
    return (int)(p - out);
}

int bn_from_dec(const char *s, int nw, bn_t *out)
{
    bn_zero(out, nw);
    while (*s == ' ' || *s == '+') s++;
    if (!*s) return 0;
    for (; *s; s++) {
        if (*s < '0' || *s > '9') return 0;
        bn_muladd(out, nw, 10, (uint64_t)(*s - '0'));
    }
    return 1;
}

/* a += v (single limb with carry propagation) */
void bn_add_u64(bn_t *a, int nw, uint64_t v)
{
    unsigned __int128 carry = v;
    for (int i = 0; i < nw && carry; i++) {
        unsigned __int128 cur = (unsigned __int128)a->w[i] + carry;
        a->w[i] = (uint64_t)cur;
        carry   = cur >> 64;
    }
}

/* -1 / 0 / +1 */
int bn_cmp(const bn_t *a, const bn_t *b, int nw)
{
    for (int i = nw - 1; i >= 0; i--)
        if (a->w[i] != b->w[i]) return a->w[i] < b->w[i] ? -1 : 1;
    return 0;
}

int tp_is_pal_str(const char *s, int len)
{
    for (int i = 0, j = len - 1; i < j; i++, j--)
        if (s[i] != s[j]) return 0;
    return 1;
}

int tp_is_pal_base(const bn_t *a, int nw, int base)
{
    bn_t t = *a;
    char buf[TP_MAX_BITS + 4];
    int n = 0;
    if (bn_is_zero(&t, nw)) return 1;               /* 0 is a palindrome  */
    while (!bn_is_zero(&t, nw)) {
        uint64_t d = bn_divmod_small(&t, nw, (uint64_t)base);
        buf[n++] = "0123456789abcdef"[d];
        if (n >= (int)sizeof(buf)) return 0;
    }
    for (int i = 0, j = n - 1; i < j; i++, j--)
        if (buf[i] != buf[j]) return 0;
    return 1;
}
