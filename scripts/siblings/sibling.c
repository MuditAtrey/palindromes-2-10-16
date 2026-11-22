/* Sibling problems in C (GMP): numbers palindromic in base 2, base 16 and base b, bit length B = 4L.
 *
 * Same algorithm as sibling_search.py, about two orders of magnitude faster, so the structure hunt can go to
 * much longer lengths. Exhaustive in the class B = 0 mod 4 (hex digits from {0,6,9,F}, mirrored).
 * Sound pruning only: the top m base-b digits are shared by every completion (interval [lo, hi]), the low m
 * digits are their mirror, so N mod b^m is known, hence N mod 2^(v m) with v = v2(b); the outer hex digits fix
 * N mod 2^(4d). Disagreement means no completion. b must be even (powers of 2 are the positive controls).
 *
 * Build: cc -O2 -o sibling sibling.c -I/opt/homebrew/include -L/opt/homebrew/lib -lgmp   (Linux: -lgmp)
 * Usage: ./sibling b L        prints "b L nodes nhits hit1 hit2 ..." (hits in decimal)
 */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXL 1024
static int B_, L_, half_, v2b, maxdig;
static mpz_t h[MAXL], suf[MAXL + 1], P[8 * MAXL], X[MAXL + 1], hi, q1, q2, top, R, tmp;
static unsigned long long nodes = 0;
static int nhits = 0;
static char *hitstr[64];

static int ndig(const mpz_t x)          /* number of base-b digits: smallest k with b^k > x */
{
    int a = 0, z = maxdig;
    while (a < z) { int mid = (a + z) / 2; if (mpz_cmp(P[mid], x) > 0) z = mid; else a = mid + 1; }
    return a;
}
static int is_pal_str(const char *s) { size_t n = strlen(s); for (size_t i = 0; i < n / 2; i++) if (s[i] != s[n - 1 - i]) return 0; return 1; }
static int is_pal_base(const mpz_t x, int base)
{
    char *s = mpz_get_str(NULL, base, x); int r = is_pal_str(s); free(s); return r;
}

static void dfs(int d)
{
    if (d == half_) {
        {
            int ok;
            if (B_ <= 62) ok = is_pal_base(X[d], B_);
            else {
                mpz_set(tmp, X[d]); int nd = 0; unsigned long dg[8 * MAXL];
                while (mpz_sgn(tmp)) dg[nd++] = mpz_fdiv_q_ui(tmp, tmp, (unsigned long)B_);
                ok = 1; for (int i = 0; i < nd / 2; i++) if (dg[i] != dg[nd - 1 - i]) ok = 0;
            }
            if (ok && is_pal_base(X[d], 2) && is_pal_base(X[d], 16)) { if (nhits < 64) hitstr[nhits] = mpz_get_str(NULL, 10, X[d]); nhits++; }
        }
        return;
    }
    static const int first[2] = {9, 15}, all[4] = {0, 6, 9, 15};
    const int *al = d == 0 ? first : all; const int na = d == 0 ? 2 : 4;
    for (int i = 0; i < na; i++) {
        nodes++;
        mpz_set(X[d + 1], X[d]); mpz_addmul_ui(X[d + 1], h[d], (unsigned long)al[i]);
        mpz_add(hi, X[d + 1], suf[d + 1]);
        int D = ndig(X[d + 1]);
        if (D == ndig(hi)) {
            int a = 0, z = D;                            /* largest p with equal top-p digits */
            while (a < z) {
                int mid = (a + z + 1) / 2;
                mpz_tdiv_q(q1, X[d + 1], P[D - mid]); mpz_tdiv_q(q2, hi, P[D - mid]);
                if (mpz_cmp(q1, q2) == 0) a = mid; else z = mid - 1;
            }
            int m = a < D / 2 ? a : D / 2;
            if (m > 0) {
                mpz_tdiv_q(top, X[d + 1], P[D - m]);    /* exactly m digits */
                mpz_set_ui(R, 0);
                for (int k = 0; k < m; k++) {           /* R = digits of top reversed: digit k from the top becomes digit k from the bottom */
                    unsigned long dk = mpz_fdiv_q_ui(top, top, (unsigned long)B_);   /* least significant = digit m-1-k from the top */
                    mpz_mul_ui(R, R, (unsigned long)B_); mpz_add_ui(R, R, dk);
                }
                long k2 = (long)v2b * m; if (k2 > 4L * (d + 1)) k2 = 4L * (d + 1);
                mpz_sub(tmp, R, X[d + 1]);
                if (!mpz_divisible_2exp_p(tmp, (mp_bitcnt_t)k2)) continue;
            }
        }
        dfs(d + 1);
    }
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s b L\n", argv[0]); return 2; }
    B_ = atoi(argv[1]); L_ = atoi(argv[2]);
    if (B_ < 4 || (B_ & 1) || L_ < 1 || L_ >= MAXL) { fprintf(stderr, "b must be even and >= 4\n"); return 2; }
    v2b = __builtin_ctz((unsigned)B_);
    half_ = (L_ + 1) / 2;
    mpz_inits(hi, q1, q2, top, R, tmp, NULL);
    for (int j = 0; j < half_; j++) {
        mpz_init(h[j]); mpz_ui_pow_ui(h[j], 16, (unsigned long)j);
        if (j != L_ - 1 - j) { mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, 16, (unsigned long)(L_ - 1 - j)); mpz_add(h[j], h[j], t); mpz_clear(t); }
    }
    for (int j = 0; j <= half_; j++) mpz_init(suf[j]);
    for (int j = half_ - 1; j >= 0; j--) { mpz_set(suf[j], suf[j + 1]); mpz_addmul_ui(suf[j], h[j], 15); }
    maxdig = 0; mpz_init_set_ui(P[0], 1);
    mpz_t lim; mpz_init(lim); mpz_ui_pow_ui(lim, 2, (unsigned long)(4 * L_ + 8));
    while (mpz_cmp(P[maxdig], lim) <= 0 && maxdig < 8 * MAXL - 1) { mpz_init(P[maxdig + 1]); mpz_mul_ui(P[maxdig + 1], P[maxdig], (unsigned long)B_); maxdig++; }
    for (int j = 0; j <= half_; j++) mpz_init(X[j]);
    dfs(0);
    printf("%d %d %llu %d", B_, L_, nodes, nhits);
    for (int i = 0; i < nhits && i < 64; i++) printf(" %s", hitstr[i]);
    printf("\n");
    return 0;
}
