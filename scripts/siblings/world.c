/* Sibling worlds: numbers palindromic in base 2, in a power-of-2 base r and in an even base b.
 *
 * Generalises sibling.c (r = 16, digits {0,6,9,F}) to any world given by a radix r = 2^k, a digit alphabet A and a set T of
 * allowed top digits. A number with exactly k*M bits that is a palindrome in base 2 and in base r has every base-r digit equal
 * to its own k-bit reversal, so its base-r digits come from A = {self-reversing k-bit digits}, mirrored. Worlds used:
 *   r = 16, A = {0,6,9,15}, T = {9,15}   bit length 4M      (the class of the real problem, with b = 10)
 *   r = 16, A = {0,2,5,7},  T = {5,7}    bit length 4M - 1  (the other large class of the real problem)
 *   r = 8,  A = {0,2,5,7},  T = {5,7}    bit length 3M      (palindromes in bases 2 and 8)
 *   r = 4,  A = {0,3},      T = {3}      bit length 2M      (palindromes in bases 2 and 4)
 * Sound pruning only (same rule as sibling.c): if every completion shares its top m base-b digits, the low m base-b digits are their
 * mirror, so N mod 2^(v m) is known (v = v2(b)); the outer d base-r digits fix N mod 2^(k d). Disagreement means no completion.
 * Every leaf is checked in full (palindrome in 2, r and b), so a wrong alphabet could only lose solutions, never invent them.
 *
 * Build: cc -O2 -o world world.c -I/opt/homebrew/include -L/opt/homebrew/lib -lgmp
 * Usage: ./world r A T b M      e.g. ./world 16 0,6,9,15 9,15 10 30
 * Output: "nodes nhits hit1 hit2 ..." (hits in decimal, at most 64 listed)
 */
#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXM 1100
#define MAXP 4400
static int R_, K_, B_, M_, half_, v2b, maxdig, nA, nT, A[16], T[16];
static mpz_t h[MAXM], suf[MAXM + 1], P[MAXP], X[MAXM + 1], hi, q1, q2, top, Rv, tmp;
static unsigned long long nodes = 0;
static int nhits = 0;
static char *hitstr[64];

static int ndig(const mpz_t x)
{
    int a = 0, z = maxdig;
    while (a < z) { int mid = (a + z) / 2; if (mpz_cmp(P[mid], x) > 0) z = mid; else a = mid + 1; }
    return a;
}
static int is_pal_base(const mpz_t x, int base)
{
    static unsigned long dg[MAXP];
    int nd = 0; mpz_set(tmp, x);
    while (mpz_sgn(tmp)) dg[nd++] = mpz_fdiv_q_ui(tmp, tmp, (unsigned long)base);
    for (int i = 0; i < nd / 2; i++) if (dg[i] != dg[nd - 1 - i]) return 0;
    return 1;
}
static int parse(const char *s, int *out)
{
    int n = 0; char *buf = strdup(s), *tok = strtok(buf, ",");
    while (tok && n < 16) { out[n++] = atoi(tok); tok = strtok(NULL, ","); }
    free(buf); return n;
}

static void dfs(int d)
{
    if (d == half_) {
        if (is_pal_base(X[d], B_) && is_pal_base(X[d], 2) && is_pal_base(X[d], R_)) {
            if (nhits < 64) hitstr[nhits] = mpz_get_str(NULL, 10, X[d]);
            nhits++;
        }
        return;
    }
    const int *al = d == 0 ? T : A; const int na = d == 0 ? nT : nA;
    for (int i = 0; i < na; i++) {
        nodes++;
        mpz_set(X[d + 1], X[d]); mpz_addmul_ui(X[d + 1], h[d], (unsigned long)al[i]);
        mpz_add(hi, X[d + 1], suf[d + 1]);
        int D = ndig(X[d + 1]);
        if (D == ndig(hi)) {
            int a = 0, z = D;
            while (a < z) {
                int mid = (a + z + 1) / 2;
                mpz_tdiv_q(q1, X[d + 1], P[D - mid]); mpz_tdiv_q(q2, hi, P[D - mid]);
                if (mpz_cmp(q1, q2) == 0) a = mid; else z = mid - 1;
            }
            int m = a < D / 2 ? a : D / 2;
            if (m > 0) {
                mpz_tdiv_q(top, X[d + 1], P[D - m]);
                mpz_set_ui(Rv, 0);
                for (int k = 0; k < m; k++) {
                    unsigned long dk = mpz_fdiv_q_ui(top, top, (unsigned long)B_);
                    mpz_mul_ui(Rv, Rv, (unsigned long)B_); mpz_add_ui(Rv, Rv, dk);
                }
                long k2 = (long)v2b * m; if (k2 > (long)K_ * (d + 1)) k2 = (long)K_ * (d + 1);
                mpz_sub(tmp, Rv, X[d + 1]);
                if (!mpz_divisible_2exp_p(tmp, (mp_bitcnt_t)k2)) continue;
            }
        }
        dfs(d + 1);
    }
}

int main(int argc, char **argv)
{
    if (argc < 6) { fprintf(stderr, "usage: %s r A T b M\n", argv[0]); return 2; }
    R_ = atoi(argv[1]); nA = parse(argv[2], A); nT = parse(argv[3], T); B_ = atoi(argv[4]); M_ = atoi(argv[5]);
    if (R_ < 4 || (R_ & (R_ - 1)) || B_ < 4 || (B_ & 1) || M_ < 1 || M_ >= MAXM) { fprintf(stderr, "bad arguments\n"); return 2; }
    K_ = __builtin_ctz((unsigned)R_); v2b = __builtin_ctz((unsigned)B_);
    int maxa = 0; for (int i = 0; i < nA; i++) if (A[i] > maxa) maxa = A[i];
    half_ = (M_ + 1) / 2;
    mpz_inits(hi, q1, q2, top, Rv, tmp, NULL);
    for (int j = 0; j < half_; j++) {
        mpz_init(h[j]); mpz_ui_pow_ui(h[j], (unsigned long)R_, (unsigned long)j);
        if (j != M_ - 1 - j) { mpz_t t; mpz_init(t); mpz_ui_pow_ui(t, (unsigned long)R_, (unsigned long)(M_ - 1 - j)); mpz_add(h[j], h[j], t); mpz_clear(t); }
    }
    for (int j = 0; j <= half_; j++) mpz_init(suf[j]);
    for (int j = half_ - 1; j >= 0; j--) { mpz_set(suf[j], suf[j + 1]); mpz_addmul_ui(suf[j], h[j], (unsigned long)maxa); }
    maxdig = 0; mpz_init_set_ui(P[0], 1);
    mpz_t lim; mpz_init(lim); mpz_ui_pow_ui(lim, 2, (unsigned long)(K_ * M_ + 8));
    while (mpz_cmp(P[maxdig], lim) <= 0 && maxdig < MAXP - 1) { mpz_init(P[maxdig + 1]); mpz_mul_ui(P[maxdig + 1], P[maxdig], (unsigned long)B_); maxdig++; }
    for (int j = 0; j <= half_; j++) mpz_init(X[j]);
    dfs(0);
    printf("%llu %d", nodes, nhits);
    for (int i = 0; i < nhits && i < 64; i++) printf(" %s", hitstr[i]);
    printf("\n");
    return 0;
}
