/* Kernel lattice { y in Z^n : sum y_j w_j = 0 mod 5^K } for the inner weights w_j = 16^j + 16^-j (j = 0..n-1, w_0 = 1)
 * versus random weights: LLL, then the shortest basis vector and the Gram-Schmidt log-norm profile, which governs
 * enumeration cost. Same profile = the inner lattices are generic; a different one = exploitable structure.
 * Build: cc -O2 gs_profile.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -lm -o gs_profile */
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static void profile(const char *tag, fmpz *w, int n, const fmpz_t M)
{
    fmpz_mat_t A; fmpz_mat_init(A, n, n);
    fmpz_t inv, t; fmpz_init(inv); fmpz_init(t);
    fmpz_invmod(inv, w + 0, M);
    fmpz_set(fmpz_mat_entry(A, 0, 0), M);
    for (int j = 1; j < n; j++) { fmpz_mul(t, w + j, inv); fmpz_neg(t, t); fmpz_mod(t, t, M); fmpz_set(fmpz_mat_entry(A, j, 0), t); fmpz_one(fmpz_mat_entry(A, j, j)); }
    fmpz_lll_t fl; fmpz_lll_context_init_default(fl); fmpz_lll(A, NULL, fl);
    double bs[128][128], n2[128];
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) bs[i][k] = fmpz_get_d(fmpz_mat_entry(A, i, k));
        for (int j = 0; j < i; j++) { double a = 0; for (int k = 0; k < n; k++) a += fmpz_get_d(fmpz_mat_entry(A, i, k)) * bs[j][k]; double mu = a / n2[j]; for (int k = 0; k < n; k++) bs[i][k] -= mu * bs[j][k]; }
        double a = 0; for (int k = 0; k < n; k++) a += bs[i][k] * bs[i][k]; n2[i] = a;
    }
    double l1 = 0, lmin = 1e300; for (int i = 0; i < n; i++) { double a = 0, s = 0; for (int k = 0; k < n; k++) { double v = fmpz_get_d(fmpz_mat_entry(A, i, k)); a += v * v; s += fabs(v); } if (a < lmin) { lmin = a; l1 = s; } }
    const double gh = sqrt(n / (2 * M_PI * M_E)) * pow(fmpz_get_d(M), 1.0 / n);
    printf("%-7s shortest |b| %.2f (L1 %.0f), Gaussian heuristic %.2f, ratio %.3f | GS log2 norms:", tag, sqrt(lmin), l1, gh, sqrt(lmin) / gh);
    for (int i = 0; i < n; i += (n > 16 ? n / 8 : 1)) printf(" %.2f", 0.5 * log2(n2[i]));
    printf(" %.2f\n", 0.5 * log2(n2[n - 1]));
    fmpz_mat_clear(A); fmpz_clear(inv); fmpz_clear(t);
}
int main(int argc, char **argv)
{
    const int n = atoi(argv[1]), K = atoi(argv[2]), trials = argc > 3 ? atoi(argv[3]) : 3;
    fmpz_t M, t, s; fmpz_init(M); fmpz_init(t); fmpz_init(s); fmpz_ui_pow_ui(M, 5, K);
    fmpz *w = _fmpz_vec_init(n);
    fmpz_t i16; fmpz_init(i16); fmpz_set_ui(t, 16); fmpz_invmod(i16, t, M);
    for (int j = 0; j < n; j++) { fmpz_set_ui(t, 16); fmpz_powm_ui(t, t, j, M); fmpz_powm_ui(s, i16, j, M); fmpz_add(w + j, t, s); fmpz_mod(w + j, w + j, M); }
    fmpz_set_ui(w + 0, 1);
    printf("n=%d K=%d (5^K = 2^%.1f, density n*2/log2 M = %.2f)\n", n, K, K * log2(5.0), 2.0 * n / (K * log2(5.0)));
    profile("inner", w, n, M);
    flint_rand_t st; flint_rand_init(st);
    for (int r = 0; r < trials; r++) { for (int j = 1; j < n; j++) fmpz_randm(w + j, st, M); fmpz_set_ui(w + 0, 1); profile("random", w, n, M); }
    return 0;
}
