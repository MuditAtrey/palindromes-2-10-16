/* Print an LLL-reduced basis of { y in Z^n : sum y_j w_j = 0 mod 5^K } for the inner weights
 * w_0 = 1 (centre), w_j = 16^j + 16^-j (mod 5^K), j = 1..n-1, or random weights (mode 1, w_0 = 1).
 * Output: n lines of n integers.
 * Build: cc -O2 lll_basis.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -o lll_basis */
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_lll.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv)
{
    const int n = atoi(argv[1]), K = atoi(argv[2]), mode = argc > 3 ? atoi(argv[3]) : 0;
    fmpz_t M, t, s, i16, inv; fmpz_init(M); fmpz_init(t); fmpz_init(s); fmpz_init(i16); fmpz_init(inv); fmpz_ui_pow_ui(M, 5, K);
    fmpz *w = _fmpz_vec_init(n);
    fmpz_set_ui(t, 16); fmpz_invmod(i16, t, M);
    flint_rand_t st; flint_rand_init(st);
    for (int j = 0; j < n; j++) {
        if (mode) fmpz_randm(w + j, st, M);
        else { fmpz_set_ui(t, 16); fmpz_powm_ui(t, t, j, M); fmpz_powm_ui(s, i16, j, M); fmpz_add(w + j, t, s); fmpz_mod(w + j, w + j, M); }
    }
    fmpz_set_ui(w + 0, 1);
    fmpz_mat_t A; fmpz_mat_init(A, n, n);
    fmpz_set(fmpz_mat_entry(A, 0, 0), M);
    for (int j = 1; j < n; j++) { fmpz_neg(t, w + j); fmpz_mod(t, t, M); fmpz_set(fmpz_mat_entry(A, j, 0), t); fmpz_one(fmpz_mat_entry(A, j, j)); }
    fmpz_lll_t fl; fmpz_lll_context_init_default(fl); fmpz_lll(A, NULL, fl);
    for (int j = 0; j < n; j++) { fmpz_print(w + j); printf(j + 1 < n ? " " : "\n"); }
    for (int i = 0; i < n; i++) { for (int k = 0; k < n; k++) { fmpz_print(fmpz_mat_entry(A, i, k)); printf(k + 1 < n ? " " : "\n"); } }
    return 0;
}
