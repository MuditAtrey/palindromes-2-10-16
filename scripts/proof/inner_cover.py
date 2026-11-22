"""Section 2 of docs/WALL.md: how many inner digit pairs make the inner sums cover every residue mod 5^m?
For random lengths L and random starting depths k, take the n positions j = k .. k+n-1 with weights W_j = 16^j + 16^(L-1-j) mod 5^m,
digits {0,6,9,15} (class B = 4L) or {0,2,5,7} (class B = 4L-1), and find the smallest n for which every residue is reachable."""
import random, sys
def min_cover(L, k, m, A):
    q = 5 ** m; S = {0}; n = 0
    while len(S) < q:
        j = k + n; W = (pow(16, j, q) + pow(16, L - 1 - j, q)) % q
        S = {(s + a * W) % q for s in S for a in A}; n += 1
        if 2 * (k + n) > L: return None
    return n
random.seed(1)
for A in ((0, 6, 9, 15), (0, 2, 5, 7)):
    for m in range(1, 9):
        worst = 0
        for _ in range(60):
            L = random.randint(400, 4000); k = random.randint(1, 30)
            worst = max(worst, min_cover(L, k, m, A))
        print(f"digits {A} mod 5^{m}: inner sums cover every residue once n >= {worst} (worst of 60 random L, k); 4^n vs 5^m: n_min = {m*1.161:.1f}")
