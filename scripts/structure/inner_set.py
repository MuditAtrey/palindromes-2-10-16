"""Is the set of reachable inner sums structured or random?

Inner sum at a node: S = sum_tau e_tau * (16^(c+tau) + 16^(c-tau)) mod 5^K over the innermost pairs tau = 1..m
(plus the centre digit, tau = 0, counted once), digits e from the alphabet ({0,6,9,15} for B = 4L, {0,2,5,7} for
B = 4L-1). The factor 16^c is a unit, so it is dropped. Compare the number of distinct residues with the same count for
random weights (same alphabet): many fewer distinct values = structured coincidences that a fast test could use.
"""
import sys, numpy as np
def sums(W, alph, M):
    s = np.zeros(1, dtype=np.int64)
    for w in W:
        s = np.concatenate([(s + a * w) % M for a in alph])
    return s
def distinct(s): return np.unique(s).size
if __name__ == "__main__":
    m = int(sys.argv[1]); alph = [0, 6, 9, 15] if (len(sys.argv) < 3 or sys.argv[2] == "0") else [0, 2, 5, 7]
    rng = np.random.default_rng(1)
    n = 4 ** (m + 1)
    print(f"m={m} (+centre): {m+1} digits, 4^{m+1} = {n} digit vectors, alphabet {alph}")
    for K in range(6, 27):
        M = 5 ** K
        if M > 2**62: break
        W = [1] + [(pow(16, t, M) + pow(16, -t, M)) % M for t in range(1, m + 1)]
        d = distinct(sums(W, alph, M))
        Wr = [int(x) for x in rng.integers(1, M, size=m + 1)]
        dr = distinct(sums(Wr, alph, M))
        exp = M * (1 - np.exp(-n / M))
        print(f"K={K:2d}  5^K/4^u={M/n:10.3g}  distinct: inner {d:10d}  random weights {dr:10d}  ideal {exp:12.0f}   inner/random {d/dr:.4f}")
