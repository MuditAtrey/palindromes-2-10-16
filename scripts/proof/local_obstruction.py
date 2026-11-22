"""Is there a congruence obstruction? For a modulus m (coprime to 10), compute the exact set of residues mod m of
  - all base-2/16 palindromes with B bits (B = 4L: hex digits {0,6,9,F}, top 9 or F; B = 4L-1: hex digits {0,2,5,7}, top 5 or 7), and
  - all D-digit decimal palindromes, for each D that numbers with B bits can have,
and test whether the two sets are disjoint. A disjoint pair would prove that no solution has that (B, D): a proof by congruences.
"""
import sys, math
import numpy as np
def shift_or(S, a, m):
    a %= m
    return S | np.roll(S, a)
def sumset(weights_alph, m, lead):
    S = np.zeros(m, bool)
    for a in lead[1]: S[(a * lead[0]) % m] = True
    for w, alph in weights_alph:
        T = np.zeros(m, bool)
        for a in alph: T |= np.roll(S, (a * w) % m)
        S = T
    return S
def hex_res(L, cls, m):
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7)
    top = (9, 15) if cls == 0 else (5, 7)
    pw = lambda j: (pow(16, j, m) + pow(16, L - 1 - j, m)) % m if j != L - 1 - j else pow(16, j, m)
    return sumset([(pw(j), A) for j in range(1, (L + 1) // 2)], m, (pw(0), top))
def dec_res(D, m):
    pw = lambda i: (pow(10, i, m) + pow(10, D - 1 - i, m)) % m if i != D - 1 - i else pow(10, i, m)
    return sumset([(pw(i), range(10)) for i in range(1, (D + 1) // 2)], m, (pw(0), range(1, 10)))
if __name__ == "__main__":
    M = int(sys.argv[1]); L0, L1 = int(sys.argv[2]), int(sys.argv[3])
    blocked = []; tested = 0
    for m in range(3, M + 1):
        if math.gcd(m, 10) != 1: continue
        for L in range(L0, L1 + 1):
            for cls in (0, 1):
                B = 4 * L - cls
                Dlo, Dhi = len(str(2 ** (B - 1))), len(str(2 ** B - 1))
                H = hex_res(L, cls, m)
                for D in range(Dlo, Dhi + 1):
                    tested += 1
                    if not (H & dec_res(D, m)).any(): blocked.append((m, B, D))
    print(f"moduli up to {M} coprime to 10, hex lengths {L0}-{L1}: {tested} (m, B, D) cases, {len(blocked)} with no common residue")
    for b in blocked[:20]: print("  blocked:", b)
