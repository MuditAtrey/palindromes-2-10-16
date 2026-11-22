"""Null model for the J-ladder: split the weak-palindrome condition into its 2-adic and 5-adic halves.

For every base-2/16 palindrome with B = 4L bits (hex digits from {0,6,9,F}, mirrored), and each J:
  R_J   = the value the last J decimal digits must have (the first J digits reversed)
  A_J   = #{N : N = R_J mod 2^J}    (2-adic half; coupled to the top bits by the binary palindrome)
  F_J   = #{N : N = R_J mod 5^J}    (5-adic half; the part the engine must test globally)
  W_J   = #{N : both}               (= weak palindrome, last J digits mirror first J)
If the 5-adic half behaves like chance, F_J / total ~ 5^-J and W_J ~ A_J * 5^-J. A deviation in the 5-adic half would be a real lead;
a deviation explained by A_J is the known 2-adic coupling. Exhaustive over all candidates, so counts are exact.

Usage: python3 j_null.py L0 L1 JMAX [procs]
"""
import itertools, math, sys
from multiprocessing import Pool
sys.set_int_max_str_digits(0)
JMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 8
if __name__ == '__main__':
    L0, L1 = int(sys.argv[1]), int(sys.argv[2]); PROCS = int(sys.argv[4]) if len(sys.argv) > 4 else 8
def work(args):
    L, first = args
    half = (L + 1) // 2
    h = [(16 ** j + 16 ** (L - 1 - j)) if j != L - 1 - j else 16 ** j for j in range(half)]
    A = [0] * (JMAX + 1); F = [0] * (JMAX + 1); W = [0] * (JMAX + 1); tot = 0
    P2 = [1 << J for J in range(JMAX + 1)]; P5 = [5 ** J for J in range(JMAX + 1)]
    base = first * h[0]
    for rest in itertools.product((0, 6, 9, 15), repeat=half - 1):
        N = base
        for e, hj in zip(rest, h[1:]): N += e * hj
        s = str(N); tot += 1
        for J in range(1, JMAX + 1):
            R = int(s[J - 1::-1])
            a = (N - R) % P2[J] == 0; f = (N - R) % P5[J] == 0
            A[J] += a; F[J] += f; W[J] += a and f
    return L, tot, A, F, W
if __name__ == '__main__':
    jobs = [(L, f) for L in range(L0, L1 + 1) for f in (9, 15)]
    tot = 0; A = [0] * (JMAX + 1); F = [0] * (JMAX + 1); W = [0] * (JMAX + 1)
    with Pool(PROCS) as pool:
        for L, t, a, f, w in pool.imap_unordered(work, jobs):
            tot += t
            for J in range(JMAX + 1): A[J] += a[J]; F[J] += f[J]; W[J] += w[J]
    print(f"all {tot:,d} base-2/16 palindromes with {4*L0}..{4*L1} bits (B = 0 mod 4)")
    print(f"{'J':>3} {'2-adic A/tot':>13} {'x 2^J':>8} {'5-adic F/tot':>13} {'x 5^J':>8} {'z(5-adic)':>10} {'weak W':>10} {'W/(A 5^-J)':>11}")
    for J in range(1, JMAX + 1):
        pa, pf = A[J] / tot, F[J] / tot
        z = (F[J] - tot * 5 ** -J) / math.sqrt(tot * 5 ** -J * (1 - 5 ** -J))
        ratio = W[J] / (A[J] * 5 ** -J) if A[J] else float('nan')
        print(f"{J:>3} {pa:>13.6f} {pa * 2**J:>8.3f} {pf:>13.8f} {pf * 5**J:>8.3f} {z:>+10.2f} {W[J]:>10} {ratio:>11.3f}")
    