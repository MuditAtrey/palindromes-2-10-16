"""Positive controls for the exact pruning rule behind Theorem 2 (and the per-length empty proofs).
Rule: for a subtree with all values in [X, X + rest], if every value shares the leading decimal digits P (len m), a decimal palindrome
N satisfies N = rev(P) (mod 10^m), so (rev(P) - X) must be 0 mod 2^k, k = min(m, number of low bits already fixed).
(1) Planted: hex candidate N0 (class 4L or 4L-1) and offset c with N0 + c a decimal palindrome; the shifted rule
    (prefix of [X+c, X+rest+c], test rev(P) - c - X) must keep N0's own branch alive at every depth.
(2) OEIS A007632 (palindromes in bases 2 and 10): along each number's own bit pairs the rule must never fire."""
import random, sys
def prefix(lo, hi):
    a, b = str(lo), str(hi)
    if len(a) != len(b): return ""
    m = 0
    while m < len(a) and a[m] == b[m]: m += 1
    return a[:m]
def path_survives(digits, pw, maxd, lowbits_per_level, c=0):
    """digits: outer-first digit choices of N; pw[i]: weight of digit i (both mirrored positions); returns (ok, depth_failed)"""
    X = 0
    for i, t in enumerate(digits):
        X += t * pw[i]
        rest = maxd * sum(pw[j] for j in range(i + 1, len(pw)))
        P = prefix(X + c, X + rest + c)
        m = len(P); k = min(m, lowbits_per_level * (i + 1))
        if k and (int(P[::-1]) - c - X) % (1 << k): return False, i
    return True, None
rnd = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
bad = tot = 0
for trial in range(400):
    cls = trial % 2; L = rnd.randrange(12, 90)
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7); TOP = (9, 15) if cls == 0 else (5, 7)
    h = (L + 1) // 2
    dig = [rnd.choice(TOP)] + [rnd.choice(A) for _ in range(h - 1)]
    pw = [16 ** j + 16 ** (L - 1 - j) if j != L - 1 - j else 16 ** j for j in range(h)]
    N0 = sum(d * w for d, w in zip(dig, pw))
    D = len(str(N0)) + rnd.randrange(2)
    while True:
        hh = [rnd.randrange(1, 10)] + [rnd.randrange(10) for _ in range((D + 1) // 2 - 1)]
        P0 = int("".join(map(str, hh + hh[: D // 2][::-1])))
        if P0 > N0: break
    ok, f = path_survives(dig, pw, max(A), 4, P0 - N0); tot += 1; bad += (not ok)
    if not ok: print(f"  PLANTED KILLED: class {cls} L {L} depth {f}")
print(f"planted controls (hex classes, L 12-89): {tot - bad}/{tot} planted branches survive every depth")
F = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../engine/data/oeis/b007632.txt")
vals = [int(l.split()[1]) for l in open(F) if l.strip() and not l.startswith('#')]
bad = tot = 0
for N in vals:
    if N < 10: continue
    B = N.bit_length(); h = (B + 1) // 2
    bits = [(N >> (B - 1 - j)) & 1 for j in range(h)]          # outer first: top bit, next, ...
    pw = [2 ** (B - 1 - j) + 2 ** j if j != B - 1 - j else 2 ** j for j in range(h)]
    ok, f = path_survives(bits, pw, 1, 1); tot += 1; bad += (not ok)
    if not ok: print(f"  A007632 KILLED: {N} at depth {f}")
print(f"A007632 controls: {tot - bad}/{tot} known base-2/10 double palindromes survive along their own path")
