# Positive control: known base-2/base-10 double palindromes (OEIS A007632).  Along the number's own
# outer bits, the exact integer-equation certificate must FIND it (never call the node dead).
import sys, time, random
from p_int import certify, decstate
F = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../engine/data/oeis/b007632.txt")
vals = [int(l.split()[1]) for l in open(F) if l.strip() and not l.startswith('#')]
def eq(N0, B, dbits, bits=None):
    b = [(N0 >> i) & 1 for i in range(B)] if bits is None else bits
    h = lambda i: (1 << i) + (1 << (B - 1 - i)) if i != B - 1 - i else 1 << i
    X = sum(b[i] * h(i) for i in range(dbits))
    inner = list(range(dbits, (B + 1) // 2))
    hi = X + sum(h(i) for i in inner)
    st = decstate(X, hi)
    if st is None: return None
    D, s = st; m = len(s)
    Y = sum(int(s[i]) * (10 ** (D - 1 - i) + 10 ** i) for i in range(m))
    coefs = [h(i) for i in inner]; ranges = [(0, 1)] * len(inner)
    for i in range(m, (D + 1) // 2):
        coefs.append(-((10 ** i + 10 ** (D - 1 - i)) if i != D - 1 - i else 10 ** i)); ranges.append((0, 9))
    return coefs, ranges, Y - X, m
lim = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
big = [v for v in vals if v.bit_length() >= int(sys.argv[2])][: int(sys.argv[3])] if len(sys.argv) > 3 else vals[-6:]
for N0 in big:
    B = N0.bit_length(); assert bin(N0)[2:] == bin(N0)[2:][::-1] and str(N0) == str(N0)[::-1]
    for dbits in (B // 8, B // 5, B // 4, B // 3):
        e = eq(N0, B, dbits)
        if e is None: print(f"  B={B} dbits={dbits}: D ambiguous"); continue
        coefs, ranges, T, m = e
        t0 = time.time(); ok, nn, real = certify(coefs, ranges, T, lim)
        hit = any(True for _ in real)
        # mutate one outer bit pair (keeping the ends 1): must almost always be dead
        b = [(N0 >> i) & 1 for i in range(B)]; i = random.Random(B * 7 + dbits).randrange(1, dbits); b[i] ^= 1; b[B - 1 - i] ^= 1
        e2 = eq(N0, B, dbits, b); res2 = "D ambiguous"
        if e2:
            ok2, nn2, real2 = certify(*e2[:3], lim); res2 = ("FOUND" if real2 else ("dead" if ok2 else "budget")) + f" ({nn2} nodes)"
        print(f"  B={B} dbits={dbits:3d} vars={len(coefs):3d} m={m:2d}: own path {'FOUND' if hit else ('DEAD (BAD)' if ok else 'budget')} ({nn} nodes, {time.time()-t0:.1f}s) | mutated path {res2}", flush=True)
