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
