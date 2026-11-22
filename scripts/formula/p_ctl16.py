# Control with a non-power-of-2 second base: palindromic in base 16 (any digits) and base 10 (OEIS A029731).
# Hex digit = x1 + 4*x2, x1,x2 in 0..3.  Along each number's own outer hex digits, the certificate must find it.
import itertools
from p_int import certify, decstate
sols = []
for L in range(6, 11):
    half = (L + 1) // 2
    for first in range(1, 16):
        for rest in itertools.product(range(16), repeat=half - 1):
            e = (first,) + rest
            N = 0
            for x in list(e) + list(e[: L // 2][::-1]): N = N * 16 + x
            s = str(N)
            if s == s[::-1]: sols.append((L, e, N))
print(len(sols), "base 16/10 double palindromes with 6..10 hex digits, e.g.", [x[2] for x in sols[:6]])
found = bad = bud = 0
for L, e, N in sols:
    h = lambda j: (16 ** j + 16 ** (L - 1 - j)) if j != L - 1 - j else 16 ** j
    for d in range(1, (L + 1) // 2):
        X = sum(e[j] * h(j) for j in range(d)); inner = list(range(d, (L + 1) // 2))
        hi = X + sum(15 * h(j) for j in inner); st = decstate(X, hi)
        if st is None: continue
        D, s = st; m = len(s)
        Y = sum(int(s[i]) * (10 ** (D - 1 - i) + 10 ** i) for i in range(m))
        coefs, ranges = [], []
        for j in inner: coefs += [h(j), 4 * h(j)]; ranges += [(0, 3), (0, 3)]
        for i in range(m, (D + 1) // 2):
            coefs.append(-((10 ** i + 10 ** (D - 1 - i)) if i != D - 1 - i else 10 ** i)); ranges.append((0, 9))
        ok, nn, real = certify(coefs, ranges, Y - X, 400000)
        if real: found += 1
        elif ok: bad += 1; print("  MISSED", L, e, d, N)
        else: bud += 1
print(f"(number, depth) pairs: found {found}, wrongly dead {bad}, budget {bud}")
