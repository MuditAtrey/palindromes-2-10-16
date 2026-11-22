# Positive control for the triple form: palindromic in base 2, base 16 (digits 0,6,9,F) and base b.
# Small ones exist for several b; along each one's own hex path the exact certificate must FIND it.
import itertools, sys
from p_int import certify
def digs(n, b):
    s = []
    while n: s.append(n % b); n //= b
    return s[::-1]
def state(lo, hi, b):
    a, c = digs(lo, b), digs(hi, b)
    if len(a) != len(c): return None
    D = len(a); p = 0
    while p < D and a[p] == c[p]: p += 1
    return D, a[:min(p, D // 2)]
def equation(L, outer, b):
    h = lambda j: (16 ** j + 16 ** (L - 1 - j)) if j != L - 1 - j else 16 ** j
    d = len(outer); X = sum(e * h(j) for j, e in enumerate(outer))
    inner = list(range(d, (L + 1) // 2)); hi = X + sum(15 * h(j) for j in inner)
    st = state(X, hi, b)
    if st is None: return None
    D, s = st; m = len(s)
    Y = sum(s[i] * (b ** (D - 1 - i) + b ** i) for i in range(m))
    coefs, ranges = [], []
    for j in inner: coefs += [6 * h(j), 9 * h(j)]; ranges += [(0, 1), (0, 1)]
    for i in range(m, (D + 1) // 2):
        coefs.append(-((b ** i + b ** (D - 1 - i)) if i != D - 1 - i else b ** i)); ranges.append((0, b - 1))
    return coefs, ranges, Y - X
bad = found = total = 0; amb = 0
for b in (3, 4, 5, 7, 9, 11, 13, 19):
    sols = []
    for L in range(6, 15):
        half = (L + 1) // 2
        for first in (9, 15):
            for rest in itertools.product((0, 6, 9, 15), repeat=half - 1):
                e = (first,) + rest
                hx = list(e) + list(e[: L // 2][::-1])
                N = int("".join("0123456789ABCDEF"[x] for x in hx), 16)
                ds = digs(N, b)
                if ds == ds[::-1] and bin(N)[2:] == bin(N)[2:][::-1]: sols.append((L, e, N))
    for L, e, N in sols:
        for d in range(1, (L + 1) // 2):
            eqn = equation(L, list(e[:d]), b)
            if eqn is None: amb += 1; continue
            # (the certificate sees only the outer digits; inner solution must be found)
            # scale ranges for base b: weights need integer 18/(hi-lo); use lcm-friendly bases only
            ok, nn, real = certify(*eqn, 400000)
            total += 1
            if real: found += 1
            elif ok: bad += 1; print("  MISSED", b, L, e, d, N)
    print(f"base {b}: {len(sols)} triple palindromes with L=6..14", flush=True)
print(f"checked {total} (number, depth) pairs: found {found}, wrongly dead {bad}, budget {total-found-bad}, ambiguous length {amb}")
