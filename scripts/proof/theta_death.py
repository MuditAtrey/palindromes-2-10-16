"""Rigorous theta-intervals with no triple palindromes, for all large L at once.
A candidate with hex length L is N = 16^(L-1) * T with T in [t0, t0+1) given by its hex digits (top digit t0). Its leading decimal
digits are those of T * 10^theta, theta = {(L-1) log10 16}, up to the normalisation to [1,10). The 2-adic check of the outside-in tree
compares the reversed shared decimal prefix with the bottom hex digits sum_{j<i} t_j 16^j (mod 2^k), which does not depend on L.
For a theta-interval J the tree is run with interval arithmetic (the shared prefix of the smallest and largest value over J and over all
completions); if every branch dies, then for every L >= L_min with theta_L in J there is no triple palindrome with that bit length.
Bisection finds such intervals; prints their total measure (a rigorous lower bound for the density of empty lengths).
Usage: python3 theta_death.py CLASS DEPTH BISECT_LEVELS   (CLASS 0: B = 4L, digits {0,6,9,F}; 1: B = 4L-1, digits {0,2,5,7})"""
import sys
from decimal import Decimal, getcontext, ROUND_FLOOR, ROUND_CEILING
getcontext().prec = 60
cls, K, LEV = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])) if __name__ == "__main__" else (0, 6, 12)
A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7); TOP = (9, 15) if cls == 0 else (5, 7); MAXD = max(A)
LN10 = Decimal(10).ln()
def pow10(t): return (Decimal(t) * LN10).exp()
SIXTEEN = Decimal(16)
def shared_prefix(lo, hi):
    # lo, hi positive Decimals; normalise lo to [1,10); if hi crosses a power of 10 the length is ambiguous -> no shared prefix
    # widen both ends far beyond the arithmetic error (exp/ln at 60 digits) and round outward, so a shared digit is never created by rounding
    lo = lo * (1 - Decimal(1).scaleb(-50)); hi = hi * (1 + Decimal(1).scaleb(-50))
    e = lo.adjusted(); lo_n = lo.scaleb(-e); hi_n = hi.scaleb(-e)
    if lo_n < 1 or hi_n >= 10: return ""
    a = format(lo_n.quantize(Decimal(1).scaleb(-40), rounding=ROUND_FLOOR), 'f').replace('.', '')
    b = format(hi_n.quantize(Decimal(1).scaleb(-40), rounding=ROUND_CEILING), 'f').replace('.', '')
    m = 0
    while m < min(len(a), len(b)) - 5 and a[m] == b[m]: m += 1
    return a[:m]
def dies(th_lo, th_hi):
    f_lo, f_hi = pow10(th_lo), pow10(th_hi)
    level = [(Decimal(0), 0)]          # (top value as hex fraction t0.t1t2..., bottom value sum t_j 16^j)
    for i in range(K):
        tail = Decimal(MAXD) * sum(SIXTEEN ** (-j) for j in range(i + 1, 200))   # largest contribution of the remaining digits
        nxt = []
        for top, bot in level:
            for t in (TOP if i == 0 else A):
                T = top + Decimal(t) * SIXTEEN ** (-i); B2 = bot + t * 16 ** i
                p = shared_prefix(T * f_lo, (T + tail) * f_hi + Decimal(1).scaleb(-45))
                m = len(p); k = min(m, 4 * (i + 1))
                if k:
                    R = int(p[::-1])
                    if (R - B2) % (1 << k): continue
                nxt.append((T, B2))
        level = nxt
        if not level: return True
    return False
def configure(c, k):
    global cls, K, A, TOP, MAXD
    cls, K = c, k
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7); TOP = (9, 15) if cls == 0 else (5, 7); MAXD = max(A)
good = []
def rec(a, b, lev):
    if dies(a, b): good.append((a, b)); return
    if lev == 0: return
    mid = (a + b) / 2; rec(a, mid, lev - 1); rec(mid, b, lev - 1)
if __name__ == "__main__":
  rec(Decimal(0), Decimal(1), LEV)
  meas = sum(b - a for a, b in good)
  print(f"class {'B=4L' if cls == 0 else 'B=4L-1'}, depth {K}: theta-intervals proved empty: {len(good)}, total measure {float(meas):.4f}")
  merged = []
  for a, b in sorted(good):
      if merged and merged[-1][1] >= a: merged[-1] = (merged[-1][0], max(b, merged[-1][1]))
      else: merged.append((a, b))
  print("  intervals:", [(round(float(a), 5), round(float(b), 5)) for a, b in merged[:30]])
