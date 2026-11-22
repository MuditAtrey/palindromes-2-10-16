"""Independent re-implementation of Theorem 2's interval proof (different code, different arithmetic).
mpmath interval arithmetic (outward rounding) for 10^theta and T * 10^theta; decimal prefix digits by integer floors of the interval ends;
its own bisection. Compares its proved theta-set with theta_death.py's."""
import sys
from fractions import Fraction
from mpmath import iv, mp
iv.dps = 50; mp.dps = 60
def run(cls, K, LEV):
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7); TOP = (9, 15) if cls == 0 else (5, 7); MAXD = max(A)
    def prefix(v):
        lo, hi = mp.mpf(v.a) * (1 - mp.mpf(10) ** -45), mp.mpf(v.b) * (1 + mp.mpf(10) ** -45)   # outward, plus margin for mp rounding below
        if lo <= 0: return 0, 0
        e = int(iv.mpf(lo).a.__floor__()) if False else None
        import math
        el = int(mp.floor(mp.log10(lo))); eh = int(mp.floor(mp.log10(hi)))
        if el != eh: return 0, 0
        # largest m with floor(lo * 10^(m-1-el)) == floor(hi * 10^(m-1-el))
        best = (0, 0)
        for m in range(1, 30):
            a = int(mp.floor(lo * mp.mpf(10) ** (m - 1 - el))); b = int(mp.floor(hi * mp.mpf(10) ** (m - 1 - el)))
            if a != b: break
            best = (m, a)
        return best
    def dies(ta, tb):
        tlo = (iv.mpf(ta.numerator) / ta.denominator).a; thi = (iv.mpf(tb.numerator) / tb.denominator).b
        f = iv.mpf(10) ** iv.mpf([tlo, thi])
        level = [(Fraction(0), 0)]
        for i in range(K):
            tail = Fraction(MAXD, 15) * Fraction(1, 16 ** i)          # MAXD * sum_{j > i} 16^-j, exact
            nxt = []
            for T, bot in level:
                for t in (TOP if i == 0 else A):
                    T2 = T + Fraction(t, 16 ** i); b2 = bot + t * 16 ** i
                    v = iv.mpf([T2.numerator, T2.numerator]) / T2.denominator
                    w = iv.mpf([(T2 + tail).numerator, (T2 + tail).numerator]) / (T2 + tail).denominator
                    val = iv.mpf([v.a, w.b]) * f
                    m, P = prefix(val)
                    k = min(m, 4 * (i + 1))
                    if k and (int(str(P)[::-1]) - b2) % (1 << k): continue
                    nxt.append((T2, b2))
            level = nxt
            if not level: return True
        return False
    good = []
    def rec(a, b, lev):
        if dies(a, b): good.append((a, b)); return
        if lev == 0: return
        mid = (a + b) / 2; rec(a, mid, lev - 1); rec(mid, b, lev - 1)
    rec(Fraction(0), Fraction(1), LEV)
    return good
if __name__ == "__main__":
    import theta_death as TD
    from decimal import Decimal
    for cls in (0, 1):
        g_ind = run(cls, 6, 12)
        TD.configure(cls, 6); TD.good.clear(); TD.rec(Decimal(0), Decimal(1), 12)
        g_td = [(Fraction(str(a)), Fraction(str(b))) for a, b in TD.good]
        S1 = set(g_ind); S2 = set(g_td)
        m1 = float(sum(b - a for a, b in g_ind)); m2 = float(sum(b - a for a, b in g_td))
        print(f"class {cls}: independent proof measure {m1:.4f} ({len(g_ind)} intervals); theta_death {m2:.4f} ({len(g_td)}); "
              f"identical interval sets: {S1 == S2}; only in theta_death: {len(S2 - S1)}; only in independent: {len(S1 - S2)}")
