"""How far can size + 2-adic information go? Fraction of lengths whose exact tree dies by depth k, as k grows.

Same exact rule as early_death.py (shared decimal prefix of the subtree, mirrored, compared 2-adically with the bottom hex digits),
but faster (prefixes by integer division, no full decimal strings) and recording the survivor count at every level, to see
whether the fraction of dead lengths tends to 1 (then almost all lengths are empty by bounded information) or to a limit q < 1
(a supercritical branching process: then the 5-adic half is needed for a positive share of lengths, at depth growing with L).

Usage: python3 death_curve.py CLASS L0 L1 STEP KMAX CAP PROCS OUT.jsonl
"""
import sys, json, bisect
from multiprocessing import Pool

def curve(args):
    L, cls, kmax, cap = args
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7)
    top = (9, 15) if cls == 0 else (5, 7)
    half = (L + 1) // 2
    pw = [16 ** j + 16 ** (L - 1 - j) if j != L - 1 - j else 16 ** j for j in range(half)]
    rest = [0] * (half + 1)
    for j in range(half - 1, -1, -1): rest[j] = rest[j + 1] + max(A) * pw[j]
    P10 = [1]
    while P10[-1] <= 16 ** L * 2: P10.append(P10[-1] * 10)
    level = [0]; counts = []
    kk_max = min(kmax, half)
    for i in range(kk_max):
        nxt = []
        for O in level:
            for a in (top if i == 0 else A):
                X = O + a * pw[i]; hi = X + rest[i + 1]
                D = bisect.bisect_right(P10, X)
                if D != bisect.bisect_right(P10, hi): nxt.append(X); continue
                p = min(D, 4 * (i + 1) + 3)                 # never need more than 4(i+1) shared digits
                q1 = X // P10[D - p]; q2 = hi // P10[D - p]
                s1, s2 = str(q1), str(q2)
                m = 0
                while m < p and s1[m] == s2[m]: m += 1
                m = min(m, D // 2)
                kk = min(m, 4 * (i + 1))
                if kk:
                    R = int(s1[:m][::-1])
                    if (R - X) % (1 << kk): continue
                nxt.append(X)
        level = nxt; counts.append(len(level))
        if not level or len(level) > cap: break
    return {"L": L, "cls": cls, "counts": counts}

if __name__ == "__main__":
    cls, L0, L1, step, kmax, cap, procs, out = sys.argv[1:9]
    cls, L0, L1, step, kmax, cap, procs = map(int, (cls, L0, L1, step, kmax, cap, procs))
    jobs = [(L, cls, kmax, cap) for L in range(L0, L1 + 1, step)]
    with Pool(procs) as pool, open(out, "w") as f:
        for r in pool.imap_unordered(curve, jobs, chunksize=1):
            f.write(json.dumps(r) + "\n"); f.flush()
