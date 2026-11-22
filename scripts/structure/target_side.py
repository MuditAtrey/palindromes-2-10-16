"""Decimal-side structure test: do the tree's targets behave like random residues with respect to the inner reachable set?

For hex length L (B = 4L, digits {0,6,9,F}, top digit 9 or F) the outer d digit pairs are enumerated with the engine's pruning
(at every level: shared decimal prefix of the subtree, low digits = prefix reversed, 2-adic check on the bits already fixed).
Each surviving node at depth d gives a 5-adic target t = R - O (mod 5^m) for the inner sum, m = its shared prefix length.
The exact inner reachable set S (mod 5^m) is computed, and the number of targets that land in S is compared with what random
residues would give (sum of |S| / 5^m over targets; Poisson). Pooled over many lengths and depths.
Usage: python3 target_side.py L0 L1 PAIRS   (PAIRS = number of inner digit pairs, 4^PAIRS inner vectors)
"""
import sys, math, numpy as np
A = (0, 6, 9, 15)
L0, L1, P = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
M26 = 5 ** 26
tot_obs = 0; tot_exp = 0.0; tot_t = 0; tot_dist = 0
for L in range(L0, L1 + 1, 2):
    d = L // 2 - P
    if d < 2: continue
    pw = [16 ** i + 16 ** (L - 1 - i) for i in range(L // 2)]
    s = np.zeros(1, dtype=np.int64)
    for i in range(d, L // 2):
        wm = pw[i] % M26
        s = np.concatenate([(s + (a * wm) % M26) % M26 for a in A])
    Scache = {}
    tg = []
    def rec(i, O):
        span = 15 * sum(16 ** j for j in range(i, L - i))
        lo, hi = str(O), str(O + span)
        m = 0
        if len(lo) == len(hi):
            while m < len(lo) and lo[m] == hi[m]: m += 1
        k = min(m, 4 * i)
        if k > 0:
            R = int(lo[:m][::-1]) if m else 0
            if (R - O) % (2 ** k): return
        if i == d:
            if 0 < m <= 26: tg.append((m, (int(lo[:m][::-1]) - O) % (5 ** m)))
            return
        for a in ((9, 15) if i == 0 else A):
            rec(i + 1, O + a * pw[i])
    rec(0, 0)
    obs = 0; exp = 0.0
    for m, t in tg:
        if m not in Scache:
            Sm = np.unique(s % (5 ** m)); Scache[m] = (set(Sm.tolist()) if Sm.size < 3_000_000 else Sm, Sm.size / 5 ** m)
        S, frac = Scache[m]
        obs += (t in S) if isinstance(S, set) else bool(np.isin(t, S)); exp += frac
    dist = len(set(tg))
    print(f"L={L:3d} d={d:2d}: {len(tg):6d} targets ({dist} distinct), m {min((m for m,_ in tg), default=0)}-{max((m for m,_ in tg), default=0)}: in S {obs}, random expectation {exp:.1f}", flush=True)
    tot_obs += obs; tot_exp += exp; tot_t += len(tg); tot_dist += dist
z = (tot_obs - tot_exp) / math.sqrt(tot_exp) if tot_exp > 0 else float('nan')
print(f"pooled: {tot_t} targets ({tot_dist} distinct), {tot_obs} land in the reachable set, random expectation {tot_exp:.1f}, z = {z:+.2f}")
