# Windows from the exact integer equation.  For y in the dual of the kernel lattice (weighted coordinates z in [-9,9]^n):
#   <y,z> = <y,aff> mod 1  and  |<y,z>| <= 9 |y|_1 ,  so if 18|y|_1 < 1 most cosets are excluded.
import sys, random, math
import numpy as np
from p_int import node_equation, survivors
from p_bkz import setup, bkz
from p_scale import dive
L = int(sys.argv[1]); ds = [int(x) for x in sys.argv[2].split(',')]; S = int(sys.argv[3]); beta = int(sys.argv[4])
if len(sys.argv) > 5: print("python tree sizes (2-adic prune only):", [len(v) for v in survivors(L, int(sys.argv[5])).values()])
rnd = random.Random(L + 1)
for d in ds:
    dead = 0; got = 0; best = []
    for _ in range(S):
        g = dive(L, d, rnd)
        if g is None: continue
        coefs, ranges, T, m, D = node_equation(L, g)
        ker, aff, n = setup(coefs, ranges, T)
        Kb = bkz(ker, beta, 2) if beta > 2 else ker
        K = np.array(Kb, dtype=float)                       # (n-1) x n
        Y = np.linalg.solve(K @ K.T, K)                     # dual basis rows: <Y_i, K_j> = delta_ij
        a = np.array(aff[:n], dtype=float)
        # also try small combinations of dual basis vectors
        cand = [Y[i] for i in range(len(Y))]
        for i in range(len(Y)):
            for j in range(i + 1, len(Y)):
                cand.append(Y[i] + Y[j]); cand.append(Y[i] - Y[j])
        killed = False; bl = 1e9
        for y in cand:
            l1 = 9 * np.abs(y).sum()
            bl = min(bl, 2 * l1)
            if 2 * l1 >= 1: continue
            f = (y @ a) % 1.0
            # need an integer k with |f + k| <= l1  ->  f <= l1 or f >= 1 - l1
            if l1 < f < 1 - l1: killed = True; break
        got += 1; dead += killed; best.append(bl)
    print(f"L={L} d={d} vars={n}: windows kill {dead}/{got}; shortest window span {np.median(best):.3f} (needs < 1)", flush=True)
