"""Ball-bounded vs box-exact (LP-tight) enumeration of the inner lattice coset.

Question: the engine's exact test enumerates the ball that contains the digit box. How many enumeration nodes would an
enumeration need that prunes exactly against the box (a partial assignment is kept only if the real relaxation can still reach
the box)? If that count stays small at higher densities, the exact test could fire much earlier in the tree.

Inner variables: n digit positions, alphabet {0,2,3,5} (hex digits {0,6,9,F} = 3*{0,2,3,5}), weights w_j mod 5^K.
Usage: python3 box_enum.py LLLB n K targets [mode]   (mode 1 = random weights)
"""
import sys, subprocess, math, random
import numpy as np
from scipy.optimize import linprog
LLLB, n, K, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
mode = sys.argv[5] if len(sys.argv) > 5 else "0"
out = subprocess.run([LLLB, str(n), str(K), mode], capture_output=True, text=True).stdout.split("\n")
w = [int(x) for x in out[0].split()]; Bm = np.array([[int(x) for x in out[1 + i].split()] for i in range(n)], dtype=float)
M = 5 ** K; ALPH = {0, 2, 3, 5}; LO, HI = 0.0, 5.0
# Gram-Schmidt
Bs = np.zeros_like(Bm); mu = np.zeros((n, n)); n2 = np.zeros(n)
for i in range(n):
    v = Bm[i].copy()
    for j in range(i):
        mu[i, j] = Bm[i] @ Bs[j] / n2[j]; v -= mu[i, j] * Bs[j]
    Bs[i] = v; n2[i] = v @ v
R2 = n * 2.5 ** 2
Cm = np.linalg.inv(Bm.T)
def ball_nodes(x0):
    # lattice points y with |x0 + y - c| <= R: target vector v = c - x0 in basis coordinates
    c = np.full(n, 2.5); v = c - x0
    coef = np.zeros(n)                       # coefficients of v along Gram-Schmidt vectors
    tv = v.copy()
    for i in range(n - 1, -1, -1): pass
    # project: v = sum_i a_i b_i (real); solve
    a = np.linalg.solve(Bm.T, v)
    cnt = 0; z = np.zeros(n)
    def rec(i, rho):
        nonlocal cnt
        ci = a[i] + sum((a[j] - z[j]) * mu[j, i] for j in range(i + 1, n))
        r = math.sqrt(max(R2 - rho, 0) / n2[i])
        for zi in range(math.ceil(ci - r), math.floor(ci + r) + 1):
            cnt += 1
            d = (zi - ci) ** 2 * n2[i]
            if rho + d > R2: continue
            z[i] = zi
            if i > 0: rec(i - 1, rho + d)
        z[i] = 0
    rec(n - 1, 0.0)
    return cnt
def box_nodes(x0):
    cnt = 0; hits = 0; z = [0] * n
    def rng(k):
        # min/max of z_k = C_k (x - x0) over x in the box with z_i = C_i (x - x0) fixed for i > k (stable: x is bounded)
        if k == n - 1: Aeq = None; beq = None
        else:
            Aeq = Cm[k + 1:]; beq = np.array([z[i] for i in range(k + 1, n)], dtype=float) + Cm[k + 1:] @ x0
        res = []
        for sgn in (1, -1):
            r = linprog(sgn * Cm[k], A_eq=Aeq, b_eq=beq, bounds=[(LO, HI)] * n, method="highs")
            if r.status == 2: return None
            if r.status != 0: raise SystemExit(f"LP failed: {r.message}")
            res.append(sgn * r.fun - Cm[k] @ x0)
        return res[0], res[1]
    def rec(k):
        nonlocal cnt, hits
        r = rng(k)
        if r is None: return
        for zk in range(math.ceil(r[0] - 1e-9), math.floor(r[1] + 1e-9) + 1):
            cnt += 1; z[k] = zk
            if k == 0:
                x = x0 + sum(z[i] * Bm[i] for i in range(n))
                if all(round(v) in ALPH and abs(v - round(v)) < 1e-6 for v in x): hits += 1
            else: rec(k - 1)
        z[k] = 0
    rec(n - 1)
    return cnt, hits
rnd = random.Random(5)
dens = 2 * n / (K * math.log2(5))
bn = []; xn = []; hs = 0
for _ in range(T):
    t = rnd.randrange(M); x0 = np.zeros(n); x0[0] = t
    bn.append(ball_nodes(x0)); c, h = box_nodes(x0); xn.append(c); hs += h
print(f"n={n:2d} K={K:2d} density {dens:.2f}{' (random weights)' if mode=='1' else ''}: ball enum nodes mean {np.mean(bn):9.1f}  box-exact (LP) nodes mean {np.mean(xn):8.1f}  ratio {np.mean(bn)/max(np.mean(xn),1e-9):7.1f}x  hits {hs}/{T}")
