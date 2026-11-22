"""Cheap box bound by weak duality: per enumeration level k, a few precomputed multipliers mu give sound bounds
  min z_k >= mu.b + 5 * sum(min(0, v)),   max z_k <= mu.b + 5 * sum(max(0, v)),   v = C_k - mu A_eq,
where A_eq fixes the already chosen coordinates and b their values. Multipliers are learned from LP duals on training targets.
Compares enumeration nodes: ball (Gram-Schmidt) / learned-dual box bound (m multipliers per level, plus the least-squares one)
/ exact LP. Usage: python3 dual_bound.py LLLB n K targets train m
"""
import sys, subprocess, math, random
import numpy as np
from scipy.optimize import linprog
LLLB, n, K, T, TR, MM = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
out = subprocess.run([LLLB, str(n), str(K), "0"], capture_output=True, text=True).stdout.split("\n")
Bm = np.array([[int(x) for x in out[1 + i].split()] for i in range(n)], dtype=float)
M = 5 ** K; LO, HI = 0.0, 5.0; ALPH = {0, 2, 3, 5}
Cm = np.linalg.inv(Bm.T)
def lp(k, z, x0):
    Aeq = Cm[k + 1:] if k < n - 1 else None
    beq = (np.array(z[k + 1:], dtype=float) + Cm[k + 1:] @ x0) if k < n - 1 else None
    out = []
    for sgn in (1, -1):
        r = linprog(sgn * Cm[k], A_eq=Aeq, b_eq=beq, bounds=[(LO, HI)] * n, method="highs")
        if r.status == 2: return None
        mu = r.eqlin.marginals if k < n - 1 else np.zeros(0)
        out.append((sgn * r.fun - Cm[k] @ x0, sgn * np.array(mu)))
    return out
# multipliers per level: least-squares projection (always) + learned LP duals
mus = {k: [] for k in range(n)}
for k in range(n - 1):
    A = Cm[k + 1:]
    mu_ls = np.linalg.lstsq(A.T, Cm[k], rcond=None)[0]
    mus[k].append(mu_ls)
mus[n - 1].append(np.zeros(0))
def learn(x0):
    z = [0] * n
    def rec(k, depth_budget):
        r = lp(k, z, x0)
        if r is None: return
        for (val, mu) in r:
            if k < n - 1 and len(mus[k]) < MM + 1 and not any(np.allclose(mu, m2, atol=1e-7) for m2 in mus[k]): mus[k].append(mu)
        lo, hi = r[0][0], r[1][0]
        for zk in range(math.ceil(lo - 1e-9), math.floor(hi + 1e-9) + 1):
            z[k] = zk
            if k > 0 and depth_budget[0] > 0: depth_budget[0] -= 1; rec(k - 1, depth_budget)
        z[k] = 0
    rec(n - 1, [400])
rnd = random.Random(11)
for _ in range(TR):
    x0 = np.zeros(n); x0[0] = rnd.randrange(M); learn(x0)
pre = {}
for k in range(n):
    lst = []
    for mu in mus[k]:
        v = Cm[k] - (mu @ Cm[k + 1:] if k < n - 1 else 0)
        lst.append((mu, HI * np.minimum(v, 0).sum() + LO * np.maximum(v, 0).sum(), HI * np.maximum(v, 0).sum() + LO * np.minimum(v, 0).sum()))
    pre[k] = lst
def dual_nodes(x0):
    cnt = 0; z = np.zeros(n); hits = 0
    Cx0 = Cm @ x0
    def rec(k):
        nonlocal cnt, hits
        b = z[k + 1:] + Cx0[k + 1:]
        lo, hi = -1e18, 1e18
        for mu, cl, ch in pre[k]:
            s = mu @ b if k < n - 1 else 0.0
            lo = max(lo, s + cl - Cx0[k]); hi = min(hi, s + ch - Cx0[k])
        for zk in range(math.ceil(lo - 1e-7), math.floor(hi + 1e-7) + 1):
            cnt += 1; z[k] = zk
            if k == 0:
                x = x0 + Bm.T @ z
                if all(abs(v - round(v)) < 1e-6 and round(v) in ALPH for v in x): hits += 1
            else: rec(k - 1)
        z[k] = 0
    rec(n - 1)
    return cnt, hits
def ball_nodes(x0, shape=False):
    """Fincke-Pohst in the ball around the box centre; shape=True adds the engine's box-shape test
    (the projection p of x - c onto span(b*_i..) of a box point satisfies |p|^2 <= sum 2.5 |p_k|)."""
    Bs = np.zeros_like(Bm); mu = np.zeros((n, n)); n2 = np.zeros(n)
    for i in range(n):
        v = Bm[i].copy()
        for j in range(i): mu[i, j] = Bm[i] @ Bs[j] / n2[j]; v -= mu[i, j] * Bs[j]
        Bs[i] = v; n2[i] = v @ v
    a = np.linalg.solve(Bm.T, np.full(n, 2.5) - x0); R2 = n * 6.25; cnt = 0; z = np.zeros(n)
    def rec(i, rho, pu):
        nonlocal cnt
        ci = a[i] + sum((a[j] - z[j]) * mu[j, i] for j in range(i + 1, n)); r = math.sqrt(max(R2 - rho, 0) / n2[i])
        for zi in range(math.ceil(ci - r), math.floor(ci + r) + 1):
            cnt += 1; d = (zi - ci) ** 2 * n2[i]
            if rho + d > R2: continue
            pv = pu + (zi - ci) * Bs[i]
            if shape and pv @ pv > 2.5 * np.abs(pv).sum() * (1 + 1e-9) + 1e-6: continue
            z[i] = zi
            if i > 0: rec(i - 1, rho + d, pv)
        z[i] = 0
    rec(n - 1, 0.0, np.zeros(n)); return cnt
rnd = random.Random(99); dn = []; bn = []; en = []; hs = 0
for _ in range(T):
    x0 = np.zeros(n); x0[0] = rnd.randrange(M)
    c, h = dual_nodes(x0); dn.append(c); hs += h; bn.append(ball_nodes(x0)); en.append(ball_nodes(x0, True))
avgm = np.mean([len(pre[k]) for k in range(n)])
print(f"n={n} K={K} density {2*n/(K*math.log2(5)):.2f}: ball {np.mean(bn):10.1f}  engine-style ball+shape {np.mean(en):9.1f}  learned-dual box bound ({avgm:.1f} multipliers/level) {np.mean(dn):9.1f}  vs engine-style {np.mean(en)/np.mean(dn):6.1f}x  hits {hs}/{T}")
