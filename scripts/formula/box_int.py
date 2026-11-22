"""Box-exact enumeration of the exact integer equation at a tree node (all information: 2-adic, 5-adic and size).

  sum_v c_v x_v = T,  x in a box (inner hex bits a,b in {0,1}, inner decimal digits in 0..9).
Solutions x = x_p + sum_i z_i k_i over an LLL-reduced kernel basis (computed in coordinates scaled to equal half-widths).
Enumeration over z from the last basis vector down; at every level the range of z_k comes from two LPs over the real box with
the equation and the already fixed z's as equalities (exact pruning against the box, no ball). Counts LP-enumeration nodes and
compares with the ball enumeration of p_int.certify.
Usage: python3 box_int.py L d1,d2,.. samples [lp_node_cap]
"""
import sys, math, random, time
import numpy as np
from scipy.optimize import linprog
from p_int import lll, node_equation, certify
from p_scale import dive

def kernel_and_particular(coefs, ranges, T):
    n = len(coefs)
    w = [18 // (hi - lo) for lo, hi in ranges]
    K = 1 << (max(abs(c) for c in coefs).bit_length() + 64)
    rows = []
    for v in range(n):
        r = [0] * (n + 2); r[v] = w[v]; r[n] = K * coefs[v]; rows.append(r)
    rows.append([0] * n + [-K * T, 1])
    Rd = lll(rows)
    Z = [list(r) for r in Rd if r[n] == 0]           # spans {(w x, t) : c.x = t T}: the kernel plus the coset generator
    if len(Z) != n: return None
    # integer elimination on the t column: one row keeps gcd(t) (must be 1), every other row gets t = 0
    while sum(1 for r in Z if r[n + 1] != 0) > 1:
        nz = [r for r in Z if r[n + 1] != 0]; nz.sort(key=lambda r: abs(r[n + 1])); a = nz[0]
        for r in nz[1:]:
            q = r[n + 1] // a[n + 1]
            for k in range(n + 2): r[k] -= q * a[k]
    piv = [r for r in Z if r[n + 1] != 0]
    if len(piv) != 1 or abs(piv[0][n + 1]) != 1: return None
    p = piv[0]; s = p[n + 1]
    ker = [r[:n] + [0] for r in Z if r[n + 1] == 0]
    ker = [r[:n] for r in lll([r[:n] for r in ker])]  # re-reduce the complete kernel
    if any(p[v] % w[v] for v in range(n)): return None
    xp = [s * p[v] // w[v] for v in range(n)]
    kx = []
    for r in ker:
        if any(r[v] % w[v] for v in range(n)): return None
        kx.append([r[v] // w[v] for v in range(n)])
    return np.array(kx, dtype=float), np.array(xp, dtype=float), w

def box_count(coefs, ranges, T, cap):
    kp = kernel_and_particular(coefs, ranges, T)
    if kp is None: return None
    Km, xp, w = kp
    m, n = Km.shape
    lo = np.array([a for a, b in ranges], float); hi = np.array([b for a, b in ranges], float)
    G = Km @ Km.T; C = np.linalg.solve(G, Km)        # z = C (x - xp) on the affine space
    cvec = np.array(coefs, float); cnorm = cvec / np.abs(cvec).max()
    Tn = float(T) / np.abs(cvec).max()
    cnt = 0; z = [0] * m; sols = 0; over = False
    def rng(k):
        rows = [cnorm]; rhs = [Tn]
        for i in range(k + 1, m): rows.append(C[i]); rhs.append(z[i] + C[i] @ xp)
        Aeq = np.array(rows); beq = np.array(rhs)
        out = []
        for sgn in (1, -1):
            r = linprog(sgn * C[k], A_eq=Aeq, b_eq=beq, bounds=list(zip(lo, hi)), method="highs")
            if r.status == 2: return None
            if r.status != 0: raise RuntimeError(r.message)
            out.append(sgn * r.fun - C[k] @ xp)
        return out
    def rec(k):
        nonlocal cnt, sols, over
        if over: return
        r = rng(k)
        if r is None: return
        for zk in range(math.ceil(r[0] - 1e-7), math.floor(r[1] + 1e-7) + 1):
            cnt += 1
            if cnt > cap: over = True; return
            z[k] = zk
            if k == 0:
                x = xp + np.array(z, float) @ Km
                if np.all(np.abs(x - np.round(x)) < 1e-6) and np.all(x >= lo - 1e-9) and np.all(x <= hi + 1e-9): sols += 1
            else: rec(k - 1)
        z[k] = 0
    rec(m - 1)
    return cnt, sols, over

if __name__ == "__main__":
    L = int(sys.argv[1]); ds = [int(x) for x in sys.argv[2].split(',')]; S = int(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
    rnd = random.Random(L)
    for d in ds:
        bx = []; bl = []; dead_box = 0; n_ok = 0; t0 = time.time(); nv = 0; overs = 0
        for _ in range(S):
            g = dive(L, d, rnd)
            if g is None: continue
            coefs, ranges, T, mm, D = node_equation(L, g); nv = len(coefs)
            res = box_count(coefs, ranges, T, cap)
            if res is None: continue
            c, s, over = res
            okb, nb, real = certify(coefs, ranges, T, 300000)
            n_ok += 1; bx.append(c); bl.append(nb if okb else float('nan')); overs += over
            if not over and s == 0: dead_box += 1
        print(f"L={L} B={4*L} d={d:2d} vars={nv:3d}: box-exact nodes mean {np.mean(bx):8.1f} max {max(bx) if bx else 0} (cap hits {overs}), dead {dead_box}/{n_ok} | ball nodes mean {np.nanmean(bl) if bl else 0:10.1f}  ({(time.time()-t0)/max(n_ok,1):.1f} s/node)", flush=True)
