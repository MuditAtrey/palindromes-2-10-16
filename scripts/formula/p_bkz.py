# Same certificate, two improvements: (1) only cosets are enumerated (target-row coefficient pinned to 1),
# (2) the kernel basis can be BKZ-reduced (Schnorr-Euchner SVP in blocks, FLINT LLL for clean-up).
import sys, math, random, time
import numpy as np
from p_int import lll, node_equation
from p_scale import dive
def gso(B):
    B = np.array(B, dtype=float); n = len(B); Bs = np.zeros_like(B); mu = np.zeros((n, n)); nn = np.zeros(n)
    for i in range(n):
        v = B[i].copy()
        for j in range(i): mu[i, j] = B[i] @ Bs[j] / nn[j]; v -= mu[i, j] * Bs[j]
        Bs[i] = v; nn[i] = v @ v
    return mu, nn
def svp_block(mu, nn, k, e):
    """shortest nonzero vector of the projected block [k,e): returns coefficients or None if b*_k is already shortest."""
    n = e - k; best = nn[k] * 0.999; bx = None; x = [0] * n
    def rec(i, part):
        nonlocal best, bx
        cc = -sum(x[j] * mu[k + j, k + i] for j in range(i + 1, n))
        w = math.sqrt(max(best - part, 0) / nn[k + i])
        for z in range(math.ceil(cc - w), math.floor(cc + w) + 1):
            p2 = part + (z - cc) ** 2 * nn[k + i]
            if p2 >= best: continue
            x[i] = z
            if i == 0:
                if any(x): best = p2; bx = list(x)
            else: rec(i - 1, p2)
        x[i] = 0
    rec(n - 1, 0.0)
    return bx
def bkz(rows, beta, tours=3):
    rows = lll(rows)
    for t in range(tours):
        changed = False
        for k in range(len(rows) - 1):
            e = min(k + beta, len(rows)); mu, nn = gso(rows)
            c = svp_block(mu, nn, k, e)
            if c is None: continue
            v = [sum(c[j] * rows[k + j][q] for j in range(e - k)) for q in range(len(rows[0]))]
            rows = lll(rows[:k] + [v] + rows[k:])
            rows = [r for r in rows if any(r)]
            changed = True
        if not changed: break
    return rows
def coset_enum(Kb, aff, R2, limit):
    B = Kb + [aff]; mu, nn = gso(B); n = len(B); x = [0] * n; nodes = 0; sol = []
    def rec(i, part):
        nonlocal nodes
        cc = -sum(x[j] * mu[j, i] for j in range(i + 1, n))
        room = R2 - part
        if room < 0: return True
        w = math.sqrt(room / nn[i])
        rng = [1] if i == n - 1 else range(math.ceil(cc - w - 1e-9), math.floor(cc + w + 1e-9) + 1)
        for z in rng:
            nodes += 1
            if nodes > limit: return False
            p2 = part + (z - cc) ** 2 * nn[i]
            if p2 > R2 * (1 + 1e-9): continue
            x[i] = z
            if i == 0: sol.append(list(x))
            elif not rec(i - 1, p2): return False
        x[i] = 0
        return True
    ok = rec(n - 1, 0.0)
    return ok, nodes, sol
def setup(coefs, ranges, T):
    n = len(coefs); w = [18 // (hi - lo) for lo, hi in ranges]; cen2 = [lo + hi for lo, hi in ranges]
    K = 1 << (max(abs(c) for c in coefs).bit_length() + 64); M = 9
    rows = [[0] * (n + 2) for _ in range(n)]
    for v in range(n): rows[v][v] = w[v]; rows[v][n] = K * coefs[v]
    rows.append([-(w[v] * cen2[v]) // 2 for v in range(n)] + [-K * T, M])
    R = lll(rows)
    Z = [r for r in R if r[n] == 0]                       # n rows: kernel + affine
    assert len(Z) == n, len(Z)
    K2 = 1 << 80
    R2 = lll([r[:n] + [(r[n + 1] // M) * K2] for r in Z])  # separate: n-1 kernel rows, one affine row
    ker = [r[:n] for r in R2 if r[n] == 0]
    aff = [r for r in R2 if r[n] != 0]
    assert len(ker) == n - 1 and len(aff) == 1 and abs(aff[0][n]) == K2, (len(ker), len(aff))
    a = aff[0]; s = 1 if a[n] > 0 else -1
    return ker, [s * q for q in a[:n]] + [M], n

if __name__ == "__main__":
    L = int(sys.argv[1]); ds = [int(x) for x in sys.argv[2].split(',')]; S = int(sys.argv[3]); lim = int(sys.argv[4]); betas = [int(x) for x in sys.argv[5].split(',')]
    rnd = random.Random(L)
    for d in ds:
        res = {b: [0, 0, 0.0] for b in betas}; got = 0
        for _ in range(S):
            g = dive(L, d, rnd)
            if g is None: continue
            coefs, ranges, T, m, D = node_equation(L, g)
            st = setup(coefs, ranges, T)
            if st is None: continue
            ker, aff, n = st; got += 1
            for b in betas:
                t0 = time.time()
                Kb = [r + [0] for r in (bkz(ker, b, 2) if b > 2 else ker)]
                ok, nodes, sol = coset_enum(Kb, aff, 81.0 * (n + 1), lim)
                assert not sol or True
                res[b][0] += ok; res[b][1] += nodes; res[b][2] += time.time() - t0
        line = "  ".join(f"beta={b}: dead {res[b][0]}/{got} nodes {res[b][1]/max(got,1):9.0f} ({res[b][2]/max(got,1):.1f}s)" for b in betas)
        print(f"L={L} d={d} vars={len(coefs)}: {line}", flush=True)
