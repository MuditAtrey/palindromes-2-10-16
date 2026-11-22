# Exact integer-equation lattice at a tree node.
#   X + sum_j e_j h_j  =  Y + sum_i d_i t_i      (hex inner digits e_j in {0,6,9,15} = 6a+9b, decimal inner digits d_i in 0..9)
# One equation over Z.  Solutions are lattice points of a coset; complete Fincke-Pohst enumeration in the
# ball that contains the whole box proves a node dead.  Weights make every coordinate's half-width 9.
import os, sys, subprocess, math, random, time
import numpy as np
LLL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lll")

def lll(rows):
    fn = f"m{os.getpid()}.txt"; out = fn + ".o"
    r, c = len(rows), len(rows[0])
    with open(fn, "w") as f:
        f.write(f"{r} {c}  " + "  ".join(" ".join(map(str, x)) for x in rows) + "\n")
    subprocess.check_call([LLL, fn, out])
    v = list(map(int, open(out).read().split()))[2:]
    os.remove(fn); os.remove(out)
    return [v[i*c:(i+1)*c] for i in range(r)]

def enum(B, R2, limit):
    """all integer combos of rows of B (floats ok: small after LLL) with norm^2 <= R2; Schnorr-Euchner."""
    B = np.array(B, dtype=float); n = B.shape[0]
    Bs = np.zeros_like(B); mu = np.zeros((n, n)); nn = np.zeros(n)
    for i in range(n):
        v = B[i].copy()
        for j in range(i):
            mu[i, j] = B[i] @ Bs[j] / nn[j]; v -= mu[i, j] * Bs[j]
        Bs[i] = v; nn[i] = v @ v
    sols = []; nodes = 0
    x = [0] * n; c = [0.0] * n; l = [0.0] * (n + 1)
    def rec(i, part):
        nonlocal nodes
        if nodes > limit: return False
        cc = -sum(x[j] * mu[j, i] for j in range(i + 1, n))
        room = R2 - part
        if room < 0: return True
        w = math.sqrt(room / nn[i])
        for z in range(math.ceil(cc - w - 1e-9), math.floor(cc + w + 1e-9) + 1):
            nodes += 1
            if nodes > limit: return False
            p2 = part + (z - cc) ** 2 * nn[i]
            if p2 > R2 * (1 + 1e-9): continue
            x[i] = z
            if i == 0:
                if any(x): sols.append(list(x))
            elif not rec(i - 1, p2): return False
        x[i] = 0
        return True
    ok = rec(n - 1, 0.0)
    return ok, nodes, sols

def certify(coefs, ranges, T, limit=200000):
    """is there x in box with sum coefs*x == T ?  returns (complete, nodes, solutions)."""
    n = len(coefs)
    w = [18 // (hi - lo) for lo, hi in ranges]           # half-width 9 everywhere
    cen2 = [ (lo + hi) for lo, hi in ranges ]            # 2*centre
    K = 1 << (max(abs(c) for c in coefs).bit_length() + 64)
    M = 9
    rows = []
    for v in range(n):
        r = [0] * (n + 2); r[v] = w[v]; r[n] = K * coefs[v]; rows.append(r)
    r = [-(w[v] * cen2[v]) // 2 for v in range(n)] + [-K * T, M]; rows.append(r)
    for v in range(n):
        assert (w[v] * cen2[v]) % 2 == 0
    Rd = lll(rows)
    # drop rows touching the K column (they are huge); keep the rest
    Bk = [r for r in Rd if r[n] == 0]
    if len(Bk) < len(Rd) - 1:
        return True, 0, []   # should not happen
    R2 = 81.0 * (n + 1)
    ok, nodes, sols = enum([r[:n] + [r[n + 1]] for r in Bk], R2, limit)
    real = []
    for z in sols:
        vec = [sum(z[i] * Bk[i][k] for i in range(len(Bk))) for k in range(n + 2)]
        if abs(vec[n + 1]) != M or vec[n] != 0: continue
        s = 1 if vec[n + 1] == M else -1
        xs = []
        good = True
        for v in range(n):
            num = s * vec[v] + (w[v] * cen2[v]) // 2
            if num % w[v]: good = False; break
            xv = num // w[v]
            if not (ranges[v][0] <= xv <= ranges[v][1]): good = False; break
            xs.append(xv)
        if good and sum(c * x for c, x in zip(coefs, xs)) == T: real.append(xs)
    return ok, nodes, real

# ---------------- the triple problem: hex palindrome {0,6,9,F}, decimal palindrome -------------------
def hexnode(L, digs):
    d = len(digs)
    h = lambda j: (16 ** j + 16 ** (L - 1 - j)) if j != L - 1 - j else 16 ** j
    X = sum(e * h(j) for j, e in enumerate(digs))
    inner = [j for j in range(d, (L + 1) // 2)]
    hi = X + sum(15 * h(j) for j in inner)
    return X, hi, inner, h

def decstate(lo, hi):
    a, b = str(lo), str(hi)
    if len(a) != len(b): return None
    D = len(a); p = 0
    while p < D and a[p] == b[p]: p += 1
    m = min(p, D // 2)
    return D, a[:m]

def survivors(L, depth):
    """outside-in tree for B = 4L with the 2-adic prune; returns nodes at each depth."""
    out = {0: [[]]}
    for d in range(1, depth + 1):
        nxt = []
        for digs in out[d - 1]:
            for e in ((9, 15) if d == 1 else (0, 6, 9, 15)):
                g = digs + [e]
                X, hi, inner, h = hexnode(L, g)
                st = decstate(X, hi)
                if st is None: continue
                D, s = st; m = len(s)
                R = sum(int(s[i]) * 10 ** i for i in range(m))
                k = min(m, 4 * d)
                if (R - X) % (1 << k): continue
                nxt.append(g)
        out[d] = nxt
    return out

def node_equation(L, digs):
    X, hi, inner, h = hexnode(L, digs)
    D, s = decstate(X, hi); m = len(s)
    Y = sum(int(s[i]) * (10 ** (D - 1 - i) + 10 ** i) for i in range(m))
    coefs, ranges = [], []
    for j in inner:
        coefs += [6 * h(j), 9 * h(j)]; ranges += [(0, 1), (0, 1)]
    for i in range(m, (D + 1) // 2):
        t = (10 ** i + 10 ** (D - 1 - i)) if i != D - 1 - i else 10 ** i
        coefs.append(-t); ranges.append((0, 9))
    return coefs, ranges, Y - X, m, D

if __name__ == "__main__":
    L = int(sys.argv[1]); dmax = int(sys.argv[2]); S = int(sys.argv[3]); lim = int(sys.argv[4]) if len(sys.argv) > 4 else 200000
    t0 = time.time(); tree = survivors(L, dmax)
    print(f"L={L} B={4*L}: tree sizes", [len(tree[d]) for d in range(dmax + 1)], f"({time.time()-t0:.1f}s)")
    rnd = random.Random(1)
    for d in range(1, dmax + 1):
        nodes = tree[d]
        if not nodes: break
        samp = rnd.sample(nodes, min(S, len(nodes)))
        dead = 0; tot = 0; mx = 0; nvars = 0; t1 = time.time(); found = 0
        for g in samp:
            coefs, ranges, T, m, D = node_equation(L, g)
            ok, nn, real = certify(coefs, ranges, T, lim)
            tot += nn; mx = max(mx, nn); nvars = len(coefs)
            if ok and not real: dead += 1
            if real: found += 1
        print(f"  d={d:2d} vars={nvars:3d} m={m:2d}: dead {dead}/{len(samp)}  enum nodes mean {tot/len(samp):9.0f} max {mx:8d}  found {found}  ({(time.time()-t1)/len(samp)*1000:.0f} ms/node)", flush=True)
