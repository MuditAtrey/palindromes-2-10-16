# Firing point of the exact integer-equation certificate vs length: random dives through the
# 2-adic-pruned tree to depth d, then certify.  Prints, per depth, fraction proved dead and enumeration cost.
import sys, random, time
from p_int import certify, hexnode, decstate, node_equation
def ok_child(L, g):
    X, hi, inner, h = hexnode(L, g); st = decstate(X, hi)
    if st is None: return False
    D, s = st; m = len(s); R = sum(int(s[i]) * 10 ** i for i in range(m)); k = min(m, 4 * len(g))
    return (R - X) % (1 << k) == 0
def dive(L, d, rnd, tries=2000):
    for _ in range(tries):
        g = []
        while len(g) < d:
            ch = [e for e in ((9, 15) if not g else (0, 6, 9, 15)) if ok_child(L, g + [e])]
            if not ch: break
            g.append(rnd.choice(ch))
        if len(g) == d: return g
    return None
if __name__ == "__main__":
    L = int(sys.argv[1]); ds = [int(x) for x in sys.argv[2].split(',')]; S = int(sys.argv[3]); lim = int(sys.argv[4])
    rnd = random.Random(L)
    for d in ds:
        dead = bud = 0; tot = 0; t0 = time.time(); nv = 0
        for _ in range(S):
            g = dive(L, d, rnd)
            if g is None: continue
            coefs, ranges, T, m, D = node_equation(L, g); nv = len(coefs)
            ok, nn, real = certify(coefs, ranges, T, lim)
            assert not real
            tot += nn
            if ok: dead += 1
            else: bud += 1
        n = dead + bud
        print(f"L={L} B={4*L} d={d:2d} vars={nv:3d}: dead {dead}/{n}, budget-hit {bud}, mean enum nodes {tot/max(n,1):9.0f} ({(time.time()-t0)/max(n,1):.2f} s/node)", flush=True)
    