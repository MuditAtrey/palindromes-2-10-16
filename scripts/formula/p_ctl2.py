# Positive control on the coset enumeration (+ optional BKZ): known base-2/10 double palindromes must be found.
import sys, time
import p_ctl_lib as C
from p_bkz import setup, coset_enum, bkz
def decode(sol, Kb, aff, ranges, coefs, T):
    B = Kb + [aff]; n = len(coefs)
    vec = [sum(sol[i] * B[i][k] for i in range(len(B))) for k in range(n + 1)]
    w = [18 // (r[1] - r[0]) for r in ranges]; xs = []
    for v in range(n):
        num = vec[v] + (w[v] * (ranges[v][0] + ranges[v][1])) // 2
        if num % w[v]: return None
        x = num // w[v]
        if not ranges[v][0] <= x <= ranges[v][1]: return None
        xs.append(x)
    return xs if sum(c * x for c, x in zip(coefs, xs)) == T else None
beta = int(sys.argv[1]); lim = int(sys.argv[2])
good = bad = bud = 0
for N0 in C.vals[-8:]:
    B = N0.bit_length()
    for dbits in (B // 4, B // 3, int(B / 2.6)):
        e = C.eq(N0, B, dbits)
        if e is None: continue
        coefs, ranges, T, m = e
        ker, aff, n = setup(coefs, ranges, T)
        Kb = [r + [0] for r in (bkz(ker, beta, 2) if beta > 2 else ker)]
        t0 = time.time(); ok, nodes, sols = coset_enum(Kb, aff, 81.0 * (n + 1), lim)
        real = [x for x in (decode(s, Kb, aff, ranges, coefs, T) for s in sols) if x]
        own = [(N0 >> i) & 1 for i in range(dbits, (B + 1) // 2)]
        hit = any(x[:len(own)] == own for x in real)
        if hit: good += 1
        elif ok: bad += 1
        else: bud += 1
        print(f"  B={B} dbits={dbits} vars={n}: {'FOUND' if hit else ('DEAD (BAD)' if ok else 'budget')} nodes {nodes} ({time.time()-t0:.1f}s), {len(real)} box solutions", flush=True)
print(f"beta={beta}: found {good}, wrongly dead {bad}, budget {bud}")
