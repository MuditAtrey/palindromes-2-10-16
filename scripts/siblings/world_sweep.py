"""Wide sibling sweep: hunt for a mechanism that makes numbers palindromic in three bases more often than chance.

A world is (r, A, T, b): numbers palindromic in base 2, in r = 2^k and in an even base b (not a power of 2), in one bit-length class.
Under the random model every world has about 2^(-c B) solutions per length, so large solutions should not exist. Any world with
more large solutions than its expectation, or a solution far beyond where the model allows one, is a mechanism to study.

Usage:
  python3 world_sweep.py check                 completeness: tree vs brute force over all binary palindromes up to 26 bits, plus
                                               the b = r positive controls (every candidate must survive) and sibling.c agreement
  python3 world_sweep.py sweep OUT PROCS CAP   run every world, lengths increasing, until one length takes CAP seconds (or 600 bits)
  python3 world_sweep.py report OUT            observed vs expected, every large solution with its surprise
"""
import sys, os, json, math, time, subprocess, itertools
from multiprocessing import Pool
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(HERE, "world")

WORLDS = {   # name: (r, A, T, bit length as a function of M)
    "hex_even": (16, (0, 6, 9, 15), (9, 15), lambda M: 4 * M),
    "hex_odd":  (16, (0, 2, 5, 7), (5, 7), lambda M: 4 * M - 1),
    "oct":      (8, (0, 2, 5, 7), (5, 7), lambda M: 3 * M),
    "quad":     (4, (0, 3), (3,), lambda M: 2 * M),
}
BIG_BITS = 32          # "large" solution: at least this many bits (the random model is meant for large numbers)
MAX_BITS = 600

def run(world, b, M, timeout=None):
    r, A, T, _ = WORLDS[world]
    out = subprocess.run([EXE, str(r), ",".join(map(str, A)), ",".join(map(str, T)), str(b), str(M)],
                         capture_output=True, text=True, check=True, timeout=timeout).stdout.split()
    return int(out[0]), int(out[1]), [int(x) for x in out[2:]]

def is_pal(x, base):
    d = []
    while x: d.append(x % base); x //= base
    return d == d[::-1]

def expected(world, b, M):
    """random model: sum over candidates of P(a random D-digit number is a base-b palindrome) = b^-(floor(D/2)),
    candidates taken uniform in value within each top digit's range, split by base-b digit count D."""
    r, A, T, _ = WORLDS[world]
    half = (M + 1) // 2; per_top = len(A) ** (half - 1)
    S = r ** (M - 1); tot = 0.0
    for t in T:
        lo, hi = t * S, (t + 1) * S
        p = 1; D = 1
        while p * b <= lo: p *= b; D += 1          # b^(D-1) <= lo < b^D
        x = lo
        while x < hi:
            nxt = min(hi, p * b)
            tot += per_top * (nxt - x) / (hi - lo) * b ** (-(D // 2))
            x = nxt; p *= b; D += 1
    return tot

def frac_div(world, b, M, t, q):
    """exact fraction of candidates with top digit t that are divisible by q (digit DP over residues mod q)"""
    r, A, T, _ = WORLDS[world]
    half = (M + 1) // 2
    W = [(pow(r, j, q) + (pow(r, M - 1 - j, q) if j != M - 1 - j else 0)) % q for j in range(half)]
    dist = [0] * q; dist[0] = 1
    for j in range(half):
        digs = (t,) if j == 0 else A
        nd = [0] * q
        for s in range(q):
            if dist[s]:
                for e in digs: nd[(s + e * W[j]) % q] += dist[s]
        dist = nd
    return dist[0] / sum(dist)

def expected_corr(world, b, M):
    """random model corrected for the one exact divisibility link: an even-length base-b palindrome is divisible by b + 1.
    Candidates divisible by b + 1 more (or less) often than 1/(b + 1) get that factor on their even-length share."""
    r, A, T, _ = WORLDS[world]
    half = (M + 1) // 2; per_top = len(A) ** (half - 1)
    S = r ** (M - 1); tot = 0.0
    for t in T:
        corr = (b + 1) * frac_div(world, b, M, t, b + 1)
        lo, hi = t * S, (t + 1) * S
        p = 1; D = 1
        while p * b <= lo: p *= b; D += 1
        x = lo
        while x < hi:
            nxt = min(hi, p * b)
            tot += per_top * (nxt - x) / (hi - lo) * b ** (-(D // 2)) * (corr if D % 2 == 0 else 1.0)
            x = nxt; p *= b; D += 1
    return tot

def bases():
    return [b for b in range(6, 201, 2) if b & (b - 1)]

def check():
    bad = 0
    for world, (r, A, T, bits) in WORLDS.items():
        k = r.bit_length() - 1
        for M in range(1, 40):
            B = bits(M)
            if B > 26: break
            if B < 2: continue
            # all binary palindromes with exactly B bits that are also base-r palindromes: must equal the alphabet enumeration
            half = (B + 1) // 2; cands = []
            for x in range(1 << (half - 1), 1 << half):
                s = format(x, "b"); full = s + s[: B // 2][::-1]; N = int(full, 2)
                if is_pal(N, r): cands.append(N)
            for b in [6, 10, 12, 14, 20, 24, 54, 98, 130, r]:
                want = sorted(N for N in cands if is_pal(N, b))
                nodes, nh, hits = run(world, b, M)
                got = sorted(hits)
                if nh != len(want) or (len(want) <= 64 and got != want):
                    bad += 1; print(f"MISMATCH {world} b={b} M={M} B={B}: tree {nh} {got[:5]} brute {len(want)} {want[:5]}")
        # positive control at larger size: with b = r every candidate is a solution, so nothing may be pruned
        for M in (12, 15, 18):
            nodes, nh, _ = run(world, r, M)
            want = len(T) * len(A) ** ((M + 1) // 2 - 1)
            if nh != want: bad += 1; print(f"CONTROL FAIL {world} b=r M={M}: {nh} vs {want}")
    # agreement with the older sibling.c on the real class
    for b, L in [(10, 30), (6, 40), (14, 30), (54, 16), (20, 50)]:
        old = subprocess.run([os.path.join(HERE, "sibling"), str(b), str(L)], capture_output=True, text=True, check=True).stdout.split()
        new = run("hex_even", b, L)
        if int(old[2]) != new[0] or int(old[3]) != new[1]: bad += 1; print(f"DISAGREE with sibling.c b={b} L={L}: {old[2:4]} vs {new[:2]}")
    print("all checks pass" if not bad else f"{bad} failures")
    return bad

def sweep_world(args):
    world, b, cap = args
    bits = WORLDS[world][3]; rows = []
    M = 1
    while bits(M) <= MAX_BITS:
        t0 = time.time()
        try:
            nodes, nh, hits = run(world, b, M, timeout=cap * 8)
        except subprocess.TimeoutExpired:
            break                                  # not finished: not counted as searched
        dt = time.time() - t0
        rows.append({"world": world, "b": b, "M": M, "bits": bits(M), "nodes": nodes, "nhits": nh, "hits": [str(h) for h in hits],
                     "expected": expected(world, b, M), "sec": round(dt, 3)})
        if dt > cap: break
        M += 1
    return rows

def sweep(out, procs, cap):
    jobs = [(w, b, cap) for w in WORLDS for b in bases()]
    jobs.sort(key=lambda j: -(math.log(j[1]) / j[1]))     # rough: expensive worlds (large b, small v2) last
    t0 = time.time(); done = 0
    with open(out, "a") as f, Pool(procs) as pool:
        for rows in pool.imap_unordered(sweep_world, jobs):
            for row in rows: f.write(json.dumps(row) + "\n")
            f.flush(); done += 1
            if done % 20 == 0: print(f"{done}/{len(jobs)} worlds, {time.time()-t0:.0f} s", flush=True)
    print(f"done: {len(jobs)} worlds in {time.time()-t0:.0f} s")

def poisson_tail(k, lam):          # P(X >= k)
    if k <= 0: return 1.0
    return 1.0 - sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))

def report(out):
    rows = [json.loads(l) for l in open(out)]
    model = "raw" if "--raw" in sys.argv else "corrected"
    if model == "corrected":
        for r in rows:
            if r["bits"] >= BIG_BITS: r["expected"] = expected_corr(r["world"], r["b"], r["M"])
    print(f"random model: {model} ({'plain' if model == 'raw' else 'with the exact b+1 divisibility link of even-length palindromes'})")
    worlds = {}
    for r in rows: worlds.setdefault((r["world"], r["b"]), []).append(r)
    exp_big = obs_big = 0.0; bigs = []; reach = {}
    for key, rs in worlds.items():
        rs.sort(key=lambda r: r["M"]); reach[key] = rs[-1]["bits"]
        for r in rs:
            if r["bits"] < BIG_BITS: continue
            exp_big += r["expected"]; obs_big += r["nhits"]
            for h in r["hits"]:
                # surprise: expected number of solutions in this world at this length or longer, over the lengths searched
                tail = sum(x["expected"] for x in rs if x["bits"] >= r["bits"])
                bigs.append((tail, key, r["bits"], int(h)))
    print(f"worlds: {len(worlds)} ({', '.join(WORLDS)}; even b in 6..200, not powers of 2)")
    for w in WORLDS:
        rb = [reach[k] for k in reach if k[0] == w]
        print(f"  {w:9s}: {len(rb)} bases, searched up to {min(rb)}-{max(rb)} bits (median {sorted(rb)[len(rb)//2]})")
    print(f"solutions with >= {BIG_BITS} bits: observed {int(obs_big)}, random model expects {exp_big:.3f}, "
          f"P(>= observed) = {poisson_tail(int(obs_big), exp_big):.3g}")
    for tail, (w, b), bits, h in sorted(bigs):
        print(f"  {w:9s} b={b:3d}: {bits} bits  {h}  (hex {h:X}, base-b expectation at this size or larger {tail:.2e})")
    if bigs:
        # look-elsewhere: P(some world has a solution whose own tail expectation is <= the smallest one seen) under the random model
        s = min(bigs)[0]; lam = 0.0
        for key, rs in worlds.items():
            for r in rs:
                if r["bits"] >= BIG_BITS and sum(x["expected"] for x in rs if x["bits"] >= r["bits"]) <= s: lam += r["expected"]
        print(f"most surprising solution: tail {s:.2e} in its own world; chance of one at least this surprising somewhere among "
              f"{len(worlds)} worlds = {1 - math.exp(-lam):.3g}")
    for w in WORLDS:
        e = sum(r["expected"] for r in rows if r["world"] == w and r["bits"] >= BIG_BITS)
        o = sum(r["nhits"] for r in rows if r["world"] == w and r["bits"] >= BIG_BITS)
        print(f"  {w:9s}: observed {o}, expected {e:.3f}")
    small = sum(r["nhits"] for r in rows if r["bits"] < BIG_BITS)
    print(f"solutions below {BIG_BITS} bits (trivial and small, not used in the test): {small}")

if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "check": sys.exit(1 if check() else 0)
    if cmd == "sweep": sweep(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]))
    if cmd == "report": report(sys.argv[2])
