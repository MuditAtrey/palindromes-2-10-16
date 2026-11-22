"""Sibling problems: numbers palindromic in base 2, base 16 and base b, for many even b.

Purpose: test the random model where it can fail. Every (2, 16, b) has about 2^(-B/4) expected solutions per length,
so under the random model no sibling should have solutions beyond B ~ 40-60. A sibling with large solutions would be a
mechanism worth studying (and maybe transferring to b = 10).

Exhaustive for bit lengths B = 4L (hex digits from {0,6,9,F}, mirrored), using only provably sound pruning:
the low m base-b digits are the mirror of the top m (fixed by the outer hex digits), so N mod b^m is known,
hence N mod 2^(v*m) with v = v2(b); the outer hex digits fix N mod 2^(4d). Disagreement = no completion.
Odd b has no 2-adic link and is not handled here.

Usage: python3 sibling_search.py LMAX PROCS b1,b2,...  [--check]
  --check : compare against brute force for L <= 12 (completeness test of the tree).
"""
import sys, bisect, itertools, math, time
from multiprocessing import Pool
sys.set_int_max_str_digits(0)

def v2(b):
    v = 0
    while b % 2 == 0: b //= 2; v += 1
    return v

class Base:
    def __init__(self, b, maxbits):
        self.b = b; self.v = v2(b)
        self.P = [1]
        while self.P[-1] < (1 << (maxbits + 8)): self.P.append(self.P[-1] * b)
    def ndig(self, x):                     # number of base-b digits of x >= 1
        return bisect.bisect_right(self.P, x)
    def digits(self, x):                   # most significant first
        out = []
        while x: out.append(x % self.b); x //= self.b
        return out[::-1] or [0]
    def is_pal(self, x):
        d = self.digits(x); return d == d[::-1]
    def common_prefix(self, lo, hi):       # largest p with lo // b^(D-p) == hi // b^(D-p); requires same digit count
        D = self.ndig(lo); a, z = 0, D
        while a < z:                       # find largest p in [0, D] with equal prefixes
            mid = (a + z + 1) // 2
            if lo // self.P[D - mid] == hi // self.P[D - mid]: a = mid
            else: z = mid - 1
        return D, a

def search_length(args):
    b, L = args
    B = 4 * L; base = Base(b, B + 4)
    half = (L + 1) // 2
    h = [(16 ** j + 16 ** (L - 1 - j)) if j != L - 1 - j else 16 ** j for j in range(half)]
    suf = [0] * (half + 1)
    for j in range(half - 1, -1, -1): suf[j] = suf[j + 1] + 15 * h[j]
    hits = []; nodes = 0
    def dfs(d, X):
        nonlocal nodes
        if d == half:
            if base.is_pal(X) and format(X, "b") == format(X, "b")[::-1] and format(X, "x") == format(X, "x")[::-1]:
                hits.append(X)
            return
        for e in ((9, 15) if d == 0 else (0, 6, 9, 15)):
            nodes += 1
            X2 = X + e * h[d]
            lo, hi = X2, X2 + suf[d + 1]
            if base.ndig(lo) == base.ndig(hi):
                D, p = base.common_prefix(lo, hi)
                m = min(p, D // 2)
                if m:
                    top = lo // base.P[D - m]            # the top m digits as a number
                    ds = base.digits(top); ds = [0] * (m - len(ds)) + ds
                    R = 0
                    for i in range(m - 1, -1, -1): R = R * b + ds[i]   # low digit i = top digit i (mirror)
                    k = min(base.v * m, 4 * (d + 1))
                    if (R - X2) % (1 << k): continue                    # 2-adic disagreement: no completion
            dfs(d + 1, X2)
    dfs(0, 0)
    return b, L, nodes, hits

def search_length_c(args):
    """the same search in C (sibling.c), ~100x faster; returns the same tuple as search_length."""
    import subprocess, os
    b, L = args
    exe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sibling")
    out = subprocess.run([exe, str(b), str(L)], capture_output=True, text=True, check=True).stdout.split()
    hits = [int(x) for x in out[4:]]
    hits += [0] * (int(out[3]) - len(hits))   # count beyond the 64 listed (controls only)
    return int(out[0]), int(out[1]), int(out[2]), hits

def brute(b, L):
    base = Base(b, 4 * L + 4); out = []
    half = (L + 1) // 2
    for first in (9, 15):
        for rest in itertools.product((0, 6, 9, 15), repeat=half - 1):
            e = (first,) + rest; ds = list(e) + list(e[: L // 2][::-1]); N = 0
            for x in ds: N = N * 16 + x
            if base.is_pal(N): out.append(N)
    return sorted(out)

def expected(b, L):
    """random model: candidates times P(base-b palindrome) = b^-(floor(D/2)), D the typical digit count."""
    half = (L + 1) // 2; cand = 2 * 4 ** (half - 1)
    D = int(math.floor((4 * L - 0.5) * math.log(2) / math.log(b))) + 1
    return cand * b ** (-(D // 2))

if __name__ == "__main__":
    LMAX = int(sys.argv[1]); PROCS = int(sys.argv[2]); bases = [int(x) for x in sys.argv[3].split(",")]
    if "--check" in sys.argv:
        bad = 0
        for b in bases:
            for L in range(1, 13):
                _, _, _, hits = search_length((b, L)); br = brute(b, L)
                if sorted(hits) != br: bad += 1; print(f"MISMATCH b={b} L={L}: tree {sorted(hits)[:5]} brute {br[:5]}")
        print(f"completeness check vs brute force, L<=12, bases {bases}: {'all match' if not bad else str(bad) + ' mismatches'}")
        sys.exit(0)
    t0 = time.time()
    # per-base length limit: the tree grows about 2^(g/2) per hex digit, g = 2 - v2(b)*log_b(16) bits per level
    def lmax_for(b):
        g = 2 - v2(b) * math.log(16) / math.log(b)
        return LMAX if g <= 0 else min(LMAX, int((48 if "--c" in sys.argv else 40) / g))
    LM = {b: lmax_for(b) for b in bases}
    jobs = sorted([(b, L) for b in bases for L in range(1, LM[b] + 1)], key=lambda x: -x[1])
    res = {}
    with Pool(PROCS) as pool:
        fn = search_length_c if "--c" in sys.argv else search_length
        for b, L, nodes, hits in pool.imap_unordered(fn, jobs):
            res[(b, L)] = (nodes, hits)
    print(f"(2, 16, b) palindromes with B = 4L <= {4*LMAX} bits, exhaustive in that class ({time.time()-t0:.0f} s)")
    print(f"{'b':>4} {'max B':>6} {'nodes':>11} {'solutions':>9} {'largest (bits)':>14} {'expected L>=8':>14} {'observed L>=8':>14}  examples")
    for b in bases:
        Lm = LM[b]
        nodes = sum(res[(b, L)][0] for L in range(1, Lm + 1))
        hits = sorted(h for L in range(1, Lm + 1) for h in res[(b, L)][1])
        big = [h for h in hits if h.bit_length() > 28]
        exp8 = sum(expected(b, L) for L in range(8, Lm + 1)); obs8 = sum(len(res[(b, L)][1]) for L in range(8, Lm + 1))
        largest = max((h.bit_length() for h in hits), default=0)
        print(f"{b:>4} {4*Lm:>6} {nodes:>11,d} {len(hits):>9} {largest:>14} {exp8:>14.4f} {obs8:>14}  {hits[:4]}{' ... big: ' + str([format(x,'X') for x in big[:3]]) if big else ''}")
