"""Positive control for box-exact enumeration of the exact equation: base-2/base-10 double palindromes (OEIS A007632).
Along each number's own outer bits the enumeration must find it; along a path with one outer bit pair flipped it should be dead."""
import sys, random
sys.argv = [sys.argv[0], "1"]
from p_ctl import eq, vals
from box_int import box_count
big = [v for v in vals if v.bit_length() >= 60][:8]
found = miss = dead = alive = 0
for N0 in big:
    B = N0.bit_length()
    for dbits in (B // 8, B // 6, B // 4):
        e = eq(N0, B, dbits)
        if e is None: continue
        coefs, ranges, T, m = e
        r = box_count(coefs, ranges, T, 200000)
        if r is None: continue
        c, s, over = r
        if s >= 1: found += 1
        else: miss += 1; print(f"  MISSED: B={B} dbits={dbits} nodes {c} over={over}")
        b = [(N0 >> i) & 1 for i in range(B)]; i = random.Random(B * 7 + dbits).randrange(1, dbits); b[i] ^= 1; b[B - 1 - i] ^= 1
        e2 = eq(N0, B, dbits, b)
        if e2:
            r2 = box_count(*e2[:3], 200000)
            if r2: dead += (r2[1] == 0 and not r2[2]); alive += (r2[1] > 0)
        print(f"  B={B:3d} dbits={dbits:3d} vars={len(coefs):3d}: own path {'FOUND' if s else 'MISSED'} ({c} nodes)", flush=True)
print(f"control: {found} found, {miss} missed; mutated paths: {dead} proved dead, {alive} with a solution")
