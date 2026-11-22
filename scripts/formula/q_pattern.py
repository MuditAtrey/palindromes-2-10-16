# Is there a pattern in the known base-2/base-10 double palindromes (OEIS A007632)?
# Every test compares the real list with a null model of the same lengths.
import random, math
F = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../engine/data/oeis/b007632.txt")
vals = sorted({int(l.split()[1]) for l in open(F) if l.strip() and not l.startswith('#')})
vals = [v for v in vals if v >= 10**6]
print(len(vals), "terms >= 10^6, largest has", len(str(vals[-1])), "digits")
rnd = random.Random(7)
def mid_digits(v):
    s = str(v); D = len(s); a, b = D // 3, D - D // 3; return s[a:(D + 1) // 2] if (D + 1) // 2 > a else ""
# 1. position inside the length interval [10^(D-1), 10^D): uniform on a log scale?
pos = [ (math.log10(v) % 1.0) for v in vals ]
pos.sort(); n = len(pos); ks = max(max(abs((i + 1) / n - p), abs(i / n - p)) for i, p in enumerate(pos))
expect_log = [ (math.log10(rnd.uniform(10**(len(str(v))-1), 10**len(str(v)))) % 1.0) for v in vals ]
print(f"1. log-position: KS distance from uniform {ks:.3f} (5% critical ~{1.36/math.sqrt(n):.3f}); note: random numbers are NOT log-uniform either")
# 2. middle decimal digits: are they uniform 0..9?
cnt = [0] * 10
for v in vals:
    for ch in mid_digits(v): cnt[int(ch)] += 1
tot = sum(cnt); chi = sum((c - tot / 10) ** 2 / (tot / 10) for c in cnt)
print(f"2. middle decimal digits: counts {cnt}, chi2 = {chi:.1f} on 9 dof (5% critical 16.9)")
# 3. middle bits: balanced?
ones = zeros = 0
for v in vals:
    b = bin(v)[2:]; B = len(b); seg = b[B // 3:(B + 1) // 2]; ones += seg.count('1'); zeros += seg.count('0')
z = (ones - (ones + zeros) / 2) / math.sqrt((ones + zeros) / 4)
print(f"3. middle bits: {ones} ones vs {zeros} zeros, z = {z:.2f}")
# 4. does one term predict the next?  agreement of leading digits / middle bits between consecutive terms vs random pairs of same lengths
def agree(a, b):
    sa, sb = str(a), str(b); k = 0
    while k < min(len(sa), len(sb)) and sa[k] == sb[k]: k += 1
    return k
cons = [agree(vals[i], vals[i + 1]) for i in range(len(vals) - 1)]
null = []
for t in range(200):
    for i in range(len(vals) - 1):
        D1, D2 = len(str(vals[i])), len(str(vals[i + 1]))
        a = rnd.randrange(10**(D1-1), 10**D1); b = rnd.randrange(10**(D2-1), 10**D2)
        null.append(agree(a, b))
print(f"4. shared leading decimal digits, consecutive terms: mean {sum(cons)/len(cons):.3f} vs random pairs {sum(null)/len(null):.3f}")
# 5. bit-level similarity between consecutive terms of equal bit length
sims = []; nul2 = []
for i in range(len(vals) - 1):
    a, b = bin(vals[i])[2:], bin(vals[i + 1])[2:]
    if len(a) != len(b): continue
    B = len(a); s = sum(x == y for x, y in zip(a[: (B + 1) // 2], b[: (B + 1) // 2])) / ((B + 1) // 2); sims.append(s)
print(f"5. consecutive same-length terms: fraction of equal bits {sum(sims)/max(len(sims),1):.3f} over {len(sims)} pairs (random: 0.5 minus the forced leading 1s)")
# 6. counts per decimal length vs the random model: N_D ~ 10^{D/2} * 2^{B/2} / 10^D  where B = D log2 10 (palindromes in both, independent)
from collections import Counter
c = Counter(len(str(v)) for v in vals)
print("6. terms per decimal length:", dict(sorted(c.items())))
# 7. mod small numbers
for m in (3, 7, 9, 11, 13, 17, 101):
    r = Counter(v % m for v in vals); print(f"7. mod {m}: " + " ".join(f"{k}:{r[k]}" for k in range(m)) if m < 20 else f"7. mod {m}: {len(r)} distinct residues of {m}")
