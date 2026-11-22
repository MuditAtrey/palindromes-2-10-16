"""J-ladder: a structure detector from near-misses.

Among numbers palindromic in base 2 and 16, count those whose first J decimal digits mirror their last J
(weak palindromes, exact engine search with TP_WEAKPAL=J). Under the random model each extra mirrored digit
divides the count by 10. A persistent excess at high J would be the fingerprint of a mechanism making
decimal palindromes more likely than chance among these numbers. Counts are exact (exhaustive), so the only
noise is Poisson.

Usage: python3 j_ladder.py B0 B1 J0 J1 [threads]
"""
import math, os, re, subprocess, sys, time
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../engine/tripal")
B0, B1, J0, J1 = (int(x) for x in sys.argv[1:5]); TH = sys.argv[5] if len(sys.argv) > 5 else "8"
CLEAN = ("TP_R5", "TP_WEAKPAL", "TP_DUAL", "TP_DUAL_K", "TP_DUAL_ALL", "TP_CVP", "TP_CVP_K", "TP_CVP_ALL", "TP_CVP_PRUNE", "TP_CVP_SEED")
def count(J):
    e = {k: v for k, v in os.environ.items() if k not in CLEAN}
    # the dual certificate with modulus 5^J is exact for this property and makes high J cheap (test suite: no hit lost)
    e.update({"TP_R5": "0", "TP_WEAKPAL": str(J), "TP_DUAL": "1", "TP_DUAL_K": str(J), "TP_DUAL_ALL": "1"})
    t0 = time.time()
    p = subprocess.run([BIN, "hexdfs", "--bits-from", str(B0), "--bits-to", str(B1), "--threads", TH], capture_output=True, text=True, env=e)
    m = re.search(r"total: (\d+) nodes, \d+ prunes, \d+ leaves, (\d+) hits", p.stdout)
    return int(m.group(2)), int(m.group(1)), time.time() - t0
rows = [(J,) + count(J) for J in range(J0, J1 + 1)]
print(f"weak palindromes among base-2/16 palindromes, B = {B0}..{B1}")
print(f"{'J':>3} {'count':>10} {'ratio to J-1':>13} {'z vs ratio 10':>14} {'nodes':>14} {'s':>6}")
prev = None
for J, c, n, dt in rows:
    if prev and prev > 0 and c > 0:
        r = prev / c
        # if count(J) ~ Poisson(prev/10): z of the observed count against prev/10
        z = (c - prev / 10) / math.sqrt(prev / 10)
        print(f"{J:>3} {c:>10} {r:>13.2f} {z:>+14.2f} {n:>14,d} {dt:>6.1f}")
    else:
        print(f"{J:>3} {c:>10} {'':>13} {'':>14} {n:>14,d} {dt:>6.1f}")
    prev = c
