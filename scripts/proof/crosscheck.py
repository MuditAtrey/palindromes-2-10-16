"""Soundness cross-check of theta_death: for every L in a range whose theta_L = {(L-1) log10 16} lies in a proved interval, the exact
integer tree for that L (early_death.tree_dies, same depth) must die."""
import math, sys
from decimal import Decimal
import theta_death as TD
from early_death import tree_dies
for cls in (0, 1):
    TD.configure(cls, 6); TD.good.clear(); TD.rec(Decimal(0), Decimal(1), 12)
    G = [(float(a), float(b)) for a, b in TD.good]
    inG = ok = bad = 0
    for L in range(30, 1501):
        th = ((L - 1) * math.log10(16)) % 1
        if any(a + 1e-9 < th < b - 1e-9 for a, b in G):
            inG += 1
            d, _ = tree_dies(L, cls, 6)
            ok += d; bad += (not d)
            if not d: print("  MISMATCH: class", cls, "L", L, "theta", th)
    print(f"class {cls}: {inG} lengths in [30,1500] fall in proved intervals ({100*inG/1471:.1f}%); exact tree dies for {ok}, fails for {bad}")
