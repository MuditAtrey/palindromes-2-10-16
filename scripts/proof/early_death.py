"""Lengths where no triple palindrome can exist, proved by a shallow tree.
For B = 4L (hex digits {0,6,9,F}, top digit 9 or F) and B = 4L-1 ({0,2,5,7}, top 5 or 7), run the outside-in tree with only the
exact prefix rule (shared decimal prefix of the subtree, mirrored to the bottom, compared 2-adically with the bottom hex digits)
to depth k. Zero survivors at depth k proves that no number with B bits is a triple palindrome (every candidate lies in some subtree,
and each subtree was refuted exactly). Prints, per class, which lengths die and how that depends on theta = {L log10 16}."""
import sys, math
def tree_dies(L, cls, k):
    A = (0, 6, 9, 15) if cls == 0 else (0, 2, 5, 7)
    top = (9, 15) if cls == 0 else (5, 7)
    pw = [16 ** j + 16 ** (L - 1 - j) if j != L - 1 - j else 16 ** j for j in range((L + 1) // 2)]
    rest = [15 * sum(pw[j] for j in range(i, (L + 1) // 2)) if cls == 0 else 7 * sum(pw[j] for j in range(i, (L + 1) // 2)) for i in range((L + 1) // 2 + 1)]
    level = [0]
    for i in range(min(k, (L + 1) // 2)):
        nxt = []
        for O in level:
            for a in (top if i == 0 else A):
                X = O + a * pw[i]; hi = X + rest[i + 1]
                lo_s, hi_s = str(X), str(hi)
                if len(lo_s) != len(hi_s): nxt.append(X); continue          # decimal length not yet fixed: keep
                m = 0
                while m < len(lo_s) and lo_s[m] == hi_s[m]: m += 1
                m = min(m, len(lo_s) // 2)
                kk = min(m, 4 * (i + 1))
                R = int(lo_s[:m][::-1]) if m else 0
                if kk and (R - X) % (1 << kk): continue
                nxt.append(X)
        level = nxt
        if not level: return True, i + 1
    return False, len(level)
if __name__ == "__main__":
    L0, L1, k = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    for cls in (0, 1):
        dead = []; tot = 0
        for L in range(L0, L1 + 1):
            d, info = tree_dies(L, cls, k); tot += 1
            if d: dead.append((L, info, round((L * math.log10(16)) % 1, 4)))
        name = "B=4L" if cls == 0 else "B=4L-1"
        print(f"{name}: {len(dead)}/{tot} lengths L in [{L0},{L1}] proved empty by a depth-{k} tree ({100*len(dead)/tot:.1f}%)")
        print("   examples (L, depth of death, theta):", dead[:12])
