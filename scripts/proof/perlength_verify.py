"""Third, independent implementation of the per-length proof (exact integers, prefix by floor division, explicit bit masks),
compared with early_death.py on every class B = 4L, 4L-1 for L in [119, 400]."""
from early_death import tree_dies
def dies2(L, cls, k):
    A = [0, 6, 9, 15] if cls == 0 else [0, 2, 5, 7]; top = [9, 15] if cls == 0 else [5, 7]; mx = A[-1]
    h = (L + 1) // 2
    w = [(16 ** j + 16 ** (L - 1 - j)) if 2 * j != L - 1 else 16 ** j for j in range(h)]
    suffix = [0] * (h + 1)
    for j in range(h - 1, -1, -1): suffix[j] = suffix[j + 1] + mx * w[j]
    nodes = [0]
    for i in range(min(k, h)):
        new = []
        for X0 in nodes:
            for a in (top if i == 0 else A):
                X = X0 + a * w[i]; Y = X + suffix[i + 1]
                DX, DY = len(str(X)), len(str(Y))
                if DX != DY: new.append(X); continue
                D = DX; m = 0
                while m < D // 2 and X // 10 ** (D - m - 1) == Y // 10 ** (D - m - 1): m += 1    # leading m digits shared
                bits = min(m, 4 * (i + 1))
                if bits:
                    P = X // 10 ** (D - m)                     # the shared leading digits as an integer
                    R = int(str(P).zfill(m)[::-1])
                    if (R - X) & ((1 << bits) - 1): continue
                new.append(X)
        nodes = new
        if not nodes: return True, i + 1
    return False, len(nodes)
agree = diff = 0; dead = 0
for L in range(119, 401):
    for cls in (0, 1):
        a = tree_dies(L, cls, 10); b = dies2(L, cls, 10)
        if a[0] != b[0] or (a[0] and a[1] != b[1]): diff += 1; print("  DIFFERENT:", L, cls, a, b)
        else: agree += 1; dead += a[0]
print(f"independent per-length check, classes L = 119-400: {agree} agree, {diff} differ; {dead} proved empty by both")
