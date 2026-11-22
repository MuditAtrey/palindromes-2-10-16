"""Check of the bad-run lemma. x_s = h V_s / 5^K mod 1, V_s = 4^s + 4^-s (mod 5^K). 'Bad' = within eps of {0, 1/3, 2/3} mod 1.
Lemma: if v_5(h) = K - r with r >= 2 and eps <= 1/200, every run of consecutive bad s has length <= 1.161 r (= log_4 5 * r).
Also computes c(eps) = max phi(x) over good x, phi(x) = |cos 6 pi x cos 9 pi x|."""
import math, random
eps = 1 / 200
grid = [i / 200000 for i in range(200000)]
def dist3(x): return min(abs(x - k / 3 - m) for k in range(3) for m in (-1, 0, 1))
c = max(abs(math.cos(6 * math.pi * x) * math.cos(9 * math.pi * x)) for x in grid if dist3(x) >= eps)
print(f"c(eps=1/200) = {c:.6f}")
worst_ratio = 0; rnd = random.Random(3); n = 400
for K in (6, 9, 12, 16, 20, 24):
    M = 5 ** K; inv = pow(4, -1, M)
    hs = [rnd.randrange(1, M) for _ in range(3000)] + list(range(1, 200)) + [M - k for k in range(1, 200)] + [5 ** j for j in range(1, K)]
    maxrun_by_r = {}
    for h in hs:
        if h % M == 0: continue
        v = 0; t = h
        while t % 5 == 0: t //= 5; v += 1
        r = K - v
        if r < 2: continue
        run = best = 0; a = pow(4, 1, M); b = inv
        p4, pi = 1, 1
        for s in range(1, n + 1):
            p4 = p4 * 4 % M; pi = pi * inv % M
            x = (h * ((p4 + pi) % M) % M) / M
            if dist3(x) < eps: run += 1; best = max(best, run)
            else: run = 0
        maxrun_by_r[r] = max(maxrun_by_r.get(r, 0), best)
    for r, mr in sorted(maxrun_by_r.items()):
        worst_ratio = max(worst_ratio, mr / (math.log(5, 4) * r))
    print(f"K={K:2d}: max bad run by r: " + ", ".join(f"r={r}:{mr} (bound {math.log(5,4)*r:.1f})" for r, mr in sorted(maxrun_by_r.items()) if r in (2, 3, K // 2, K)))
print(f"largest observed run / bound = {worst_ratio:.3f} (must be <= 1)")
