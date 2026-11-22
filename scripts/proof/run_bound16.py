"""Check of Theorem 1's lemma exactly as used: phases x = h (4^s + 4^-s) / 5^K mod 1 over s in steps of 2, bad = within 1/1000 of (1/3)Z,
bound: bad runs <= 0.581 r + 1 where 5^r is the exact denominator. Random and adversarial h, both parities of s."""
import random, math
eps = 1 / 1000
def dist3(x): return min(abs(x - k / 3 - m) for k in range(3) for m in (-1, 0, 1))
rnd = random.Random(17); worst = 0.0; n = 300
for K in (8, 12, 16, 20, 24, 28):
    M = 5 ** K; inv = pow(4, -1, M)
    hs = [rnd.randrange(1, M) for _ in range(1500)] + list(range(1, 150)) + [M - k for k in range(1, 150)] + [5 ** j * u for j in range(1, K) for u in (1, 2, 3)]
    for h in hs:
        if h % M == 0: continue
        v = 0; t = h
        while t % 5 == 0: t //= 5; v += 1
        r = K - v
        if r < 2: continue
        for par in (0, 1):
            run = best = 0
            for k in range(n):
                s = par + 2 * k
                x = (h * ((pow(4, s, M) + pow(inv, s, M)) % M) % M) / M
                if dist3(x) < eps: run += 1; best = max(best, run)
                else: run = 0
            worst = max(worst, best / (0.581 * r + 1))
    print(f"K={K}: worst run / bound so far {worst:.3f}", flush=True)
print(f"largest observed bad run / (0.581 r + 1) = {worst:.3f} (must be <= 1)")
