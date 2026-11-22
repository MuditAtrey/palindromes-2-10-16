# Are known base-2/base-10 double palindromes (A007632) closer to being hex palindromes than random binary palindromes?
# Statistic: fraction of mirrored hex digit pairs that agree. Null: random binary palindromes of the same bit length.
import random, math
F = __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../../engine/data/oeis/b007632.txt")
vals = sorted({int(l.split()[1]) for l in open(F) if l.strip() and not l.startswith('#')})
vals = [v for v in vals if v.bit_length() >= 24]
def agree(N):
    h = format(N, "x"); L = len(h); pairs = L // 2
    return sum(h[i] == h[L - 1 - i] for i in range(pairs)), pairs
rnd = random.Random(11)
def rand_binpal(B):
    half = [1] + [rnd.getrandbits(1) for _ in range((B + 1) // 2 - 1)]
    bits = half + half[: B // 2][::-1]
    return int("".join(map(str, bits)), 2)
obs = 0; exp = 0.0; var = 0.0; tot_pairs = 0
for N in vals:
    a, p = agree(N); B = N.bit_length()
    sims = [agree(rand_binpal(B))[0] for _ in range(400)]
    mu = sum(sims) / len(sims); sd2 = sum((x - mu) ** 2 for x in sims) / (len(sims) - 1)
    obs += a; exp += mu; var += sd2; tot_pairs += p
z = (obs - exp) / math.sqrt(var)
print(f"{len(vals)} double palindromes with >= 24 bits, {tot_pairs} mirrored hex pairs")
print(f"agreeing pairs: observed {obs}, random binary palindromes of the same lengths {exp:.1f} +- {math.sqrt(var):.1f}  (z = {z:+.2f})")
best = sorted(vals, key=lambda N: -agree(N)[0] / agree(N)[1])[:5]
print("closest to hex palindromes:", [(format(N, 'X'), f'{agree(N)[0]}/{agree(N)[1]}') for N in best])
