"""docs/WALL.md section 3: every hex palindrome N of hex length L satisfies N = (-4)^(L-1) * (sum of its hex digits) (mod 25).
(-4 is the 5-adic square root of 16 that is 1 mod 5, so 16^((L-1)/2) = (-4)^(L-1), and 16^x + 16^-x = 2 mod 25 for every x.)"""
import random
random.seed(7); bad = 0; n = 0
for _ in range(5000):
    L = random.randint(1, 120); half = (L + 1) // 2
    for A in ((0, 6, 9, 15), (0, 2, 5, 7), tuple(range(16))):
        e = [random.choice(A) for _ in range(half)]; a = e + e[: L // 2][::-1]
        N = sum(d * 16 ** p for p, d in enumerate(a)); n += 1
        if (N - pow(-4, L - 1, 25) * sum(a)) % 25: bad += 1
print(f"N = (-4)^(L-1) * digit sum (mod 25): {n} hex palindromes, {bad} failures")
