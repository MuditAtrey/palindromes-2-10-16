"""Count numbers with B bits that are palindromes in base 2 and base 16, by class B mod 4 (exact, via binary palindromes)."""
import math
for B in range(5, 45):
    h = (B + 1) // 2; cnt = 0
    for top in range(1 << (h - 1), 1 << h):
        s = bin(top)[2:]; n = int(s + s[: B // 2][::-1], 2)
        x = format(n, 'x')
        if x == x[::-1]: cnt += 1
    print(f"B={B:2d} (B mod 4 = {B % 4}): {cnt:8d} base-2/16 palindromes   2^(B/4) = {2 ** (B / 4):10.1f}", flush=True) if B % 4 in (1, 2) or B > 36 else None
