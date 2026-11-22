"""Rigorous small-denominator case of the exponential-sum lemma.
F(h) = prod_s phi(h V_s / 5^K), phi(x) = |cos(6 pi x) cos(9 pi x)| (hex digit 6a + 9b), V_s = 4^s + 4^-s (5-adic unit times the
pair weight). If v_5(h) = K - r (r >= 1) then h V_s / 5^K = h' V_s / 5^r mod 1 with 5 not dividing h', and V_s mod 5^r is periodic
in s with period P_r = ord(4 mod 5^r). Over any P_r consecutive s the factors multiply to rho(r, h') < 1 (no factor is 1: x = a/5^r
with 5 not dividing a is never 0, 1/3, 2/3 mod 1). Hence |F(h)| <= rho_max(r)^floor(n / P_r): exponential decay for all such h when
n >> P_r. This prints P_r and rho_max(r) (worst h') and the per-factor worst case."""
import math
def phi(x): return abs(math.cos(6 * math.pi * x) * math.cos(9 * math.pi * x))
for r in range(1, 6):
    M = 5 ** r
    P = next(k for k in range(1, 4 * M) if pow(4, k, M) == 1)
    inv = pow(4, -1, M)
    V = [(pow(4, s, M) + pow(inv, s, M)) % M for s in range(P)]
    worst = 0.0; worst_factor = 0.0
    for hp in range(1, M):
        if hp % 5 == 0: continue
        logp = 0.0
        for v in V:
            f = phi((hp * v % M) / M); worst_factor = max(worst_factor, f); logp += math.log(f)
        worst = max(worst, math.exp(logp))
    print(f"r={r}: period P={P:5d}, worst product over one period rho_max={worst:.3e} (per factor {worst**(1/P):.4f}); worst single factor {worst_factor:.6f}")
