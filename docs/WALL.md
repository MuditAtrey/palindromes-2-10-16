# The powers-of-5 wall: what any proof for every length must get past

Written 2026-11-22 during the autonomous attack on the wall. Everything here is either proved in a few lines or measured; each item says which.

## 1. Bounded information has a hard ceiling (measured)

Information available at depth k of the outside-in tree (k outer hex digit pairs fixed, any length L >= 2k + 2):
- size: the leading decimal digits shared by the whole subtree, a function of theta_L = {(L - 1) log10 16} and the top digits;
- 2-adic: those digits mirrored to the bottom, compared with the bottom hex digits (Theorem 2's rule).

`scripts/proof/death_curve.py` runs this exact tree to depth 40 for every L in [100, 4100] (both classes, 8,002 lengths;
`results/wall/death_curve_c0.jsonl`, `_c1.jsonl`). Cumulative share of lengths proved empty by depth k:

| depth k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 11 and beyond |
|---|---|---|---|---|---|---|---|---|---|
| B = 4L | 6.5% | 9.9% | 11.4% | 12.3% | 12.9% | 13.3% | 13.5% | 13.6% | **13.65%** |
| B = 4L - 1 | 3.2% | 6.2% | 8.2% | 9.4% | 10.2% | 10.5% | 10.8% | 10.8% | **10.87%** |

No length dies after depth 11; every other length passes 3,000 surviving branches by depth 12-25 and keeps growing about x1.7 per
level. This is a supercritical branching process: it dies out with a fixed probability (13.65% / 10.87%) and otherwise never.
Theorem 2's rigorous 10.3% / 6.9% (depth 6, theta-intervals of width 2^-12) is close to this ceiling.

Consequence for the Baker route (prove the surviving thetas are "special", then use linear forms in logarithms to keep theta_L away from
them): dead. The surviving set has measure about 0.86 / 0.89, not a thin exceptional set.

## 2. Adding 5-adic (or any coprime modular) information at bounded depth adds nothing (proved)

At depth k the leading m ~ 1.2k decimal digits are known, so N mod 5^m is known. N = O + I with O the outer part (fixed) and
I = sum over the n = L/2 - k inner positions of e_j W_j. Each W_j = 16^j + 16^(L-1-j) is 2 mod 5, so the digit 6 (class 4L) or 2
(class 4L-1) times W_j is a unit u_j mod 5^m. Proof that n >= 5^m - 1 inner positions reach every residue: if a set S of residues
satisfies S + {0, u} = S for a unit u, then S is closed under adding u, which generates Z/5^m, so S is everything; hence each new
position (digit 0 or the unit multiple) enlarges S until it is everything. In practice far fewer are needed: `scripts/proof/inner_cover.py`
finds full coverage once n >= 1.2 m + 3 (m <= 8, worst of 60 random L and depths, both classes). For every fixed depth k and every L
large enough, the 5-adic condition at depth k removes nothing. The same argument
applies to every modulus coprime to 10 (mod 11, 17, 3: the only ones with forced residues, `checksums.py`). So the 2-adic + size ceiling of
section 1 is the ceiling for all bounded-depth information.

**Barrier statement.** A proof that covers every length must use the decimal mirror at a depth that grows with L: the 5-adic test only
removes branches once the number of inner digits n satisfies 4^n < 5^m, i.e. at depth about L / 4.79 (measured in the engine: STATUS 39-40,
99.8% of the nodes that reach that level die there). Beyond the measured ceiling, no argument of the form "look at the outer digits only"
can work, for any number of outer digits that does not grow with L.

## 3. What the 5-adic condition says, exactly (proved, verified)

Write the hex palindrome as N = sum_p a_p 16^p over all L positions (a_p = a_(L-1-p)), c = (L - 1)/2, and lambda = log(16) in the 5-adic
numbers (v5(lambda) = 1). The 5-adic square root of 16 that is 1 mod 5 is -4, so 16^c = (-4)^(L-1), and by symmetry

    N = (-4)^(L-1) * sum_i lambda^(2i) / (2i)! * mu_(2i),     mu_(2i) = sum_p a_p (p - c)^(2i)   (even central moments)

with v5(lambda^(2i)/(2i)!) = 2i - v5((2i)!). The first case:

    N = (-4)^(L-1) * (sum of all hex digits)   (mod 25)

(checked on 15,000 random hex palindromes, any digits, 0 failures: `scripts/proof/moment_identity.py`). With the decimal mirror (last two decimal digits = first two
reversed), the hex digit sum of any triple palindrome is fixed mod 25 by L and its first two decimal digits. Mod 5^K the condition is a
linear condition on the first ~0.67K even central moments of the digit sequence. These moments are global: every digit, including the
innermost, contributes to them, which is the same barrier as section 2 seen from the other side. (The moments are small integers,
mu_(2i) < 15 L^(2i+1), but the number of distinct moment vectors needed exceeds 5^K, so this gives no compression either.)

## 4. Where that leaves a proof for every length

- Every exact test that is uniform in L (size, 2-adic, congruences) is captured by sections 1 and 2 and stops at 13.65% / 10.87% of lengths.
- The remaining lengths need, for each L, a statement of the form "none of the ~2^(0.17 L) surviving prefixes has its 5-adic target in the
  inner reachable set", where both the targets and the reachable sets measure as random (STATUS 39, Fourier spectrum flat, RESEARCH_PLAN B2).
- The 372-world sibling sweep (RESEARCH_PLAN, last section) shows no mechanism anywhere: the random model predicts every sibling exactly.
- So a proof for every length would have to explain a random-looking event for each L. That is the shape of Erdős's ternary conjecture
  (2^n has a digit 2 in base 3 for n > 8) and the 82000 problem, both open. No tool in number theory today proves statements of this shape.

## 5. The centre does not help either (measured)

The innermost positions have weights W_j = 4^(L-1) (4^s + 4^-s) with s = L - 1 - 2j small, so the inner reachable set at depth d is
4^(L-1) C_n with C_n = {sum over the n innermost s of e_s (4^s + 4^-s)}: the same set for every length of the same parity. This is real
structure (the engine's cross-length lattice cache uses it). But the targets t(p) 4^-(L-1) that must land in C_n measure as random
relative to it: `scripts/structure/target_side.py`, 42-54 hex digits, 11-13 inner pairs, the number of targets that land matches the random
expectation at every depth (e.g. 30 vs 36.2, 100 vs 109.3, 28 vs 31.5, 15 vs 15.2, 5 vs 5.2). Pooled over two runs (L 36-48 and
40-50, 12,255 targets): 263 land vs 268.9 expected (z = -0.36) and 195 vs 214.3 (z = -1.32). No excess, so no hidden mechanism
favouring hits, and no deficit strong enough to suggest one forbidding them.

## 6. The picture in one sentence

A triple palindrome is an integer whose real size (leading digits, base 16 and 10), 2-adic part (trailing bits) and 5-adic part (trailing
decimal digits mod 5^K) all agree with mirrored digit patterns; size and the 2-adic part can be matched uniformly in L and give the
13.65% / 10.87% ceiling, while the 5-adic part is tied to the size only through the fact that N is an integer, i.e. at full precision,
which is exactly where every known method, mine and the literature's, stops.
