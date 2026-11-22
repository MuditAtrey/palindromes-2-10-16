# Proof research: is there a number above 9 that is a palindrome in bases 2, 10 and 16?

Started 2026-11-21. Goal: a proof of non-existence (or of existence), so that further search is only done if a proof says it can pay.

## 1. The problem in proof form

A triple palindrome N with B bits lies in one of four classes by B mod 4. For B = 4L its hex digits are in {0,6,9,F} (each hex digit must
be its own 4-bit reversal); for B = 4L-1 they are in {0,2,5,7}; the classes B = 4L-2, 4L-3 are forced into tiny families. So N is a base-16
palindrome with L/2 free digit pairs from a 4-letter alphabet (about 2^L candidates of size 2^B) that must also be a decimal palindrome
(D/2 free digits, about 10^(D/2) = 2^(B/2) numbers of size 2^B). In the language of digit-restricted sets: a set of "dimension" 1/4 in base 16
intersected with a set of dimension 1/2 in base 10, two multiplicatively independent bases. Dimensions sum to 3/4 < 1, so the random model
predicts finitely many solutions (about 2^(-B/4) per length; 1.6e-35 expected beyond 2^472).

## 2. What is known (literature, checked)

- Bérczes and Ziegler, *On simultaneous palindromes* (arXiv 1403.0787, 2014): finiteness for bases (g, h) with h | g type relations via
  Baker's method; for bases 2 and 10 they only guess finiteness. Nothing proven for 2 and 10.
- Goins (Integers 9, 2009, arXiv 0909.5452): finitely many (exactly 203) d-digit base-10 palindromes that are d-digit palindromes in another
  base: same digit count only, not this case.
- Bašić (IJNT 8, 2012; paper read): for every K there exist K bases, constructed from the number (b_i = n/a_i - 1), in which a number is a
  d-digit palindrome. The bases are chosen by the proof; nothing about fixed bases 2, 10, 16.
- Senge-Straus (effective by Stewart): numbers with boundedly many nonzero digits in two multiplicatively independent bases are finite.
  The candidates here have about half their digits nonzero: not applicable, except to sparse sub-families.
- Burrell and Yu, *Digit expansions of numbers in different bases* (arXiv 1905.00832): "digit-special" integers below N number at most
  C(eps) N^eps; dimension-zero statements. A density bound, not finiteness.
- Analogous open problems of the same type (restricted digits in independent bases, dimension sum below 1):
  - Erdős: 2^n has a digit 2 in base 3 for n > 8. Open; verified to n ~ 5.9e21 (Saye 2022); Lagarias (JLMS 2009) studies a weak form and
    Hausdorff dimension of 3-adic intersections; recent preprint arXiv 2610.03789 (no proof).
  - 0, 1, 82000 are the only numbers using digits 0/1 in bases 2, 3, 4, 5 (OEIS A146025): open; conjecture recorded by Burrell and Yu.
- Furstenberg intersection conjecture, proved by Shmerkin and Wu (Annals 2019): dimension of the intersection of x p and x q invariant sets is
  at most max(0, sum - 1). It is about real numbers and dimensions; for integers it gives counting at scale, not finiteness. The search tree's
  growth equals this bound exactly (notebook 5.12).

## 3. Results of this program so far

- **No congruence obstruction exists** (`scripts/proof/local_obstruction.py`). For every modulus m <= 2000 coprime to 10, every bit length
  (hex lengths 20-40, and 30-34 for m up to 2000) and every compatible decimal length: 19,132 cases, the residues of candidates and of decimal
  palindromes always intersect. Reason: decimal palindromes are restricted mod m only through 11 (even length), the hex candidates only
  through 3 and 17; these never conflict. So non-existence cannot be proved by congruences; any proof must be global (size and
  Diophantine approximation, or dynamics), as in the analogous open problems.

## 4. Proof strategies, with an honest feasibility note

Non-existence (or finiteness):
1. **Density bound N^eps for triple palindromes** (adapt Burrell-Yu; Fourier analysis of digit-restricted sets). Achievable in principle,
   publishable, but proves rarity, not absence.
2. **Finiteness via a Furstenberg/Shmerkin-Wu type argument for integers.** Would settle a problem of the same type as Erdős's ternary
   conjecture and the 82000 problem; no known method. Research-level open.
3. **Effective finiteness with an explicit bound, then search to that bound.** Needs a Baker-type lower bound for linear forms with about
   B/4 free digits; known methods only handle a bounded number of terms. Out of reach as stated; reachable for sparse sub-families (proves
   those empty for all sizes).
4. **Conditional proofs.** Under a quantitative equidistribution hypothesis for digits of digit-restricted numbers in a second base (the
   kind assumed in random models), finiteness follows from counting; worth writing precisely, since it states exactly what is missing.

Existence: no non-constructive route is available, because the expected number of solutions is far below 1; existence could only be
shown by exhibiting N, i.e. by search.

## 4b. Reading of Burrell-Yu (arXiv 1905.00832) and where it stops for palindromes

Their Lemma 4.2: integers whose base-p and base-q digits avoid given digits (sets invariant under x p and x q) with dimension sum below
1 number at most C N^eps. Method: renormalise n by p^m to a point on a line through the origin whose slope is q^{m log p / log q mod 1},
so the integers become points of a line meeting the product of two invariant Cantor sets; Wu's uniform Furstenberg theorem bounds the
covering number of such intersections by 2^(eps m) uniformly in the line. Requires digit sets defined by the same rule at every position.
Palindromes are not of that kind: the condition ties the top digits to the reversed bottom digits, and reversal does not commute with
x p or x q. So the theorem does not transfer. No counting bound for palindromes in two bases exists in the literature (searched:
Banks-Hart-Sakata, Banks-Shparlinski, Col, Tuxanidy-Panario; the 2024 survey in arXiv 2307.16637 calls exponential-sum bounds over
palindromes "rather limited").

## 4c. A theorem that looks provable: a power-saving upper bound

Claim (to prove): the number of triple palindromes N <= X is O(X^(c + eps)) with c well below the trivial 1/4 (c ~ 0.105 from the
argument below, possibly ~0.042 with more work).

Sketch. Fix a bit length B = 4L (B = 4L - 1 is the same with digits {0,2,5,7}). Fix d outer digit pairs (4^d choices). They fix the top
~1.2 d decimal digits of N (up to a bounded carry ambiguity), hence by the decimal mirror the residue of N mod 5^K with K ~ 1.2 d, as a
target t. The inner digits (u = L - 2d bits, 4^(L/2 - d) choices) must then give an inner sum = t (mod 5^K). If the inner sums are
equidistributed mod 5^K with a power-saving error, every fibre has size 4^(L/2-d) 5^(-K) (1 + o(1)), and
   #solutions <= 4^d * 4^(L/2-d) * 5^(-K) * (1 + o(1)) = 2^L * 5^(-1.2 d) (1 + o(1)).
Equidistribution needs the dense regime, roughly 5^K <= 2^(u(1 - delta)), i.e. 2.79 d <= (L - 2d)(1 - delta), d up to ~L/4.8. Then the count
per length is 2^(L - 2.79 L/4.8) = 2^(0.42 L) = X^(0.105). Using the 2-adic half too (only ~2^(0.8 d) of the 4^d prefixes survive it) would
give 2^(0.8 d) * 2^(L - 2d) * 5^(-K) at d ~ L/4.8, i.e. X^(0.042) (the perfect-oracle count), but that needs a rigorous count of 2-adic
survivors, which is exactly a Furstenberg-type statement for reversal-coupled sets (not available).

Key lemma needed: for h not divisible by 5^K, with W_tau = 16^tau + 16^(-tau) (times a unit),
   | prod_{tau in inner} cos(6 pi h W_tau / 5^K) cos(9 pi h W_tau / 5^K) | <= 4^(-(L/2 - d) delta')
(the digit sum's Fourier transform factors because each hex digit is 6a + 9b). The factor equals 1 only when h W_tau / 5^K is near 0, 1/3
or 2/3 mod 1, which a 5-power denominator avoids; the work is a uniform lower bound on how many tau are far from those points, for every h,
including h divisible by high powers of 5 (then W_tau mod 5^r has period 5^(r-1) in tau, and W_tau = 2 mod 5 always). Numerically (B2) the
largest coefficients over all h are 0.17 (10 positions, K = 9) and 0.077 (12 positions, K = 10), decaying as expected.

Value: a first theorem on this problem (and, with bases 2 and 10 alone, on double palindromes) - rarity with an explicit exponent. Not
finiteness: finiteness needs the count per length to go to 0, i.e. c < 0 per length, which the information budget never gives (the
expected count per length is 2^(-B/4), below any power bound only if the 5-adic and 2-adic information can be used all the way, the
open part).

## 4d. Result: Theorem 1 (docs/THEOREM_1.md)

Proved: T(X) <= C X^(1/4) 2^(-0.04 sqrt(log_16 X)) for the number T(X) of triple palindromes up to X, below the trivial X^(1/4).
Correction to 4c: the X^0.105 sketch needed every Fourier coefficient below 5^(-K) at density near 1, which is impossible (even random weights
have a largest coefficient near 2^(-n)); the sound route gives the sub-power saving above. Key lemma: bad stretches of the phases
h (4^s + 4^-s) / 5^K are at most 0.581 r + 1 long (x16 dynamics, separation of 5^r from 1/765), checked numerically. Small-denominator
frequencies have explicit exponential decay (`scripts/proof/small_denominators.py`). The obstruction to a power saving, and to finiteness,
is named exactly: bounded bad runs for every h (a statement about short segments of 4^s mod 5^r with r >> log n).

## 5. Next steps (in order)

1. (done) Burrell-Yu read: method needs shift-invariant digit sets; palindromes are reversal-coupled, so it does not transfer (4b).
2. Prove the key exponential-sum lemma of 4c (uniform over h, including high powers of 5), then the X^(0.105) bound; test the lemma
   numerically at larger sizes first.
3. Write the conditional finiteness theorem (strategy 4) precisely.
4. Baker-type exclusion of every sparse family for all sizes (strategy 3, partial).
5. Keep a list of every place where a proof attempt needs a statement that is equivalent to an open problem, so that the obstruction is
   named exactly.
