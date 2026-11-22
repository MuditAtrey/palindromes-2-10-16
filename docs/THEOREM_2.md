# Theorem 2: no triple palindrome exists for a positive proportion of all bit lengths

Status: proved 2026-11-21 by exact computation plus a uniformity argument; cross-checked (`scripts/proof/`). This is a non-existence
result (not a rarity bound): for every bit length in an explicit infinite set, there is no triple palindrome of that length at all.

## Statement

Let theta_L = {(L - 1) log10 16} (fractional part). There are explicit finite unions of intervals G0, G1 in [0, 1) with measures
|G0| >= 0.1028 and |G1| >= 0.0686 (depth 6; deeper trees enlarge them) such that for every L >= 12:
- if theta_L is in G0, no number with exactly 4L bits is a palindrome in bases 2, 10 and 16;
- if theta_L is in G1, no number with exactly 4L - 1 bits is.

Since log10 16 is irrational, theta_L is equidistributed, so these hold for a set of L of natural density |G0| (resp. |G1|): at least
10.3% (resp. 6.9%) of all lengths, out to infinity. Lengths far beyond any search are covered.

Concrete consequence beyond the searched range: exact per-length trees (no uniformity needed) prove 70 of the 564 bit lengths
between 473 and 1600 (classes 4L and 4L - 1) empty, e.g. B = 475, 527, 531, 535, 567, 583, ..., 1600
(`results/empty_bit_lengths_473_1600.txt`).

## Proof

1. Any number with 4L bits that is a palindrome in bases 2 and 16 has hex digits in {0,6,9,F} forming a palindrome, top digit 9 or F
   (each hex digit equals its own 4-bit reversal); with 4L - 1 bits, digits {0,2,5,7}, top digit 5 or 7.
2. Outside-in tree: fix the top i digits t_0 ... t_(i-1) (and by the palindrome the bottom i). All completions lie in an interval
   [X, X + R_i] with R_i the largest contribution of the remaining digits. If every number of that interval has the same first m decimal
   digits P (m at most half the decimal length), a decimal palindrome among them satisfies N = rev(P) (mod 10^m), hence
   N = rev(P) (mod 2^m). The remaining digits change N only by multiples of 16^i at the bottom (and far above at the top), so
   N = sum_{j<i} t_j 16^j (mod 2^(4i)). If rev(P) and sum_{j<i} t_j 16^j differ mod 2^min(m, 4i), the whole subtree contains no
   triple palindrome. This rule is exact.
3. Uniformity in L: N = 16^(L-1) T with T = t_0 . t_1 t_2 ... (base 16) in [t_0, t_0 + 1), and 16^(L-1) = 10^(E + theta_L) with E an
   integer. So the leading decimal digits of the subtree interval are those of [T 10^theta_L, (T + tail) 10^theta_L]: they depend on L only
   through theta_L. The bottom value sum_{j<i} t_j 16^j does not depend on L at all (the mirrored top terms are multiples of 16^(L-1-j),
   far above bit 4i once L >= 2i + 2). Hence the depth-k tree's survival is a function of theta_L alone.
4. Interval arithmetic: for a theta-interval J, the shared decimal prefix of [T_min 10^min J, (T + tail) 10^max J] is a prefix shared for
   every theta in J; applying step 2 with it is valid for every theta in J. Bisection of [0,1) into intervals of width 2^-12 and a depth-6
   tree (60-digit decimal arithmetic, 5 guard digits) give the intervals of G0 and G1 (`scripts/proof/theta_death.py`).
5. Equidistribution of theta_L (Weyl) gives the density statement.

## Checks

- Cross-check (`scripts/proof/crosscheck.py`): every L in [30, 1500] with theta_L in G0 or G1 (141 and 99 lengths) was run through the exact
  integer tree for that specific L (`early_death.py`); all 240 die. Their frequencies (9.6%, 6.7%) match the measures.
- Consistency with the full search: lengths below 2^472 in these sets were also found empty by the engine (as all were).

## Limits (stated exactly)

This method uses only the size of N and the 2-adic half of the decimal mirror, which depend on L only through theta_L. Outside the death
set the tree keeps branching (about x1.73 per level), like a supercritical branching process, which dies out only with a fixed probability.
So this method proves a fixed fraction of lengths empty (deeper trees approach that fraction), never all of them. The 5-adic half of the
mirror bites only when few inner digits remain, i.e. at depths proportional to L, so it cannot be made uniform in L in this way. A proof
for all lengths needs a new idea at exactly that point (see PROOF_RESEARCH.md).

## Adversarial checks (2026-11-21)

1. **Rounding hole found and closed.** The first version read decimal digits from 60-digit arithmetic rounded to nearest; a carry from
   rounding at digit 40 could in principle create a false shared prefix. Fixed: both ends widened by a relative 1e-50 and rounded outward
   (floor / ceiling) before digits are read. The proved sets did not change (63 and 69 intervals, measures 0.1028 and 0.0686).
2. **Independent re-implementation** (`scripts/proof/theta_verify.py`): different code, mpmath interval arithmetic with outward rounding,
   digits by integer floors, its own bisection. It proves exactly the same interval sets (identical, none extra on either side).
3. **Positive controls of the pruning rule** (`scripts/proof/controls.py`): 800 planted solutions in both hex classes (L = 12-89, a
   decimal-palindrome offset c, the shifted rule) survive at every depth; all 207 known base-2/10 double palindromes (OEIS A007632) survive
   along their own branches.
4. **The controls can fail** (mutation test): with a deliberately wrong rule (one bit too many) they kill 399/400 planted solutions and
   199/207 A007632 terms; and 400/400 random non-solution branches are killed by the correct rule.
5. **Cross-check per length**: all 240 lengths in [30, 1500] whose theta falls in a proved interval die in the exact integer tree for that L.
6. **Separate C engine** on the claimed-empty bit lengths (`tripal hexdfs --bits-from B --bits-to B`, no lattice tests): B = 475 (hex
   length 119) finishes with 10 nodes and 0 hits. The engine's per-level tables stop at 64 levels (hex length 128), so it aborts on the
   claims from B = 527 up; it was not modified for this check.
8. **Third independent implementation** of the per-length proof (`scripts/proof/perlength_verify.py`: exact integers, prefix by floor
   division, bit masks): agrees with `early_death.py` on all 564 classes B = 4L, 4L-1 for L in [119, 400], the same 70 proved empty at the
   same depths.
9. **Digit-set characterisation** checked against brute force over all binary palindromes (B = 20-40, both classes): the sets of base-2/16
   palindromes equal the {0,6,9,F} / {0,2,5,7} enumerations exactly (not just in size).
7. **Line-by-line review** of the proof's assumptions: the mirror relation holds for every m <= D (not only D/2); bottom independence needs
   L >= 2i + 2 (L >= 12 at depth 6, as stated); undecided digits change N only by multiples of 16^(i+1); the tail bound covers all remaining
   digits; reversed prefixes with leading zeros parse correctly.
