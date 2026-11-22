# What I tried, what worked, and what failed

A record of every serious approach to the question "is there a number above 9 that is a palindrome in bases 2, 10 and 16?",
including the ones that failed, so nobody has to repeat them. Each item says whether it was measured, proved, or reasoned
from a known theorem, and points to the script or document with the details.

The short version: the search is fast because of a handful of exact tests on the powers of 5; no method I found, from any
field, gets past the obstacle described in [WALL.md](WALL.md), and the question is open.

## 1. Making the search faster

The engine went from 13.4 hours to 17 seconds on one length (92 hex digits) and from 1387 s to about 110 s on another
(114 hex digits) through the steps below. Every step was checked with weak-palindrome and planted-solution controls
before it was trusted.

| Idea | Outcome |
|---|---|
| Build candidates as mirrored hex digit pairs from the outside in, with the 2-adic test (decimal prefix mirrored against the known low bits) | **Works.** The base of everything. |
| 5-adic reach tables (subset sums of undecided digits modulo 5^J) | **Works.** Each extra power of 5 removes about 1.9x of the tree; 5^13 is the practical memory limit (4.7 GB). |
| Digit-sum and mod 11 sieves | **Small win** (1.2 to 1.5x). |
| Carried decimal representation, chunk kernel, prefetching | **2 to 3x** on the constant factor. |
| Bloom filter of reachable residues deeper in the tree | **Works**, about 3x fewer nodes; no false negatives by construction. |
| Dual-lattice window certificate (short dual vectors from LLL give windows a reachable residue must fall in) | **Works, and is exact with integer arithmetic.** Fires at about 0.079 B undecided bits versus 32 for the Bloom filter. Used for the certified 116-118 hex digit runs. |
| Exact lattice enumeration (closest-vector style) inside the digit box | **Fastest** configuration, but steered by floating point, so used only for speed, never for certification. |
| Lattice reduction alone to find solutions directly (embedding, no tree) | **Fails at size.** Finds planted solutions up to about 144 bits, not at 192; a BKZ cost estimate at 460 bits is 2^60 against 2^40 for the tree. |
| GPU port of the core arithmetic | Measured 3.2x on the inner arithmetic; **not built** into the certified engine. |
| Probabilistic pruning (accept some risk of missing a solution) | **No gain.** Planted solutions were found at 0.55 to 0.9 times the rate of exact runs per unit of compute. |

The tree still grows like roughly 2^(0.07 to 0.10 per bit) of the target size. A perfect 5-adic reachability oracle would
give about 2^(0.042 B); nothing I built or found reaches it.

## 2. Looking for structure that would make solutions likely or impossible

| Idea | Outcome |
|---|---|
| Random model: treat decimal and binary digits as independent | **Fits everything I could check.** It predicts the 212 known base-2/10 double palindromes (212 observed, 206.5 expected, chi^2 per degree of freedom 0.96) and the counts in sibling problems. Under it, the expected number of triple palindromes above 2^472 is about 10^-35. |
| Sibling problems: palindromes in bases 2, r = 2^k and an even base b (372 combinations: r = 16, 8, 4; b from 6 to 200) | **No mechanism.** 20 large solutions observed against 20.4 expected (P = 0.56). See RESEARCH_PLAN.md, last section; `scripts/siblings/`. |
| Machine learning on near-misses (neural nets, 1.2 million near-palindromes) | **No signal** beyond the 2-adic rule until almost every digit is known. A trap worth knowing: samples from a timed-out search are biased by DFS order and faked an AUC of 0.96. |
| Fourier analysis of the digit measure modulo 5^9 and 5^10 | **Flat spectrum**, the same as random weights (`scripts/structure/fourier_decay.py`). |
| Forced residues (cyclotomic checksums), all moduli below 20,000 | Only 11 for even-length decimal palindromes (used), 17 and 3 for even-length hex and binary; they constrain the decimal side only (`scripts/structure/checksums.py`). |
| Closed formulas and structured families (two-block, periodic, sparse, windowed digit patterns; 1.5 billion candidates) | **None below 10^5000.** By the Senge-Straus theorem (effective form by Stewart), only finitely many numbers have few nonzero digits in two independent bases, so no sparse formula can give infinitely many solutions. |
| Whether the middle digits favour or avoid solutions | **Random rate.** The reachable sets of the innermost digits are the same for every length, but the targets land in them exactly as often as chance predicts (263 vs 268.9 and 195 vs 214.3 over 12,255 targets; `scripts/structure/target_side.py`). |

## 3. Looking for a proof that covers every length

| Idea | Outcome |
|---|---|
| Rarity bound | **Proved:** at most X^(1/4) 2^(-c sqrt(log X)) triple palindromes below X ([THEOREM_1.md](THEOREM_1.md)). True, but it cannot show there are none. |
| Empty lengths | **Proved:** no triple palindrome for a positive proportion of bit lengths; 70 lengths between 473 and 1600 bits shown empty ([THEOREM_2.md](THEOREM_2.md)). Covers about 10% of lengths, not all. |
| How far "outer digits only" arguments can go | **Measured ceiling:** exactly 13.65% of lengths (B = 4L) and 10.87% (B = 4L-1), reached by depth 11, over 8,002 lengths ([WALL.md](WALL.md), section 1). |
| Adding 5-adic or any coprime modular information at bounded depth | **Proved useless:** the free inner digits reach every residue modulo 5^m once there are about 1.2m + 3 of them (WALL.md, section 2). |
| Baker's method (linear forms in logarithms) on the surviving lengths | **Dead:** the surviving lengths have measure about 0.86, not a thin exceptional set. |
| Burrell-Yu style digit results | **Do not transfer** to palindromes ([PROOF_RESEARCH.md](PROOF_RESEARCH.md), 4b). |
| Decision diagrams, cross-branch learning (clause learning), automata | **Closed:** every node's 5-adic target is distinct (61,867 nodes checked); by Cobham's theorem, reading base 16 and base 10 together needs unbounded memory. |
| Coppersmith's method (degree 1 and 2), belief propagation, transfer operators, dissection, representation trees, dequantised backtracking | **Closed**, each with a measurement or a theorem (RESEARCH_PLAN.md, "Results so far"). |
| A new identity | **Proved and checked:** every hex palindrome N with L digits satisfies N = (-4)^(L-1) x (sum of its hex digits) mod 25 (WALL.md, section 3). It restates the 5-adic condition through digit moments; it does not break the obstacle. |

## 4. Why it stops here

A solution needs its size (leading digits), its 2-adic part (trailing bits) and its 5-adic part (trailing decimal digits
mod 5^K) to agree with mirrored digit patterns. Size and the 2-adic part can be matched for all lengths at once; the 5-adic
part is tied to the size only through the fact that N is an integer, which shows up only when almost every digit is known.
That is the same shape as Erdős's conjecture that 2^n has a digit 2 in base 3 for n > 8, and the "82000 problem" (numbers
whose base 2, 3, 4 and 5 digits are all 0 or 1), both open for decades. [WALL.md](WALL.md) states this precisely.

What would change the picture: an exact test that rejects most 5-adic targets at a depth that does not grow with the
length, or a reason (not a model) why the 5-adic targets of mirrored digit patterns avoid the reachable sets. Contributions
and challenges are welcome as GitHub issues.
