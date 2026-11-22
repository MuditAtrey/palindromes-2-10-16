# Research plan: a fundamentally faster way to find (or rule out) a base-2/10/16 palindrome

Goal: a method 100x-1000x faster than the current engine, on this Mac, no spend. The current engine is an outside-in tree
(cost about 2^(0.07-0.10 B)) whose only expensive part is a 5-adic reachability test (a generic-looking modular subset sum).

What is already settled by measurement (STATUS 38-41, notebook 5.12-5.13): the 2-adic half of the tree grows at exactly the
Furstenberg intersection dimension (Shmerkin-Wu), so it is optimal for the information it uses; the inner reachable sets and lattices
look random; the targets look random; more modular information does not help the enumeration; box-exact enumeration is a constant
factor; four reformulations conserve the difficulty. So a new method must come from somewhere outside the lattice / subset-sum frame,
or find structure that those tests could not see.

Each track below names: the field, the idea, why it might beat the barrier, a first experiment small enough to run here, and the
result that would kill it. Tracks are ordered by (chance of a big win) x (cheapness of the first test), not by field.

---

## A. Discrete mathematics and logic

**A1. Decision diagrams over carries (BDD / ZDD / automata).** The whole problem is a carry computation: binary digits in, decimal
digits out, both reversed. A zero-suppressed decision diagram builds the set of all partial solutions with sharing: two different
prefixes that reach the same carry state merge, which the tree never does. If carry states collide often (they are bounded by the
width of the undecided window), the diagram can be exponentially smaller than the tree.
First test: build the ZDD of base-2/16 palindromes whose last J decimal digits mirror the first J, for growing B, and compare its node
count with the tree's node count. Kill: the diagram is not smaller than the tree by a growing factor.

**A2. Conflict-driven clause learning (SAT).** CDCL learns clauses from dead branches and reuses them in other branches; the tree
never transfers anything between branches. If dead subtrees die for shared reasons (the same few carry bits), learning removes them
wholesale. First test: hand-written CDCL is out of scope; instead, measure whether dead nodes share small "explanations": for each
dead node, the smallest subset of decided digits that already makes the 5-adic target unreachable. If explanations are short and repeat,
learning pays. Kill: explanations are as long as the full prefix.

**A3. Cobham / Büchi arithmetic and automatic sequences.** Sets definable in base 2 and in base 10 at once are ultimately periodic
(Cobham); palindromes are not automatic, but the "mirror check of the last J digits" is. Look for a transducer formulation in which the
two reversals become a single synchronized two-tape automaton over (base 16, base 10) digit pairs read from both ends. First test:
literature and a small construction for base pairs (2, 4) where the answer is known. Kill: the state count grows like the tree.

## B. Analysis, dynamics and fractal geometry

**B1. Rotation structure of the leading digits.** The top decimal digits of N are a function of the hex prefix through the irrational
rotation by log10(16). Leading digits of a hex-digit-restricted set are a Sturmian-like coding. The 2-adic check then asks for a
match between a rotation coding and a 2-adic (odometer) coding. Hochman's inverse theorems and the Shmerkin-Wu method give not just a
dimension bound but a description of where intersections concentrate. First test: plot survivors of the 2-adic tree on the circle
(fractional part of log10 N) and in Z_2: look for clustering beyond the dimension count. Kill: survivors equidistribute.

**B2. Fourier decay of digit-restricted measures.** The {0,6,9,F} Cantor measure in base 16 has Fourier transform given by a Riesz
product; its decay along 5-adic characters controls how evenly the inner sums spread mod 5^K. If decay is slow at specific frequencies,
those frequencies are short dual vectors in disguise and give strong windows at shallower levels. First test: compute the Riesz product
at the characters e(k x / 5^K) for the inner digit set and find the largest values. Kill: flat spectrum (consistent with the measured
randomness).

**B3. Transfer operators and exact counting.** Counts of digit-restricted numbers with a modular condition are given by powers of a
transfer matrix. Exact counting cannot separate 0 from 1 on its own, but combined with the tree it can count solutions in a whole
subtree at once and prune subtrees whose count is provably 0 by a parity or modular argument. First test: transfer matrix for "hex
digits in {0,6,9,F}, value = t mod 5^j" and its spectrum. Kill: no subtree-level certificate cheaper than the reach tables.

## C. Number theory

**C1. p-adic analysis of N as a function of its digits.** 16^t = exp(t log 16) 5-adically; N mod 5^K is a polynomial in the digit
positions of degree about 0.67K (notebook section 3). Strassmann / Weierstrass preparation bounds the number of zeros of p-adic power
series. Look for a version that bounds or locates digit vectors in terms of their moments. First test: for small K, compute the variety
of moment vectors reachable by digit vectors and see if it has algebraic structure (lies on few hypersurfaces). Kill: generic.

**C2. Linear forms in logarithms and S-unit equations.** Bérczes-Ziegler used Baker's method for structured candidates. A palindrome
in base 10 that is a palindrome in base 16 gives N = a 10^n + rev10(a) = b 16^m + rev16(b): two-term S-unit-like equations when a and b
are short. Extend to "N has a long run of zeros in either base" families and prove them empty to huge sizes, shrinking the search space
where solutions could hide. First test: which families does the structured search (32.1) not cover that Baker handles cheaply. Kill:
families are a negligible fraction of the space.

**C3. Coppersmith / Howgrave-Graham small roots.** The inner part Z of N is a small unknown satisfying a modular equation
Z = t (mod 5^K) together with "Z has restricted hex digits". Coppersmith finds small roots of modular polynomials; the digit restriction
can be written as polynomial equations in the digit variables (x(x-6)(x-9)(x-15) = 0). Multivariate Coppersmith with lattice
reduction sometimes beats plain lattice enumeration when unknowns are tied by polynomial relations. First test: small instance, compare
with the CVP test. Kill: same or worse density threshold.

**C4. Checksums from cyclotomic structure.** 16 has order 4 mod 257 and 10 has order 4 mod 101; both give exact constraints on
alternating digit sums. Find all primes p where both 16 and 10 have small multiplicative order (p | 16^a - 1 and p | 10^b - 1 with small
a, b): each gives a global constraint that mixes both representations. First test: enumerate such p, measure the bits each contributes
inside the lattice test (like mod 11 did). Kill: total information small and the tests stay compute-limited.

## D. Algorithms and complexity

**D1. Dissection and multi-list algorithms (Dinur-Dunkelman-Keller-Shamir).** Dissection beats meet-in-the-middle on problems that
split into many layers with small interfaces. The palindrome problem is layered (digit pairs) with interfaces = carries. First test: model
the problem as a k-layer composition and compute the dissection exponent. Kill: exponent not below 0.071B.

**D2. Representation technique with structure (BCJ / BJMM on digit alphabets).** Generic representation algorithms give 2^(0.29 n) per
instance; with the {0,6,9,F} alphabet split as 6a+9b, representations multiply. Find whether many targets can share the top levels of the
representation tree. Kill: per-query cost above the tree's per-node cost at the information limit.

**D3. Belief propagation / survey propagation.** Message passing solves random constraint problems near their threshold far faster than
search. This problem is a sparse random-like system (digit variables, carry constraints). Use BP marginals to order branching or to prune
(only as a discovery tool unless a certificate follows). First test: BP on the carry factor graph for planted instances; does it find
the planted solution? Kill: BP does not converge or gives uniform marginals.

**D4. Quantum-inspired classical sampling.** Montanaro backtracking has no classical analog, but "dequantized" sampling sometimes
does. Low priority; literature check only.

## E. Learning and search

**E1. Learned certificates.** Train a model to predict which nodes die and why, then extract the learned rule as a provable test
(e.g. a short dual vector or a modular invariant). First test: correlate node death with simple features of the target and outer digits.
Kill: no feature beyond the known ones.

**E2. Program search for new exact tests (FunSearch style).** Generate candidate pruning predicates, keep only those that never kill a
planted or weak-control solution and are provably sound, select for speed.

## F. Reframing the goal

**F1. Probabilistic certification.** If the target is to find the number (not to prove its absence), a search that covers each length
with known probability p at much lower cost may reach far larger sizes. Requires recall measured with planted solutions.

**F2. Where would a solution be most likely?** Counting with the measured constant factors (the 1.3x low-residue excess, class sizes per
bit length, the B = 4L-1 class dominating) to rank lengths by expected number of solutions per unit of compute, and search those first.

---

## Results so far (2026-11-21)

- **A1 decision diagrams: closed.** Diagrams beat the tree only if branches reach the same state. Across 61,867 nodes at the levels that
  matter every 5-adic target is distinct, and apart from the 5-adic residue a node's future depends on the exact decimal remainder of its top
  part. The only merging available is by 5-adic residue, which is the engine's old reach tables.
- **C4 cyclotomic checksums: closed** (`scripts/structure/checksums.py`). Over all moduli below 20,000 and lengths 20-39, the only forced
  residues are 11 for even-length decimal palindromes (already used, section 38) and 17 for even-length hex palindromes (and 3 for even
  binary length). 17 and 3 constrain only the decimal digits, so they would help a decimal-side tree, not the hex tree.
- **B2 Fourier decay: closed** (`scripts/structure/fourier_decay.py`). Over every frequency mod 5^9 and 5^10, the inner digit measure has
  the same spectrum as random weights (largest coefficients 0.17 vs 0.12 and 0.08 vs 0.10; 752 vs 726 and 168 vs 160 frequencies above 0.02).
  No hidden short dual vectors.
- **C3 Coppersmith: closed** (`scripts/formula/coppersmith.py`). Lattices of polynomials vanishing mod 5^K at every root, kept when the sum of
  their coefficients is below 5^K (then they vanish over the integers). 12-bit instances, unreachable random targets: degree 1 certifies 9/12 at
  density 0.65 and none above (the window tests again); degree 2, which adds the bit relations x^2 = x, finds more usable polynomials
  (3.6 vs 0.9 per target) but certifies 0 at every density. The polynomial structure adds no pruning power.
- **D1 dissection: closed (by the known result).** Dissection algorithms (Dinur, Dunkelman, Keller, Shamir) improve the memory of layered
  meet-in-the-middle problems, not the time below meet-in-the-middle; memory is not the constraint here, and the layers are linked through
  the global 5-adic condition rather than small interfaces.
- **D3 belief propagation: closed** (`scripts/structure/bp_inner.py`). The inner test as base-5 column additions with carries, exact factor
  messages, planted instances: marginals predict 50.0% of planted bits at 12 bits (density 0.65) and 61.5% at 16 bits (density 0.69);
  BP-guided decimation solves 0/10 and 1/6. A faint signal, never a certificate, too weak to steer the search at real sizes.
- **C2 Baker families: closed (by counting).** Families with few nonzero digits or few blocks hold polynomially many candidates per length
  against about 2^L in all; the structured search already covered ~2^29 of them to 10^5000. Baker's method could make those exclusions
  permanent, but it cannot reach the generic candidates, so it cannot speed up finding a solution.
- **B3 transfer operators: closed (reduces to reach tables).** An exact subtree count needs the full 5-adic residue as state; with a small
  modulus 5^J the count is zero only on small subtrees, which is the engine's former 5-adic reach filter.
- **E2 program search for exact tests: closed (subsumed).** A discovered test must be provably sound; the provable families are linear windows,
  modular invariants, lattices and polynomial certificates, all measured in tracks above or in STATUS 27-41.
- **F2 where to search: done.** Random model (validated on the 212 base-2/10 double palindromes and 25 sibling bases), measured cost per
  length: expected solutions per hour of compute fall monotonically with length (about 3e-27 at Hex 100, 1e-35 at Hex 120, 7e-40 at Hex 130,
  2e-48 at Hex 150), so searching lengths in increasing order is already the best order. Expected number beyond 2^472 under the model: 1.6e-35.
- **A2 cross-branch learning: closed (by measurement).** At the level where the exact test fires, 99.8% of nodes die (STATUS 39-40 depth
  histograms), and a node's death depends on its 5-adic target, which depends on every decided digit and the decimal prefix; the targets are
  all distinct (A1). A learned clause would have to name the whole prefix, so nothing transfers between branches.
- **A3 automata / Cobham: closed (by theorem).** Reading base-16 and base-10 digits in step from both ends needs an unbounded buffer, because
  digit boundaries in multiplicatively independent bases never realign (Cobham); the buffer is the carry, i.e. the decimal remainder of A1.
- **B1 clustering of survivors: closed (by measurement).** Survivor counts already equal the Furstenberg intersection bound (notebook 5.12), so
  clustering cannot reduce how many survivors must be visited; it could only help if it correlated with 5-adic reachability, and the
  decimal-side test (STATUS 39) shows targets fall in the reachable set at exactly the random rate.
- **C1 p-adic moment structure: closed (by measurement).** The inner sum mod 5^K is a linear function of the digit moments; the reachable
  sets have the size and lattice shape of random weights (STATUS 39) and no Fourier structure (B2).
- **D2 shared representation trees: closed (by derivation).** Representation algorithms cost about 2^(0.29 u) per instance; at the information
  limit (u ~ 0.58 L) that is 2^(0.17 L) per node on top of 2^(0.17 L) nodes, 2^(0.084 B) in all, worse than the engine.
- **D4 dequantized backtracking: closed (by theorem).** Montanaro's quadratic speedup comes from quantum walks; the classical analogues
  (tree-size estimation) estimate how big a tree is, they do not decide whether it holds a solution.
- **E1 learned certificates: closed (by earlier measurement).** Neural predictors on 1.2M near-palindromes found no signal beyond the 2-adic
  rule (STATUS 28).
- **F1 probabilistic coverage: closed (by earlier measurement).** Planted-solution calibration (STATUS 36-37): pruned runs find solutions at
  0.55-0.9x the rate of exact runs per unit of compute; no gain.

## Status

All 17 tracks are closed. None gives a method beyond constant factors. The constant factors that are real and measured: mod 11 and extra
prefix digits (1.2-1.5x on the integer-only certificate, built), box-exact bounds in the enumeration (estimated 1.5-3x, not built), and the
GPU (3.2x on the core arithmetic, not built).

## Order of work

1. A1 decision diagrams over carries (could change the exponent; cheap to test).
2. A2 explanations of dead nodes (tests whether cross-branch learning could pay).
3. B1 survivors on the circle and in Z_2 (tests for clustering the dimension count cannot see).
4. C4 cyclotomic checksums (cheap, extends the mod 11 result).
5. B2 Fourier decay of the digit measure along 5-adic characters.
6. D1 dissection exponent.
7. C3 Coppersmith on a small instance.
8. D3 belief propagation on planted instances.
9. The rest as earlier results suggest.

Every track records its result, positive or negative, in STATUS.md.

## Wide sibling sweep (2026-11-22): looking for a mechanism

The only experiment that could overturn the random model: if some mechanism makes numbers palindromic in three bases more often than
chance, it should show up in problems built the same way. `scripts/siblings/world.c` + `world_sweep.py` search 372 worlds exhaustively:
palindromes in base 2, in r (16 with digits {0,6,9,F}, 16 with {0,2,5,7}, 8 with {0,2,5,7}, 4 with {0,3}) and in every even b from 6 to 200
that is not a power of 2, each length complete until one length takes 20 s (or 600 bits). 43.6 billion nodes, 56 minutes on 8 cores.
Checks first: brute force over all binary palindromes up to 26 bits (297 cases, 110 solutions, all match), b = r positive controls (no
candidate pruned), identical node and hit counts to the older sibling.c on the real class.

Random model, corrected for the one exact divisibility link (an even-length base-b palindrome is divisible by b + 1): calibrated at 12-31 bits,
where solutions are plentiful, it tracks the counts in every family.

Result (`results/world_sweep_report.txt`): **20 solutions with 32 bits or more, against 20.4 expected** (hex 3 / 1.6, hex odd 0 / 1.25, octal
13 / 15.2, base 4 4 / 2.3). The largest is 75 bits (octal world, b = 146), with chance 0.29 of something at least that surprising
somewhere among 372 worlds. No large solutions share a pattern. The base-54 oddity 0xF6F99F6F is ordinary at this scale.

Reading: no sign of a mechanism in any of 372 sibling problems, including the base-4 version and the two classes of the real problem with other
third bases. The random model predicts these worlds exactly, which is the strongest evidence so far that it also holds for b = 10, where it
predicts 1.6e-35 solutions beyond the searched 2^472.
