# Method

This describes how the search works and why each pruning step is sound. The implementation is in [engine/](engine)
(see its [README](engine/README.md)), and matches the cryptographic commitment described in [COMMITMENT.md](COMMITMENT.md).

## 1. Candidates

A number that is a palindrome in base 2 is odd. If its bit length B is a multiple of 4, it is also a palindrome in base 16
exactly when every hex digit is its own 4-bit reversal, that is one of 0, 6, 9 or F, and the hex digits are mirrored.
So a length with L hex digits has about 4^(L/2) = 2^(B/4) candidates. Bit lengths that are not multiples of 4 force
stronger coupling between the binary and hex mirrors and leave fewer candidates. All four classes are searched.

## 2. The search tree

Hex digits are chosen in mirrored pairs from the outside in: the leading digit together with the last digit, then the
next pair, and so on. After d pairs are fixed:

- the low 4d bits of N are known exactly (the undecided digits only contribute multiples of 16^d), and
- N lies in a known interval [lo, hi], whose endpoints share their first m decimal digits.

If N is a decimal palindrome, its last m decimal digits are those first m digits reversed. That fixes N modulo 10^m,
which splits into a condition modulo 2^m and a condition modulo 5^m.

## 3. Pruning, every rule exact

| Test | What it checks | Why it is sound |
|---|---|---|
| 2-adic | N mod 2^min(m, 4d) from the mirrored decimal digits must equal the known low bits | both sides are exact functions of the fixed digits |
| 5-adic reach tables | the undecided hex digits must be able to produce the required residue mod 5^J | exhaustive table of reachable residues for each depth (a subset-sum table) |
| Dual-lattice interval certificate | for an integer lambda, every reachable residue r satisfies lambda*r mod 5^K inside a computable interval; outside it, no completion exists | one line of integer arithmetic; integer-only, no floating point |
| Lattice enumeration (used only for speed, not for the certified run) | complete enumeration of the coset of the reachability lattice inside the ball containing the digit box | exact in principle; its empty-search verdict relies on floating-point margins, so it is excluded from the certified configuration |
| Leaf check | full verification in bases 2, 10 and 16 | direct |

The certified configuration uses only the 2-adic test, the dual-lattice certificate, and the leaf check
(integer arithmetic throughout). The 5-adic tables were used for lengths 1 to 115.

## 4. How the search is checked

A pruning rule that is wrong would hide a solution, and a search that finds nothing proves nothing unless it is shown to
find solutions where they exist. Four kinds of positive control are used:

1. **Brute force and known sequences.** Exhaustive agreement with brute force below 2^32, and exact reproduction of all terms
   of OEIS A007632 (palindromes in bases 2 and 10) below 2^48 when the base-16 condition is dropped.
2. **Weak palindromes.** Replace "decimal palindrome" with "the last J decimal digits mirror the first J". Thousands of such
   numbers exist; every pruning rule must find exactly the same set as the unpruned search (tested for J = 3 to 9).
3. **Planted solutions at full scale.** Search for N such that N + c is a decimal palindrome. For c = 0 this is the real problem
   (identical node counts). For a planted instance, pick a candidate N0 and a decimal palindrome P0 and set c = P0 - N0; the
   search must report P0. Plants are found at every length tested, up to and including the certified configuration.
4. **Internal self-checks.** Optional modes recompute every node the slow way, verify that genuine subset sums are never
   rejected by the certificates, and compare thread counts and work splits.

## 5. Why it is hard

Base 16 is a power of 2, so the binary and hexadecimal conditions are compatible and local. Base 10 = 2 x 5 shares the 2 but
not the 5, and 16 is a unit modulo every power of 5: every undecided hex digit changes N mod 5^k. The 5-adic half of the
decimal information can therefore only be used through a reachability test over all undecided digits, which is where all
the effort goes. The cost grows roughly like 2^(0.07 to 0.10 per bit) of the target size.
