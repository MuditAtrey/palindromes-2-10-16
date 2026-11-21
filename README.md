# Numbers palindromic in bases 2, 10 and 16

**Result.** Apart from 0, 1, 3, 5, 7 and 9, there is no integer below **2^472** (about 1.22 x 10^142) that reads the same
forwards and backwards in binary, decimal and hexadecimal at the same time.

This is an exhaustive computer search: every one of the roughly 2^118 candidate numbers below 2^472 that is a palindrome in
base 2 and base 16 was either checked directly or ruled out by an exact mathematical test. No probabilistic or heuristic
step is involved, and every pruning test in the certified configuration uses integer arithmetic only.

Whether any larger such number exists is open.

## What exactly was done

| Range | Lengths | Configuration | Nodes visited | Hits |
|---|---|---|---|---|
| N < 2^460 | 1 to 115 hex digits | 2-adic test + exact 5-adic reachability tables | 2.01 x 10^12 | 1, 3, 5, 7, 9 only |
| 2^460 <= N < 2^472 | 116 to 118 hex digits | 2-adic test + dual-lattice interval certificate | 6.72 x 10^10 | none |
| **Total** | **1 to 118** | **integer-only throughout** | **2.07 x 10^12** | **1, 3, 5, 7, 9 only** |

Every length, with its node count, run time, configuration and hits, is in [results/per_length.csv](results/per_length.csv).
The search ran on a single 10-core Apple M4 machine. The last length (118 hex digits, 52.5 billion nodes) took 13 minutes.

## How I know the search would have found a solution

A search that finds nothing is only meaningful if it is shown to find solutions where they exist. The same engine and
configuration were checked against:

1. **Brute force and known sequences.** Exact agreement with brute force below 2^32, and exact reproduction of every term of
   OEIS [A007632](https://oeis.org/A007632) (palindromes in bases 2 and 10) below 2^48 when the base-16 condition is dropped.
2. **Weak palindromes.** Numbers whose last J decimal digits mirror their first J exist in large numbers; every pruning test
   finds exactly the same set as the unpruned search (J = 3 to 9).
3. **Planted solutions at full size.** The engine can search for N such that N + c is a decimal palindrome. With c = 0 this is
   the real problem (identical results). Choosing a candidate N0 and a decimal palindrome P0 and setting c = P0 - N0 plants a
   guaranteed solution, which the search must report. Results: 91 of 91 plants found with the fast configuration at 16 to 100 hex digits (numbers up to 400 bits), and 7 of 7 with the exact integer-only configuration used for 116 to 118 hex digits, at 100, 105, 110, 114, 116, 117 and 118 hex digits.
4. **Self-checks** that recompute every node independently, verify that genuine subset sums are never rejected by the
   certificates, and compare runs across thread counts and work splits. The engine's test suite has 47 tests.

## Proof of priority (commit-reveal)

The engine source, its tests and all result logs are fingerprinted (SHA-256) in [commitment/manifest.txt](commitment/manifest.txt),
and that file is timestamped with [OpenTimestamps](https://opentimestamps.org) in
[commitment/manifest.txt.ots](commitment/manifest.txt.ots), anchored in the Bitcoin blockchain. This proves the exact code and
results existed on the stated date and cannot be changed afterwards. See [COMMITMENT.md](COMMITMENT.md) for how to verify.

The source code is not public yet and will be released later; it will hash to the values in the manifest.
Researchers who want to check the search before then can ask for access.

## Method

Candidates are built digit by digit from both ends, keeping all three bases consistent and abandoning a branch as soon as
no completion can work. The difficult part is that 16 shares no factor with 5, so half of the decimal information can only be
used through a reachability test over all undecided digits. [METHOD.md](METHOD.md) describes each test and why it is sound.

## Other findings

- **The random model fits.** The 212 known base-2/10 double palindromes match a model where decimal and binary digits are
  independent (212 observed against 206.5 expected). Under that model the expected number of base-2/10/16 triple palindromes
  above 2^472 is far below one.
- **Sibling problems behave the same way.** For 25 other bases b (even, from 6 to 60), numbers palindromic in bases 2, 16 and b
  were searched exhaustively to between 144 and 320 bits. In every case the largest solution has at most 32 bits.
- **No sparse family can work.** By a theorem of Senge and Straus (made effective by Stewart), only finitely many numbers have
  few nonzero digits in two multiplicatively independent bases, so no formula built from sparse digit patterns produces
  infinitely many solutions. A search of structured families (two-block, periodic, sparse and windowed digit patterns, 1.5 billion
  candidates) found none below 10^5000.

## Prior work

- Bérczes and Ziegler, *On simultaneous palindromes* (2014), study numbers palindromic in two bases and list the base-2/10
  double palindromes below 10^18.
- [mtschoen/multibase-palindrome-search](https://github.com/mtschoen/multibase-palindrome-search) reports the same problem
  searched to 2^255 (November 2026).
- Base-2/10 double palindromes are listed in OEIS A007632; base-10/16 in A029731.

## Contact

Questions, scrutiny and requests for access are welcome through GitHub issues on this repository.
