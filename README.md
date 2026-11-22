# Numbers palindromic in bases 2, 10 and 16

**Result.** Apart from 0, 1, 3, 5, 7 and 9, there is no integer below **2^472** (about 1.22 x 10^142) that reads the same
forwards and backwards in binary, decimal and hexadecimal at the same time.

This is an exhaustive computer search: every one of the roughly 9 x 10^35 (about 2^119) candidate numbers below 2^472 that is a palindrome in
base 2 and base 16 was either checked directly or ruled out by an exact mathematical test. No probabilistic or heuristic
step is involved, and every pruning test in the certified configuration uses integer arithmetic only.

Whether any larger such number exists is open. This repository has the engine, the logs, the checks, and an honest account
of every approach that did not settle it.

## Contents

| Where | What |
|---|---|
| [engine/](engine) | `tripal`, the search engine (C), its 47-test suite, and how to build and run it |
| [results/](results) | every search log; [per_length.csv](results/per_length.csv) has one row per hex length |
| [METHOD.md](METHOD.md) | how the search works and why each pruning step is sound |
| [docs/WHAT_I_TRIED.md](docs/WHAT_I_TRIED.md) | everything I tried, what worked, what failed, and why |
| [docs/WALL.md](docs/WALL.md) | the obstacle any proof for every length must get past |
| [docs/THEOREM_1.md](docs/THEOREM_1.md), [docs/THEOREM_2.md](docs/THEOREM_2.md) | two proved partial results (rarity; empty lengths) |
| [docs/RESEARCH_PLAN.md](docs/RESEARCH_PLAN.md), [docs/PROOF_RESEARCH.md](docs/PROOF_RESEARCH.md) | the research log behind them |
| [scripts/](scripts) | sweeps, sharding for several machines, planted controls, and the research scripts the documents cite |
| [commitment/](commitment), [COMMITMENT.md](COMMITMENT.md), [verify/](verify) | the timestamped commitment and a script that checks this release against it |
| [docs/GETTING_STARTED.md](docs/GETTING_STARTED.md) | start here if you are new: setup, first runs, reading order, glossary, how to contribute |
| [CREDITS.md](CREDITS.md) | everyone whose work this builds on, with links |

## What exactly was done

| Range | Lengths | Configuration | Nodes visited | Hits |
|---|---|---|---|---|
| N < 2^460 | 1 to 115 hex digits | 2-adic test + exact 5-adic reachability tables | 2.01 x 10^12 | 1, 3, 5, 7, 9 only |
| 2^460 <= N < 2^472 | 116 to 118 hex digits | 2-adic test + dual-lattice interval certificate | 6.72 x 10^10 | none |
| **Total** | **1 to 118** | **integer-only throughout** | **2.07 x 10^12** | **1, 3, 5, 7, 9 only** |

The search ran on a single 10-core Apple M4 machine. The last length (118 hex digits, 52.5 billion nodes) took 13 minutes.

## Try it

```bash
cd engine && make && make test        # build, then the 47-test suite (about 5 minutes)
./tripal verify 585                   # 585 is a palindrome in bases 2 and 10, not 16
./tripal hexdfs --hex-from 1 --hex-to 40 --threads 10      # re-run every length below 2^160
```

[engine/README.md](engine/README.md) has the certified configuration, every switch, and how to split a length across
machines.

## How I know the search would have found a solution

A search that finds nothing is only meaningful if it is shown to find solutions where they exist. The same engine and
configuration were checked against:

1. **Brute force and known sequences.** Exact agreement with brute force below 2^32, and exact reproduction of every term of
   OEIS [A007632](https://oeis.org/A007632) (palindromes in bases 2 and 10) below 2^48 when the base-16 condition is dropped.
2. **Weak palindromes.** Numbers whose last J decimal digits mirror their first J exist in large numbers; every pruning test
   finds exactly the same set as the unpruned search (J = 3 to 9).
3. **Planted solutions at full size.** The engine can search for N such that N + c is a decimal palindrome. With c = 0 this is
   the real problem (identical results). Choosing a candidate N0 and a decimal palindrome P0 and setting c = P0 - N0 plants a
   guaranteed solution, which the search must report. Results: 91 of 91 plants found with the fast configuration at 16 to
   100 hex digits (numbers up to 400 bits), and 7 of 7 with the exact integer-only configuration used for 116 to 118 hex
   digits ([results/planted_controls_integer_only_config.txt](results/planted_controls_integer_only_config.txt)).
4. **Self-checks** that recompute every node independently, verify that genuine subset sums are never rejected by the
   certificates, and compare runs across thread counts and work splits.

## Proof of priority

The engine source, its tests and all result logs were fingerprinted (SHA-256) in [commitment/manifest.txt](commitment/manifest.txt)
and timestamped with [OpenTimestamps](https://opentimestamps.org) on 2026-11-21, before the code was public. The released
files match it byte for byte:

```bash
git checkout v1.0-as-committed && python3 verify/verify_commitment.py
```

Later commits only add comments and documentation; the compiled engine is identical. See [COMMITMENT.md](COMMITMENT.md).

## Other findings

- **The random model fits.** The 212 known base-2/10 double palindromes match a model where decimal and binary digits are
  independent (212 observed against 206.5 expected). Under that model the expected number of base-2/10/16 triple palindromes
  above 2^472 is far below one.
- **Sibling problems behave the same way.** Across 372 related problems (palindromes in base 2, a power of 2, and an even base
  b from 6 to 200), large solutions occur exactly as often as the random model predicts (20 observed, 20.4 expected).
- **No sparse family can work.** By a theorem of Senge and Straus (made effective by Stewart), only finitely many numbers have
  few nonzero digits in two multiplicatively independent bases, so no formula built from sparse digit patterns produces
  infinitely many solutions. A search of structured families (1.5 billion candidates) found none below 10^5000.
- **A proof for every length is out of reach of known methods.** Every argument that looks only at a bounded number of outer
  digits covers at most 13.65% of lengths; the rest needs control of a random-looking 5-adic event, the same shape as
  Erdős's ternary-digit problem. See [docs/WALL.md](docs/WALL.md).

## Prior and related work

- Bérczes and Ziegler, *On simultaneous palindromes* (2014), study numbers palindromic in two bases and list the base-2/10
  double palindromes below 10^18.
- [mtschoen/multibase-palindrome-search](https://github.com/mtschoen/multibase-palindrome-search) searched the same problem
  to 2^255 (November 2026); [nevdelap/pp](https://github.com/nevdelap/pp) searched the prime version to 256 bits.
- OEIS: [A007632](https://oeis.org/A007632) (bases 2 and 10), [A029731](https://oeis.org/A029731) (bases 10 and 16),
  [A340559](https://oeis.org/A340559) (bases 2 and 16).
- The question was posed publicly by Rick Regan, [In Search of Decimal/Binary/Hexadecimal Palindromes](https://www.exploringbinary.com/in-search-of-decimal-binary-hexadecimal-palindromes/)
  (Exploring Binary, 2010).

The full list of people, papers, data and software this builds on, with links, is in [CREDITS.md](CREDITS.md).

## Citing and license

If you use the engine, the results or the methods, please cite this repository; [CITATION.cff](CITATION.cff) has the details
(GitHub shows a "Cite this repository" button). Code is under the [MIT License](LICENSE); documents, results and data are
under [CC BY 4.0](LICENSE-docs), which requires credit. The bundled OEIS data file is CC BY-SA 4.0 (see CREDITS.md).

Questions, scrutiny and challenges are welcome as GitHub issues.
