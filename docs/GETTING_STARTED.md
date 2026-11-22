# Getting started

This page is for someone arriving fresh: what the problem is, how to set up, what to run first, what to read in what order,
the vocabulary, and how to contribute without breaking the one thing that matters (never missing a solution).

## 1. The problem in two minutes

Find an integer N > 9 that reads the same forwards and backwards in base 2, base 10 and base 16. The only known examples
are 0, 1, 3, 5, 7 and 9. Rick Regan asked publicly in 2010 whether one exists
([Exploring Binary](https://www.exploringbinary.com/in-search-of-decimal-binary-hexadecimal-palindromes/)); Matt Schoen's
search reached 2^255 ([mtschoen/multibase-palindrome-search](https://github.com/mtschoen/multibase-palindrome-search)); this
repository reaches 2^472. Nobody knows whether a larger one exists, and no proof either way is known.

Why it is hard, in one paragraph: being a palindrome in bases 2 and 16 is easy to arrange (16 is a power of 2), and it
already forces the hex digits to come from {0, 6, 9, F} (or {0, 2, 5, 7} for bit lengths one less than a multiple of 4).
Base 10 = 2 x 5 shares the 2 with them but not the 5, and every undecided hex digit changes N modulo every power of 5. So
the decimal condition can only be used through a test over all undecided digits at once. [METHOD.md](../METHOD.md) has the
details; [WALL.md](WALL.md) explains why this blocks proofs too.

## 2. Setup

You need a C compiler, make, Python 3, and FLINT (which brings GMP and MPFR).

```bash
# macOS (Homebrew)
brew install flint
# Debian / Ubuntu
sudo apt-get install -y build-essential python3 libflint-dev libgmp-dev libmpfr-dev
```

The research scripts in `scripts/` also use NumPy, SciPy and mpmath:

```bash
python3 -m pip install numpy scipy mpmath
```

On a fresh Linux machine, `bash scripts/cloud/setup_linux.sh` installs everything, builds the engine and runs the checks.

## 3. First runs (about ten minutes in total)

```bash
cd engine
make                                  # builds ./tripal
make test                             # 47 tests, about 5 minutes on a 10-core machine; all must pass
./tripal verify 585 9 0x6996          # checks numbers in bases 2, 10 and 16
./tripal hexdfs --hex-from 1 --hex-to 40 --threads 10      # every length below 2^160: about 4 s, 0.7 GB
```

The `hexdfs` output has one row per bit length (nodes visited, prunes, time) and a `HIT n` line for every solution; you
should see exactly 1, 3, 5, 7 and 9. Then try the positive control that shows the search can find a solution when one
exists:

```bash
cd ..
python3 scripts/discovery/plant.py control 30 31 2     # plants 4 solutions at 120-124 bits; all must be FOUND
```

## 4. Reproducing the record

| Range | Command | Cost on a 10-core Apple M4 |
|---|---|---|
| lengths 1 to 115 hex digits | `python3 scripts/frontier_sweep_v18.py --from 1 --to 115 --threads 10` | about 8.2 hours in total; the longest length (115 hex digits) took 2.2 hours. Memory 1.1 GB (5^12 tables) or 4.7 GB (5^13, chosen automatically with 12 GB of RAM or more) |
| 116 to 118 hex digits, integer-only certificate | `cd scripts/cloud && bash certify_rerun.sh` | 54 s, 168 s and 13 minutes |

[results/per_length.csv](../results/per_length.csv) lists the node count of every length; your counts should match exactly
(node counts are deterministic). To split one length across machines, run `python3 scripts/cloud/run_shard.py HEX HEX K N` on
machine K of N and merge with `merge_shards.py`.

To check that this release is the code that was timestamped before publication, see [COMMITMENT.md](../COMMITMENT.md).

## 5. Reading order

1. [README.md](../README.md): the result and how it was checked.
2. [METHOD.md](../METHOD.md): the search and why each pruning step is sound.
3. [engine/README.md](../engine/README.md), then the long comment at the top of `engine/src/search.c`: how the code is laid out.
4. [WHAT_I_TRIED.md](WHAT_I_TRIED.md): what worked, what failed, and why. Read this before trying an idea.
5. [WALL.md](WALL.md): the obstacle for any proof covering every length.
6. [THEOREM_1.md](THEOREM_1.md) and [THEOREM_2.md](THEOREM_2.md): the two proved partial results.
7. [RESEARCH_PLAN.md](RESEARCH_PLAN.md) and [PROOF_RESEARCH.md](PROOF_RESEARCH.md): the full research log.
8. [CREDITS.md](../CREDITS.md): the papers and prior work, with links.

## 6. Glossary

| Term | Meaning |
|---|---|
| B, L | bit length of N, and hex length L = ceil(B / 4). "Hex 118" means bit lengths 469 to 472. |
| site, orbit | a site is one bit position; bits tied together by the binary and hex mirrors form an orbit, which is one yes/no decision |
| classes | B = 4L (hex digits {0,6,9,F}), B = 4L - 1 (digits {0,2,5,7}); B = 4L - 2 and 4L - 3 have far fewer candidates |
| node | one partial assignment in the search tree; `hexdfs` decides one mirrored pair of hex digits per level |
| lo, hi | the smallest and largest number a node can still complete to |
| 2-adic test | lo and hi share their first m decimal digits, so N must end with them reversed; compare with N's known low bits |
| 5-adic target | the same mirrored digits fix N modulo 5^m; the target is the residue the undecided digits must supply |
| reach set, R5 tables | all residues modulo 5^J that the undecided digits can produce (subset sums); a target outside it kills the node |
| dual-lattice certificate | an exact integer test that a target is unreachable, using short vectors of a lattice (no tables) |
| weak palindrome | a number whose last J decimal digits mirror its first J; thousands exist, so they test the pruning (`TP_WEAKPAL=J`) |
| planted solution | `TP_OFFSET=c` makes N0 a guaranteed solution; the search must report it |
| random model | decimal and binary digits treated as independent; it predicts every count I could check |

## 7. Contributing

Ways to help, roughly from easiest to hardest:

- **Reproduce** a range on your hardware and report node counts and times (they must match exactly).
- **Break it**: find a decimal palindrome the pruning would reject, or a flaw in the enumeration argument. Matt Schoen's open
  [challenge issues](https://github.com/mtschoen/multibase-palindrome-search/issues) apply here too.
- **Speed**: the GPU port measured 3.2x on the core arithmetic and is not built; the lattice test still uses floating point
  to steer and could be made exact.
- **Theory**: an exact test that rejects most 5-adic targets at a depth that does not grow with the length would change the
  picture (WALL.md, section 4).

The rule for any change to the search: a new pruning test must leave the hit set of the weak-palindrome controls unchanged
and must find every planted solution. Run `make test` and `python3 scripts/discovery/plant.py control ...` before and after;
the test suite already includes a deliberately unsound pruner (`TP_BADPRUNE`) to show that these checks catch one.

Questions and results are welcome as GitHub issues.
