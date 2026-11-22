# tripal: the search engine

`tripal` searches for integers that are palindromes in bases 2, 10 and 16 at the same time. It is the engine behind the
result in the [top-level README](../README.md): nothing above 9 exists below 2^472.

The C sources here are the exact code that produced the certified runs (tag `v1.0-as-committed` holds them byte for byte;
later commits only add comments, and the compiled object files are identical before and after).

## Build

Needs a C11 compiler and, for the dual-lattice certificate used above 2^460, [FLINT](https://flintlib.org) (with GMP and
MPFR). The Makefile finds FLINT in `/opt/homebrew`, `/usr` or `/usr/local`; without it everything except `--dual` / `TP_DUAL`
and `TP_CVP` still builds.

```bash
brew install flint          # macOS; on Debian/Ubuntu: apt install libflint-dev
make                        # portable -O3 build, produces ./tripal
make native                 # tuned for the host CPU
```

## Check it before trusting it

```bash
make test                   # full suite: 47 tests, about 5 minutes on a 10-core M4
./tripal selftest           # brute force, OEIS A007632 ground truth, independent reference enumerator
./tripal verify 585 0x249   # check any number in bases 2, 10 and 16
```

The suite includes the positive controls that matter for a search whose answer is "nothing": weak palindromes (last J
decimal digits mirror the first J, thousands exist) must be found identically with every pruning rule on and off, and
planted solutions at full size must be found. See section 4 of [METHOD.md](../METHOD.md).

## Run

```bash
./tripal hexdfs --hex-from 1 --hex-to 40 --threads 10          # every bit length 1..160, 5-adic tables (J=12)
TP_DUAL=1 ./tripal hexdfs --hex-from 116 --hex-to 116 --r5 0 --threads 10   # the certified configuration above 2^460
```

Each line of output is one bit length: nodes visited, time, and `HIT n` lines for every solution (only 1, 3, 5, 7, 9 exist
below 2^472). `scripts/frontier_sweep_v18.py` runs a range and writes one JSON line per hex length;
`scripts/cloud/` splits a length into shards for several machines (`--shard K/N`) and merges the results.

Other commands: `search` (the original bit-by-bit orbit DFS, kept for cross-checks), `bench`, `estimate`, `orbits`,
`ref`, `brute`. Run `./tripal` with no arguments for the list.

## Switches

Search features are switched with environment variables, read once per bit length in `hexplan_prepare` (src/search.c).

| Variable | Effect |
|---|---|
| `--r5 J` / `TP_R5=J` | 5-adic reach tables modulo 5^J (default 12, or 13 on machines with 12 GB or more; 0 = off) |
| `TP_DUAL=1` | dual-lattice window certificate (exact, integer arithmetic; used for the certified 116-118 hex digit runs) |
| `TP_CVP=1` | exact lattice enumeration test (fast; uses floating point to steer, so not used for certification) |
| `--bloom U:b` / `TP_BLOOM` | deep Bloom reach filter at U undecided orbits (no false negatives) |
| `TP_SIEVE=1` | leading-two-digit / digit-sum sieve |
| `TP_NOMOD11=1`, `TP_NOFASTK=1`, `TP_NOCARRY9=1` | turn off the mod 11 leaf sieve, the chunk kernel, the carried decimal |

Test-only switches (they exist to prove the pruning sound, never for real runs):

| Variable | Effect |
|---|---|
| `TP_WEAKPAL=J` | accept N whose last J decimal digits mirror the first J (thousands of positive cases) |
| `TP_OFFSET=c` | search for N with N + c a decimal palindrome; c = P0 - N0 plants a known solution |
| `TP_DECCHK`, `TP_R5CHK`, `TP_DUALCHK`, `TP_CVPCHK`, `TP_BLOOMCHK`, `TP_FASTCHK`, `TP_XM` | recompute values the slow way at every node, or sample real completions, and abort on any disagreement |
| `TP_NOPRUNE=1` | no pruning at all (reference node counts) |
| `TP_BADPRUNE=1` | a deliberately unsound pruner, to show the tests catch one |

## Files

| File | What it does |
|---|---|
| `src/orbits.c` | which bits are tied together by the binary and hex mirrors; one free decision per orbit |
| `src/search.c` | the search: `hexdfs` (production kernel), the 2-adic test, the 5-adic reach tables, the lattice certificates, threading |
| `src/main.c` | command line, self-test, verify, brute-force referee |
| `src/bigint.c` | fixed-width integers and base conversion |
| `src/checkpoint.c` | crash-safe logs for long `search` runs |
| `include/tripal.h` | shared types; `include/viable_outer.h` is a generated table of possible outer digit pairs |
| `tests/run_tests.py` | the 47-test suite (`make test`) |
| `verify/check.py` | independent Python check of a result file |
| `data/oeis/b007632.txt` | OEIS A007632, palindromes in bases 2 and 10, the ground truth for the self-test |

Each C function has a comment saying what it does; the long comment at the top of `src/search.c` explains the whole
pruning argument and how the file is organised.
