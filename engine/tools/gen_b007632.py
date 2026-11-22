#!/usr/bin/env python3
"""Regenerate / verify the OEIS A007632 ground truth without network access.

Enumerates decimal palindromes by their leading half, keeps the ones that are
also base-2 palindromes, and compares against data/oeis/b007632.txt.

    python3 tools/gen_b007632.py --digits 13          print the terms found
    python3 tools/gen_b007632.py --digits 13 --verify  check them against the file

Terms with more than ~13 decimal digits are expensive to reach this way; the
checked-in file gets those from the OEIS b-file, which is complete below 10^64.
"""

from __future__ import annotations

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def decimal_palindromes(max_digits: int):
    """Yield every decimal palindrome with at most max_digits digits, ascending."""
    yield 0
    for d in range(1, max_digits + 1):
        half = (d + 1) // 2
        for lead in range(10 ** (half - 1), 10 ** half):
            s = str(lead)
            if d % 2 == 0:
                n = int(s + s[::-1])
            else:
                n = int(s + s[-2::-1])
            yield n


def is_bin_pal(n: int) -> bool:
    b = bin(n)[2:]
    return b == b[::-1]


def load_file(path: str) -> list[int]:
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            out.append(int(line.split()[1]))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--digits", type=int, default=13)
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--file", default=os.path.join(ROOT, "data", "oeis", "b007632.txt"))
    a = ap.parse_args()

    found = [n for n in decimal_palindromes(a.digits) if is_bin_pal(n)]
    print(f"decimal palindromes with <= {a.digits} digits that are also base-2 "
          f"palindromes: {len(found)}", file=sys.stderr)

    for n in found:
        print(n)

    if a.verify:
        ref = load_file(a.file)
        limit = 10 ** a.digits
        ref_in_range = [n for n in ref if n < limit]
        ok = ref_in_range == found
        print(f"\nfile has {len(ref)} terms, {len(ref_in_range)} below 1e{a.digits}; "
              f"generated {len(found)}", file=sys.stderr)
        print("MATCH" if ok else "MISMATCH", file=sys.stderr)
        if not ok:
            print(f"  only in file: {sorted(set(ref_in_range) - set(found))[:5]}", file=sys.stderr)
            print(f"  only in gen : {sorted(set(found) - set(ref_in_range))[:5]}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
