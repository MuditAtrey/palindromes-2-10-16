#!/usr/bin/env python3
"""Independent verifier for tripal results.

Nothing here shares code with the C engine: palindrome tests use Python's own
arbitrary-precision integers and string formatting, so a bug in the C search
cannot hide behind a matching bug in its verifier.

Usage
-----
  python3 verify/check.py 12345 0xdeadbeef ...      verify given numbers
  python3 verify/check.py --log runs/x/hits.log     verify every recorded hit
  python3 verify/check.py --enumerate-bits 44       exhaustive cross-check
  python3 verify/check.py --selftest                sanity-check this script
"""

from __future__ import annotations

import argparse
import sys

BASES = (2, 10, 16)


def digits(n: int, base: int) -> str:
    """Representation of n in `base`, lowercase, no prefix."""
    if n < 0:
        raise ValueError("negative input")
    if n == 0:
        return "0"
    out = []
    while n:
        n, r = divmod(n, base)
        out.append("0123456789abcdef"[r])
    return "".join(reversed(out))


def is_palindrome(n: int, base: int) -> bool:
    s = digits(n, base)
    return s == s[::-1]


def check(n: int) -> dict:
    return {
        "n": n,
        "dec": str(n),
        "hex": "0x%x" % n,
        "reps": {b: digits(n, b) for b in BASES},
        "pal": {b: is_palindrome(n, b) for b in BASES},
    }


def is_triple(n: int) -> bool:
    return all(is_palindrome(n, b) for b in BASES)


def enumerate_binary_palindromes(max_bits: int):
    """Every number whose base-2 representation is a palindrome, with at most
    `max_bits` bits.  Built from the free half, not from the orbit machinery."""
    for bits in range(1, max_bits + 1):
        half = (bits + 1) // 2
        # bit 0 must be 1 (it is the mirror of the top bit), so the free half
        # is odd; that also guarantees the value really has `bits` bits.
        for hi in range(1, 1 << half, 2):
            v = hi
            for i in range(half, bits):
                v |= ((hi >> (bits - 1 - i)) & 1) << i
            yield v


def triples_up_to_bits(max_bits: int):
    return [0] + [v for v in enumerate_binary_palindromes(max_bits) if is_triple(v)]


# ----------------------------------------------------------------------

def cmd_numbers(nums) -> int:
    bad = 0
    for n in nums:
        r = check(n)
        print(f"n     = {r['dec']}")
        print(f"hex   = {r['hex']}")
        print(f"base2 {r['reps'][2]}")
        print(f"base10 {r['reps'][10]}")
        print(f"base16 {r['reps'][16]}")
        ok = all(r["pal"].values())
        print("=> " + ("TRIPLE PALINDROME" if ok else "not a triple palindrome"))
        print()
        if not ok:
            bad += 1
    return 1 if bad else 0


def cmd_log(path: str) -> int:
    """hits.log lines look like:  <bits> <decimal> 0x<hex>"""
    n_ok = n_bad = 0
    seen = set()
    try:
        lines = open(path, "r", encoding="utf-8").read().splitlines()
    except OSError as exc:
        print(f"cannot read {path}: {exc}", file=sys.stderr)
        return 2

    for line in lines:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            print(f"MALFORMED: {line!r}")
            n_bad += 1
            continue
        bits = int(parts[0])
        dec = parts[1]
        n = int(dec, 10)
        ok = is_triple(n) and n.bit_length() == bits
        if len(parts) > 2:
            ok = ok and int(parts[2], 16) == n
        if n in seen:
            print(f"DUPLICATE: {n}")
            ok = False
        seen.add(n)
        print(f"{'PASS' if ok else 'FAIL'}  B={bits}  n={n}  0x{n:x}")
        n_ok += ok
        n_bad += not ok

    print(f"\n{n_ok} verified, {n_bad} rejected, {len(seen)} distinct")
    return 1 if n_bad else 0


def cmd_enumerate(max_bits: int) -> int:
    res = triples_up_to_bits(max_bits)
    for v in res:
        print(f"{v.bit_length() or 0} {v} 0x{v:x}")
    print(f"\nexhaustive enumeration of all binary palindromes with B <= {max_bits}: "
          f"{len(res)} triple palindromes (including 0)", file=sys.stderr)
    return 0


def cmd_selftest() -> int:
    fails = 0

    def expect(cond, msg):
        nonlocal fails
        print(("PASS  " if cond else "FAIL  ") + msg)
        if not cond:
            fails += 1

    expect(is_palindrome(585, 2) and is_palindrome(585, 10),
           "585 is palindromic in bases 2 and 10")
    expect(not is_palindrome(585, 16), "585 (0x249) is not palindromic in base 16")
    expect(all(is_palindrome(9, b) for b in BASES), "9 is a triple palindrome")
    expect(not is_palindrome(33, 16), "33 (0x21) is not palindromic in base 16")
    expect(digits(255, 16) == "ff", "digits(255,16) == ff")
    expect(digits(0, 2) == "0", "digits(0,2) == 0")

    # every A007632 term must be a base-2 and base-10 palindrome
    terms = []
    for path in ("data/oeis/b007632.txt", "../data/oeis/b007632.txt"):
        try:
            for line in open(path, encoding="utf-8"):
                if line.startswith("#") or not line.strip():
                    continue
                terms.append(int(line.split()[1]))
            break
        except OSError:
            continue
    expect(len(terms) == 213, f"A007632 b-file loaded ({len(terms)} terms)")
    if terms:
        expect(all(is_palindrome(v, 2) and is_palindrome(v, 10) for v in terms),
               "every A007632 term is palindromic in bases 2 and 10")
        expect(all(terms[i] < terms[i + 1] for i in range(len(terms) - 1)),
               "A007632 terms strictly increasing")

    expect(triples_up_to_bits(24) == [0, 1, 3, 5, 7, 9],
           "exhaustive search to B=24 gives [0,1,3,5,7,9]")

    print(f"\n{'FAILED' if fails else 'PASSED'}")
    return 1 if fails else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="independent tripal verifier")
    ap.add_argument("numbers", nargs="*", help="decimal or 0x-prefixed numbers")
    ap.add_argument("--log", metavar="FILE", help="verify a hits.log")
    ap.add_argument("--enumerate-bits", type=int, metavar="N",
                    help="exhaustive enumeration of binary palindromes up to N bits")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        return cmd_selftest()
    if a.enumerate_bits:
        return cmd_enumerate(a.enumerate_bits)
    if a.log:
        return cmd_log(a.log)
    if a.numbers:
        nums = [int(s, 0) if s.lower().startswith("0x") else int(s, 10)
                for s in a.numbers]
        return cmd_numbers(nums)
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
