#!/usr/bin/env python3
"""Check the released files against the timestamped commitment.

commitment/manifest.txt was timestamped with OpenTimestamps on 2026-11-21, before any code was public. It lists the SHA-256 of
every engine source file, the test suite, the scripts and the result logs as they were in the private research repository.
This script recomputes those hashes for the files in this repository and reports, for each one:

  OK         the released file is byte-identical to the committed one
  MISMATCH   the file exists but differs (this would mean the release does not match the commitment)
  REDACTED   the committed file is deliberately not published byte-for-byte (see COMMITMENT.md)

The private repository kept the engine in a folder named agent_submission_v18/tripal; here it lives in engine/. Nothing else
is renamed. Run it from the repository root, on the tag v1.0-as-committed (later commits add comments to the source, which
changes the bytes but not the code):

    git checkout v1.0-as-committed && python3 verify/verify_commitment.py
"""
import hashlib, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ROOT, "commitment", "manifest.txt")
RENAMES = {"agent_submission_v18/tripal/": "engine/"}
# Published only as copies with the machine's host name removed (COMMITMENT.md, "Note on the public result files").
REDACTED = {"results/certify_hex116_118_integer_only.jsonl"}
MANIFEST_SHA256 = "efff8c5d592c48a79739cd427ee66ff5e54e1fd4df48fe0f90e8057a90b6fed2"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    bad = 0
    got = sha256(MANIFEST)
    print(f"manifest.txt sha256 {got} {'OK' if got == MANIFEST_SHA256 else 'MISMATCH'}")
    bad += got != MANIFEST_SHA256
    for line in open(MANIFEST, encoding="utf-8"):
        parts = line.split()
        if len(parts) != 2 or len(parts[0]) != 64:
            continue                      # header lines, not file entries
        want, committed_path = parts
        path = committed_path
        for old, new in RENAMES.items():
            if path.startswith(old):
                path = new + path[len(old):]
        if committed_path in REDACTED:
            print(f"REDACTED  {path}")
            continue
        full = os.path.join(ROOT, path)
        if not os.path.exists(full):
            print(f"MISSING   {path}"); bad += 1; continue
        ok = sha256(full) == want
        print(f"{'OK      ' if ok else 'MISMATCH'}  {path}")
        bad += not ok
    print("all published files match the commitment" if not bad else f"{bad} problem(s)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
