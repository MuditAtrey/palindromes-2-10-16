#!/usr/bin/env bash
# One-time setup of a fresh Ubuntu/Debian machine (cloud VM or spare box) for sharded search.
# Usage: bash setup_linux.sh   (then: python3 run_shard.py ...)
set -euo pipefail
sudo apt-get update -y
sudo apt-get install -y build-essential python3 libflint-dev libgmp-dev libmpfr-dev git
cd "$(dirname "$0")/../../engine"
make clean >/dev/null 2>&1 || true
make
# self-checks before any shard is trusted: brute force, OEIS A007632, re-verification
./tripal selftest
python3 tests/run_tests.py
echo "machine ready: $(nproc) cores"
