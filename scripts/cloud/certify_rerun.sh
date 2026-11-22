#!/usr/bin/env bash
# Re-run Hex 116-118 with the integer-only certificate (--dual, no floating-point enumeration), so that the
# frontier past Hex 115 can be cited as proven. Shardable: pass K N to split across machines.
# Usage: bash certify_rerun.sh [K N] [threads]
set -euo pipefail
K=${1:-0}; N=${2:-1}; T=${3:-$(sysctl -n hw.ncpu 2>/dev/null || nproc)}
cd "$(dirname "$0")"
python3 run_shard.py 116 118 "$K" "$N" --mode dual --threads "$T" --out certify_116_118.jsonl
python3 merge_shards.py certify_116_118.jsonl
