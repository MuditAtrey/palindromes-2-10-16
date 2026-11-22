#!/usr/bin/env bash
# Checkpointed sweep of a hex-digit range. Safe to interrupt and re-run: every
# finished work unit is recorded in <workdir>/done.log and is skipped next time.
#
#   scripts/sweep.sh --from-hex 33 --to-hex 38 --threads 10
#
set -euo pipefail

FROM=33; TO=38; THREADS=10; RATE=5.5e6; WORKDIR=""; YES=0; BUDGET=0
while [ $# -gt 0 ]; do
  case "$1" in
    --from-hex) FROM="$2"; shift 2;;
    --to-hex)   TO="$2";   shift 2;;
    --threads)  THREADS="$2"; shift 2;;
    --rate)     RATE="$2"; shift 2;;
    --node-budget) BUDGET="$2"; shift 2;;
    --workdir)  WORKDIR="$2"; shift 2;;
    --yes|-y)   YES=1; shift;;
    *) echo "unknown option $1" >&2; exit 2;;
  esac
done
[ -n "$WORKDIR" ] || WORKDIR="runs/hex${FROM}-${TO}"

cd "$(dirname "$0")/.."

if [ ! -x ./tripal ]; then
  echo "building..."
  (make native 2>/dev/null || make) >/dev/null
fi

echo "=================================================================="
echo " triple-palindrome sweep: hex digits $FROM..$TO  (B = $((4*(FROM-1)+1))..$((4*TO)))"
echo " threads $THREADS   workdir $WORKDIR"
echo "=================================================================="
./tripal estimate --hex-from "$FROM" --hex-to "$TO" --threads "$THREADS" --rate "$RATE"

if [ "$YES" -ne 1 ]; then
  printf '\nproceed? [y/N] '
  read -r ans
  case "$ans" in y|Y) ;; *) echo "aborted."; exit 0;; esac
fi

mkdir -p "$WORKDIR"
LOG="$WORKDIR/run.log"
mkdir -p "$WORKDIR"

ARGS=(search --hex-from "$FROM" --hex-to "$TO" --threads "$THREADS" --workdir "$WORKDIR")
[ "$BUDGET" != "0" ] && ARGS+=(--node-budget "$BUDGET")

# keep the machine awake on macOS
RUNNER=()
if command -v caffeinate >/dev/null 2>&1; then RUNNER=(caffeinate -is); fi

echo "--- starting: ${RUNNER[*]:-} ./tripal ${ARGS[*]}"
( "${RUNNER[@]}" ./tripal "${ARGS[@]}" 2>&1 ) | tee -a "$LOG"
rc=${PIPESTATUS[0]}

echo
echo "--- run finished (exit $rc); log: $LOG"

# ---- verify every hit independently -------------------------------------
HITS="$WORKDIR/hits.log"
ok=0; bad=0
if [ -f "$HITS" ] && [ -s "$HITS" ]; then
  echo "--- independently verifying hits from $HITS"
  if out=$(python3 verify/check.py --log "$HITS" 2>&1); then
    echo "$out" | sed 's/^/    /'
    ok=$(echo "$out" | tail -1 | sed -n 's/^\([0-9]*\) verified.*/\1/p')
    ok=${ok:-0}
  else
    echo "$out" | sed 's/^/    /'
    bad=1
  fi
else
  echo "--- no hits found in this range"
fi

SUMMARY="hex $FROM..$TO done: $ok verified, $bad failed"
echo "--- $SUMMARY"

if command -v osascript >/dev/null 2>&1; then
  osascript -e "display notification \"$SUMMARY\" with title \"tripal sweep\"" || true
fi

[ "$bad" -eq 0 ] || exit 1
