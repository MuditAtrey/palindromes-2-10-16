#!/usr/bin/env python3
"""Sequential per-hex-length sweep with the v18 engine (5-adic reach filter).

  python3 scripts/frontier_sweep_v18.py --from 1 --to 110 [--r5 12] [--threads 10]
                                         [--sieve] [--out runs/v18_sweep]

Each hex length h covers bit lengths 4h-3 .. 4h.  The raw engine output of every
length is kept (<out>/hex<h>.txt) and a summary line is appended to
<out>/summary.jsonl, so the sweep is resumable (finished lengths are skipped) and
auditable.  Any HIT line is flagged loudly; a nonzero return code stops the sweep.
"""
import argparse, json, os, re, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "agent_submission_v18", "tripal", "tripal")

ap = argparse.ArgumentParser()
ap.add_argument("--from", dest="h0", type=int, default=1)
ap.add_argument("--to", dest="h1", type=int, required=True)
ap.add_argument("--r5", type=int, default=12)
ap.add_argument("--threads", type=int, default=10)
ap.add_argument("--sieve", action="store_true")
ap.add_argument("--out", default=os.path.join(ROOT, "runs", "v18_sweep"))
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
summ = os.path.join(a.out, "summary.jsonl")
done = set()
if os.path.exists(summ):
    for line in open(summ):
        try: done.add(json.loads(line)["hex"])
        except Exception: pass

for h in range(a.h0, a.h1 + 1):
    if h in done:
        continue
    b0, b1 = max(1, 4 * h - 3), 4 * h
    cmd = [BIN, "hexdfs", "--bits-from", str(b0), "--bits-to", str(b1),
           "--threads", str(a.threads), "--r5", str(a.r5)]
    if a.sieve:
        cmd.append("--sieve")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    open(os.path.join(a.out, f"hex{h}.txt"), "w").write(p.stdout + p.stderr)
    m = re.search(r"total: (\d+) nodes, (\d+) prunes, (\d+) leaves, (\d+) hits", p.stdout)
    hits = [l for l in p.stdout.splitlines() if l.startswith("HIT ")]
    rec = dict(hex=h, bits=[b0, b1], rc=p.returncode, seconds=round(dt, 2),
               nodes=int(m.group(1)) if m else None, hits=int(m.group(4)) if m else None,
               hit_lines=hits, r5=a.r5, sieve=a.sieve, threads=a.threads,
               engine="agent_submission_v18+" + ("dual+" if os.environ.get("TP_DUAL") else "") + ("cvp" if os.environ.get("TP_CVP") else "fastk"))
    with open(summ, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"hex {h:3d} bits {b0}-{b1}: {rec['nodes']} nodes, {rec['hits']} hits, {dt:.1f}s"
          + ("   <<<<<< HIT >>>>>>" if (rec['hits'] or 0) > 0 and h > 3 else ""), flush=True)
    if p.returncode != 0 or m is None:
        print("engine failed; stopping", file=sys.stderr); sys.exit(1)
