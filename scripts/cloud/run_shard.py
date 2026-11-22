"""Run one shard of the exhaustive search and write a self-describing result record.

Every machine runs the same command with its own K (0..N-1); the shards partition the tree exactly
(test suite: 'shards partition the tree exactly'). Results are JSON lines, merged by merge_shards.py.

Usage: python3 run_shard.py HEX_FROM HEX_TO K N [--threads T] [--out results.jsonl] [--mode dual|stack]
  mode dual  : --dual only (integer-only certificate; use for coverage you want to cite as proven)
  mode stack : --dual --cvp (fastest; the enumeration's empty-search verdict relies on floating-point margins)
"""
import argparse, json, os, platform, re, socket, subprocess, time, hashlib
ap = argparse.ArgumentParser()
ap.add_argument("hex_from", type=int); ap.add_argument("hex_to", type=int); ap.add_argument("k", type=int); ap.add_argument("n", type=int)
ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
ap.add_argument("--out", default="shard_results.jsonl")
ap.add_argument("--mode", choices=("dual", "stack"), default="dual")
a = ap.parse_args()
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../agent_submission_v18/tripal")
binp = os.path.join(root, "tripal")
sha = hashlib.sha256(open(binp, "rb").read()).hexdigest()[:16]
src = hashlib.sha256(open(os.path.join(root, "src/search.c"), "rb").read()).hexdigest()[:16]
env = dict(os.environ, TP_DUAL="1")
if a.mode == "stack": env["TP_CVP"] = "1"
for L in range(a.hex_from, a.hex_to + 1):
    t0 = time.time()
    p = subprocess.run([binp, "hexdfs", "--hex-from", str(L), "--hex-to", str(L), "--threads", str(a.threads),
                        "--r5", "0", "--shard", f"{a.k}/{a.n}"], capture_output=True, text=True, env=env)
    m = re.search(r"total: (\d+) nodes, (\d+) prunes, (\d+) leaves, (\d+) hits", p.stdout)
    hits = re.findall(r"HIT\s+(\S+)", p.stdout)
    rec = {"hex": L, "shard": a.k, "of": a.n, "mode": a.mode, "ok": p.returncode == 0 and m is not None,
           "nodes": int(m.group(1)) if m else None, "hits": int(m.group(4)) if m else None, "hit_values": hits,
           "seconds": round(time.time() - t0, 1), "host": socket.gethostname(), "cpu": platform.machine(),
           "threads": a.threads, "binary_sha": sha, "source_sha": src, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    with open(a.out, "a") as f: f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
    if hits: print("!!! HIT FOUND: verify with  ./tripal verify " + " ".join(hits), flush=True)
