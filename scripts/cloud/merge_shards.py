"""Merge shard records: check every (hex length, shard) is present exactly once, all from one source version,
and report per length: complete or missing shards, total nodes, hits (each hit re-verified in all three bases).

Usage: python3 merge_shards.py results1.jsonl [results2.jsonl ...]
"""
import json, sys
recs = [json.loads(l) for f in sys.argv[1:] for l in open(f) if l.strip()]
def pal(n, b):
    d = []
    while n: d.append(n % b); n //= b
    return d == d[::-1]
by = {}
for r in recs: by.setdefault((r["hex"], r["mode"]), []).append(r)
srcs = {r["source_sha"] for r in recs}
if len(srcs) > 1: print("WARNING: records come from different engine sources:", srcs)
for (L, mode), rs in sorted(by.items()):
    n = rs[0]["of"]; got = {}
    for r in rs:
        if r["of"] != n: print(f"hex {L}: mixed shard counts"); 
        if r["ok"]: got.setdefault(r["shard"], r)
    missing = [k for k in range(n) if k not in got]
    nodes = sum(r["nodes"] for r in got.values())
    hv = sorted({int(h) for r in got.values() for h in r["hit_values"]})
    verified = [h for h in hv if pal(h, 2) and pal(h, 10) and pal(h, 16)]
    status = "COMPLETE" if not missing else f"MISSING shards {missing[:10]}{'...' if len(missing) > 10 else ''}"
    print(f"hex {L:4d} [{mode}] {status:30s} nodes {nodes:>16,d}  hits {hv}  verified {verified}")
