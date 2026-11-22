"""Planted solutions at full scale (engine offset mode, TP_OFFSET=c).

The engine can search for base-2/16 palindromes N such that N + c is a decimal palindrome. With c = 0 that is the real
problem; for a planted instance we pick a candidate N0 (hex digits from {0,6,9,F}, mirrored, B = 4L) and a decimal
palindrome P0 > N0, and set c = P0 - N0. The search must then report P0. This gives:
  * a positive control for the whole engine stack at real lengths (exact mode must find every plant), and
  * a measurement of discovery mode: which plants a pruned/seeded run still finds, and how fast it is.

Usage:
  python3 plant.py control  L0 L1 PER_L [ENGINE_ENV...]      exact stack must find every plant
  python3 plant.py discover L PLANTS q SEEDS [ENGINE_ENV...]   recall and time of pruned runs vs exact, same plants
ENGINE_ENV items look like TP_DUAL=1 (defaults: TP_DUAL=1 TP_CVP=1, --r5 0).
"""
import os, random, re, subprocess, sys, time
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../agent_submission_v18/tripal/tripal")
CLEAN = [k for k in os.environ if k.startswith("TP_")]

def make_plant(L, rnd):
    half = (L + 1) // 2
    e = [rnd.choice((9, 15))] + [rnd.choice((0, 6, 9, 15)) for _ in range(half - 1)]
    hx = e + e[: L // 2][::-1]
    N0 = 0
    for x in hx: N0 = N0 * 16 + x
    D = len(str(N0))
    while True:
        dd = D if rnd.random() < 0.8 else D + 1
        h = [rnd.randrange(1, 10)] + [rnd.randrange(10) for _ in range((dd + 1) // 2 - 1)]
        P0 = int("".join(map(str, h + h[: dd // 2][::-1])))
        if P0 > N0: return N0, P0, P0 - N0

def run(L, c, extra, threads=10):
    e = {k: v for k, v in os.environ.items() if k not in CLEAN}
    e.update({"TP_DUAL": "1", "TP_CVP": "1", "TP_OFFSET": str(c)})
    for kv in extra:
        k, v = kv.split("=", 1)
        if v == "-": e.pop(k, None)          # e.g. TP_CVP=- : run without the lattice enumeration
        else: e[k] = v
    t0 = time.time()
    p = subprocess.run([BIN, "hexdfs", "--hex-from", str(L), "--hex-to", str(L), "--threads", str(threads), "--r5", "0"],
                       capture_output=True, text=True, env=e)
    if p.returncode != 0: raise SystemExit(p.stderr[-2000:])
    hits = {int(h) for h in re.findall(r"HIT\s+(\d+)", p.stdout)}
    m = re.search(r"total: (\d+) nodes", p.stdout)
    return hits, int(m.group(1)) if m else -1, time.time() - t0

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "control":
        L0, L1, per = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); extra = sys.argv[5:]
        rnd = random.Random(2026); ok = bad = 0
        for L in range(L0, L1 + 1):
            for _ in range(per):
                N0, P0, c = make_plant(L, rnd)
                hits, n, dt = run(L, c, extra)
                found = P0 in hits
                ok += found; bad += not found
                print(f"  hex {L:3d}: plant N0=0x{N0:X} P0={P0}  {'FOUND' if found else 'MISSED'}  ({n:,d} nodes, {dt:.1f} s, {len(hits)} hits)", flush=True)
        print(f"control: {ok} found, {bad} missed  (exact mode must miss none)")
    elif mode == "discover":
        L, n_pl, q, seeds = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]); extra = sys.argv[6:]
        rnd = random.Random(7 + L); te = tp = 0.0; fe = fp1 = fpu = 0
        for i in range(n_pl):
            N0, P0, c = make_plant(L, rnd)
            he, _, dte = run(L, c, extra); te += dte; fe += P0 in he
            found_any = False; f1 = False
            for s in range(1, seeds + 1):
                hp, _, dtp = run(L, c, extra + [f"TP_CVP_PRUNE={q}", f"TP_CVP_SEED={s}"]); tp += dtp
                f = P0 in hp
                if s == 1: fp1 += f; f1 = f
                found_any |= f
            fpu += found_any
            print(f"  plant {i+1}: exact {'found' if P0 in he else 'MISSED'} {dte:.1f} s | pruned q={q}: seed1 {'found' if f1 else 'missed'}, any of {seeds} {'found' if found_any else 'missed'}", flush=True)
        print(f"hex {L}, {n_pl} plants: exact found {fe}/{n_pl} in {te:.0f} s total;"
              f" pruned q={q}: 1 seed {fp1}/{n_pl}, {seeds} seeds {fpu}/{n_pl}, time {tp/seeds:.0f} s per seed-pass "
              f"(speed {te/(tp/seeds):.2f}x per pass, hits per compute {(fpu/max(fe,1))/(tp/te):.2f}x of exact)")
