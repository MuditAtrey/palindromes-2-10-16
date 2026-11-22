#!/usr/bin/env python3
"""tripal test suite.

  python3 tests/run_tests.py            full suite
  python3 tests/run_tests.py --quick    skip the slower checks
  python3 tests/run_tests.py -v         show subprocess output
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "tripal")
TMP = os.path.join(ROOT, "runs", "_test")
VERBOSE = False

results: list[tuple[bool, str, str]] = []


# ----------------------------------------------------------------------

def run(args, timeout=1800, check=False, env=None):
    if VERBOSE:
        print("    $", " ".join(args))
    p = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                       timeout=timeout, env=env)
    if check and p.returncode != 0:
        print(p.stdout)
        print(p.stderr, file=sys.stderr)
        raise RuntimeError(f"command failed: {' '.join(args)}")
    return p


def record(ok, name, detail=""):
    results.append((ok, name, detail))
    print(("  \033[32mPASS\033[0m  " if ok else "  \033[31mFAIL\033[0m  ") + name
          + (f"   [{detail}]" if detail else ""))


KV = re.compile(r"(\w+)=([-\w.]+)")


def parse_blines(text):
    """Parse the per-bit-length report lines into dicts."""
    out = []
    for line in text.splitlines():
        if line.startswith("B=") and "nodes=" in line:
            d = {k: v for k, v in KV.findall(line)}
            out.append(d)
    return out


def to_i(d, k):
    return int(d[k])


# ----------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------

def t_selftest(quick):
    """The C engine's own suite: brute force, A007632, re-verification,
    exhaustive binary-palindrome enumeration."""
    ref = 40 if quick else 46
    p = run([BIN, "selftest", "--brute-bits", "26", "--dual-bits", "40",
             "--ref-bits", str(ref)], timeout=1800)
    ok = "selftest PASSED" in p.stdout
    record(ok, "engine selftest (brute force + OEIS A007632 + re-verification)",
           f"ref_bits={ref}")
    if not ok:
        print(p.stdout)
        print(p.stderr)


def t_a007632_count():
    """Reproduce a long prefix of OEIS A007632 in the dual-base mode."""
    p = run([BIN, "selftest", "--brute-bits", "24", "--dual-bits", "48",
             "--ref-bits", "32"], timeout=900)
    m = re.search(r"OEIS terms in range : (\d+)", p.stdout)
    m2 = re.search(r"search      found   : (\d+)", p.stdout)
    ok = bool(m and m2) and m.group(1) == m2.group(1)
    record(ok, "A007632 reproduced exactly for all terms < 2^48",
           f"{m.group(1) if m else '?'} terms")


def t_ref_vs_python(quick):
    """C reference enumerator vs the independent Python enumerator."""
    bits = 40 if not quick else 32
    p = run([BIN, "ref", "--bits-from", "1", "--bits-to", str(bits)])
    c_vals = sorted(int(l.split()[1]) for l in p.stdout.splitlines() if l.strip())
    sys.path.insert(0, ROOT)
    from verify.check import triples_up_to_bits
    py_vals = sorted(v for v in triples_up_to_bits(bits) if v > 0)
    record(c_vals == py_vals, f"C reference enumeration == Python enumeration, B<={bits}",
           f"{len(py_vals)} values")


def t_determinism():
    """Node counts must not depend on thread count or work-unit split."""
    ref = None
    ok = True
    details = []
    for threads, units in ((1, 1), (1, 64), (2, 256), (4, 4096)):
        p = run([BIN, "search", "--bits-from", "184", "--bits-to", "184",
                 "--threads", str(threads), "--units", str(units),
                 "--no-checkpoint"])
        b = parse_blines(p.stderr + p.stdout)
        if not b:
            ok = False
            break
        sig = (to_i(b[0], "nodes"), to_i(b[0], "prunes"), to_i(b[0], "hits"))
        details.append(f"T{threads}/U{units}:{sig[0]}")
        if ref is None:
            ref = sig
        elif sig != ref:
            ok = False
    record(ok, "results identical across (threads, split depth)",
           " ".join(details))


def t_checkpoint_decomposition():
    """Splitting a run into checkpointed pieces must reproduce the same total."""
    shutil.rmtree(TMP, ignore_errors=True)
    base = [BIN, "search", "--bits-from", "200", "--bits-to", "200",
            "--threads", "1", "--units", "4096"]

    p = run(base + ["--workdir", os.path.join(TMP, "whole")])
    whole = parse_blines(p.stderr + p.stdout)[0]
    total_nodes = to_i(whole, "nodes")
    enum = to_i(whole, "enum")

    pieces, enum_sum = [], 0
    wd = os.path.join(TMP, "pieces")
    for _ in range(4):
        p = run(base + ["--workdir", wd, "--limit-units", "12"])
        b = parse_blines(p.stderr + p.stdout)[0]
        pieces.append(to_i(b, "nodes"))
        enum_sum += to_i(b, "enum")
        if to_i(b, "pend") == 0:
            break

    # each piece re-runs the split enumeration, so its cost is counted per piece
    expected = total_nodes + enum_sum - enum
    ok = sum(pieces) == expected
    record(ok, "checkpoint pieces sum to the whole run",
           f"{sum(pieces)} == {expected} over {len(pieces)} pieces")


def t_kill_and_resume():
    """SIGKILL mid-run, then resume; the run must complete and skip done units."""
    wd = os.path.join(TMP, "kill")
    shutil.rmtree(wd, ignore_errors=True)
    args = [BIN, "search", "--bits-from", "200", "--bits-to", "200",
            "--threads", "1", "--units", "4096", "--workdir", wd]

    proc = subprocess.Popen(args, cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    time.sleep(0.7)
    proc.send_signal(signal.SIGKILL)
    proc.wait()

    done_log = os.path.join(wd, "done.log")
    n_before = sum(1 for _ in open(done_log)) if os.path.exists(done_log) else 0

    p = run(args, timeout=600)
    b = parse_blines(p.stderr + p.stdout)
    ok = (n_before > 0 and b and to_i(b[0], "skip") > 0
          and "COMPLETE" in p.stdout)
    record(ok, "SIGKILL mid-run then resume to completion",
           f"{n_before} units survived the kill")

    # a third run must have nothing left to do
    p = run(args, timeout=600)
    b = parse_blines(p.stderr + p.stdout)
    ok2 = b and to_i(b[0], "units") == 0 and to_i(b[0], "pend") == 0
    record(ok2, "a fully checkpointed range does zero work on rerun")


def t_hits_log_and_verifier():
    """Every recorded hit must survive the independent Python verifier."""
    wd = os.path.join(TMP, "hits")
    shutil.rmtree(wd, ignore_errors=True)
    run([BIN, "search", "--bits-from", "1", "--bits-to", "40",
         "--threads", "2", "--units", "64", "--workdir", wd])
    hits = os.path.join(wd, "hits.log")
    if not os.path.exists(hits):
        record(False, "hits.log written", "missing")
        return
    p = run([sys.executable, os.path.join(ROOT, "verify", "check.py"),
             "--log", hits])
    ok = p.returncode == 0 and " 0 rejected" in p.stdout
    record(ok, "independent Python verifier accepts every hit",
           p.stdout.strip().splitlines()[-1])


def t_search_finds_known_solutions():
    """All triple palindromes below 2^32 must be exactly {0,1,3,5,7,9}."""
    p = run([BIN, "search", "--bits-from", "1", "--bits-to", "32",
             "--no-checkpoint"])
    found = sorted(int(m) for m in re.findall(r"HIT\s+B=\d+\s+n=(\d+)",
                                                p.stdout + p.stderr))
    expect = [1, 3, 5, 7, 9]
    record(found == expect, "search finds exactly {1,3,5,7,9} below 2^32",
           str(found))


def t_python_verifier_selftest():
    p = run([sys.executable, os.path.join(ROOT, "verify", "check.py"),
             "--selftest"])
    record(p.returncode == 0 and "PASSED" in p.stdout,
           "python verifier self-test")


# ----------------------------------------------------------------------

def t_hexdfs_brute_parity():
    """Outside-in hex-digit DFS must agree with brute force below 2^32."""
    b = run([BIN, "brute", "--bits-from", "1", "--bits-to", "32"], timeout=300)
    h = run([BIN, "hexdfs", "--hex-from", "1", "--hex-to", "8"], timeout=300)
    bs = {int(x) for x in b.stdout.split() if x.isdigit()}
    hs = {int(m.group(1)) for m in re.finditer(r"HIT (\d+)", h.stdout)}
    want = {1, 3, 5, 7, 9}
    ok = hs == want and want <= bs
    record(ok, "hexdfs agrees with brute force below 2^32",
           f"hexdfs={sorted(hs)}")


def t_hexdfs_node_accounting():
    """Lock the outside-in tree size; it must reproduce upstream's shape."""
    e = dict(os.environ); e["TP_R5"] = "0"      # the reference tree is the unfiltered one
    p = run([BIN, "hexdfs", "--hex-from", "31", "--hex-to", "38"], timeout=300, env=e)
    m = re.search(r"total: (\d+) nodes", p.stdout)
    n = int(m.group(1)) if m else -1
    ok = n == 1728434 and "HIT" not in p.stdout
    record(ok, "hexdfs hex 31-38 tree reproduces 1,728,434 nodes, 0 hits",
           f"{n} nodes")


def t_pruner_soundness_positive_control():
    """The single-endpoint pruner must agree with the two-conversion reference
    at every node -- and the check that says so must be able to fail.

    TP_XM recomputes the shared decimal prefix and the prune verdict the slow
    way (both endpoints converted, prefix compared digitwise) and counts any
    disagreement.  TP_BADPRUNE deliberately claims two extra shared digits, i.e.
    an unsound pruner; the same check must then report errors.  Without that
    second half this test would be vacuous -- which is exactly the bug that the
    old differential test had.
    """
    def xm(env, bits):
        e = dict(os.environ); e.update(env)
        p = subprocess.run([BIN, "hexdfs", "--bits-from", str(bits),
                            "--bits-to", str(bits), "--threads", "1"],
                           capture_output=True, text=True, env=e, timeout=300)
        m = re.search(r"\[xm\] (\d+) calls, m bad (\d+), verdict bad (\d+)"
                      r" \(we-prune/ref-keeps (\d+), ref-prunes/we-keep (\d+)\)", p.stdout + p.stderr)
        n = re.search(r"total: (\d+) nodes", p.stdout)
        return (tuple(int(g) for g in m.groups()) if m else None,
                int(n.group(1)) if n else -1)

    bad = []
    for bits in (121, 123, 124, 128):        # 1-pad, 1-pad, 0-pad, 0-pad
        (calls, mbad, vbad, fwd, back), nodes = xm({"TP_XM": "1"}, bits)
        if not calls or mbad or vbad or fwd:
            bad.append(f"B={bits} clean run: calls={calls} m_bad={mbad} "
                       f"verdict_bad={vbad} we-prune/ref-keeps={fwd}")
    record(not bad, "pruner agrees with the two-conversion reference at every node",
           "; ".join(bad))

    bad = []
    for bits in (124, 128):
        (calls, mbad, vbad, fwd, back), nodes = xm({"TP_XM": "1", "TP_BADPRUNE": "1"}, bits)
        if not calls or not fwd:
            bad.append(f"B={bits}: the deliberately unsound pruner was NOT caught "
                       f"(calls={calls} we-prune/ref-keeps={fwd}) -- check is vacuous")
    record(not bad, "positive control: a deliberately unsound pruner IS caught",
           "; ".join(bad))


def t_sharding_exact_and_hits_once():
    """Shards must partition the tree exactly and must not each report the hits
    that the shared serial prefix finds.  Both were broken: --shard was ignored
    entirely on one thread (the single-thread entry point never split), and every
    shard printed the same HIT lines."""
    def nodes_and_hits(args):
        p = run([BIN, "hexdfs", "--bits-from", "4", "--bits-to", "28"] + args, timeout=300)
        n = re.search(r"total: (\d+) nodes", p.stdout)
        return (int(n.group(1)) if n else -1,
                len(re.findall(r"^HIT ", p.stdout, re.M)),
                int(re.search(r"total: \d+ nodes, \d+ prunes, \d+ leaves, (\d+) hits",
                              p.stdout).group(1)) if re.search(r"total: .* (\d+) hits", p.stdout) else -1)

    unsplit, uhits, ucnt = nodes_and_hits(["--threads", "2"])
    tot = 0; hits = 0; cnt = 0
    for k in range(3):
        n, h, c = nodes_and_hits(["--threads", "1", "--cut", "2", "--shard", f"{k}/3"])
        tot += n; hits += h; cnt += c
    record(unsplit > 0 and tot == unsplit,
           "shards partition the tree exactly (3 shards sum to the unsplit run)",
           f"{tot} vs {unsplit}")
    record(hits == uhits and cnt == ucnt,
           "a hit is reported once across shards, not once per shard",
           f"HIT lines {hits} vs {uhits}, counter {cnt} vs {ucnt}")


def t_leaf_sieve_sound():
    """The leaf sieve skips the decimal test for numbers that cannot be
    palindromes (a decimal palindrome with an even digit count is divisible by
    11).  It must never change the hit list, and the tree it walks must be
    untouched: only the work done at the leaves may differ."""
    def run_it(env):
        e = dict(os.environ); e.update(env)
        p = subprocess.run([BIN, "hexdfs", "--bits-from", "1", "--bits-to", "80",
                            "--threads", "2"], capture_output=True, text=True,
                           env=e, timeout=300)
        hits = re.findall(r"^HIT (\d+)", p.stdout, re.M)
        m = re.search(r"total: (\d+) nodes", p.stdout)
        return hits, int(m.group(1)) if m else -1

    on, non = run_it({})
    off, noff = run_it({"TP_NOMOD11": "1"})
    record(on == off and len(on) == 5 and non == noff,
           "leaf sieve changes no hit and no node (5 known hits, same tree)",
           f"on={on} off={off} nodes {non} vs {noff}")


def _has_flint():
    return any(os.path.exists(d + "/flint/fmpz_lll.h") for d in ("/opt/homebrew/include", "/usr/include", "/usr/local/include"))


def _hexdfs(env_extra, b0, b1, threads=4, timeout=900):
    e = dict(os.environ)
    for k in ("TP_R5", "TP_SIEVE", "TP_WEAKPAL", "TP_NOMOD11", "TP_R5CHK", "TP_DECCHK", "TP_BLOOM", "TP_BLOOM_K", "TP_BLOOMCHK", "TP_DUAL", "TP_DUAL_K", "TP_DUAL_ALL", "TP_DUALCHK", "TP_CVP", "TP_CVP_K", "TP_CVP_ALL", "TP_CVP_BUDGET", "TP_CVPCHK", "TP_NOSHARE", "TP_CVP_SEED", "TP_CVP_PRUNE", "TP_SPEC", "TP_DUAL_POOL", "TP_OFFSET"):
        e.pop(k, None)
    e.update(env_extra)
    p = subprocess.run([BIN, "hexdfs", "--bits-from", str(b0), "--bits-to", str(b1),
                        "--threads", str(threads)], cwd=ROOT, capture_output=True,
                       text=True, env=e, timeout=timeout)
    m = re.search(r"total: (\d+) nodes, (\d+) prunes, (\d+) leaves, (\d+) hits", p.stdout)
    return (int(m.group(1)), int(m.group(4)), p) if m else (-1, -1, p)


def t_weak_palindrome_positive_control():
    """Every real search finds zero hits above 9, so 'no hit lost' proves nothing.
    This builds a weaker property that MANY numbers have -- the last J decimal
    digits mirror the first J -- and has the engine count them with each pruning
    rule off and then on.  The rules are sound for that property exactly when they
    are sound for full palindromes (they only ever use the first J digits), so the
    counts must be identical and must be nonzero.  A pruner that over-prunes loses
    hits here, where it would be invisible on the real problem."""
    base = {"TP_WEAKPAL": "5", "TP_NOMOD11": "1"}
    n0, h0, _ = _hexdfs(base, 60, 105)
    record(h0 > 1000, "weak-palindrome control is non-vacuous", f"{h0} hits, {n0} nodes")
    for name, extra in (("5-adic reach filter, J=5", {"TP_R5": "5"}),
                        ("5-adic reach filter, J=3", {"TP_R5": "3"}),
                        ("digit-sum / leading-2 sieve (--sieve)", {"TP_SIEVE": "1"}),
                        ("reach filter + sieve", {"TP_R5": "5", "TP_SIEVE": "1"})):
        e = dict(base); e.update(extra)
        n1, h1, _ = _hexdfs(e, 60, 105)
        record(h1 == h0 and 0 < n1 < n0,
               f"weak control: {name} loses no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")


def t_reach_filter_selfcheck():
    """Draw random real completions of every surviving node, build the actual
    number with bigint arithmetic, and require (a) its residue mod 5^J to be in
    the reach set and (b) its first J decimal digits to equal those the filter
    read off the node's lower endpoint."""
    n, h, p = _hexdfs({"TP_R5": "9", "TP_R5CHK": "1"}, 100, 200, threads=1)
    m = re.findall(r"\[r5chk\] (\d+) sampled completions, (\d+) NOT in reach set; "
                   r"(\d+) prefix samples, (\d+) prefix MISMATCH", p.stderr)
    last = m[-1] if m else None
    record(bool(last) and int(last[0]) > 100000 and last[1] == "0" and last[3] == "0",
           "reach filter: sampled real completions are all inside the reach set",
           f"{last[0]} samples, {last[1]} outside; {last[2]} prefix samples, {last[3]} mismatch" if last else "no output")


def t_deep_filter_loses_no_hit():
    """The Bloom deep reach filter must lose no weak-palindrome hit (K = J so the
    test genuinely fires), must actually prune, and its carried residue, target and
    sampled members must all check out."""
    e0 = {"TP_R5": "0", "TP_WEAKPAL": "6"}
    n0, h0, _ = _hexdfs(e0, 80, 124)
    e1 = dict(e0, TP_BLOOM="10:12", TP_BLOOM_K="6")
    n1, h1, _ = _hexdfs(e1, 80, 124)
    record(h0 > 1000 and h0 == h1 and n1 < n0 / 3,
           "weak control: deep reach filter (--bloom) loses no hit",
           f"hits {h1} vs {h0}, nodes {n1} vs {n0}")
    n, h, p = _hexdfs({"TP_R5": "10", "TP_BLOOM": "22:12", "TP_BLOOMCHK": "1"}, 230, 260, threads=2)
    bad = "BLOOMCHK" in p.stderr or p.returncode != 0
    m = re.findall(r"\[bloomchk\] d=\d+: (\d+) of (\d+) real sums missing", p.stderr)
    record(bool(m) and not bad and all(x[0] == "0" for x in m),
           "deep reach filter: carried residue, target and sampled members all correct",
           f"{len(m)} filters checked, nodes {n}")


def t_dual_certificate_loses_no_hit():
    """The dual-lattice certificate must lose no weak-palindrome hit (K = J so it
    genuinely fires), must prune, and must never reject a real subset sum."""
    if not _has_flint():
        record(True, "dual certificate: skipped (no FLINT)", "")
        return
    for J, b0, b1 in ((6, 80, 124), (7, 80, 124), (8, 90, 124)):
        e0 = {"TP_R5": "0", "TP_WEAKPAL": str(J)}
        n0, h0, _ = _hexdfs(e0, b0, b1)
        e1 = dict(e0, TP_DUAL="1", TP_DUAL_K=str(J), TP_DUAL_ALL="1")
        n1, h1, _ = _hexdfs(e1, b0, b1)
        record(h0 > 10 and h0 == h1 and n1 < n0 / 2,
               f"weak control J={J}: dual certificate loses no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")
    for J, b0, b1 in ((6, 80, 124), (8, 90, 124)):
        e0 = {"TP_R5": "0", "TP_WEAKPAL": str(J)}
        n0, h0, _ = _hexdfs(e0, b0, b1)
        e1 = dict(e0, TP_DUAL="1", TP_DUAL_K=str(J), TP_DUAL_ALL="1", TP_DUAL_POOL="3000")
        n1, h1, _ = _hexdfs(e1, b0, b1)
        record(h0 > 10 and h0 == h1 and n1 < n0 / 2,
               f"weak control J={J}: learned windows lose no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")
    n, h, p = _hexdfs({"TP_R5": "10", "TP_DUAL": "1", "TP_DUALCHK": "1", "TP_DUAL_POOL": "3000"}, 228, 260, threads=2)
    m = re.findall(r"\[dualchk\] d=\d+: (\d+) of (\d+) real sums rejected", p.stderr)
    record(bool(m) and p.returncode == 0 and all(x[0] == "0" for x in m),
           "learned windows: real sums never rejected",
           f"{len(m)} levels checked, nodes {n}")
    sp = {"TP_SPEC": "1", "TP_SPEC_ITERS": "30", "TP_SPEC_M": "300", "TP_SPEC_TRIP": "5000", "TP_DUAL_EXTRA": "2", "TP_DUAL_LEVELS": "10"}
    for J, b0, b1 in ((6, 96, 112), (8, 100, 116)):
        e0 = {"TP_R5": "0", "TP_WEAKPAL": str(J)}
        n0, h0, _ = _hexdfs(e0, b0, b1)
        e1 = dict(e0, TP_DUAL="1", TP_DUAL_K=str(J), TP_DUAL_ALL="1", **sp)
        n1, h1, _ = _hexdfs(e1, b0, b1)
        record(h0 > 0 and h0 == h1 and n1 < n0,
               f"weak control J={J}: Fourier certificate loses no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")
    n, h, p = _hexdfs(dict({"TP_R5": "10", "TP_DUAL": "1", "TP_DUALCHK": "1"}, **sp), 228, 244, threads=2)
    m = re.findall(r"\[dualchk\] d=\d+: (\d+) of (\d+) real sums rejected", p.stderr)
    record(bool(m) and p.returncode == 0 and all(x[0] == "0" for x in m),
           "Fourier certificate: real sums never rejected",
           f"{len(m)} levels checked, nodes {n}")
    n, h, p = _hexdfs({"TP_R5": "10", "TP_DUAL": "1", "TP_DUALCHK": "1"}, 228, 260, threads=2)
    m = re.findall(r"\[dualchk\] d=\d+: (\d+) of (\d+) real sums rejected", p.stderr)
    record(bool(m) and "DUALCHK" not in p.stderr and p.returncode == 0 and all(x[0] == "0" for x in m),
           "dual certificate: real sums never rejected, targets match the decimal string",
           f"{len(m)} levels checked, nodes {n}")


def t_lattice_enumeration_loses_no_hit():
    """Exact lattice enumeration (--cvp) must lose no weak-palindrome hit, must prune,
    and must find every genuine digit combination."""
    if not _has_flint():
        record(True, "lattice enumeration: skipped (no FLINT)", "")
        return
    for J, b0, b1 in ((6, 80, 124), (7, 80, 124), (8, 90, 124)):
        e0 = {"TP_R5": "0", "TP_WEAKPAL": str(J)}
        n0, h0, _ = _hexdfs(e0, b0, b1)
        e1 = dict(e0, TP_CVP="1", TP_CVP_K=str(J), TP_CVP_ALL="1", TP_CVP_BUDGET="2000")
        n1, h1, _ = _hexdfs(e1, b0, b1)
        record(h0 > 10 and h0 == h1 and n1 < n0,
               f"weak control J={J}: lattice enumeration loses no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")
    n, h, p = _hexdfs({"TP_R5": "0", "TP_CVP": "1", "TP_CVPCHK": "1", "TP_CVP_BUDGET": "2000"}, 300, 336, threads=2)
    m = re.findall(r"\[cvpchk\] d=\d+: (\d+) of (\d+) genuine sums rejected", p.stderr)
    record(bool(m) and p.returncode == 0 and all(x[0] == "0" for x in m),
           "lattice enumeration: every genuine sum is found",
           f"{len(m)} levels checked, nodes {n}")


def t_flags_reach_threaded_path():
    """The threaded and sharded entry points once ignored every TP_* switch except
    TP_NOPRUNE, so 'sieve off' runs on several threads silently ran the sieve."""
    a, ha, _ = _hexdfs({"TP_WEAKPAL": "5", "TP_NOMOD11": "1"}, 60, 90, threads=1)
    b, hb, _ = _hexdfs({"TP_WEAKPAL": "5", "TP_NOMOD11": "1"}, 60, 90, threads=4)
    record(a == b and ha == hb and ha > 0,
           "TP_* switches behave identically on 1 and 4 threads",
           f"nodes {a}/{b}, hits {ha}/{hb}")


def t_mod11_residue_carried_exactly():
    """The carried N mod 11 must equal the true residue at every leaf (the first
    leaf sieve double-added a term and wrongly rejected 7.5% of multiples of 11)."""
    n, h, p = _hexdfs({"TP_DECCHK": "1"}, 1, 200, threads=1)
    record(p.returncode == 0 and "MOD11 carried residue wrong" not in p.stderr and h == 5,
           "mod-11 leaf sieve: carried residue equals the true residue", f"hits {h}")



def t_sibling_sharing_is_exact():
    """Children reuse the parent's 256-bit products (a(base - add) = a base - a add mod N) and the
    enumeration runs only on survivors of the cheap tests: node counts must be identical to the
    per-child reference path (TP_NOSHARE)."""
    if not _has_flint():
        record(True, "sibling sharing: skipped (no FLINT)", "")
        return
    base = {"TP_R5": "0", "TP_DUAL": "1", "TP_CVP": "1"}
    n1, h1, _ = _hexdfs(base, 300, 330)
    n2, h2, _ = _hexdfs(dict(base, TP_NOSHARE="1"), 300, 330)
    record(n1 == n2 and h1 == h2 and n1 > 0, "sibling sharing gives identical nodes to the per-child path",
           f"nodes {n1} vs {n2}")
    e0 = {"TP_R5": "0", "TP_WEAKPAL": "7"}
    n0, h0, _ = _hexdfs(e0, 80, 124)
    e1 = dict(e0, TP_DUAL="1", TP_DUAL_K="7", TP_DUAL_ALL="1", TP_CVP="1", TP_CVP_K="7", TP_CVP_ALL="1", TP_CVP_BUDGET="2000")
    n1, h1, _ = _hexdfs(e1, 80, 124)
    record(h0 > 10 and h0 == h1, "weak control J=7: shared dual + lattice path loses no hit", f"hits {h1} vs {h0}, nodes {n1} vs {n0}")


def t_seeded_basis_is_exact():
    """A rerandomized lattice basis (TP_CVP_SEED, discovery mode) describes the same lattice:
    without pruning, the enumeration must still lose no weak-palindrome hit."""
    if not _has_flint():
        record(True, "seeded basis: skipped (no FLINT)", "")
        return
    e0 = {"TP_R5": "0", "TP_WEAKPAL": "7"}
    n0, h0, _ = _hexdfs(e0, 80, 124)
    for seed in ("1", "2"):
        e1 = dict(e0, TP_CVP="1", TP_CVP_K="7", TP_CVP_ALL="1", TP_CVP_BUDGET="2000", TP_CVP_SEED=seed)
        n1, h1, _ = _hexdfs(e1, 80, 124)
        record(h0 > 10 and h0 == h1, f"weak control J=7: lattice enumeration on seeded basis {seed} loses no hit",
               f"hits {h1} vs {h0}, nodes {n1} vs {n0}")


def t_planted_offset_solutions_found():
    """Offset mode (TP_OFFSET=c: N + c must be a decimal palindrome) with c = 0 equals the real search, and planted
    solutions (N0 a base-2/16 palindrome, P0 a decimal palindrome, c = P0 - N0) are found by the full stack on several
    thread counts. This is the engine's only full-scale positive control with real solutions."""
    if not _has_flint():
        record(True, "planted offsets: skipped (no FLINT)", "")
        return
    import random
    n0, h0, _ = _hexdfs({"TP_R5": "8", "TP_DUAL": "1", "TP_CVP": "1", "TP_NOMOD11": "1"}, 200, 300)
    n1, h1, _ = _hexdfs({"TP_R5": "8", "TP_DUAL": "1", "TP_CVP": "1", "TP_OFFSET": "0"}, 200, 300)
    record(n0 == n1 and h0 == h1, "offset mode with c = 0 equals the real search", f"nodes {n1} vs {n0}")
    rnd = random.Random(31); found = 0; tot = 0
    for L in (20, 30, 38, 44, 52):
        half = (L + 1) // 2
        e = [rnd.choice((9, 15))] + [rnd.choice((0, 6, 9, 15)) for _ in range(half - 1)]
        N0 = 0
        for x in e + e[: L // 2][::-1]: N0 = N0 * 16 + x
        D = len(str(N0))
        while True:
            h = [rnd.randrange(1, 10)] + [rnd.randrange(10) for _ in range((D + 1) // 2 - 1)]
            P0 = int("".join(map(str, h + h[: D // 2][::-1])))
            if P0 > N0: break
        for th in (1, 4):
            tot += 1
            _, _, p = _hexdfs({"TP_R5": "8", "TP_DUAL": "1", "TP_CVP": "1", "TP_OFFSET": str(P0 - N0)}, 4 * L, 4 * L, threads=th)
            found += str(P0) in re.findall(r"HIT\s+(\d+)", p.stdout)
    record(found == tot, "planted solutions found by the full stack (1 and 4 threads)", f"{found} of {tot}")

def main():
    global VERBOSE
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--no-build", action="store_true")
    a = ap.parse_args()
    VERBOSE = a.verbose

    os.chdir(ROOT)
    os.makedirs(os.path.join(ROOT, "runs"), exist_ok=True)

    if not a.no_build:
        print("building ...")
        p = run(["make", "-j"], timeout=600)
        if p.returncode != 0:
            print(p.stdout, p.stderr)
            return 1

    print("\n=== tripal test suite ===\n")
    t0 = time.time()
    for fn in (t_python_verifier_selftest,
               t_search_finds_known_solutions,
               t_selftest,
               t_a007632_count,
               t_ref_vs_python,
               t_determinism,
               t_hits_log_and_verifier,
               t_checkpoint_decomposition,
               t_kill_and_resume,
               t_hexdfs_brute_parity,
               t_hexdfs_node_accounting,
               t_pruner_soundness_positive_control,
               t_sharding_exact_and_hits_once,
               t_leaf_sieve_sound,
               t_weak_palindrome_positive_control,
               t_reach_filter_selfcheck,
               t_deep_filter_loses_no_hit,
               t_dual_certificate_loses_no_hit,
               t_lattice_enumeration_loses_no_hit,
               t_sibling_sharing_is_exact,
               t_seeded_basis_is_exact,
               t_planted_offset_solutions_found,
               t_flags_reach_threaded_path,
               t_mod11_residue_carried_exactly):
        try:
            if fn in (t_selftest, t_ref_vs_python):
                fn(a.quick)
            else:
                fn()
        except Exception as exc:            # noqa: BLE001
            record(False, fn.__name__, f"exception: {exc}")

    npass = sum(1 for ok, _, _ in results if ok)
    nfail = len(results) - npass
    print(f"\n{npass} passed, {nfail} failed in {time.time() - t0:.1f}s")
    print("RESULT:", "\033[32mALL TESTS PASSED\033[0m" if not nfail
          else "\033[31mFAILURES\033[0m")
    return 0 if not nfail else 1


if __name__ == "__main__":
    sys.exit(main())
