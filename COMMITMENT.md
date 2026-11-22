# Commitment and timestamp

Before any code was public, the exact contents of the search engine, its test suite and all result logs were fixed by
[commitment/manifest.txt](commitment/manifest.txt), whose existence on 2026-11-21 is proven by the OpenTimestamps proof
[commitment/manifest.txt.ots](commitment/manifest.txt.ots) (anchored in the Bitcoin blockchain). The code is now released,
and the released files can be checked against that commitment.

## What the manifest fixes

- the git commit of the private research repository that held the engine and results,
- the SHA-256 of `git archive --format=tar <commit>` (the complete private source tree at that commit),
- the SHA-256 of each engine source file, the test suite, the scripts, the sweep and certification logs, and the
  planted-solution controls.

SHA-256 of the manifest itself: `efff8c5d592c48a79739cd427ee66ff5e54e1fd4df48fe0f90e8057a90b6fed2`

## How to verify

1. **The timestamp.** Install the OpenTimestamps client (`pip install opentimestamps-client`) and run
   `ots verify commitment/manifest.txt.ots`. It reports the Bitcoin block whose time bounds the manifest's existence.
2. **The released code and logs.** The tag `v1.0-as-committed` holds them exactly as committed:

   ```bash
   git checkout v1.0-as-committed
   python3 verify/verify_commitment.py
   ```

   The script recomputes the SHA-256 of every file listed in the manifest and prints OK for each one that matches. The
   private repository kept the engine in a folder named `agent_submission_v18/tripal`; here it is `engine/`, and the script
   maps the names. Nothing else is renamed.
3. **The later, commented engine.** Commits after the tag add comments to the C sources and documentation; they do not
   change the code. You can confirm this by compiling both versions with the same flags and comparing the object files,
   which are byte-identical.

The archive hash of the whole private tree cannot be checked from this repository, because that tree also held unrelated
working files that are not published. The per-file hashes cover everything the result depends on.

## Note on the public result files

One committed log, `results/certify_hex116_118_integer_only.jsonl`, is published as a copy with the machine's host name
removed, so its bytes differ from the committed original and the verify script reports it as REDACTED. Its numbers (nodes,
times, hits) are unchanged and also appear in [results/per_length.csv](results/per_length.csv), which matches the
commitment exactly. All other result files are byte-identical to the committed ones.
