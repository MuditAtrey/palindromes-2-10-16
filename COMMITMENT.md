# Commitment and timestamp

The search engine, its test suite and all result logs are not public yet. Their exact contents are fixed by
[commitment/manifest.txt](commitment/manifest.txt), whose existence on 2026-11-21 is proven by the OpenTimestamps proof
[commitment/manifest.txt.ots](commitment/manifest.txt.ots) (anchored in the Bitcoin blockchain).

## What the manifest fixes

- the git commit of the private repository that holds the engine and results,
- the SHA-256 of `git archive --format=tar <commit>` (the complete source tree at that commit),
- the SHA-256 of each engine source file, the test suite, the sweep and certification logs, and the planted-solution controls.

SHA-256 of the manifest itself: `efff8c5d592c48a79739cd427ee66ff5e54e1fd4df48fe0f90e8057a90b6fed2`

## How to verify

1. **The timestamp.** Install the OpenTimestamps client (`pip install opentimestamps-client`) and run
   `ots verify commitment/manifest.txt.ots`. It reports the Bitcoin block whose time bounds the manifest's existence.
   (Bitcoin confirmation takes a few hours after stamping; until the proof is upgraded, `ots info` shows the calendar
   attestations it is waiting on.)
2. **The code, once released.** In the released repository, run
   `git archive --format=tar <commit> | shasum -a 256` with the commit named in the manifest. The output must equal the
   archive hash in the manifest. Each individual file can also be checked against its listed hash.

## Note on the public result files

The files in [results/](results) are copies of the committed logs with the machine's host name removed and the engine label
shortened, so their bytes differ from the committed originals. The originals are the ones fixed by the manifest.
