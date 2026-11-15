# Portable wave replay amendment

This package implements repairs for all four findings in the
[rows 29/30 review](../../notes/wave_review_20261115/REVIEW.md), preserving the
original sealed evidence and failed byte replays. Read
[AMENDMENT.md](AMENDMENT.md) for the exact changes and their scope, and
[POLICY.json](POLICY.json) for the frozen numerical policy.

The amended adapter passes on Windows/Python 3.13.2/NumPy 2.5.1 and
Ubuntu WSL/Python 3.12.3/NumPy 2.3.5. Each environment executes **1,009,734
native calls across overlapping lineage**, with every original scientific
assertion enabled. The saved outputs also agree across these environments
under the amended policy. These are repeated histories, not independent
scientific experiments.

| Finding | Implemented repair |
| --- | --- |
| WR-01: CRLF breaks custody | Extract and check raw pinned Git blobs; scoped `.gitattributes` protects future checkouts; provided bundles must have exact bytes and file inventory |
| WR-02: floating serialization breaks byte replay | Exact discrete/structural checks; `1e-12` absolute plus `1e-12` relative float budget; independent archive/hash binding; byte mismatches still reported |
| WR-03: optimization removes assertions | Reject optimized entry; sanitize child environments; confirm child optimization flag; explicit adapter failures and retained process logs |
| WR-04: roundoff changes orbital reference | Fix internal reference at port 0 and positive spatial ground sign; independently reconstruct amended reference readouts from sealed records |

The adapter compares **4,621,142 exact discrete scalar values** across the
free, orbital and unified record archives per environment. Floating readouts
are compared separately. The orbital change is a declared diagnostic decoder
amendment, requiring disclosure in the author's candidate/pin proposal.

## Evidence and regression checks

- [Final Windows replay](results/windows_final/REPLAY.json)
- [Final Linux replay](results/linux_final/REPLAY.json)
- [Final cross-environment verification](results/FINAL_VERIFICATION.json)
- [Windows regression receipt](results/windows_final_tests/TESTS.json) and
  [Linux regression receipt](results/linux_final_tests/TESTS.json)
- [Input freeze](REGISTRATION.json), made before the final executions
- [Development attempts](attempts/): the initial Windows package-discovery
  failure and the subsequently discovered provided-bundle inventory gap

The final suite has 21 tests in each environment, including optimized entry,
assertion-enabled children, corrupted custody, extra unsealed files, altered
records with updated hashes, bad embedded hashes, nonfinite/structural
changes, wrong orbital readouts and retained timeout failures. A harmless
recompression control passes. Old successful development runs in
`results/windows`, `results/linux_numpy235` and `results/CROSS_ENVIRONMENT.json`
bind the previous input freeze retained in `attempts/02_source_inventory`.
Use the **final** receipts above for the delivered code.

The root `MANIFEST.sha256` seals this package, including all attempts and logs.
The review package retains its original separate seal.

## Run

Requires Git, Python 3.10+ and NumPy. From the repository root:

```sh
python -B -X utf8 shared/unified_wave_replay_amendment_20261115/replay.py --custody-only
python -B -X utf8 shared/unified_wave_replay_amendment_20261115/replay.py --output /path/to/new-run
python -B -X utf8 shared/unified_wave_replay_amendment_20261115/run_tests.py --run /path/to/new-run --output /path/to/new-test-logs
python -B -X utf8 shared/unified_wave_replay_amendment_20261115/verify_results.py /path/to/new-run --output /path/to/new-verification.json
```

Every output destination must be new. Keep new runs outside this sealed
package. To compare environments, supply both run directories to
`verify_results.py`. With an exact extracted copy of the original bundle,
`replay.py --bundle /path/to/bundle` can replay without Git; it rejects a
CRLF-converted working tree. No command rewrites the historical bundle.

## Submission boundary

This repairs the local replay path under a disclosed numerical policy.
It **does not pass the original exact-byte G4 rule** or supply outside-group
G6 review. The policy and orbital diagnostic amendment need adjudication,
along with the full written free-wave proof/source dictionary. The author
still selects the complete candidate and published pin; prospective platform
registration/execution and maintainer decisions follow. The
[updated submission handoff](NEXT_SUBMISSION_TASK.md) makes those steps concrete.

The repository-wide gate checker still has the pre-existing missing-workspace
and CRLF-manifest failures recorded in `results/REPOSITORY_CHECK.json`.
No row status, candidate registry or official mark is changed by this repair.
