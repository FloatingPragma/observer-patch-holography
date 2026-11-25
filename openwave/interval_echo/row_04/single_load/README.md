# OpenWave M10.5 candidate packet

**Owner copy prepared for adjudication. Formal filing status remains HOLD.** Proposed source `69ad9f5eea6fa9b2553f08c4e49e5303b83114c9` is a direct child of the adopted M10.4 pin and has local tag `interval-echo-m10.5-one-load-proposal`. The particle under test is the bounded, gateway-prepared one-load excitation on blank memories. The home observer reads a port. This is an explicit change from the two-load state with a home-seam detector. The official task-selected replay has no anchor at the checked OpenWave base. The complete development grid was known before registration. [Candidate declaration](../../../../model/candidate_laws/interval_echo/openwave/row_04_compact/m10_5/CANDIDATE.md) · [Task and precise draw](../../../../model/candidate_laws/interval_echo/openwave/row_04_compact/m10_5/TASK.md)

| Deliverable | What it contains |
| --- | --- |
| [OpenWave patch](openwave.patch) · [Method note](METHOD.md) | Eight model-local files against `bab6818b4a278ffef083c1e86c12134e5ff14ab5`. Includes the canonical amendment, roadmap, briefing, task, method note, standalone replay and compact data. |
| [Tagged source bundle](source.bundle) | Proposed source and both refs. Its sole prerequisite is the public M10.4 commit `09f062f77f7f540f2783f128ae6314c5d8b818ea`. |
| [Draft PR and proposed cell](PR.md) | The requested decisions, scope, prior information and cell wording. |
| [All outcomes](all_settings.npz) · [Summary](survey.json) · [Historical pilot](pilot.json) | Native probabilities for 20 homes × 12 starts × 2 directions × 24 exposures × 4 arms × 14 bins. Array key `probabilities`. Directions −1, +1. Exposures 15 through 38. |
| [Audit](REVIEW.md) · [File hashes](SHA256SUMS) | Scientific checks, clean replay, public dependencies and outstanding qualifications on this version. |

## Replay

Use a clean RER checkout at the proposed pin with NumPy, SciPy and SymPy. The existing OPH unit tests also use pytest. An existing RER clone with the public M10.4 commit can import the local proposal with `git fetch /absolute/path/to/source.bundle refs/heads/candidate/openwave-m10.5-one-load:refs/heads/review/m105`. Create a detached worktree at the proposed pin to keep other work intact. In that source tree use these commands.

```sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_compact.m10_5.survey
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_compact.m10_5.experiment --pilot --record pilot-replay.json
python -B -m pytest -q model/candidate_laws/interval_echo/openwave/row_04_compact/m10_5
```

To replay the platform side, apply the patch to the stated OpenWave base in a clean checkout. Give `openwave/xperiments/m10_oph/research/scripts/m10_5_one_load.py --board` the RER file `model/executable_model/carriers/FEDERATION_LEVEL0.npz`. The script reads the board and the compact pilot receipt, reconstructs all four arms, and compares all 56 bins. The board remains in the author's repository.

## Official confirmation

The trusted OpenWave checkout supplies the first official main first-parent commit carrying the exact task bytes after the stated base. Use `--status /path/to/openwave` on the source experiment to inspect it. When present, `--confirm /path/to/openwave --record held.json` executes both implementations at its settings, checks every gate, and keeps an incomplete or rejected receipt if an eligible run fails. Output paths are exclusive. The source bundle and patch are local. Owner release authorization applies to the exact outgoing versions under the project's filing guidelines. [Filing rules](../../../../docs/OPENWAVE_SUBMISSION_GUIDELINES.md#9-filing-gates)
