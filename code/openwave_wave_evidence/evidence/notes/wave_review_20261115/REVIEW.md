# Rows 29/30: proof-and-replay review

**Verdict: the checked free-wave mathematics and finite scientific checks are
supported; the current bundle is not ready for reliable in-platform replay.**

Four reproducible packaging/verification findings should be resolved before
filing. One involves a genuinely different numerical readout, not merely
last-digit floating-point differences. No blocking error was found in the
checked conditional free Gaussian proof chain.

Reviewer: GPT-6 Astra Pro (OpenAI), user-requested review on 2026-11-15.
This session is not an appointed outside-group reviewer and awards no G6,
candidate pin or OpenWave cell. The original evidence remains unchanged.

## Subject and custody

- Tracking commit: `60519f4d50cd4b16f7eeaad1a2d1e92544a037dd`.
- Candidate: `MFR-unified-native-field-v1`, free preparation `c=0`.
- Bundle (`B` below): `shared/unified_wave_submission_20261115`.
- Outer manifest SHA-256:
  `1f8ac01dc06b2180ee565a4c7aa926a921220b6232758838d9ee102c7dd644a6`.
- Verified: 420 outer manifest entries, 409 copied evidence pins, all six
  original packet/review manifests, and every extracted Git blob ID.
- Original native endpoint SHA-256:
  `dd7c87bcda73f55990690effd0633020ff05b9131205d2314793329f8e0c811d`.

Read [PROOF_REVIEW.md](PROOF_REVIEW.md) for the argument-by-argument assessment
and [SUBMISSION_TASK_DRAFT.md](SUBMISSION_TASK_DRAFT.md) for the proposed next
platform task. The latter is a draft, not a registration.

## Reproduction results

The publisher recorded macOS/ARM64, Python 3.13.7 and NumPy 2.3.5.
This review used exact committed source bytes in isolated scratch copies:

| Check | Windows x86-64, Python 3.13.2 / NumPy 2.5.1 | Ubuntu WSL x86-64, Python 3.12.3 / NumPy 2.3.5 |
| --- | --- | --- |
| Committed-byte custody | Pass | Pass |
| Original top-level replay, assertions enabled | **Exit 1**, at the free replay | **Exit 1**, at the free replay |
| Four unmodified scientific leaf procedures, assertions enabled | All exit 0 | All exit 0 |
| Original byte-equality requirement | Fails | Fails |
| Original replay with inherited `PYTHONOPTIMIZE=1` | **Exit 0 despite failed byte comparisons** | **Exit 0 despite failed byte comparisons** |
| Deliberately wrong origin digest, ordinary custody | Rejected, exit 1 | Rejected, exit 1 |
| Same wrong origin digest, optimized custody | **Accepted, exit 0** | **Accepted, exit 0** |

The normal leaf executions reproduce the scientific assertions for 169,364
free, 347,533 orbital-ancestor and 492,837 unified native calls: **1,009,734
calls in overlapping lineage**, plus the boost analysis without native calls.
They do not pass the original wrapper's byte gate.

Raw logs, original/generated artifact digests and structural/numerical
comparisons are in:

- [Windows receipt](results/windows/REPRODUCTION.json)
- [Linux receipt](results/linux_numpy235/REPRODUCTION.json)
- [Independent finite review](results/INDEPENDENT.json)
- [Profile/readout diagnosis](results/PROFILE_DIAGNOSTICS.json)

The diagnostic `1e-10` comparison is labelled as such. It has not replaced,
relaxed or retroactively passed an original scientific or reproduction gate.

## Findings

### WR-01 — CRLF checkout breaks sealed custody

**Priority: reproduction setup blocker.**

In this clean Windows checkout, **398 of 409** evidence leaves differ from
their pins. Every one differs only by CRLF versus LF. The original command
fails immediately on `evidence/boost/CELL.md`.

The committed Git objects are intact and pass all manifests. Even `git archive`
can apply local line-ending conversion; the reviewer therefore disables that
conversion for extraction and checks every blob ID independently.

**Required repair:** distribute an immutable byte-preserving archive or document
and enforce a fresh checkout with line-ending conversion disabled. Protect
sealed files with a scoped attributes policy in the publication repository.
Do not regenerate historical manifests from a converted checkout.

### WR-02 — Floating serialization makes the advertised replay nonportable

**Priority: in-platform replay blocker.**

Locations: `B/replay.py:57–69`,
`B/evidence/free/scripts/task_common.py:37–41`, and
`B/evidence/orbital/scripts/task_orbital.py:20–23`.

Both environments pass source custody and every normal scientific leaf
assertion, but the free wrapper fails its receipt/archive byte comparison.
For the free packet the maximum numeric receipt discrepancy is approximately
`1.04e-14`; differences in its decompressed floating readouts are at most
`5.33e-15`. All discrete free history values agree exactly.

Likewise the decompressed orbital and unified archives retain exactly equal
discrete records. Their floating fields differ slightly, changing the embedded
raw-data hashes and compressed bytes. Pinning NumPy to the recorded 2.3.5
does not restore byte equality on Linux.

**Required repair:** freeze a portable reproduction policy that compares source
bytes and discrete records exactly, applies justified finite numerical
tolerances to derived floats, rejects nonfinite values and structural changes,
and checks embedded raw hashes against the corresponding regenerated archive.
Alternatively provide a genuinely reproducible complete numeric environment.
Keep these failed byte comparisons as the packaging amendment's ancestry.
WR-04 must also be addressed; simply increasing a tolerance is inadequate.

### WR-03 — Optimized Python disables verification

**Priority: verification blocker.**

Locations: `B/replay.py:32–35,63–69,79`; the free/orbital entry points and
their scientific scripts also use removable assertions.

With `PYTHONOPTIMIZE=1`, the unchanged top-level replay exits **0** while its
JSON explicitly reports `byte_replay: false`. The same environment is inherited
by subprocesses, disabling their scientific assertions as well.

A separate disposable-copy probe changes only
`ORIGINS.json.free.manifest_sha256` to 64 zeros. Ordinary custody rejects it;
optimized custody prints `"custody": "pass"`. The top-level explicit copied-leaf
check still runs, but it does not protect this disabled origin check.

**Required repair:** reject optimized execution or re-execute with optimization
disabled, remove inherited `PYTHONOPTIMIZE`, and enforce critical acceptance
conditions with explicit exceptions. Every child must have assertions enabled.
Retain the wrong-origin and failed-byte-replay cases as regression checks.

### WR-04 — An orbital dependency selects the readout by floating-point noise

**Priority: deterministic-observable blocker for the current unified bundle.**

Location:
`B/evidence/orbital/evidence/openwave/row_31_orbital_levels/work/2026-11-14_codex_production_field/probe.py:72–79,248,271–280`.

The internal reference is selected with:

```python
u = proj[:, np.argmax(np.diag(proj))]
```

Vertex transitivity makes all twelve exact diagonal entries equal to the
band's rank divided by twelve. Numerical eigensolver noise therefore chooses
the reference port. Independently reconstructing the stored readouts identifies:

| Positive band | Published reference port | Windows replay | Linux replay |
| --- | ---: | ---: | ---: |
| `5-sqrt(5)` | 7 | 7 | 7 |
| `6` | 5 | 7 | 7 |
| `5+sqrt(5)` | 6 | 10 | 2 |

Those inferred ports reproduce every respective four-channel ground readout
within `2.45e-15`. For example, a stored ground readout `0.8059537648223112`
becomes about `0.156520121200273` in both replay environments: a changed
observable, not a harmless last digit.

Spatial ground vectors also have arbitrary global signs. Their rank-one
projectors agree within `2.78e-15`; the sign differences should not be confused
with the internal reference change.

This affects the orbital receipt required by the unified wrapper. The
rows-29/30 free source uses stable full-band projector readouts; this finding
does not refute its KG dispersion or Maxwell quotient.

**Required repair:** choose the reference port by an explicitly fixed,
source-defined convention before the run, or publish complete band-projector
observables. Canonicalize the spatial ground-vector sign where a signed vector
is reported. Regenerate affected diagnostic readouts under a documented
amendment and preserve the original records and failed comparisons.

## Independent scientific checks

The reviewer code imports no producer result routines. It reads the pinned
primitive in isolation only to compare six endpoint cases against independent
floor/ceiling averaging.

| Checked object | Result |
| --- | --- |
| Geometric graph and native internal spectrum | Reconstructed from icosahedral coordinates; 30 seams; eigenvalue discrepancy below `4.45e-15` |
| Free branch census | All 37,908 cases agree with independent endpoint arithmetic |
| Complete actual free tape | All 131,072 currents, seeded seam/coin draws, occupations and clock records agree |
| Field and band readout | Local curl discrepancy below `8.89e-16`; spectral-projector readout discrepancy below `3.34e-15` |
| Fifteen finite population modes | Covariance discrepancy below `2.49e-14`; cycle-hit discrepancy below `4.18e-14` |
| Moment mechanism | Separate finite-state recursion at 64, 256, 1024 and 4096 events corroborates convergence of second/fourth moments; direct covariance sum agrees |
| KG error and speed bounds | 2,000 sampled parameter cases pass; this corroborates the separately reviewed analytic inequalities |
| Clock recurrence and costs | Twenty direct killed-chain solves agree with the finite-cycle formula; multiplier discrepancy below `3.06e-16` |
| Maxwell quotient | Positive rank-two Gram factorization checked in three directions; wrong reflection-sign control has a negative eigenvalue |
| Native changed-law population | Exact rational F12 contrast `1/4`; contact-adjusted ratios `3/4` and `1` |

The 100-million-event population calculations are not actual tapes or empirical
precision estimates. The full free-field claim relies on the written limits.
The independent script does not re-enumerate the 384 full-vector spot checks;
the normal original free leaf executes those.

## Clause readiness and next action

| Clause | Evidence assessment | Remaining step |
| --- | --- | --- |
| Maxwell recovered | Coherent conditional source-to-Gaussian-to-field proof, with actual finite record checks | Resolve replay findings; full outside-group proof/source-admissibility review |
| Transverse waves at common `c` | Positive rank-two limit and shared native symbol; finite contacts retained | Adjudicate the stated limiting observation/clock dictionary in-platform |
| Emergent KG matter waves | Native positive internal bands and controlled mass-shell limit; common ruler/clock | Same external review and platform execution; empirical species are a separate claim |

Recommended order:

1. Repair the byte-preserving publication/replay adapter, optimization handling
   and unstable reference selection under a documented amendment.
2. Confirm the repaired adapter on a clean second environment, including the
   recorded negative verification probes.
3. Obtain full outside-group review of this source and its limiting dictionary.
4. Have the author choose the candidate and publish the pin-change proposal,
   explicitly addressing its relationship to static-source PR #723.
5. Register and execute the two wave-row tasks prospectively, then request
   separate maintainer decisions.

The review narrows the immediate work to specific reproducibility and
adjudication tasks. It neither establishes the full OPH hypothesis nor finds
a general impossibility for observer-first constructions.

## Review-package checks

`verify_review.py` passes 100 consistency/tool checks, including exact bindings
to captured subprocess output and executed reviewer source snapshots,
nonfinite/structural comparator negatives and reference-port reconstruction.
The review's own seal is checked separately.

`git diff --check` passes. The repository-wide `tools/gate.py check` exits 1
with 335 workspace violations: 78 missing paths in the absent sibling
workspaces and 257 changed manifest hashes in this CRLF checkout. Its full
output is retained in `results/repository_gate.stdout.txt`. The gate's inputs
(`STATUS.json`, candidate registry, DAG and original manifests) were not changed
by this review, so these are not newly introduced row-gate failures.
