# V1 Protocol Erratum

Date: 2026-07-16

This erratum supersedes the historical V1 scoring procedure. It does not
rewrite the frozen v2 target or the historical runtime artifact.

Immutable records:

- `falsification/frozen_targets/hadronic_closure_target_2026-07-14_v2.json`
  has SHA-256
  `8dcefe93f124e21295b8e7cd85f2524ab17c65ad003c4fdac6205213d03ac6b2`.
- `runtime/ward_projected_payload_bracket_current.json` has SHA-256
  `395b1b6bc3d53ec40222e03004cced920236da3ed3b3e46cd407e71ab4852a15`.

## Historical V1 Disposition

Historical V1 cannot be promoted. Seven independent defects force that
disposition:

1. `run_bracket.py` contained target constants and a comparison tolerance.
2. The directing session had read target-side canon values.
3. The emitter evaluated the internal root
   `P = 1.630972095858897...`, while v2 declares
   `P = 1.630968209403959...`.
4. The old scoring text compared an emitted `Delta_source` total near 7 to 10
   with a residual coordinate near 0.041. Its S containment used the intended
   S coordinate, but the artifact remained a sampled grid envelope rather than
   a point payload or interval certificate.
5. V2 assigned one scalar source target to CL-1 and CL-2, although CL-2 adds
   `alpha_U(P)` and therefore requires a different source total.
6. The purported hadronic target silently treated the open electroweak
   remainder as zero, conflating `S_QEW_effective` with `S_hadronic`.
7. The present source-method choices were made after target access. Removing
   target numerals from the emitter prevents direct leakage but cannot make
   that method selection prospectively blind.

Numerical containment does not cure source contamination, coordinate
conflation, or an evaluation-point mismatch. The old artifact remains an
exploratory engineering record.

## Corrective Machine Contract

`falsification/frozen_targets/hadronic_closure_target_2026-07-16_v3.json`
is an inactive erratum scaffold that specifies the corrected coordinate and
map shape. It is not a retroactive or activatable target and is explicitly marked
`corrective_contract_not_yet_externally_timestamped_not_scorable`.

It corrects more than the total-versus-residual type error:

- `Delta_source_total` and
  `Delta_source_residual_vs_implemented` are different coordinate kinds. A
  residual is diagnostic and cannot substitute for a total.
- `S_QEW_effective = (Delta_hadronic + Delta_EW) / Delta_quark_naive` and
  `S_hadronic = Delta_hadronic / Delta_quark_naive` are different coordinate
  kinds while `Delta_EW = 0` remains unproved.
- CL-1 and CL-2 have independent completed maps. CL-2 adds `alpha_U(P)`, so
  one scalar target cannot close both rows.
- Point values in v3 are diagnostics only. An eligible payload must be a
  target-blind function or a certified interval enclosure over the registered
  P basin, with numerical and theory errors and a derivative or Lipschitz
  bound.
- A sampled grid's extrema are not an interval certificate.

V3 permits no promotion or falsification and must never be populated or
activated in place. An eligible future run needs a new immutable successor
whose complete method and dependency cone were frozen before genuinely
withheld data existed, or a documented clean-room team with no target access.
The successor must also preregister confidence, combined uncertainty, maximum
eligible prediction width, and three-way `COMPATIBLE`/`FAIL`/`INCONCLUSIVE`
rules. Historical V1 predates the correction and remains permanently
non-promotable.

## Repaired Process Boundary

The source process and comparison process are separate:

1. Place only source code and allowed source-side inputs in the emitter
   environment. Exclude `falsification/`, prior score files, and target-bearing
   prose.
2. Run `run_bracket.py` at its source-derived internal Stage-5 P. The program
   has no P override and rejects a measurement-located or target-supplied P.
   Its output is an uncertified singleton diagnostic, not an eligible closure
   payload.
3. Seal both the file hash and its embedded
   `content_sha256` before comparison.
4. Transfer the sealed artifact to a separate scorer environment.
5. Run `score_bracket.py` with the sealed artifact and corrective v3. Today it
   returns `NOT_EVALUABLE` with a nonzero status because v3 is inactive. Even
   under a future activated contract, the current singleton grid fails the
   certification gate.
6. Do not convert point-diagnostic containment into a verdict. A future scorer
   must solve and interval-certify CL-1 and CL-2 independently.

Example source emission:

```bash
python3 run_bracket.py --output /sealed/source-bracket.json
```

Example separate scoring step:

```bash
python3 score_bracket.py \
  /sealed/source-bracket.json \
  ../../../../falsification/frozen_targets/hadronic_closure_target_2026-07-16_v3.json \
  --output /sealed/source-bracket-score.json
```

The expected current result is `NOT_EVALUABLE`; there is no successful scalar
score. Future eligibility requires a newly frozen target version, blind or
clean-room method-selection chronology, an independently verified external
anchor preceding payload work, a clean target-free dependency cone, the full
receipt set named by v3, a certified P-domain source object, a closed
`Delta_EW` gate, a preregistered three-way decision policy, and separate
completed-map interval solves. The current scorer intentionally lacks that
solver; its receipt objects are schema-checked bindings, not proof that an
external signature, timestamp, or interval certificate was independently
verified.
