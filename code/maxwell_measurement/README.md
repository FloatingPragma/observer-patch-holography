# Serial Maxwell measurement adapter

This package prepares and checks captures of one fixed Maxwell **emulator**.
Its bounded observer-like patches have 42 local potential registers, ports,
baseline and response records, destructive pair-average probes, and feedback
that restores the measured registers. Two gauge executions each contain three
slices, 60 probes per slice, 180 feedback cycles and 585 semantic events.

No apparatus, physical raw data, frozen experiment or observed result is shipped.
The adapter does not identify encoded values with electromagnetic fields in
nature. It does not close a physical realization or calibration obligation.

## Draft and freeze

From the repository root:

```sh
python3 code/maxwell_measurement/contract.py --output /tmp/maxwell-preregistration.json
```

The command writes a `DRAFT`, refuses to overwrite an existing contract, and
prints its hash. It pins the adapter producer, verifier and adversarial tests;
a changed implementation invalidates a frozen contract. Supply an actual
study/run identity, device/firmware identity,
capture kind, calibration references, signed gain and offset intervals, raw
uncertainty, and a separately calibrated hardware clock. The 42 oriented channel
maps must describe the apparatus. Negative gain requires negative orientation;
an edge channel then reverses its physical endpoint order. Raw units are `V`
or `raw_count`, and transformed quantities are dimensionless model potentials.

Each channel calibration contains at least two distinct external reference
samples, recorded before the freeze. A reference has `id`, `raw`, `model_value`,
`captured_utc` and `source: "external_reference"`. Its reference-certificate hash
is retained. No episode sample, including the held-out third slice, may enter
calibration. A hardware operator must inspect the actual references; declaring
their hashes is not a verification of a calibration laboratory.

One unknown gain and offset must explain **all** readings on a channel:
external references, the three controls, both gauge executions and every
slice, including the held-out slice. The raw-error displacement may differ
between readings. For supplied raw value `r_i` and required model value `y_i`,
the contract is `y_i = g (r_i + delta_i) + b`, with one `(g,b)` in the frozen
gain/offset box and `|delta_i| <= raw_error` for each reading. Calibration drift
between readings is not part of this contract.

The verifier checks this common feasibility with exact rational arithmetic.
For the fixed gain sign `s` and raw-error bound `e`, each reading gives
`y_i - (r_i + s e) g <= b <= y_i - (r_i - s e) g`.
Include the frozen offset bounds as constant lower/upper lines. A common
offset exists exactly when every lower line is at most every upper line at
the same gain. Eliminating the offset therefore gives linear bounds on `g`;
their intersection with the frozen gain interval must be nonempty. Inclusive
endpoints preserve valid singleton and line-segment solutions, for either
gain sign. No floating-point tolerance or fitted midpoint is used.

This is an additional rejection check, not a calibration estimator. The
verifier never replaces the declared gain/offset intervals with the feasible
subset or reports a selected gain/offset pair. All coordinate, field, action
and control uncertainty still uses the full frozen box. In particular,
episode and held-out readings cannot tighten their own acceptance bounds.

The supplied design ceilings bound coordinate uncertainty by `1/1000`, field
uncertainty by `1/20`, and field/source-action uncertainty by `1/2` in the
declared model units. They may be tightened. Disabled-feedback controls must
resolve at least `1/4` of model potential. These are acceptance-design choices,
not measured device specifications. Freeze the completed file before capture,
retain its bytes externally, and pass its digest separately to the verifier.
Changing a threshold after capture requires a new experiment.

## Raw capture interface

The closed capture schema is `oph.maxwell_measurement.capture.v1`. Its fields
are `study_id`, `run_id`, `capture_kind`, `device_id`, `firmware_sha256`,
`captured_utc`, `controls`, and `executions`, in addition to `schema`.
UTC timestamps use ISO spelling ending in `Z`. Every numerical measurement,
interval endpoint and raw device timestamp is a canonical rational string,
such as `"-3/2"`, with at most 200 characters; JSON floats are rejected.
Exact arithmetic may produce longer output endpoints. Gauge comparisons use
the retained rational intervals before report serialization, so the external
input limit does not restrict internally computed field or action precision.
Output numerator and denominator integers use exact `Decimal(int)` fixed
formatting, preserving the caller's decimal context and the process-wide
integer-string digit limit even when accumulated results exceed that limit.

The three ordered controls are `sham_no_probe`, `pair_mean`, and
`feedback_disabled`. Each has `name`, `start_tick`, `end_tick`, and four
`samples`, whose roles are `before_left`, `before_right`, `after_left`,
`after_right`. Inputs are model potentials `-1/2` and `3/2` on `slot0` and
`slot12`. The sham preserves both values. The mean and disabled-feedback
controls end at `1/2` on both ports. The latter must visibly fail restoration.
This is a calibrated held-out control, not a calibration-independent ratio:
different channel gains do not cancel.

The two executions have `gauge: false` and `gauge: true`, in that order. Each
has 585 `events`. An event has `trace`, `start_tick`, `end_tick`, and `samples`.
Its `trace` uses the serial parent's exact operation/argument/read/writer/write
schema. It must represent actual device consumption. Copying the software
parent is not physical evidence. Every logical read needs the register,
writer identity and consumed value version; claimed parents cannot invent
reads. An independent witness must bind this trace to the apparatus capture.

Every sample has exactly `role`, `channel`, `raw`, `tick`, and a unique `nonce`.
Samples occur strictly inside their event windows, in the order below; event
windows and the preceding controls cannot overlap. Clock uncertainty enters
these inequalities. The required roles are:

| Operation | Samples, in order |
| --- | --- |
| `seed`, `advance` | `loaded_0` through `loaded_41` |
| `baseline` | `input`, `retained_baseline` |
| `probe` | `before_left`, `before_right`, `after_left`, `after_right` |
| `response` | `response`, `retained_response` |
| `feedback` | `restored_left`, `restored_right` |
| `decode` | `decoded_0` through `decoded_41` |
| Other operations | Empty sample array; full trace retained |

`slot0` through `slot11` encode node potentials; `slot12` through `slot41`
encode the thirty canonically oriented edge potentials. Baseline and response
records use their parent-node channel's transfer. Capture adapters must measure
those retained records under that transfer, or declare an apparatus where this
encoding is actually implemented. Channel assignments cannot be inferred from
captured values. The complete population has 3,540 raw samples.

The verifier reconstructs potentials from the raw baselines/responses using
`A = 2 response - baseline`. It checks the independently captured public decode,
both endpoint reconstructions, held-out slice 2, Maxwell residuals and the
field/source action. Exact rational interval arithmetic propagates calibration,
offset and raw-error bounds into `E = -2 ΔA - D phi` and `B = C A`, followed by
the action. The immutable parent supplies the action step `h = 1/2` and charged
paths. The verifier freshly executes its independent proof on every call.
These currents are declared source registers, not independently measured
physical matter currents. Agreement of a submitted trace with the parent
cannot by itself attest that hardware consumed those values; that separate
fact is part of the externally witnessed scope.

Clock calibration interval width may be at most `1/1000` of its positive lower
endpoint. Physical-device decode durations are reported in seconds under the declared
clock calibration. They are not identified with the model action parameter,
the feedback count, a Lorentz clock or signal travel time. No common-clock or
physical spacetime claim follows from monotone captures.

## Verification and external trust

```sh
python3 code/maxwell_measurement/verify.py \
  --preregistration /path/to/frozen-preregistration.json \
  --capture /path/to/capture.json \
  --preregistration-sha256 OPERATOR_RETAINED_PREREGISTRATION_HASH \
  --capture-sha256 OPERATOR_RETAINED_CAPTURE_HASH
```

Synthetic captures return `SYNTHETIC_CONFORMANCE_ONLY`. Hardware captures
without external attestation return `CAPTURE_CONFORMS_UNATTESTED` with exit
code 1. Malformed or nonconforming captures return `INVALID` with exit code 2.

For independently witnessed data, also supply `--operator-attestation` and
`--operator-attestation-sha256`. The operator retains this hash independently
of the claimant. The closed `oph.maxwell_measurement.operator_attestation.v1`
packet names the same study/run, device/firmware, preregistration/capture hashes
and UTC dates. It also names a distinct `operator_identity` and
`independent_witness_identity`; either may name a person or an organization.
It does not require three institutions. Its `scope` must explicitly attest
all these predicates:

- `physical_device_and_raw_capture_witnessed`
- `actual_read_write_consumption_witnessed`
- `complete_run_and_control_population`
- `firmware_and_device_identity_checked`
- `calibration_references_checked`
- `calibration_independent_of_episode_and_heldout_data`
- `preregistration_precedes_capture`
- `clock_reference_and_capture_timestamps_checked`

The result then says `CAPTURE_CONFORMS_RELATIVE_TO_OPERATOR_ATTESTATION`.
This is conditional on the operator's externally chosen trust declarations;
the adapter does not authenticate an organization, verify a witness signature,
or detect a coalition that supplies false external attestations. It never
certifies physical truth or a complete physical episode. The existing class-H
hardware evidence verifier is a separate, unchanged policy boundary.
A synthetic capture is ineligible for attestation even if a trust file is
supplied. No automatic change is made to an instrument or observation ledger.

## Controls

```sh
python3 -m pytest -q code/maxwell_measurement/test_measurement.py
```

Tests create explicitly synthetic captures in temporary directories. They
exercise signed transfers, exact interval field/action reconstruction, strict
schemas and fresh parent replay. Coherently rehashed mutations cover omitted
channels/events/probes, changed held-out coordinates, forged writers/parents,
extra synthetic reads, missing feedback, wrong current, calibration leakage,
uncertainty inflation, clock substitution, nonce reuse and capture reordering.
No synthetic capture is committed as a physical run or valid observed receipt.

Shared-calibration controls reproduce the false acceptance present at
`efa92897`: individually compatible readings could require different gains.
For example, with gain `[1, 1 + 10^-9]`, offset `-7/10` and zero raw error,
`(raw, model) = (7/10, 0)` forces gain 1, while
`(17/10, 1 + 17/10^10)` forces gain `1 + 10^-9`. Each individual interval
contains its target, but the pair is impossible. Full 3,540-sample controls
also put incompatible readings in different slices or gauges, while all
unchanged readings and references admit either endpoint. Feasible companions
must still pass. An independent exact boundary-intersection oracle checks the
feasibility solver against the original raw-error equations, including
negative gains, degenerate boxes and pairwise-compatible but jointly
incompatible triples.

An additional valid-input regression at `efa92897` replaces only six channel
raw-error bounds by `1/1000000007`, `1/1000000009`, `1/1000000021`,
`1/1000000033`, `1/1000000087`, and `1/1000000093`. Each 12-character input
tightens the original `1/1000000000` bound, and the unchanged singleton gains
and offsets explain every reading with zero displacement. The exact action
endpoints have 277 and 275 characters. Reapplying the external 200-character
input limit during the final gauge comparison would return `INVALID`, so
the verifier retains its exact intervals for both gauge comparisons.
API and CLI controls compare all decoded coordinates, fields, action endpoints
and device durations against an independent original-input scalar center/radius
calculation. A separate case crosses the field-endpoint length boundary;
another uses full 200-character inputs whose action integers exceed Python's
default 4,300-digit formatting ceiling. An independent base-10 chunk formatter
checks those outputs, including under restricted decimal contexts. Overlong
external calibration and capture values must be rejected.

The verifier rejects inconsistent synthetic or hardware-emulator submissions
and admits conforming captures with long exact intermediate results.
Contracts frozen against other adapter bytes fail the adapter source pin.
A future run needs a contract frozen against these bytes, not a replacement
hash attached to another experiment. The live source-current inventory
records this adapter as
downstream-contaminated with no shipped capture. The immutable serial parent,
voltage decoder, paper claims and physical-realization obligations are
unaffected by this verifier, which belongs to
[evidence audit #1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033).
