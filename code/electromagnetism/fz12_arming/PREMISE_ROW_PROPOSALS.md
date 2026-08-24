# FZ-12 arming premise row proposals

This document proposes premise-register rows for the seven arming premises of
the frozen dispersion row FZ-12. It is a proposal pending registration: this
document registers nothing by itself, amends no register, and changes no frozen
byte. The rows below become premises only through the owner and integrator
acting on `tracking/premise_register.json` (schema `oph.premise_register.v3`).
The id values `PR-70` through `PR-76` are OWNER/INTEGRATOR SLOTS: the final ids
are assigned at registration against the register state at that time. The
highest registered id at the time of writing is `PR-69`.

## Source of the clauses

The comparison protocol of frozen row FZ-12 (`claims/frozen_prediction_register.json`,
row id `FZ-12`; mirrored in `docs/FROZEN_PREDICTION_LADDER.md`) declares
physical comparison ineligible and unarmed, and names the premises an eligible
comparison requires: a source-derived homogeneous position action with the
complete edge orbit as its sole direct support, equal source weights,
continuous field and same-operator sector attachment, cofinal gluing, finite
scale, coherent frame and boost transport, frozen nuisance and coverage rules,
and a post-custody dataset-specific contract; a null verdict additionally
requires the source-derived positive lower bound owned by issue #664. Those
premises exist in the corpus only as prose inside the frozen row. The Lean
module `Lean/Screen/DispersionArmingInterface.lean` types them; this document
proposes their register rows. The post-custody dataset-specific contract is
typed in the same module as an abstract token and is not proposed as a
premise row, because it is dataset specific and cannot be registered in
advance of a dataset.

Each `type` value is drawn from the register's six-value enum:
`empirical_import`, `external_mathematics`, `numerical_input`,
`representation_choice`, `selection_rule`, `structural_rule`. Each
`disposition` is `remove` (the premise is expected to be discharged by a
source theorem, and the row tracks that obligation) or `axiomatize` (the
premise is a declared choice recorded as such). Consuming lanes are the
FZ-12 arming lanes 733, 739, and 742.

## Proposed rows

### Slot PR-70 (OWNER/INTEGRATOR SLOT): source homogeneous edge position action

```json
{
  "id": "PR-70",
  "name": "source homogeneous edge position action",
  "type": "selection_rule",
  "disposition": "remove",
  "statement": "The physical position action of the compared regime is homogeneous, its sole direct support is the complete thirty-direction edge orbit of the committed seam current, and its source weights are equal across the thirty directions. The committed corpus proves the control counterpart only: the seam-derived control symbol has complete edge support with one common squared chord and equal control weight (source_seam_edge30_control_certificate), and no committed theorem attaches that control symbol to a physical action.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "Lean/Screen/SeamCurrentEdge30Moment.lean",
    "Lean/Screen/SeamCurrentHomogeneousAction.lean",
    "claims/frozen_prediction_register.json"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "Lean/Screen/SeamCurrentEdge30Moment.lean": "conditional_consumer",
    "Lean/Screen/SeamCurrentHomogeneousAction.lean": "conditional_consumer",
    "claims/frozen_prediction_register.json": "statement"
  },
  "notes": "Removal needs a source theorem selecting the homogeneous direct-only equal-weight edge action as the physical position action. The interface fields homogeneousPositionAction and equalSourceWeights type the clause; prop_fields_carry_no_separation proves their content is external input."
}
```

Disposition justification: `remove`, because the finite control packet is
committed and the clause is exactly the shape a source selection theorem
would discharge.

### Slot PR-71 (OWNER/INTEGRATOR SLOT): continuous field and same-operator sector attachment

```json
{
  "id": "PR-71",
  "name": "continuous field same-operator sector attachment",
  "type": "selection_rule",
  "disposition": "remove",
  "statement": "The compared observable lives in a continuous field sector, and the comparison reads the same operator sector that the committed source theorem fixes for the seam-derived symbol. No committed theorem selects the physical sector or identifies the committed symbol with a laboratory operator.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean",
    "code/a5_fingerprint/carrier_class_dispersion_certificate.py",
    "claims/frozen_prediction_register.json"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean": "conditional_consumer",
    "code/a5_fingerprint/carrier_class_dispersion_certificate.py": "conditional_consumer",
    "claims/frozen_prediction_register.json": "statement"
  },
  "notes": "Removal needs a source attachment theorem for the continuous field sector and operator. The interface field sectorAttachment types the clause as a required-to-hold proposition whose content is external input."
}
```

Disposition justification: `remove`, because sector selection is a source
obligation the frozen row itself names as unproved.

### Slot PR-72 (OWNER/INTEGRATOR SLOT): cofinal gluing into the compared regime

```json
{
  "id": "PR-72",
  "name": "cofinal gluing into the compared regime",
  "type": "structural_rule",
  "disposition": "remove",
  "statement": "The finite committed patch embeds cofinally into the compared regime, so that the committed edge-branch symbol governs propagation at the compared scales. No committed theorem supplies this gluing for the FZ-12 comparison path.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "claims/frozen_prediction_register.json"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "claims/frozen_prediction_register.json": "statement"
  },
  "notes": "Removal needs a machine-checked cofinal gluing statement covering this clause. The interface field cofinalGluing types the clause; its content is external input."
}
```

Disposition justification: `remove`, because gluing is a structural source
obligation rather than a declared convention.

### Slot PR-73 (OWNER/INTEGRATOR SLOT): declared finite carrier scale

```json
{
  "id": "PR-73",
  "name": "declared finite carrier scale",
  "type": "numerical_input",
  "disposition": "axiomatize",
  "statement": "One declared positive carrier scale fills the committed coefficient slot, giving C4 = -a^2/20, B0 = a^4/840, and B6 = -a^4/12600 at that scale. The committed manifold ratios are scale free: they hold at every positive scale (edgeScale_ray), no positive scale is singled out (committed_ray_selects_no_scale), and the scale-boundary module proves the metric coefficient counterfamily, so the value is external declared input.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "Lean/Screen/PrimitivePortScaleBoundary.lean",
    "claims/frozen_prediction_register.json"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "Lean/Screen/PrimitivePortScaleBoundary.lean": "no_go",
    "claims/frozen_prediction_register.json": "statement"
  },
  "notes": "Declared choice recorded as declared. A later source scale theorem would convert this row to disposition remove; the committed corpus proves nonforcing of the scale, not a value. No scale number is proposed here."
}
```

Disposition justification: `axiomatize`, because the committed corpus proves
the ratios are scale free, so the value is a declared numerical input until a
source scale theorem exists.

### Slot PR-74 (OWNER/INTEGRATOR SLOT): coherent frame and boost transport

```json
{
  "id": "PR-74",
  "name": "coherent frame and boost transport",
  "type": "structural_rule",
  "disposition": "remove",
  "statement": "One coherent frame carries the comparison, and a boost transport connects that frame to the source frame of the committed symbol, so that the intrinsic coefficients are read in a resolved frame. No committed theorem supplies frame resolution or boost transport for the FZ-12 comparison path.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "claims/frozen_prediction_register.json"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "claims/frozen_prediction_register.json": "statement"
  },
  "notes": "Removal needs a transport theorem on the completed carrier connecting source and comparison frames. The interface field frameBoostTransport types the clause; its content is external input. The frozen decision rule marks unresolved-frame outcomes INCONCLUSIVE."
}
```

Disposition justification: `remove`, because frame transport is a structural
source obligation the frozen decision rule depends on.

### Slot PR-75 (OWNER/INTEGRATOR SLOT): frozen nuisance and coverage readout rules

```json
{
  "id": "PR-75",
  "name": "frozen nuisance and coverage readout rules",
  "type": "representation_choice",
  "disposition": "axiomatize",
  "statement": "The readout of any eligible comparison isolates nuisance parameters and calibrates joint coverage under rules fixed before exposure, as the frozen decision rule requires for every verdict class. These are declared analysis conventions of the comparison protocol, not source-derivable structure.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "claims/frozen_prediction_register.json",
    "docs/FROZEN_PREDICTION_LADDER.md"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "claims/frozen_prediction_register.json": "statement",
    "docs/FROZEN_PREDICTION_LADDER.md": "statement"
  },
  "notes": "Declared choice recorded as declared. The interface field nuisanceCoverageRules types the clause as a required-to-hold proposition whose content is external input. No coverage rule content is proposed here."
}
```

Disposition justification: `axiomatize`, because readout, nuisance isolation,
and coverage calibration are declared conventions of the comparison protocol.

### Slot PR-76 (OWNER/INTEGRATOR SLOT): source-derived positive exclusivity lower bound

```json
{
  "id": "PR-76",
  "name": "source-derived positive exclusivity lower bound",
  "type": "numerical_input",
  "disposition": "remove",
  "statement": "A pre-exposure same-action, same-sector source theorem fixes a positive lower bound a_min for the admissible scales in the comparison units, and the preregistered sensitivity and remainder contract covers every admissible scale at or above a_min. The frozen append-only clarification makes a null verdict conditional on exactly this obligation, owned by issue #664. The committed corpus supplies no such bound: every positive value is realized as the lower bound of a stipulated interface inhabitant (committed_corpus_selects_no_lower_bound), so a stipulated value is inadmissible for a null verdict.",
  "consuming_lanes": [733, 739, 742],
  "evidence": [
    "Lean/Screen/DispersionArmingInterface.lean",
    "claims/frozen_prediction_register.json",
    "docs/FROZEN_PREDICTION_LADDER.md"
  ],
  "evidence_roles": {
    "Lean/Screen/DispersionArmingInterface.lean": "statement",
    "claims/frozen_prediction_register.json": "statement",
    "docs/FROZEN_PREDICTION_LADDER.md": "statement"
  },
  "notes": "Removal needs the source lower-bound theorem the frozen clarification names. Disposition remove rather than axiomatize because the protocol itself rules out a declared value for null verdicts; the interface fields exclusivityLowerBound and lowerBoundSourceDerived type the obligation."
}
```

Disposition justification: `remove`, because the frozen protocol requires this
bound to be a source theorem and forbids substituting a declared value in a
null verdict.

## Boundary

This document and the Lean module it cites type the arming path of frozen row
FZ-12 and propose register rows; nothing more. The frozen bytes of FZ-12 are
immutable and untouched. No arming occurred and none is claimed: no premise is
discharged, no comparison is scored, no frozen record is created, and no
comparison data, bound value, or grain-scale number enters this document or
the Lean module. Whether any proposed row can be discharged from source is
open. Registration, id assignment, and any edit to
`tracking/premise_register.json` are owner and integrator actions outside the
scope of this document.
