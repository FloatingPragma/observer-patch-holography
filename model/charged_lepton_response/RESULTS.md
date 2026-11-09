# Charged-field readout results

The configured coherent program has a full-CAR realization and a controlled free Dirac limit.
All lattice valleys and quantum archives are retained. The interacting chiral limit and measured-mass attachment are unproved.

The phase below is calibrated to the muon/electron ratio. Tau is excluded from that calibration.

| Quantity | Action result | Experimental reference | Role |
| --- | ---: | ---: | --- |
| muon/electron | 206.768282700000 | 206.7682827 ± 0.0000046 | Calibration |
| tau/electron | 3477.441636283849 | 3477.23 ± 0.23 | Conditional comparison |

Reference: [2022 CODATA adjustment](https://physics.nist.gov/cuu/pdf/all.pdf).
The held-out tau ratio lies inside this reported one-standard-uncertainty interval.
No joint statistical significance or actual interaction-correction bound is claimed.

The exact rational phase, coefficient intervals and input ancestry are in [results.json](results.json).
The [calibration proof](CALIBRATED_SANITY.md) and [dependency map](PROOF.md) state the physical requirements.
The [source-phase audit](SOURCE_PHASE_SELECTION.md) distinguishes C/E duties from orientation reversal and the conditional pixel law.
Physically identified particles: **none**. Source-selected physical phase: **unproved**.

## Conditional C-duty spectrum

The complete memory construction realizes the supplied pixel-dependent law.
Its four coordinate branches give the following ratios. Each lightest coefficient is normalized to one.

| Declared branch | Middle/light | Heavy/light |
| --- | ---: | ---: |
| forward | 206.768300595181 | 3477.365536596023 |
| measured_endpoint | 206.768282739123 | 3477.365262268936 |
| mssm_beta4_proxy | 206.768297922262 | 3477.365495531309 |
| sm_color_beta7 | 206.782621310398 | 3477.644439943646 |

Only the measured-alpha endpoint lies inside both quoted individual ratio intervals.
The measured-alpha endpoint inherits a measured fine-structure input. Other branches retain their specified coordinates or conditional census and scale law.
The [capacity proof](P_DEPENDENCE.md) gives interval bounds and input ancestry.

The [candidate map](CANDIDATE_MAP.md) covers eighteen target channels and preserves the separate eleven-entry protected-object inventory.
Its [machine-readable map](candidate_map.json) assigns charged spectral modes as hypotheses; physical masses and defect-to-species links are unassigned.
