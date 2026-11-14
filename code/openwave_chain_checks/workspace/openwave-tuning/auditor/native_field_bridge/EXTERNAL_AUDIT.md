# Force and field source audit

The reviewed N20 solution correctly distinguishes the native accumulated
Green field from a response law. Its preparation, creation prior and
event/sweep clocks remain explicit. In particular its Theorem 5 uses a
supplied recorded pair creation. A fixed creation edge retains a cycle
component, so a homogeneous creation ensemble cannot be substituted for one
configured preparation. The new feedback construction in this directory
adds an executed response through a declared scheduler. It does not remove
N20's coupling or preparation inputs from the unchanged native law.

N20's small-momentum Laplacian coefficient is 1/8. For port density twelve
per unit cell the inverse-Laplacian coefficient is
(1/12)/(4*pi*(1/8))=1/(6*pi). This is a static Green normalization. The native
diffusion coefficients are 1/(16M) per event and nu/16 for a per-seam Poisson
rate nu. The number 1/sqrt(8) belongs to the separate coherent H_m clock.

## Literal recorded events versus recycling jumps

`NATIVE_JUMP_RATES.json` reexecutes the existing independent exact witness
`code/openwave_audit/native_jump_records.py` against the current
`N3D_cubic_routing/n3dlib.py::kraus_bare` function. The witness reconstructs
the local matrices and imports no producer runner or acceptance case.
All five checks pass, including the exact superoperator identity and actual
Kraus comparison.

For a unit load on a degree-six port and a per-seam rate r, the instantaneous
rates are:

| Recorded event | Rate |
| --- | ---: |
| All scheduled native events | M r |
| Actual load transfer, C or C/E mixture | 3 r |
| Nonzero literal C/E-mixture outcome | 9 r / 2 |
| Jump of the equivalent recycling unravelling | 9 r |

Literal C alone has nonzero-current rate 3r. The E outcomes include nonzero
values without a load transfer. The recycling operators give the same
unconditional dissipator, while changing the conditional record law.
Therefore H_eff=H-i(9r/2)P_port is the no-jump Hamiltonian of the recycling
unravelling. Its norm does not give the probability of no scheduled native
event, no nonzero C/E outcome, or no load transfer. A coherent state with
port weight w has initial recycling hazard 9rw and literal C/E nonzero
hazard (9r/2)w. The unconditional N3D Lindblad propagation is unaffected by
this distinction.

The N3D `commits` and `lindblad` artifacts use H_eff and labels such as
`undetected_norm`, `final_fraction_undetected` and `mean_free_path_cells`.
Those conditional labels must name the recycling instrument or be replaced
by a calculation for the literal instrument. An unconditional-generator
identity alone does not perform that replacement. The separate earlier
massive zero-mode numerical fix does not address this record distinction.

## Current worker evidence and common-source contracts

The latest read communication entries are Fable F009 and Codex C027, with
the joined row-one executor handed over in auditor A022. Fable's N31
directory `issues/N31_orbital_3d` is empty at this audit. G3/WQ-11 is queued
in `GAPS.md`; no G3 or WQ-11 receipt is present. Neither label is executed
evidence. These facts should be refreshed at integration if producers land
new files.

The TG package reproduction receipt records 5, 3 and 6 passing tests for
`openwave_autonomous_coulomb`, `openwave_common_electromagnetism` and
`openwave_current_coulomb`. Their owning proofs retain their supplied
Hamiltonians and domains. The first gives an autonomous eighteen-state
invariant physical sector of a specialized native gauge Hamiltonian and
an exact current-acceleration contribution. The current-Coulomb result
distinguishes a measured quadratic observable from a replaced energy law.
The common-electromagnetism construction adds explicit full integer rotors
and a quantum reference. These are useful constructive precedents; they
do not identify their Hamiltonians, states or clocks with the literal
native C process or the new feedback scheduler.

The common-source theorem at `model/particle_common_contract/PROOF.md`
requires actual source, preparation, reference and clock compatibility.
This review applies that requirement without demanding architectural
uniqueness. The declared feedback scheduler is a distinct finite
construction with a positive executed response. Its result does not
promote rows 20, 21, 29 or 31 to platform validation.
