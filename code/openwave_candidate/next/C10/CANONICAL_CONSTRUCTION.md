# One C10 construction and readout contract

The canonical entrypoint is `canonical_construction.py`. It uses `JointNativeSource` from `joint_source.py`, which composes the existing boundary source and periodic field adapter on the same state space. All completed readouts below use that one interacting Hamiltonian and its recomputed global reference. Earlier bare, uniform, character and standalone field experiments are controls, not imported measurements.

From the public bundle root:

```sh
OPENBLAS_NUM_THREADS=1 python next/C10/joint_source.py
OPENBLAS_NUM_THREADS=1 python next/C10/canonical_construction.py --write
OPENBLAS_NUM_THREADS=1 python next/C10/canonical_construction.py --check
```

The first command reproduces the complete interacting reference, all-number spectra, charged-reservoir responses and native instrument controls. The second selects one protocol and writes its complete readouts and source contract. The third recomputes the selected protocol, compares its scientific values and poststate, and verifies the source contract. No Fable acceptance flag or grid is invoked. The runtime uses existing portable CAR providers; the public export manifest binds the complete scientific bundle.

## State, Hamiltonian and reference

There are six CAR orbitals: three ports and three transit positions. Every source occupation number0..6 is retained. Six shared Z6 half-records lie in an affine Gauss fiber, with the exact integer boundary record `B=3−Nsource`. This gives384 source basis states and768 coordinates with the two-occupation reservoir. The fixed source Gauss background is `(111000)`; the reservoir quantity `nr−B` is conserved. The state labels and complete H are retained in `joint_source_states.npz`.

The one Hamiltonian is CAR Hm, periodic-half-record electric energy and one cyclic face term, at m=.23, kappa=.37 and g=.61. No uniform or character compression is used after the periodic energy is attached. The actual global all-number ground is unique and lies in N3, with energy −2.7262233938630343. Its normalized projector is the canonical reference. The separately reported fixed-N3 reference coincides here; it is a control equality, not a second selected source or an assumed transfer. All-number spectra and every charged energy response are recomputed from this interacting H.

Setting kappa=g=0 gives exactly the same384-state native CAR half-hop Hamiltonian. This is the baseline relation checked by the canonical entrypoint. It does not make bare C9 spectral gaps equal to the interacting gaps.

## Clock and selected native protocol

The source clock is the original dimensionless half-hop parameter. The selected finite protocol evolves for .37 in that parameter, then applies the complete MIX native event with probability .29; no commit has probability .71. MIX draws one of the three seams uniformly, C/E fairly and the native coin fairly, grouping only indistinguishable local outcomes. It does not measure all positions or individual records. These schedule values are supplied; no physical time conversion or derived continuum commit rate is claimed. Because the reference is an exact coherent eigenstate, it remains stationary during the coherent interval. Monitoring generally excites it and does work.

The canonical receipt keeps every native outcome and the no-commit complement, their probabilities, branch energies and the complete poststate. Standalone C and E remain component controls in the joint receipt. They are not additional simultaneous channels in the selected protocol.

## Readouts on this source

| Completed observable | Actual result |
| --- | --- |
| All source-number spectra and interacting global reference | `JOINT_SOURCE.json`: `actual_all_number_spectra`, `references.global_allN` |
| Complete dressed addition/removal/no-change response, initial reservoir labels and final number/energy eigenspaces | `JOINT_SOURCE.json`: `references.global_allN.complete_reservoir_response` |
| Native C/E/MIX monitoring work on the same reference | `JOINT_SOURCE.json`: `references.global_allN.native_instruments` |
| Charged pulse's unchanged unconditioned field energy | `JOINT_SOURCE.json`: `references.global_allN.field_energy_before_reservoir` and `field_energy_after_reservoir` |
| Selected split-step probabilities, energy work, occupations, whole-region charge and current | `CANONICAL_CONSTRUCTION.json`: `readouts` |
| Complete selected poststate | `canonical_states.npz`: `selected_protocol_poststate` |
| Source and result ancestry | `CANONICAL_SOURCE_CONTRACT.json` |

Whole-region charge here counts every occupation at port0 and the whole transit0; it is not the separate half-transit effective source observer. Every integer charge label0..6 is retained, including zero outcomes. Position and same-H current apparatus use the actual full basis. Energy measurements resolve complete degenerate eigenspaces, not arbitrary internal eigenvectors. The source-number measurements use actual occupations; no generic integer recovery from insufficient modular residues is assumed.

## Explicit boundaries

The charged apparatus is the complete CAR reservoir swap dressed by B, with both initial reservoir settings and all four outcomes retained. The supplied reservoir energy is .23 nr; extra boundary energy is explicitly zero. Apparatus work is computed for that declared source/reservoir budget, not an autonomous battery or full physical preparation cost. The pulse has no derived duration.

The selected Gauss background is not automatically charge-conjugation invariant. Modulus6 remains a finite periodic register choice. Statistics, transport phase, field coefficients, ideal spectral reference preparation, readout projectors and the split-step schedule are supplied. The reference is a finite mathematical ground, not an identified physical vacuum. No physical electron/quark/neutrino species, spin, calibrated charge unit, GeV mass, outgoing radiation, Lorentz symmetry or gravity is inferred. No official OpenWave row is promoted by this contract. The contract gives subsequent work one source to extend and one place to attach new observables.

## Continuing grouped spectral readout

The complete `(source number,energy)` readout is explicitly a grouped Lüders PVM. Every branch, including zero probability, is retained as a shared sector basis plus a small full complex subnormalized coefficient matrix in `canonical_states.npz`. Reconstruct it as `Embed(indices,basis @ coefficient @ basis†)`. Coherences within a degenerate eigenspace are preserved; no hidden rank-one measurement or forgotten basis label is substituted. The pre-readout `joint_source_states.npz` labels remain unchanged.

The receipt continues one actually populated degenerate branch from an archived joint output under the same H and selected MIX protocol, retaining its original subnormalization and complete joint native outcomes. An explicit coherent preparation in the actual N3 doublet at energy −.7613052173 rejects fine-basis measurement-and-forgetting: its poststate differs from the grouped Lüders poststate by trace distance .5. Basis rotation accompanied by the inverse coefficient rotation reconstructs the same state. This spectral PVM and coherent control are added apparatus choices, not operations claimed to be synthesized by native C/E.
