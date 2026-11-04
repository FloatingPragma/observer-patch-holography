# Formal witnesses for the conditional matter action

[PROOF.md](PROOF.md) composes the finite witnesses below with a local
field-action grammar and explicit source-sector matching conditions.
The listed Lean statements prove their typed conclusions. The native
matrix identities use exact executable algebra. The invariant-action
classification, singular-value argument and spectral-transfer estimates
are analytic arguments in the proof, rather than additional Lean claims.

## Gauge and matter representation

The input structure is
`OPH.SMStructureAdequacySurface.SMStructurePremiseData` in
[SMStructureAdequacySurface.lean](../../Lean/Screen/SMStructureAdequacySurface.lean).
It supplies barycentric vertex share, constant normalized face weights,
the compact-classification dimension menu, center restrictions, and a
nonempty chiral anomaly-free mask of the ten-component exterior grammar.
The compact-classification interpretation of the dimension menu is an
input; the integer enumeration does not construct a Lie-group isomorphism.

| Exact Lean witness | Conclusion on that input |
|---|---|
| `OPH.SMStructureAdequacySurface.SMStructurePremiseData.forced_gauge_type` | Center dimension one and simple-ideal dimensions three and eight. |
| `OPH.SMStructureComposition.smStructure_two_coarse_masks` | The selected coarse family/mask pairs are the G label with either conjugate parity mask. The coarse label does not specify a continuous bracket-family point. |
| `OPH.SMStructureComposition.smStructureComposition_receipt` | Threads the same bundle through coarse selection, dimension count, anomaly forms, parity grading and central kernel. |
| `OPH.SMStructureAdequacySurface.SMStructurePremiseData.masked_generation_count` | The selected mask has total complex dimension 15. |
| `OPH.SMStructureAdequacySurface.SMStructurePremiseData.masked_anomaly_forms` | Its four charge-weighted anomaly forms vanish. |
| `OPH.SMStructureAdequacySurface.SMStructurePremiseData.masked_z6_kernel` | The common stabilizer of the selected component weights is the tensor-action kernel. |
| `OPH.SMStructureAdequacySurface.smStructure_generation_weights` | Identifies the five even-sector weights with the committed matter weights. |
| `OPH.ExteriorSelection.witten_automatic` | Each surviving parity mask has an even number of weak doublets. |
| `OPH.SMStructureAdequacySurface.smStructure_global_form_menu` | The locally compatible central quotients form the four-subgroup menu of Z6. |

The composition statements are in
[SMStructureComposition.lean](../../Lean/Screen/SMStructureComposition.lean);
the finite 1024-mask scan is
[ExteriorSelection.lean](../../Lean/Screen/ExteriorSelection.lean).
`smStructure_loop_kernel_intertwiner` is an abstract loop/character
identification. It does not instantiate a physical line-operator algebra
or select the maximal quotient.

For a concrete representation,
[native_matter_attachment.py](../common_source/native_matter_attachment.py)
consumes the declared architecture in
[native_sm_current.py](../common_source/native_sm_current.py).
`representation()` checks all 66 current brackets after the declared
five-dimensional lift and restriction to `Lambda2(C5)+Lambda4(C5)`.
`decoded_multiplets()` computes joint charge and Casimir records;
`global_cover_contract()` computes their central kernel. These supply
the concrete representation dictionary under that architecture. Physical
Spin and the global line/loop realization are separate attachments.

## Conditional family count and Weyl conventions

The extended structure
`OPH.SMStructureAdequacySurface.SMThreeFamilyPremiseData` supplies a
single complete band from the four-band adjacency menu, faithful action,
dimension in `[3,5]`, and strict minimality in the specified seam-cost order.

| Exact Lean witness | Conclusion and required attachment |
|---|---|
| `OPH.SMStructureAdequacySurface.SMThreeFamilyPremiseData.band_identified` | The admissible minimum is `("3",3,(5,-1),1)`. The comparison rule is supplied; identification of its minimizer is proved. |
| `OPH.SMStructureAdequacySurface.SMThreeFamilyPremiseData.three_family_receipt` | The band has rank three; the selected mask has 45 components in three copies, vanishing charge-weighted anomalies and the same central action. |
| `OPH.A5FamilyBand.family_projector_refinement_natural` | The committed finite coordinate transports preserve the family projector, without asserting an arbitrary continuum refinement theorem. |
| `OPH.A5FamilyBand.factor_transports_commute` | Family transport and matter action operate on separate factors. An isometry from the physical family carrier to this factor is an input. |
| `OPH.MatterGrammarIndexBridge.matterGrammarIndexBridge_receipt` | The two ten-component grammars agree in charge, color/doublet indicators, conjugation and Weyl parity through an explicit permutation. |
| `OPH.MatterGrammarIndexBridge.selected_components_weyl_parity` | Relates the bundled mask to that index dictionary, without deriving spacetime chirality. |

The concrete projectors are in
[A5FamilyBand.lean](../../Lean/Screen/A5FamilyBand.lean), and the index map is
in [MatterGrammarIndexBridge.lean](../../Lean/Screen/MatterGrammarIndexBridge.lean).
The multiplicity follows from the complete-band conditions, rather than
three distinct Yukawa singular values or three topological windings.

[WeylYukawaConventions.lean](../../Lean/Screen/WeylYukawaConventions.lean)
provides `OPH.WeylYukawaConventions.weyl_up_neutral`,
`weyl_down_neutral`, `weyl_lepton_neutral`, `right_charge_forced`, and
`epsilon_SL2_invariant`. These fix all-left charge conjugation and the
weak epsilon contraction. `mixed_convention_not_neutral` rejects the
inconsistent convention. They do not select a spacetime spin representation.

## Same-current Yukawa tensors and normalized masses

The executable witnesses are functions of the same native current:

| Exact-algebra witness | Content |
|---|---|
| `native_matter_attachment.invariant_higgs_lines()` | Solves all twelve current-invariance equations on the matter-pair/Higgs grammar; each of three allowed charge origins has one invariant line. |
| `native_matter_attachment.action_intertwiner()` | Exact 15-row permutation, signs and tensor scales joining those lines to the field action. |
| `native_matter_attachment.kinetic_attachment()` | Both Higgs components, all 27 complex Yukawa symbols, the 90-row unitary map, three Weyl principal matrices and antisymmetric CAR pairing. |
| `native_matter_attachment.relative_spectrum_freedom()` | The invariant tensors allow independent complex family matrices; the representation fixes no singular-value hierarchy. |
| `conditional_standard_model.spectrum.mass_kernel()` | The retained 45-row action with directed conjugate Yukawa blocks and transposes at `H=(0,vEW/sqrt2)`. |
| `conditional_standard_model.spectrum.singular_spectrum()` | The complete family Gram spectrum, including zero and repeated singular values. |

The normalization agrees with
[calibrated_matter.py](../common_source/calibrated_matter.py):
`neutral_higgs_norm`, `family_yukawas`, `mass_kernel`, `pairing_kernel`,
`mass_basis_change`, and `bosonic_tree_control`. Its parameter packet
is an empirical common-scale example, not a selection theorem. Its
conventions include `gY=sqrt(3/5)*g1_GUT` and
`lambda_SM=lambda_paper/4`. The same principal symbol and full Nambu
conjugation are defined in
[chiral_candidate/kinetic.py](../common_source/chiral_candidate/kinetic.py).

## Higgs vacuum algebra and local action identities

[ElectroweakBreakingComposition.lean](../../Lean/Screen/ElectroweakBreakingComposition.lean)
defines the shell radius `vev=w`, not the conventional `vEW=sqrt2*w`.
Within namespace `OPH.ElectroweakBreakingComposition`, the witnesses are:

| Exact Lean witness | Typed conclusion |
|---|---|
| `minimizer_iff` | For positive quartic coefficient, the minima are exactly `H^dagger H=w^2`. |
| `stabilizer_exactly_one_line` | For positive gauge couplings and `w>0`, the electroweak stabilizer is one-dimensional. |
| `goldstone_exactly_three_parameters` | The zero second-variation subspace has three real gauge-orbit coordinates. |
| `gaugeMassForm_exact` | The exact infinitesimal-action norm at the supplied vacuum. Canonical real-vector mass squared is twice this quadratic form. |
| `wMass_value`, `zMass_normalized`, `photon_massless`, `tree_mass_ratio` | Finite gauge quadratic-form values, normalized neutral directions and their ratio. |
| `yukawa_mass_at_vacuum` | The invariant line becomes its Yukawa coefficient times the shell radius. |
| `electroweakBreakingComposition_receipt` | Composes these claims on `EWBreakingPremiseData`, with positive scalar/gauge parameters and a supplied Yukawa coefficient. |

The three family matrices and full native invariant tensors are supplied
by the executable witnesses, rather than inferred from the Lean receipt's
single colorless Yukawa line.

[LocalGaugeJetAction.lean](../../Lean/Screen/LocalGaugeJetAction.lean) proves
`OPH.LocalGaugeJetAction.covariant_first_jet`, `curvature_second_jet`,
`quadratic_trace_ward`, and `anticommuting_spin_pair`. These identities retain
the derivative of a local parameter and the noncommuting field order.
They are ring/trace identities, not a construction of a chiral continuum
gauge theory or a selected quantum vacuum.

## Source spectral transfer

[matching.py](matching.py) checks a source/target isometry, the entire
Hamiltonian intertwining residual, and supplied conserved-charge maps.
The Frobenius residual gives a certified operator-error upper bound.
The spectral-distance, probability-tail and resolvent conclusions follow
from the spectral theorem as proved in [PROOF.md](PROOF.md). Exact
intertwining transfers spectral measures; approximate matching alone
does not determine exact level multiplicities, spectral exhaustion or a protected defect
origin. Topological sector maps, full Fock interactions, the physical
reference/clock and continuum control are independent premises.

No theorem in this composition identifies the radial defect catalogue
with these matter multiplets or derives the empirical parameter packet.
The conditional field action and its mass matrices are defined without
asserting those identifications.
