# Horizontal particle-family dependencies

[`families.json`](families.json) gives nine left-to-right lanes covering the
eighteen elementary target labels in the
[authenticated vertical atlas](../particle_candidate_atlas/audit.json).
Every numerical mass, uncertainty, confidence level and scheme is a lookup
in its frozen [reference table](../particle_candidate_atlas/reference_values.json).
The family data do not replace those values or compute physical source masses.

The four stages in each lane are the exact internal source structure,
the attached field representation, the propagation or mass mechanism,
and the appropriate observable. The labels distinguish finite theorems,
supplied field completions and unproved physical matching. The lanes form
a complete dependency map, not a claim that all implications have been
discharged. No protected source defect is assigned to an individual species.

## The gauge theorem is more than a dimension count

The [finite-current fixed-space theorem](../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf#page=7)
starts with a faithful complete twelve-port response whose image is closed
under commutators in a unitary Lie algebra. The proper carrier rotations
must act by inner automorphisms of that same response. Their transitive
alternating-group action has a one-dimensional fixed space. Compact
reductivity, the centre bound and the exclusion of the centreless
`su(2)^4` branch then force `u(1)+su(2)+su(3)` through compact-simple
classification. The conclusion uses brackets and the inner-action
premise; the arithmetic `11=3+8` by itself supplies no Lie algebra.

The [Lean bridge](../../Lean/Screen/A2HolonomyBridge.lean) proves the
faithful-response dimension and finite fixed-space selection while
retaining the centre/classification adapters as explicit hypotheses.
It does not manufacture a source response or holonomy-fullness witness.
The JSON records preserve that division of mathematical work.

The [conditional matter theorem](../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf#page=14)
also requires the declared matrix current, spin lift and exterior menu.
The anomaly scan selects a charge-conjugate projector pair inside that
menu. Its fifteen weights have the Standard Model quantum numbers and
admit the three charged Yukawa channels with a supplied Higgs doublet.
The diagonal `Z6` kernel determines the maximal faithful local matter
image. The [global-form boundary](../../paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf#page=15)
keeps four compatible covers/quotients distinct from the physical group.

There is a positive finite interaction result as well: the
[full-Fock Gauss construction](../standard_model_source_map/GAUGE_ATTACHMENT.md)
and [quantum Higgs join](../particle_candidate_atlas/INTERACTING_CHIRAL_JOIN.md)
use the same currents in an actual finite Hamiltonian, with every local
Gauss commutator zero. Their coefficients, physical state and continuum
interpretation are separate data. The graphic must place that finite
construction between the representation dictionary and the conditional
field-theory completion, rather than suggesting the representation alone
determines an interaction strength.

The [family-band selector](../../Lean/Screen/A5FamilyBand.lean) supplies a
useful conditional shortcut. Its exact cost order is
`5-sqrt(5)<6<5+sqrt(5)`. Among the faithful bands in the stated dimension
window, the rank-three band minimizes that cost. Identifying physical
family multiplicity with that band and making the cost operational are
additional premises. A spatial rank-three readout does not automatically
count three particle generations.

## Reading the family lanes

| Lane | Conditional field content | Mass or observable boundary |
|---|---|---|
| Charged leptons | Three charged Dirac branches | Electron charged infrared threshold; muon and tau resonances and widths |
| Up quarks | Three up-type color triplets | Renormalized mass scheme; top reconstruction convention; no free isolated colored pole |
| Down quarks | Three down-type color triplets | Running mass and QCD matching; a common unit cannot fix the rejected light-quark ratios |
| Neutrinos | Neutral components plus a declared Dirac or Majorana extension | Flavor is a superposition of mass fields; absolute masses are symbolic |
| Photon | Unbroken transverse electromagnetic vector | Zero action mass; physical detector POVM and field attachment required |
| Gluons | Eight adjoint gauge directions | Zero action masses; interacting color and inclusive/hadronic observables required |
| W signs | Two charged massive vector directions | Higgs-phase quadratic masses, matched resonance mass and width |
| Z | One massive neutral vector direction | Neutral mixing and line-shape convention |
| Higgs | One radial scalar direction | Normalized potential curvature, then scalar resonance matching |

Physical Lorentz spin is distinct from internal weak isospin. The configured
CAR algebra supplies finite fermionic statistics, while its spacetime
spin-statistics interpretation uses the local Lorentzian completion.
Conjugate charges alone do not establish equal measured particle and
antiparticle masses; the conditional field theory includes the stated
CPT and vacuum hypotheses.

A massive branch costs energy above its chosen vacuum even at zero
momentum. The free relativistic formula is `E^2=c^2 p^2+m^2 c^4` once its
physical generator, clock and unit are matched. Repair effort is not
automatically that energy. A static source texture also has background and
binding energies, and its localized frequency is not automatically a
relativistic rest mass. The [compact source proof](../particle_candidate_atlas/DEFECT_FIELD_JOIN.md#which-identification-is-excluded)
states the exact translation and Pauli-matrix obstruction.

## Calibrations and their removal

Ten entries in `calibration_paths` identify the supplied data, their positive
use, a concrete derivation route and its success condition. They cover the
one-ratio phase fit, three-input attenuator, capacity law, common clock,
family/chiral structure, gauge/Higgs parameters, quark flavor, neutrino
extension, bare-to-effective matching and defect-to-species dictionary.

The [balanced one-ratio construction](../charged_lepton_response/CALIBRATED_SANITY.md)
uses muon/electron to select a phase and holds tau out of that fit. The
[three-input attenuator](../particle_candidate_atlas/CALIBRATED_TRIPLET.md)
uses all three charged central masses and therefore has no held-out mass
prediction. The [capacity law](../charged_lepton_response/P_DEPENDENCE.md)
has distinct target-informed architecture and endpoint ancestry. These
calculations must remain separately labeled in the rendering.

The physical route is to compute complete same-state correlations and
their temporal/spatial inverse kernels, kinetic forms and self-energy,
then compare the appropriate charged threshold, resonance or running
parameter with bounded errors. Selecting an arbitrary allowed correction
because it improves agreement does not perform that derivation. The
[conditional Standard Model construction](../particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md#concrete-source-to-benchmark-obligations)
states these common dependencies.

## Link and coverage checks

The content checks require nine lanes with four stages each, exactly
eighteen unique target IDs matching the frozen reference table, explicit
source-identification flags set to false, and resolvable repository paths
for every stage, gauge theorem and calibration entry. Paper links use
physical PDF page indices checked against the release PDF text:
particle §§2.2–2.5 on pages4–7, §4.5 on page14, §6.1 on page60,
§6.3 on page61, §§6.4–6.6 on page62 and §8.2 on page66;
gauge theorem3.1 on page7, §6 on pages13–15 and §7 on page16.
These are explanatory links to retained results, not additional source
premise certificates.
