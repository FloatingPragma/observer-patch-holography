# Corrected Einstein Branch: Theorem and Dependency Index

This index covers GitHub issue #578. The Lean modules mechanize exact finite
algebra and the composition *conditional on* the listed continuum,
asymptotic, and physical inputs. They do not construct an inhabited
Einstein-admissible source tower.

Toolchain: `leanprover/lean4:v4.29.1`, Mathlib `v4.29.1`. The CI build imports
`EinsteinBranch/AxiomAudit.lean`, which applies `assert_no_sorry` to every
public theorem and prints a per-theorem axiom report. The reports contain at
most Lean's standard extensionality/choice/quotient axioms (`propext`,
`Classical.choice`, `Quot.sound`).

## Theorem-to-paper map

| Lean declaration | Compact-paper anchor | Level | What is proved / what remains input |
|---|---|---|---|
| `BareConsensusTower` | `def:compact-bare-finite-consensus-reduct` | exact finite structure | Finite states/quotients, mismatch, repair, normal forms, boundary data, and typed coarse maps. |
| `bare_consensus_not_einstein_complete` | `thm:compact-bare-consensus-not-einstein-complete` | model-theoretic exact | Two geometric extensions share one nondegenerate reduct and disagree on Einstein truth. |
| `TypedArrowCommutation`, `ReadoutNaturality` | `def:einstein-admissible-realized-tower`, EB2 | explicit receipt predicates | Common-domain typed arrows and refinement diagrams are failable premises, not constructor truths. |
| `EinsteinAdmissibleTower` | `def:einstein-admissible-realized-tower` | bundled conditional structure | Packages one bare source tower, its typed readouts, boundary fibre, arrow commutation, and refinement naturality; no inhabitant is constructed. |
| `boundary_fiber_readout_composition` | EB1; D1 composition | exact finite composition | Same protected boundary plus consistent fibre uniqueness gives identical geometry/modular/event/stress/entropy/scale readouts. |
| `tomographyDirections_null` | `thm:null-tomography` | exact algebra | Each of the displayed nine directions is Minkowski-null. |
| `tomographyDecoder_design`, `tomographyDecoder_charge`, `tomographyCharge_injective` | `thm:null-tomography` | exact algebra | An explicit left inverse for the nine-direction design map, hence exact rank-nine reconstruction of the trace-free symmetric tensor. |
| `null_cone_determines`, `jacobson_step` | `thm:null-tomography`; `lem:nullreconstruct` | exact algebra | Null data leave exactly a metric multiple. |
| `directSumEntropy_eq_bulk_add_edge` | `thm:bulk-edge-central-first-law` preamble | exact finite | Corrected convention: Shannon entropy belongs to `S_bulk`; `S_edge = sum p log d`. |
| `centralVariation_eq_edgeEntropyVariation` | `thm:bulk-edge-central-first-law` (2) | exact finite | `z_alpha = log d_alpha` makes the central variation exactly the edge variation. |
| `central_edge_normalization_defect` | `prop:entropy-coefficient-countermodel` | exact finite negative control | Computes the defect for a wrong central normalization. |
| `finite_bulk_edge_central_first_law` | `thm:bulk-edge-central-first-law` | exact finite composition | Entropy differential = modular pairing, `K=2pi B+Z`, and `Z=L_edge` imply the boxed first law. These three equalities remain explicit premises. |
| `maxEnt_envelope_identity` | `thm:maxent-lagrange-stationarity` | exact finite | Normalization tangent + unit constraint tangent + Gibbs differential imply `dS/dt=lambda`. MaxEnt-family existence is not asserted. |
| `diamondKernelIntegral_eq` | `lem:smallball` | exact coefficient conditional on kernel | Evaluates the radial polynomial integral as `4pi ell^4/15`; identification with a continuum modular charge remains a premise. |
| `bulkSmallBallCoefficient` | `lem:smallball` | exact arithmetic | Multiplies the kernel by `2pi` to obtain `8pi^2 ell^4/15`. |
| `fixedVolumeAreaCoefficient` | `thm:fixed-volume-small-ball-area` | exact arithmetic | Checks the coefficient algebra; the smooth fixed-volume geometric expansion is not derived. |
| `restFrameEinsteinRelation` | `thm:einstein` | exact conditional algebra | Stationarity + continuum kernel identity + fixed-volume area identity imply the rest-frame coefficient `8pi G`. |
| `unit_timelike_determines`, `tensor_upgrade` | corrected first-variation content of `cor:einstein` | exact algebra | All future unit-timelike quadratic data determine every symmetric-tensor component. Coverage is a premise. |
| `ddiv_lam_eta`, `step_invariant_of_divergence_free` | `lem:lambda-constancy` | exact discrete composition | Ward+Bianchi divergence cancellation forces componentwise invariance of the metric residue. |
| `einstein_equation_with_constant_symm` | `lem:lambda-constancy`; `thm:absolute-einstein-base-condition` | exact algebra conditional on conservation | Null metric ambiguity + Ward + Bianchi + symmetric connectivity gives one constant per component. Continuum conservation laws are premises. |
| `strictManifest_dependency` | strict dependency audit | exact manifest logic | Every required finite/analytic/physical key has mathematical evidence. |
| `strictManifest_erase_required_fails` | artifact deletion logic; `prop:no-hidden-geometry` boundary | exact syntactic negative control | Erasing any required key fails strict validation. This does not claim semantic minimality. |
| `continuumEinstein_from_explicit_premises` | `thm:absolute-einstein-base-condition` | conditional physical composition | Ward/Bianchi constancy is fixed by VR and the coupling is translated by the independent scale premise. |
| `composedEinsteinBranch` | `thm:einstein-branch-entry-composed` | composed conditional theorem | Returns the boundary/readout, first-law, MaxEnt, small-ball, common-tail, universal-source, dependency, and Einstein certificates from explicit premises. It proves no premise truth and no tower nonemptiness. |
| `composedEinsteinAdmissibleTower` | `thm:einstein-branch-entry-composed` | bundled composed theorem | The same conditional implication with the finite/common-domain receipts supplied by an explicit `EinsteinAdmissibleTower` value. |

## Premises deliberately not proved

These appear either as fields of `ContinuumEinsteinPremises` and
`SmallBallPremises`, as `ScalingTailReceipt`, or as keys in the strict
manifest. They are not Lean axioms: they are ordinary theorem arguments.

| Premise family | Lean surface | Status boundary |
|---|---|---|
| Source-derived common-domain readouts and factorization | `TypedArrowCommutation`, `ReadoutNaturality`, manifest keys | One source tower must supply every readout; shared labels/digests alone are insufficient. |
| Continuum diamond modular kernel and local stress identification | `SmallBallPremises.continuumKernel` | Analytic/physical input. The polynomial integral only checks its coefficient. |
| Smooth fixed-volume small-geodesic-ball area expansion | `SmallBallPremises.fixedVolumeArea` | Imported continuum geometric identity, including its remainder control. |
| One common certified `o(ell^4)` family | `ScalingTailReceipt` | Full-tail statement; finite regression is insufficient. |
| Lorentzian regularity and all-timelike variation coverage | manifest keys and the hypothesis of `unit_timelike_determines` | Coverage and smoothness are not inferred from finite consensus. |
| Weak Ward identity | `ContinuumEinsteinPremises.ward` | Semiclassical conservation input. |
| Contracted Bianchi identity and metric compatibility | `ContinuumEinsteinPremises.bianchi`, manifest key | Smooth-geometric input. |
| Connected component | `ContinuumEinsteinPremises.connected` | Constants may differ between disconnected components. |
| Universal coupling (`UC`) | `ContinuumEinsteinPremises.universalCoupling` | The entropy and null channels must use the same stress source. |
| Vacuum reference (`VR`) | `ContinuumEinsteinPremises.vacuumReference` | Fixes the otherwise undetermined metric residue. |
| Independent physical scale | `ContinuumEinsteinPremises.physicalScale`, manifest stabilizer key | Identifies `kappa=8pi G`; SI conversion is separate. |
| Inhabited Einstein-admissible tower | no constructor or theorem | Open source-construction/certification problem. |

## Receipt non-vacuity and deletion boundary

The modules include explicit passing and failing finite receipts for repair,
boundary fibres, first-law data, MaxEnt tangents, scaling-family behavior,
normalization defects, connectivity, and scale stabilizers. The manifest maps
each key to an equality, inequality, limit, rank, conservation equation, or
stabilizer condition; no key maps to constant `True`.

`strictManifest_erase_required_fails` proves syntactic deletion failure. It
does not establish that every receipt is semantically independent of every
other receipt; that stronger result requires the isolated full-tower
countermodel matrix described in the compact paper.
