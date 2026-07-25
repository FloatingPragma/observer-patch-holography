# Flavor and de Sitter Two-Sided Claim Audit

## Scope and decision rule

This audit covers the flavor, de Sitter shock, and Phase 0 source memoranda
against the registered paper, code, proof, and tracking surfaces. It applies
the same burden to favorable and adverse statements.

A favorable statement binds only on the domain proved by its receipt. An
adverse statement binds only when all of the following conditions hold:

1. the tested model class is stated exactly;
2. the search is exhaustive inside that class;
3. the comparison input is independent of the solve path;
4. conventions and uncertainty assumptions match the tested observable;
5. every live escape route outside the class is named; and
6. the result is transferred to the OPH core only through a proved dependency.

Failure of any condition gives an open bridge or a diagnostic. It does not
give a theory-wide negative result.

## Audit verdict

The audited material contains several exact successes and no mathematical
contradiction with the recovered OPH core. The strongest adverse result
retrospectively rejects the complete six-assignment family of one declared
register-Clebsch implementation of down-type flavor. It does not reject the
pairing theorem, the target-free unordered weight selection under its declared
rules, the charged-lepton circulant identity, gauge reconstruction, or OPH as
a whole.

The complete theory-of-everything claim has no theory-wide falsifier discharged
by this packet. It also has no complete physical establishment: the canonical
scoreboard contains zero claims in the `physical_establishment` class and
seven undischarged physical identifications. This is a missing-bridge
classification, with no mathematical contradiction implied.

## Exact successes

| Result | Closed content | Evidence and boundary |
| --- | --- | --- |
| Phase 0 projection and public-surface discipline | The registry, matrices, dependency graph, scoreboard, public quantitative annotations, provenance checks, and negative-control scorecard form machine-checked projections. Public numerical rows resolve to a producer and claim class. | [`tools/check_public_surface_claims.py`](../tools/check_public_surface_claims.py), [`claims/claim_registry.yaml`](../claims/claim_registry.yaml), [`tracking/claims_scoreboard.md`](../tracking/claims_scoreboard.md). This establishes reproducibility and typing, with no automatic physical truth claim. |
| Flavor channel pairing | The declared exterior-matter packet has one invariant line on each allowed down-type and charged-lepton channel, with a zero-dimensional forbidden-channel control. | [`clebsch_register_pairing_selection.json`](../code/particles/runs/flavor/clebsch_register_pairing_selection.json). Independent invariant lines are compatible channels; the receipt does not equate their Yukawa coefficients. |
| Clebsch unordered multiset | Exhaustive enumeration of the declared alphabet `{1/3, 1, 3}` under measure balance and register faithfulness selects the unordered multiset `{1/3, 1, 3}` without measured masses or angles in the solve path. | [`clebsch_register_pairing_selection.json`](../code/particles/runs/flavor/clebsch_register_pairing_selection.json). The alphabet and the two selection rules are declared premises. No physical generation order follows. |
| Positive-chamber Koide identity | A positive Hermitian regular-`C3` circulant satisfies `Q = 1/3 + (2/3)(abs(b)/a)^2`, hence `Q = 2/3` exactly when `abs(b)/a = 1/sqrt(2)`. The phase drops out of `Q` and carries the two mass ratios jointly. | [`koide_circulant_identity.json`](../code/particles/runs/leptons/koide_circulant_identity.json), [`KoideCirculant.lean`](../Lean/ObserverPatchHolography/KoideCirculant.lean), [standalone paper](../extra/koide_identity_from_positive_c3_face_circulants.tex). The result is restricted to the nonnegative-eigenvalue chamber. |
| Finite tracial Koide balance | Equal rank-two event blocks in the declared finite tracial Gelfand-Naimark-Segal packet give `abs(b)/a = 1/sqrt(2)` exactly. | [`charged_koide_orientation_isometry.json`](../code/particles/runs/leptons/charged_koide_orientation_isometry.json). The event packet is conditional. A physical chiral-family attachment, phase selection, and the two ratios are open. |
| Icosahedral axis spectrum | The complete 31-axis real three-dimensional menu, including multiplicities, is exact. Its `31.717474...` degree entry reproduces `arctan(1/phi)` as a hard geometric self-test. | [`icosahedral_axis_angle_spectrum.json`](../code/particles/runs/flavor/icosahedral_axis_angle_spectrum.json). This is a strong finite-geometry result independent of the Cabibbo comparison. |
| `W5` stabilizer theorem | Threefold and fivefold fixed loci are one-dimensional and doubly degenerate. The twofold fixed locus has projective dimension two and admits simple spectrum. | [`w5_stabiliser_spectrum_bound.json`](../code/particles/runs/flavor/w5_stabiliser_spectrum_bound.json), [`W5Stabilizer.lean`](../Lean/ObserverPatchHolography/W5Stabilizer.lean). The theorem locates the dynamical selector obligation and permits a screen-derived invariant potential. |
| Pure de Sitter shock normalization | For pure de Sitter space, `mu^2 = (d-2) kappa r_c = d-2 = lambda_1(S^(d-2))`; the radius and cosmological constant cancel. | [`screen_shock_spectrum.json`](../code/geometry/runs/screen_shock_spectrum.json), [`DeSitterCapacityShock.lean`](../Lean/ObserverPatchHolography/DeSitterCapacityShock.lean), [de Sitter note](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex). Observer-mass corrections and non-pure-de-Sitter backgrounds lie outside this identity. |
| Finite entropy and capacity transfer | Maximization over sector probabilities gives `S_gen^max = log M`. Uniform admissible integer depletion gives the exact change `log(1-f)`. The positive-real interpolation is strictly decreasing. | [`screen_shock_spectrum.json`](../code/geometry/runs/screen_shock_spectrum.json), [`DeSitterCapacityShock.lean`](../Lean/ObserverPatchHolography/DeSitterCapacityShock.lean). The finite statement and its analytic interpolation are typed separately. |
| Exact Hessian split | The logarithmic sector observable has nonzero symmetric-point gradient and Hessian `(1/(n d^2))(I-(2/n)J)`. Fixed-`M` tangent curvature is positive and homogeneous curvature is negative. | [`screen_shock_spectrum.json`](../code/geometry/runs/screen_shock_spectrum.json). This exact result rules out an unconstrained interior-maximum argument and supports the separate one-sided transfer statement. |
| Port and edge spectral identity | The icosahedral line graph contains the complete port Laplacian spectrum plus eigenvalue `10` with multiplicity `18`. The low port and edge spectra are identical. | [`screen_shock_spectrum.json`](../code/geometry/runs/screen_shock_spectrum.json). This closes the port-versus-edge choice only inside the nearest-neighbour combinatorial-Laplacian class. |
| Repair-domain separation | The declared repair generator acts on one fixed-sector state space. A sector-dimension shock changes the carrier domain and requires a separately typed extension. | [de Sitter note](../extra/de_sitter_time_advance_sign_from_fixed_screen_capacity.tex). The kernel equals the intersection of the projection ranges; reduction to constants requires an additional irreducibility result. |

## Corrections to source formulations

| Source formulation | Audit result |
| --- | --- |
| The down-type ratio misses FLAG by about `24 sigma`. | The canonical FLAG 2024 inputs give `9.126` and `9.498` standard uncertainties under independent propagation, and `7.850` and `7.480` under maximal positive correlation. The covariance is unavailable, no OPH theory uncertainty is supplied, and the rejection gate is retrospective. The papers therefore report the `15.2%` and `12.8%` central discrepancies without a covariance-aware significance claim. |
| The generation ordering is discharged by a sixfold margin. | The six-assignment scan uses observed mass discrepancies to rank the assignments. The adopted assignment is uniquely least discrepant under that retrospective metric. It is target-informed and does not derive a physical generation order. |
| The distinct-assignment ratio menu contains `1` or `27`. | Distinct assignments from `{1/3, 1, 3}` give `{1/9, 1/3, 3, 9}`. Values `1`, `27`, and `1/27` do not occur in that menu. |
| The register-scale ratio is independent of loop order and threshold matching. | Common down-lane and same-scale QCD factors cancel on the declared one-loop lane. Flavor-dependent charged-lepton self-running makes the unqualified low-scale substitution inexact, and a generation-dependent threshold changes the ratio. The protected identity is restricted to the stated register-scale common-transport class. |
| The Gatto-Sartori-Tonin row is a second failure and nature independently obeys it. | The displayed `sqrt(m_d/m_s)` is a restatement of the same rejected mass ratio under an assumed texture. The lane supplies no up/down matrix pair or relative left-handed eigenbasis, so it gives no independent Cabibbo prediction or test. |
| The Particle Data Group central Koide coordinate is `0.666660511`. | The cited 2026 central masses give `0.6666644634026367`. The declared response coordinate is `0.6666644634090389`. Their proximity is a retrospective diagnostic because the response architecture is historically target-informed. |
| OPH predicts the sign and part of the measured Koide deviation. | The exact theorem gives balance under a conditional finite event packet. The MCPR deviation belongs to a historically target-informed architecture, so its agreement in sign and magnitude has diagnostic status and no prospective evidential weight. |
| Icosahedral symmetry cannot produce the Cabibbo angle. | The exact receipt excludes direct equality with an acute angle between two axes in the canonical 31-axis real three-dimensional menu. Spinorial representations, broken symmetry, source-derived dynamics, and general overlap constructions lie outside the tested class. No universal higher-order correction bound is proved. |
| The golden-ratio solar angle is a live `A5` prediction. | The exact axis menu contains `arctan(1/phi)`. The screen does not select the required residual pair, physical family attachment, or charged basis. This is an exact available angle and a candidate selector target, with no frozen physical prediction. |
| The 280-expression null grammar gives `39.3%`, `20.2%`, and `11.6%`. | The 280-expression grammar is underspecified. The fully declared 602-expression depth-one grammar gives `38.65%`, `20.65%`, and `10.45%`. These rates calibrate expression mining and grant no evidence for or against OPH. |
| Fixed capacity makes the symmetric sector point an interior area maximum. | The symmetric point has nonzero gradient. Its fixed-`M` tangent Hessian is positive. The valid maximum is the one-sided boundary at zero observer transfer under the explicit uniform-depletion law. |
| A black-hole horizon has no fixed capacity budget and is therefore a minimum. | The finite-screen receipt contains no black-hole capacity theorem. The black-hole comparison requires its own capacity dictionary, observer placement, and gravitational attachment. |
| The normalized icosahedral shock spectrum needs only a gauge premise. | The physical spectrum requires both `DS-GAUGE` and `DS-LAPLACIAN`. The first identifies the rotation triplet as exact gauge. The second identifies the physical kinetic operator with the scaled nearest-neighbour graph Laplacian. |
| The repair generator has kernel equal to the constants from projection algebra alone. | Projection algebra gives `ker L_rep = intersection_v Ran(P_v)`. Equality with the constants requires irreducibility or a separate intersection theorem. The fixed-domain separation is valid without that stronger kernel claim. |
| `mu^2 = lambda_1` supplies the physical finite-screen operator. | The equality is exact for the imported pure de Sitter shock equation. It supplies no discrete kinetic operator, gravitational coefficient, screen attachment, or source-bound Einstein branch. |

## Adverse results under the high-burden rule

| Apparent negative | Exact tested class and closure strength | Live escape routes | Verdict and OPH-core reach |
| --- | --- | --- | --- |
| Down-type register-Clebsch mismatch | The class fixes the declared alphabet and rules, exhausts all six assignments, imports the historically stipulated charged-lepton coordinate at the register scale, and uses common down-type multiplicative transport. The four distinct assignment ratios are `{1/9, 1/3, 3, 9}`. Both FLAG rows reject every assignment under the retrospective conservative experimental-only gate. | A different source-derived cross-sector coefficient relation, a different source-derived alphabet, generation-dependent threshold transport, a physical charged-lepton attachment, and a different flavor potential define other classes. A source-derived generation order would repair ancestry, although no order in the tested menu repairs the numerical mismatch. | **Conditional branch rejection.** The complete declared common-transport assignment family is retrospectively rejected. The pairing theorem, unordered multiset theorem, charged-lepton algebra, gauge reconstruction, and recovered core are unaffected. |
| Reciprocal-ray quark candidate | The stored common-scale reciprocal-ray class misses held-out Yukawa coordinates across the tested scale range. The strongest reported held-out miss is `21.56%`. | The general six-scalar quark interface, a source-derived flavor carrier, another orbit selector, and another target-clean transport law. | **Conditional candidate rejection.** The restricted reciprocal-ray class is rejected. No generic quark nonexistence theorem or recovered-core falsifier follows. |
| Direct residual-axis Cabibbo construction | All pairs in the canonical 31-axis real three-dimensional menu are enumerated exactly. The smallest nonzero acute angle is `20.905157...` degrees, while the comparison coordinate is `13.002878...` degrees. | Spinorial or higher representations, higher-order breaking, weighted overlaps, dynamical mixing, and non-axis constructions. | **Exact narrow no-go.** It reaches one direct geometric ansatz and no broader `A5` or OPH mixing claim. |
| `W5` symmetry-only mass-ratio selection | The fixed loci and spectrum degeneracies are exact for `W5 = Sym^2_0(R^3)` under the real `A5` rotation action. | A specific source-derived invariant potential, other dynamics, trivial-stabilizer points, and other representations. | **Exact selector boundary.** Residual stabilizer symmetry alone does not select the ratios. The theorem leaves the physical orbit-selection program viable. |
| Weighted-cycle neutrino candidate | The stored hand-written, target-informed template fails its NuFIT 6.1 correlated profile, and the claimed shared-basis recovery is algebraically an identity for the stored construction. | A source-derived neutral-family kernel, physical charged basis, independent selector, absolute scale, and another neutrino mechanism. | **Conditional candidate rejection.** The stored weighted-cycle candidate has no theorem or prediction status. The recovered core is independent. |
| RSCC correction ledger | On the stored retrospective metric, deleting the `w^2` terms and `delta_g` improves both the maximum residual and raw residual sum. | A source-derived loss, an independently frozen comparison, another observable, or a derived refinement law. | **Conditional model-selection rejection.** The detailed RSCC correction ledger lacks support from its own benchmark. No recovered-core theorem is tested. |
| Twelve-port Standard Model selector uniqueness | The declared twelve-port response realizes the Standard Model Lie type, while the exact coefficient classification also permits two other compact Lie types. Seven alternative deltahedral producers are undeclared and therefore unknown. | A source theorem that selects the twelve-port producer and excludes the viable same-carrier alternatives. | **Missing selector.** Standard Model structure is available on the declared carrier; uniqueness and physical source selection are open. No negative conclusion about the recovered finite witness follows. |
| De Sitter interior-maximum argument | Exact differentiation refutes an interior stationary maximum of the logarithmic sector observable at the symmetric point. | The one-sided uniform capacity-transfer family supplies a distinct exact boundary maximum. Other transfer laws require separate analysis. | **Argument correction.** The invalid Hessian interpretation falls. The finite entropy maximum, boundary transfer sign, and graph identities are exact. |
| Physical de Sitter shock sign | No adverse physical datum is present. The screen-to-gravity map requires capacity-to-area, capacity-to-observer-mass, Einstein-branch, coefficient, and physical-scale identifications. | Discharge those dictionaries on one source-bound branch, or construct another shock carrier and kinetic operator. | **Missing physical bridge.** The finite sign mechanism is conditional. It neither confirms nor falsifies the OPH gravitational continuation. |
| Discrete icosahedral shock spectrum | The graph spectrum and line-graph identity are exact. The normalized shock interpretation depends on two undischarged physical premises. | Failure of exact gauge transport, a weighted or nonlocal kinetic operator, refinement, or another carrier. | **Conditional diagnostic.** The set `{-2, 0, 1+3/sqrt(5), 1+sqrt(5)}` is not a physical prediction without both premises. |
| Fixed-sector repair gap versus a growing shock mode | The repair generator and the shock deformation have different domains. | A separately constructed dimension-changing repair or shock generator. | **Category separation.** The repair gap neither excludes nor produces the shock mode. There is no OPH contradiction. |
| Static-patch trace obstruction | Chen, Stanford, Tang, and Yang test a positive cyclic Hilbert-space trace interpretation of a Euclidean static-patch path integral with an observer worldline. Their negative-shock result obstructs that proposal under its stated assumptions. | Horizon-screen descriptions without that trace, or another operator-algebraic interpretation outside the tested assumptions. | **External conjecture obstruction.** OPH is affected only if it adopts the tested trace representation. The finite-screen calculation supplies no such trace. |
| Zero `physical_establishment` claims | The registry assigns every claim to its declared epistemic class and lists seven undischarged physical identifications. | Source-bound physical carriers, clocks, currents, horizon maps, and Einstein attachments. | **Missing establishment.** The count blocks an empirical claim of complete TOE establishment. It is not evidence of mathematical inconsistency. |

## TOE reach classification

| Category | Finding |
| --- | --- |
| Mathematical contradiction with the recovered core | None in the audited packet. |
| Exact favorable theorem or receipt | Pairing and unordered Clebsch selection, the positive-chamber Koide identity, finite tracial balance, exact axis geometry, the `W5` fixed-locus theorem, pure de Sitter normalization, finite entropy and transfer laws, the Hessian split, and the regular line-graph identity. |
| Conditional branch rejection | The complete declared register-Clebsch common-transport assignment family, the restricted reciprocal-ray quark candidate, the stored weighted-cycle neutrino candidate, and the detailed RSCC correction ledger on its stored benchmark. |
| Missing physical bridge | Physical family attachment, Yukawa-coefficient equality, generation order, the `W5` potential, capacity-to-area, capacity-to-observer-mass, source-bound Einstein realization, shock coefficient and scale, `DS-GAUGE`, and `DS-LAPLACIAN`. |
| External conjecture obstruction | The positive cyclic static-patch trace proposal tested by Chen, Stanford, Tang, and Yang. |
| Theory-wide OPH or TOE falsification | None supplied by these materials. |
| Complete physical TOE establishment | Absent under the canonical registry classification. |

The favorable results retain their full exact content without relying on the
rejected flavor branch or the open gravitational attachments. The adverse
results retain their force inside their closed classes without transfer to
independent theorem rows. This separation gives the strongest defensible case
for the successes and the strongest defensible boundary on every negative.

## Tracking precision checks

The following state-tracking details require explicit treatment:

- The de Sitter capacity-transfer registry statement should use
  `c_obs` or another capacity symbol instead of `m`; this prevents an
  unproved identification with observer mass.
- The finite integer depletion law and the positive-real analytic
  interpolation should occupy separate registry clauses, matching the
  executable receipt.
- The physical shock-sign assumption list should name the
  capacity-ledger-to-observer-mass and shock-coefficient dictionaries
  explicitly.
- `DS-GAUGE` and `DS-LAPLACIAN` should each have a tracked discharge route.
  A normalized physical spectrum has no promotion path while both premises
  lack live gates.
- Comparative TOE scorecards require a declared competitor set, scoring
  rubric, and evidence links. Without those objects, a claim that one program
  leads every competitor is outside the machine-checked claim surface.
