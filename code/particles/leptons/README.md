# Charged-Lepton Lane

This directory is the active `/particles` charged-lepton completion lane.

The live chain is:

1. read back the current ordered charged package on the local family
2. certify that the current support is exhausted
3. expose the minimal beyond-support extension scalar
4. turn the ordered coefficients into the charged excitation-gap map
5. build the log-spectrum readout
6. attach the scale/norm lane and emit the forward charged candidate

## Frozen Stage-5 Postdiction Audit

The historical Stage-5 formula has been preserved as a reproducible
postdiction, separate from the theorem-grade completion lane:

- `derive_charged_stage5_frozen_candidate.py` freezes the legacy
  `P = 1.63094` branch and emits
  `(m_e,m_mu,m_tau) = (0.510882243295, 105.635282871,
  1776.579124017) MeV`, with maximum comparison error `0.02284%`.
- `audit_charged_stage5_frozen_candidate.py` records that the same formula on
  the truncated D10 comparison probe has maximum error `0.08347%`, and that no
  charged RG/threshold map to the compared pole-mass scheme is supplied.
- `derive_charged_z3_maxent_balance_audit.py` proves the conditional block
  identity and records the counterresult: natural trace MaxEnt weights the
  real singlet and charged doublet as `1:2`. Koide balance needs a separately
  derived one-bit compensation and an ensemble-to-Yukawa bridge.
- `derive_charged_z3_phase_holonomy_no_go.py` certifies that the arithmetic
  identity `2/9` is not a family-holonomy theorem. Hypercharge acts as a common
  generation phase and does not generate the cyclic family shift.

The frozen formula does not read charged reference masses at runtime, but its
balanced carrier, phase, exponent prescription, `2^(1/6)` normalization,
determinant attachment, family assignment, and mass-scheme bridge are not
derived from current OPH axioms. It is therefore an accurate historical
postdiction, not a promoted prediction.

## Scoped Non-Identifiability And Exact Port Frontier

- `derive_charged_current_oph_moduli_independence.py` binds the relevant
  structural receipts and proves a carefully scoped result. D9 charged-channel
  admissibility, fixed D10 `P,v(P)`, and invariant 12/24 data permit the exact
  family
  `diag(exp(mu+x), exp(mu+y), exp(mu-x-y))`; they do not select its two centered
  coordinates or common affine coordinate. Mass-dependent electromagnetic
  endpoint transport is explicitly outside this scope.
- `derive_charged_a5_w5_port_frontier.py` proves
  `R^12 = W1 + W3 + W3' + W5`, constructs
  `P5=(I+A_antipode)/2-11^T/12`, and verifies
  `Q^*Q=(8/5)P5`. Only a source-emitted `W5` port record can feed the traceless
  quadrupole. The current uniform 12-port load and count-only 24-slot
  orientation data have zero `W5` projection.
- Three geometric lines require the strict discriminant gate
  `(tr Q^2)^3/2-3(tr Q^3)^2>0`. Physical charged lines additionally require a
  source-defined generation-space `A5` action and a normalized, signed
  intertwiner. Absolute masses still require `log|det Y_e|` and a declared mass
  scheme.

These are negative/frontier theorems for the current source signature, not a
claim that every future OPH extension must fail.

## Icosahedral Face-Incidence Completion Conjecture

The reviewed Pro archive is preserved unchanged at workspace
`correspondence/oph_charged_icosahedral_completion_v1.zip` (SHA-256
`4dff09dbc66368a99dbddac4e238e3dcd42906dada7d2f084ff059dda32b04ef`).

- `derive_charged_icosahedral_face_carrier_frontier.py` proves the valid new
  source-side geometry: the twenty outward-oriented faces form `A5/C3`, and
  each local three-corner fiber carries the regular `C3` action. The sixty
  face-corner flags form the regular `A5` torsor, so a physical global
  generation attachment is still additional data.
- `audit_charged_icosahedral_completion_v1.py` is compare-only. It verifies the
  archive and independently rebuilds the result, inverts the charged triple
  into `(kappa,chi,zeta)`, and evaluates coherent public and source-audit
  branches.
- Against the frozen rounded comparison fields, the submitted hybrid branch
  gives a maximum residual `0.000300 ppm`; the tau residual against the
  packet's higher-precision API central value is `-1.387289 ppm`. The
  coherent canonical public branch gives `0.424 ppm`, and the source-audit
  branch gives `84.10 ppm`.

The candidate is integrated as
`FROZEN_RETROSPECTIVE_TARGET_INFORMED_FACE_CANDIDATE_NOT_COMPLETION`.
Its three equalizers, MaxEnt compensation, phase law, determinant character,
physical attachment, coherent source tuple, and pole-scheme map remain proof
obligations. It does not replace the public theorem boundary.

The follow-up archive
`correspondence/oph_charged_face_incidence_theorem.zip` (SHA-256
`cc550b311c5f2bc416d78f65b4ec9e1bf32f10f9b73fc4135fd05de626b48d4f`)
is integrated as a conditional proof supplement:

- `derive_charged_face_incidence_conditional_theorem.py` proves that the
  explicitly declared affine repair map is a strict contraction with one
  globally attracting fixed point. It consumes no charged references at
  runtime but records the map's historically target-informed ancestry.
- The same receipt exposes the bare-to-endpoint amplitude staging gap and
  three source-multiplier witnesses that retain the displayed Jacobian and
  contraction while changing the masses.
- `audit_charged_face_incidence_theorem_bundle.py` is compare-only and
  verifies archive safety, internal hashes, byte reproduction, and exact
  numerical agreement with the canonical conditional receipt.

This closes declared-map fixed-point uniqueness and iteration stability. It
does not derive or rigidify the map from OPH, prove refinement naturality,
repair the hybrid branch tuple, or close the physical attachment and pole
scheme.

The subsequent archive
`correspondence/oph_charged_source_law_rigidity_v1.zip` (SHA-256
`792dad5c1524984c52d7582e751e1123099978164def80fad34310a34e314dff`)
is integrated as a second conditional supplement:

- `derive_charged_source_law_rigidity_conditional.py` proves the exact
  implication `CFQ1-CFQ8 => one scalar affine response law` inside a declared
  finite packet. Rank-one normalized traces reproduce the eight reciprocal
  path weights and the prior affine kernel.
- The same receipt reconstructs all eight finite graph models without an
  external graph library, but labels them mathematical models rather than
  physical response evidence. Every physical CFQ gate and evidence hash is
  unset.
- `audit_charged_source_law_rigidity_v1.py` is compare-only. It never extracts
  or executes archive code; it verifies safe members, CRC, internal hashes,
  the frozen receipt, graph/path/register agreement, and exact agreement with
  the prior fixed point and conditional mass readout.

This is packet-weight rigidity, not an OPH derivation of the packet. CFQ2-CFQ7
stipulate the register dimensions, tracial rank-one weights, clock, eight
paths, coupling monomials, signs, multiplicity one, and path exhaustion that
carry the coefficients. Existing OPH constrained MaxEnt does not automatically
give the uniform trace. Moreover, rank-one projections in a full matrix
algebra are noncentral, while completed OPH records are normally central, so a
process/readout record-algebra bridge is still required. The one-bit balance,
`2/9` phase, and `6^-14` determinant baselines remain separate conditional
attachments. The old multiplier witnesses therefore remain valid against the
weaker current OPH premises, and public mass promotion remains forbidden.

## How To Read The Active Charged-Lepton Files

The live charged scripts use the same compact derivation header:

- `Chain role`: where the file sits between the shared flavor carrier and the
  forward charged candidate
- `Mathematics`: which ordered-gap, support-rank, or spectral step it performs
- `OPH-derived inputs`: which values are inherited directly from the active
  `/particles` flavor/lepton artifacts
- `Output`: the emitted artifact and, when applicable, the exact missing scalar

The active charged completion tail is:

- `derive_charged_sector_local_current_support_obstruction_certificate.py`
- `derive_charged_sector_local_minimal_source_support_extension_emitter.py`
- `derive_lepton_excitation_gap_map.py`
- `derive_lepton_log_spectrum_readout.py`
- `build_forward_charged_leptons.py`

Current scripts:

- [`derive_charged_sector_local_current_support_obstruction_certificate.py`](derive_charged_sector_local_current_support_obstruction_certificate.py)
- [`derive_charged_sector_local_minimal_source_support_extension_emitter.py`](derive_charged_sector_local_minimal_source_support_extension_emitter.py)
- [`derive_lepton_excitation_gap_map.py`](derive_lepton_excitation_gap_map.py)
- [`derive_lepton_log_spectrum_readout.py`](derive_lepton_log_spectrum_readout.py)
- [`build_forward_charged_leptons.py`](build_forward_charged_leptons.py)
- [`test_no_koide_import.py`](test_no_koide_import.py)
- [`test_no_experiment_label_matching.py`](test_no_experiment_label_matching.py)
- [`test_channel_norm_not_fit.py`](test_channel_norm_not_fit.py)
- [`test_ratio_only_not_promoted.py`](test_ratio_only_not_promoted.py)

These scripts do **not** claim charged leptons are theorem-level. They
give the `e30` closure lane a concrete local home under `/particles`, and they
keep the ordered-shape/hierarchy problem separate from absolute-scale closure.

Current same-carrier constructive blocker:

- `eta_source_support_extension_log_per_side`

Theorem-facing absolute route:

- `oph_generation_bundle_branch_generator_splitting`
- `refinement_stable_uncentered_trace_lift`
- then the determinant-line section and `A_ch` are induced data

Operator-side exact obstruction certificate:

- `oph_generation_bundle_branch_generator_splitting_obstruction_certificate`
  This fixes the current operator-side claim boundary sharply: the displayed
  proxy ordered spectrum together with the displayed level-0 and level-1
  projector system does not force the descended-commutator clause. The
  obstruction witness is the nonzero natural centered commutator with operator
  norm `0.04861550547372144` against projector defect
  `0.06363734112184061`, so the displayed current data do not imply exact
  vanishing or uniform quadratic smallness after the central split.

The post-promotion lift slot is packaged in carrier form:

- `oph_charged_uncentered_trace_lift_cocycle_reduction`
  This keeps the single-slot frontier unchanged while making its exact content
  explicit: after centered promotion, the ambiguity is one scalar
  affine cocycle primitive `mu`, with `C_tilde_e = C_hat_e + mu I`,
  `s_det = 3 mu`, and `A_ch = mu`.

- `oph_charged_mu_physical_descent_reduction`
  This expresses the same carrier one step further under the same refinement-
  stability contract required by the lift: on theorem-grade physical
  `Y_e`, the refinement cocycle vanishes and the exact residual object is one
  physical affine scalar `mu_phys(Y_e)`.

- `oph_charged_centered_operator_mu_phys_no_go`
  This closes the false post-promotion shortcut: even a theorem-grade centered
  `C_hat_e` cannot emit `mu_phys(Y_e)` by itself, because centered operator
  data is invariant under the common-shift action that `mu_phys` breaks.

- `oph_charged_p_to_affine_anchor_reduction`
  This fixes the smallest supported bridge target for the `P`-driven absolute
  lane. A theorem-grade landing from the D10 descendants of `P` on theorem-
  grade physical `Y_e(P)` or on the charged determinant line forces
  `A_ch(P) = (1/3) * log|det(Y_e(P))|`, so the bridge target is smaller than a
  free-standing affine-scalar theorem.

  The same descent package also gives a sharper P-threaded reformulation:
  once `P -> Y_e(P)` lands on theorem-grade physical charged data and the
  post-promotion charged descent hypotheses hold, the physical affine scalar
  `mu_phys(Y_e(P))` is forced by the identity-mode equalizer and does not form
  a separate bridge-sized theorem.

- `oph_charged_determinant_character_frontier`
  This records the proposed source-facing charged frontier. The present
  `q_e = sqrt(g_e d_e)` values are template-dependent refinement-arrow defects,
  not an independently descended nonzero same-object multiplicative character;
  the identity arrow gives `q=0`, so its logarithm is not finite. Therefore the
  candidate expression `S_M = sum_e M_e^ch log q_e` does not yet define a
  physical determinant character on the realized support. The formal no-go
  isolates one additive determinant normalization
  attaching that source-side character to the physical charged determinant
  line. Equivalently, the exact smaller theorem is the sector-
  isolated trace-lift attachment
  `3 mu(r) = sum_e M_e^ch log q_e(r)` on the charged determinant channel, not
  a fresh mass ansatz and not an end-to-end charged closure.

- `oph_charged_determinant_trace_normalization_no_go`
  This packages the formal underdetermination theorem behind the new charged
  frontier wording. The current corpus fixes the populated source-side
  determinant character and the post-promotion identity-mode equalizer
  separately, and it leaves one additive determinant normalization defect
  `N_det(P) = s_det(P) - sum_e M_e^ch log q_e(P)` unresolved.

- `oph_charged_trace_lift_theorem`
  This is the fail-closed issue-546 gate. It audits the proposed central-sector
  factorization, the uncentered lift constant, and the interval for `N_det`.
  On the current corpus it emits `no_go_confirmed_new_source_needed`: `M_ch`
  is only formal, the numerical `q_psi` arrow readback is source-open and not
  stage-indexed, and the D10 matter-determinant landing is conditional. A D9
  representation-role projector does isolate a supplied charged Yukawa channel
  with exact zero quark leakage, but it neither attaches D10 nor normalizes the
  determinant line. The exact countermodel `Y_e -> exp(kappa) Y_e` preserves all
  declared antecedents and ratios while shifting `log|det Y_e|` by `3 kappa`;
  the proposed reference-stage `q` also fails to provide a finite trivialization.
  The ledger may flip only for a source-only certificate with zero leakage and
  a singleton residual interval `[0, 0]`; an interval merely containing zero
  remains non-promoting.

- `oph_charged_physical_class_affine_scalar_reduction`
  This is the strongest supported selected-surface charged theorem packaged on the
  current local surface. It does not claim a quark-style public selected-class
  charged source descent. It states that on theorem-grade physical charged
  `Y_e(P)`, once theorem-grade `C_hat_e` exists, the determinant line fixes the
  canonical physical affine scalar `mu_phys(Y_e(P)) = (1/3) * log|det(Y_e(P))|`
  and therefore the canonical uncentered lift, determinant-line, affine-anchor,
  and charged mass formulas.

- `oph_charged_physical_determinant_line_canonical_uncentered_lift`
  This is the exact conditional theorem behind `#151`: on theorem-grade
  physical charged `Y_e`, theorem-grade `C_hat_e` and the determinant line
  canonically fix
  `C_tilde_e(Y_e) = C_hat_e(Y_e) + mu_phys(Y_e) I`
  with
  `mu_phys(Y_e) = (1/3) * log|det(Y_e)|`.
  The identity-mode equalizer on a fixed physical fiber is then tautological.

- `oph_charged_mass_readout_from_affine_anchor`
  This is the exact downstream theorem beneath the source-landing frontier:
  once theorem-grade `A_ch(P)` exists on theorem-grade physical charged data,
  the charged absolute scale and charged mass triple are pure algebraic
  readouts,
  `g_e(P) = exp(A_ch(P))` and
  `m_i(P) = exp(A_ch(P) + ell_i_centered(P))`.
  The charged source theorem is the determinant trace-lift attachment / zero-
  normalization identity beneath the broader landing from the D10 calibration
  descendants of `P` to theorem-grade physical charged data or the charged
  determinant line.

- `oph_charged_determinant_trace_lift_attachment_required`
  This is the explicit missing-theorem gate for public charged-lepton masses.
  The downstream map `A_ch(P) -> m_e,m_mu,m_tau` is algebraic, but the live
  corpus does not prove `P -> A_ch(P)`. Promotion requires
  `3*A_ch(P)=sum_psi M_ch[psi]*log(q_psi(P))` with source-emitted `M_ch`,
  source-emitted same-label `q_psi(P)`, and a no-target-leak parent DAG.

Layered frontier ledger on disk:

- `oph_charged_absolute_frontier_factorization`
  This separates the current-surface missing affine object `A_ch` from the
  strictly smaller post-promotion single slot `refinement_stable_uncentered_trace_lift`.

Smaller same-carrier primitive on disk:

- `oph_charged_sector_local_support_extension_source_scalar_pair_readback`
  This collects the ordered `eta` then `sigma` invariant readbacks beneath the
  full support-extension shell without promoting either scalar value.

That is the first charged scalar that leaves the exhausted current
support and changes the hierarchy. Without it, the forward charged
artifact remains a structured gap surface rather than a promotable prediction.

The theorem-facing absolute tail lies above that same-carrier frontier.
Once the latent centered charged operator is promoted, the supported post-promotion
single slot is the refinement-stable uncentered trace lift. The determinant-line
section and affine absolute anchor `A_ch` then follow as induced data.

Additional guards:

- [`test_shared_budget_not_silently_localized.py`](test_shared_budget_not_silently_localized.py)
- [`test_channel_norm_refinement_limit.py`](test_channel_norm_refinement_limit.py)

## License And Patent Policy

This particle lepton surface is part of the OPH public repository. See the main
[LICENSE](../../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../../PATENTS.md).
