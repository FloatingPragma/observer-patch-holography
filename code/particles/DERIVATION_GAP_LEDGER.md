# Particle Derivation Gap Ledger

Generated: `2026-06-12T05:21:42Z`

Systematic claim-safe queue after the five-equation P-trunk simplification.

## P-Trunk Status

- Artifact: `P_derivation/runtime/p_closure_trunk_current.json`
- Exists: `True`
- Claim status: `compressed_candidate_trunk_not_final_particle_root`
- May feed promoted particle predictions: `False`
- Candidate P: `1.63097210492078846050203640439`
- Candidate alpha^-1: `136.994020662724205139718642793`
- Source report mode: `thomson_structured_running_asymptotic`

## Electroweak Hierarchy Certificate

- Artifact: `particles/hierarchy`
- Exists: `True`
- Claim status: `closed_local_electroweak_hierarchy_certificate`
- May feed local hierarchy claim: `True`
- Public endpoint P: `1.630968209403959324879279847782648941`
- Public endpoint alpha_U: `0.041124336195630495`
- Public endpoint v/E_star: `2.0199803239725553e-17`
- Public endpoint log10(E_star/v): `16.69465286086613`
- Source-audit P: `1.63097209569432901817967892561191884270169`
- Source-audit alpha_U: `0.04112424744557487`
- Source-audit v/E_star: `2.0198114150099223e-17`
- R_U interval: `[0.041123336195630494, 0.041125336195630496]`
- Krawczyk image: `[0.041124335718554498251685561120150817089903543894567180291985915541514084341979470, 0.041124336672706466103571496781875444719601176485634484342812958790445912366833206]`
- Derivative interval: `[-10.995767998506369371985481516028357811043316285244607, -10.985284409796418191205377821955038437448873048903355]`
- Krawczyk interior inclusion: `True`
- Forbidden DAG paths into protected targets: `[]`
- Boundary: This certificate closes the local P -> alpha_U -> v/E_star hierarchy lane. Excluded gates are the public Thomson endpoint transport, the Higgs RG naturality defect bound, theorem-grade W/Z promotion, and the full no-G clock stack for SI gravity.
- Local/global tick status: `particles/hierarchy/certificates/R_N_global_repair_tick_certificate.json` closes the global repair-tick theorem with derived round count: under the area-law counting model `F(N)=pi/rho_read^2` (a D6-consistent modeling identification of `Cap_read`), the corpus readback fixed-point equation `N_CRC=F(N_CRC)` is equivalent to the closure transport `G_N(1)=rho_star`, which forces the full-cycle multiplier `(N_CRC/pi)^(-1/2)`. The representation-to-spectrum certificate `R_m_rep_24_certificate.json` derives `m_rep=2*(8+3+1)=24` from the doubled observer-visible product-adjoint repair spectrum, giving one-tick contraction `|g_*'|=(N_CRC/pi)^(-1/48)` with no electroweak inputs.
- Local/global EW projection status: `particles/hierarchy/certificates/R_EW_tick_projection_certificate.json` closes the electroweak tick-projection bridge for the EW-refined capacity `N_CRC^EW`. The projection map `Pi_EW(P,N)=24*pi/(alpha_U(P)*log(N/pi))` is derived from the D10 transmutation law `t_tr=2*pi/(beta_EW*alpha_U)` and the global tick `|g_*'|=(N/pi)^(-1/48)` by ratio. Under the OPH local/global EW resonance condition `Pi_EW(P_star,N_CRC^EW)=4*P_star=beta_EW*P_star`, this is equivalent to the bridge residual `B_EW(P,N)=alpha_U(P)*log(N/pi)-6*pi/P=0`, with the integer factor `4=beta_EW=N_c+1` traced to the D10 transmutation theorem and `12=2*m_rep/beta_EW=48/4` traced algebraically to the global repair-tick and representation-to-spectrum theorems plus D10. The exact-capacity certificate `R_EW_global_capacity_certificate.json` supplies the unique source-side bridge fixed point `N_CRC^EW=pi*exp[6*pi/(P_star*alpha_U(P_star))]` as a `lambda=1/2` Banach contraction. The rounded `3.31e122` cosmological capacity display remains a diagnostic label. The certificate emits a machine-readable `derivation_chain`, `factor_origins`, `acceptance_criteria_status`, and scope record for the resonance target. The geometric strengthening that the certificate previously flagged as belonging to a separate theorem is now supplied by the icosahedral screen-sieve theorem `R_screen_sieve_icosahedral_certificate.json` (deriving the `1/12` port-read factor from the `A5/C5` orbit on the declared triangulated S^2 screen branch) and is wired into the umbrella closeout below.
- Local/global resonance status: `particles/hierarchy/certificates/R_local_global_hierarchy_resonance_closeout_335.json` closes the full local/global hierarchy-resonance theorem on the selected branch, composing the eight prerequisite theorems (global repair-tick, representation-to-spectrum, finite readback-resolution, joint product fixed point, icosahedral screen-sieve, EW tick-projection, EW-refined exact-capacity, RG/Higgs naturality) into the target relation `t_tr(P_star)=(P_star/12)*log(N_CRC^EW/pi)` via a machine-readable nine-step `derivation_chain`. The umbrella `factor_origins` traces every integer factor (`4=beta_EW` from D10; `24=2*dim(su(3)+su(2)+u(1))` from the representation-to-spectrum theorem; `48=2*m_rep` from the global tick; `12=|A5|/|C5|=60/5` from the icosahedral screen-sieve theorem; total curvature charge `12` from discrete Gauss-Bonnet) to a corpus theorem with an explicit source artifact. The only residual definitional residue surfaced honestly under `factor_origins.cell_entropy_factor_one_over_four` is the OPH structural identification `cell entropy = P/beta_EW = P/4`, consumed by the screen-sieve composition as input. The promotion flag `full_theorem_grade_resonance_promoted` is computed (not stipulated) from the conjunction of all dependency-status checks plus the bridge-residual checks; no promotion gates remain open in the hierarchy package on the selected branch (OPH product gauge branch + declared triangulated S^2 screen branch + selected RG/Higgs natural exact branch + product-separated joint branch).

## Bundle Execution Plan

Open work is grouped into coupled closure packets with bundle-level promotion questions.

| Bundle | Status | Gaps | Promotion question |
| --- | --- | --- | --- |
| `electroweak-root-closure-bundle` | `endpoint_package_closed_source_measure_payload_absent` | `pclosure.compressed-trunk-artifact`, `d10.ward-projected-thomson-endpoint`, `d10.source-residual-map-and-interval-certificate`, `d10.rg-matching-threshold-scheme`, `pclosure.certified-codepath-adoption` | Can one source-emitted map Delta_Th(P), with declared matching and interval bounds, certify the compressed P trunk as the particle root without importing alpha(0)? |
| `spectrum-source-bundle` | `closed_corpus_limited_source_boundaries_emitted` | `charged.determinant-normalization-transport`, `quark.selected-class-vs-global-classification`, `neutrino.pmns-status-and-absolute-rows` | Is there one OPH excitation dictionary and sector-isolated trace-lift theorem that explains the charged affine anchor, quark selected-class boundary, and neutrino PMNS comparison surface without hidden target fitting? |
| `strong-cp-closure-bundle` | `open_physical_invariant_gap` | `qcd.strong-cp-angle` | Can the exact quark/Yukawa branch be extended to the physical strong-CP invariant, including theta_QCD, bar(theta), and a vanishing theorem on the realized branch? |
| `qcd-thomson-backend-bundle` | `source_backend_boundary_empirical_policy_emitted` | `d10.ward-projected-thomson-endpoint`, `hadron.production-backend-systematics`, `hadron.empirical-ee-spectral-closure` | Can the source-only backend or the empirical e+e- spectral payload supply the Ward-projected hadronic term needed by the Thomson endpoint with row class recorded? |
| `top-codomain-bridge-bundle` | `closed_corpus_limited_codomain_no_go` | `calibration.direct-top-bridge` | Can the exact top coordinate be mapped into the auxiliary direct-top extraction codomain without using Q007TP as a calibration input? |
| `particle-root-integration-gate` | `keep_candidate_with_constructive_next_artifacts` | `pclosure.compressed-trunk-artifact`, `d10.ward-projected-thomson-endpoint`, `d10.rg-matching-threshold-scheme`, `pclosure.certified-codepath-adoption`, `charged.determinant-normalization-transport`, `quark.selected-class-vs-global-classification`, `neutrino.pmns-status-and-absolute-rows`, `calibration.direct-top-bridge`, `hadron.production-backend-systematics`, `hadron.empirical-ee-spectral-closure` | Do the returned packets jointly close the endpoint, matching, interval, and source-object requirements strongly enough to promote the compressed trunk into particle builders? |

## Bundle Packet Results

- `electroweak-root-closure-bundle`: `endpoint_package_closed_source_measure_payload_absent`. Constructive result. The admissible endpoint object is explicit and the endpoint package computes the residual inverse-alpha packet. The source-spectral reduction theorem is emitted. Delta_Th(P) must split into source lepton transport, a Ward-projected hadronic spectral density rho_had(s;P), a certified electroweak/scheme remainder, RG/matching certificates, quadrature bounds, and an interval-level fixed-point certificate. The local implementation targets are P_derivation/runtime/thomson_endpoint_contract_current.json and P_derivation/runtime/source_spectral_theorem_current.json.
- `spectrum-source-bundle`: `closed_corpus_limited_source_boundaries_emitted`. No promotion. Charged leptons are closed as a corpus-limited no-go by the end-to-end impossibility theorem: the same-family witness and conditional algebraic readout remain, but no theorem-grade A_ch(P) is emitted. Quarks remain selected-class on f_P with global classification closed as a corpus-limited no-go, and neutrino PMNS rows remain visible comparison-tension rows.
- `strong-cp-closure-bundle`: `open_physical_invariant_gap`. No promotion. The selected-class exact quark theorem fixes the running sextet and forward Yukawas on the public class f_P. The available corpus does not emit the determinant-line phase contribution, the bare theta_QCD term, or a theorem forcing the physical strong-CP phase to vanish.
- `qcd-thomson-backend-bundle`: `source_backend_boundary_empirical_policy_emitted`. Constructive result with two surfaces. The source-only primitive remains production_ward_projected_hadronic_spectral_measure_export and requires a real OPH hadron backend. The empirical surface uses a separate e+e- payload class and cannot promote the source-only theorem.
- `top-codomain-bridge-bundle`: `closed_corpus_limited_codomain_no_go`. No-go result. The exact top theorem row remains on Q007TP4. The auxiliary direct-top row Q007TP is compare-only because the available corpus emits no source-side extraction-response kernel into that codomain.
- `particle-root-integration-gate`: `keep_candidate_with_constructive_next_artifacts`. No promotion. The first wave emits constructive next artifacts, so the compressed P trunk remains candidate/audit metadata until those artifacts are populated and certified.

## Open Gates

| ID | Lane | Status | Issue | Next action |
| --- | --- | --- | --- | --- |
| `pclosure.compressed-trunk-artifact` | P closure | `closed_canonical_guarded_candidate_artifact` | #224 | Keep emitting p_closure_trunk_current.json and block live prediction promotion. |
| `d10.ward-projected-thomson-endpoint` | D10 electromagnetic endpoint | `closed_blocker_isolated_endpoint_package` | #223 | Keep the package as the closed blocker-isolation artifact for issue #223. |
| `d10.source-residual-map-and-interval-certificate` | D10 electromagnetic endpoint | `closed_blocker_isolated_source_residual_no_go` | #235 | Populate the Ward-projected source spectral measure payload, including rho_had(s;P) or an equivalent spectral primitive, matching remainder, certified quadrature bounds, and the interval certificate for the final map. |
| `d10.rg-matching-threshold-scheme` | D10 running and matching | `closed_declared_convention_contract` | #32 | Keep the declared-convention status visible in prediction surfaces and require a separate theorem before treating those conventions as OPH-derived. |
| `pclosure.certified-codepath-adoption` | P closure | `closed_guarded_codepath_adoption` | #224 | Switch live particle builders only after the source spectral measure payload emits R_Q(P), the interval certificate proves the full map, and the compressed trunk is promoted beyond guarded candidate metadata. |
| `charged.determinant-normalization-transport` | Charged leptons | `closed_corpus_limited_charged_end_to_end_no_go` | #201 | Keep charged masses suppressed on the public theorem lane. Reopen only for a theorem-grade uncentered trace lift proving 3 mu(r) = sum_e M_e^ch log q_e(r), equivalently zero normalization defect N_det(P), on the physical charged branch. |
| `quark.selected-class-vs-global-classification` | Quarks | `selected_class_closed_global_classification_no_go` | #199 | Keep every public exact-quark claim explicitly selected-class. Reopen only for a new source-emitted global public-frame classifier artifact. |
| `neutrino.pmns-status-and-absolute-rows` | Neutrinos | `theorem_rows_with_visible_comparison_tension` | #117 | Do not hide PMNS residuals behind the exact absolute-splitting rows; either prove a better branch or leave the comparison tension explicit. |
| `qcd.strong-cp-angle` | Strong CP | `open_theta_qcd_bar_theta_vanishing_gap` | #155 | Keep strong CP explicit as an open branch. Reopen only for a theorem-grade descent from exact quark/Yukawa phase data to the determinant-line phase contribution, together with a theorem fixing the topological-angle contribution and proving the physical strong-CP phase vanishes on the realized branch. |
| `calibration.direct-top-bridge` | D11/top codomain | `closed_corpus_limited_codomain_no_go` | #207 | Keep Q007TP compare-only while the theorem row remains anchored on Q007TP4. Reopen only for a concrete source-side extraction-response kernel. |
| `hadron.production-backend-systematics` | Hadrons | `source_backend_absent_empirical_policy_emitted` | #153 | Keep source-only hadron rows suppressed. Use empirical hadron closure rows only through the documented e+e- spectral payload. Promote source-only hadron rows only after a working OPH hadron backend emits production hadron output, Ward-projected spectral data, and systematics. |
| `hadron.empirical-ee-spectral-closure` | Hadrons | `policy_scaffold_emitted_dataset_absent` | n/a | Populate oph_empirical_ee_hadronic_spectral_measure from PDG, HEPData, alphaQED, or an equivalent documented compilation, then feed the empirical Thomson endpoint builder. |

## Claim Policy

- The compressed P trunk is an audit/candidate artifact until the endpoint and certificate gates close.
- Open blockers should be worked as coupled bundles.
- The particle pipeline must keep compare-only, continuation, selected-class, and theorem-grade rows mechanically distinct.
- Golden-ratio torus or resonance language requires a separate representation-to-spectrum theorem before it can serve as a derivation input.
