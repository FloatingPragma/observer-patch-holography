# Particle Derivation Gap Ledger

Generated: `2026-07-12T02:19:29Z`

Systematic claim-safe queue after the five-equation P-trunk simplification.

## P-Trunk Claim Boundary

- Artifact: `P_derivation/runtime/p_closure_trunk_current.json`
- Exists: `True`
- Claim label: `compressed_candidate_trunk_not_final_particle_root`
- May feed promoted particle predictions: `False`
- Candidate P: `1.63097210492078846050203640439`
- Candidate alpha^-1: `136.994020662724205139718642793`
- Source report mode: `thomson_structured_running_asymptotic`

## Electroweak Hierarchy Certificate

- Artifact: `particles/hierarchy`
- Exists: `True`
- Claim label: `closed_local_global_hierarchy_and_naturality_certificate`
- May feed local hierarchy claim: `True`
- May feed Higgs naturality claim: `True`
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
- Local/global resonance status: `closed_full_local_global_hierarchy_resonance`
- Full theorem-grade resonance promoted: `True`
- Remaining promotion gates: `[]`
- Exact EW bridge capacity: `3.5323546226929906511187512962330547600462096590942035604238177731136802717148740820434110040644403858228611984E+122`
- Bridge residual: `0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`
- Fixed-point residual in log capacity: `0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000`
- Source v/E_cell: `2.5797040595276223506583341903922018695490035369903290256096683452314630559578761079150341073645683725521032233E-17`
- Higgs naturality defect: `epsilon_H=0`
- Higgs naturality interval: `[0, 0]`
- Boundary: This certificate closes the selected local P -> alpha_U -> v/E_star hierarchy lane, the local/global resonance bridge, and the Higgs naturality defect epsilon_H=0. Separate non-promoted gates are the public Thomson endpoint transport, theorem-grade W/Z promotion, charged-lepton absolute masses, source-only hadron masses, Strong CP, and the full no-G clock stack for SI gravity.

## Bundle Claim Gates

Claim gates are grouped into coupled closure packets rather than a one-blocker-at-a-time queue.

| Bundle | Claim label | Gaps | Promotion question |
| --- | --- | --- | --- |
| `electroweak-root-closure-bundle` | `endpoint_package_closed_source_measure_payload_absent` | `pclosure.compressed-trunk-artifact`, `d10.ward-projected-thomson-endpoint`, `d10.source-residual-map-and-interval-certificate`, `d10.rg-matching-threshold-scheme`, `d10.same-scheme-anchor-bridge`, `pclosure.certified-codepath-adoption` | Can one source-emitted map Delta_Th(P), with declared matching and interval bounds, certify the compressed P trunk as the particle root without importing alpha(0)? |
| `spectrum-source-bundle` | `closed_corpus_limited_source_boundaries_emitted` | `charged.determinant-normalization-transport`, `quark.selected-class-vs-global-classification`, `neutrino.pmns-status-and-absolute-rows` | Is there one OPH excitation dictionary and sector-isolated trace-lift theorem that explains the charged affine anchor and quark selected-class boundary while also deriving a source-closed neutrino operator, charged basis, and mass-label rule without target fitting? |
| `strong-cp-closure-bundle` | `open_physical_invariant_gap` | `qcd.strong-cp-angle` | Can the exact quark/Yukawa branch be extended to the physical strong-CP invariant, including theta_QCD, bar(theta), and a vanishing theorem on the realized branch? |
| `qcd-thomson-backend-bundle` | `source_backend_boundary_empirical_policy_emitted` | `d10.ward-projected-thomson-endpoint`, `hadron.production-backend-systematics`, `hadron.empirical-ee-spectral-closure` | Can the source-only backend or the empirical e+e- spectral payload supply the Ward-projected hadronic term needed by the Thomson endpoint with row class recorded? |
| `top-codomain-bridge-bundle` | `closed_corpus_limited_codomain_no_go` | `calibration.direct-top-bridge` | Can the exact top coordinate be mapped into the auxiliary direct-top extraction codomain without using Q007TP as a calibration input? |
| `particle-root-integration-gate` | `keep_candidate_with_constructive_next_artifacts` | `pclosure.compressed-trunk-artifact`, `d10.ward-projected-thomson-endpoint`, `d10.rg-matching-threshold-scheme`, `pclosure.certified-codepath-adoption`, `charged.determinant-normalization-transport`, `quark.selected-class-vs-global-classification`, `neutrino.pmns-status-and-absolute-rows`, `calibration.direct-top-bridge`, `hadron.production-backend-systematics`, `hadron.empirical-ee-spectral-closure` | Do the returned packets jointly close the endpoint, matching, interval, and source-object requirements strongly enough to promote the compressed trunk into particle builders? |

## Bundle Packet Results

- `electroweak-root-closure-bundle`: `endpoint_package_closed_source_measure_payload_absent`. Constructive result. The admissible endpoint object is explicit and the endpoint package computes the residual inverse-alpha packet. The source-spectral reduction theorem is emitted. Delta_Th(P) must split into source lepton transport, a Ward-projected hadronic spectral density rho_had(s;P), a certified electroweak/scheme remainder, RG/matching certificates, quadrature bounds, and an interval-level fixed-point certificate. The local implementation targets are P_derivation/runtime/thomson_endpoint_contract_current.json and P_derivation/runtime/source_spectral_theorem_current.json.
- `spectrum-source-bundle`: `closed_corpus_limited_source_boundaries_emitted`. No promotion. Charged leptons are closed as a corpus-limited no-go by the end-to-end impossibility theorem: the same-family witness and conditional algebraic readout remain, but no theorem-grade A_ch(P) is emitted. Quarks remain selected-class on f_P with global classification closed as a corpus-limited no-go. The weighted-cycle neutrino point is a rejected target-informed candidate; no physical PMNS, ordering, Majorana, or absolute-mass row is emitted.
- `strong-cp-closure-bundle`: `open_physical_invariant_gap`. No promotion. The selected-class quark support wrapper conditionally carries the running sextet and forward Yukawas on the public class f_P, but the source sigma selector is still missing. The available corpus does not emit the determinant-line phase contribution, the bare theta_QCD term, or a theorem forcing the physical strong-CP phase to vanish.
- `qcd-thomson-backend-bundle`: `source_backend_boundary_empirical_policy_emitted`. Constructive result with two surfaces. The source-only primitive remains production_ward_projected_hadronic_spectral_measure_export and requires a real OPH hadron backend. The empirical surface uses a separate e+e- payload class and cannot promote the source-only theorem.
- `top-codomain-bridge-bundle`: `closed_corpus_limited_codomain_no_go`. No-go result. The exact top theorem row remains on Q007TP4. The auxiliary direct-top row Q007TP is compare-only because the available corpus emits no source-side extraction-response kernel into that codomain.
- `particle-root-integration-gate`: `keep_candidate_with_constructive_next_artifacts`. No promotion. The first wave emits constructive next artifacts, so the compressed P trunk remains candidate/audit metadata until those artifacts are populated and certified.

## Claim Gates

| ID | Lane | Claim label | Gate |
| --- | --- | --- | --- |
| `pclosure.compressed-trunk-artifact` | P closure | `closed_canonical_guarded_candidate_artifact` | Keep emitting p_closure_trunk_current.json and block live prediction promotion. |
| `d10.ward-projected-thomson-endpoint` | D10 electromagnetic endpoint | `closed_blocker_isolated_endpoint_package` | Keep the package as the closed blocker-isolation artifact for issue #223. |
| `d10.source-residual-map-and-interval-certificate` | D10 electromagnetic endpoint | `closed_blocker_isolated_source_residual_no_go` | Populate the Ward-projected source spectral measure payload, including rho_had(s;P) or an equivalent spectral primitive, matching remainder, certified quadrature bounds, and the interval certificate for the final map. |
| `d10.rg-matching-threshold-scheme` | D10 running and matching | `closed_declared_convention_contract` | Keep the declared-convention status visible in prediction surfaces and require a separate theorem before treating those conventions as OPH-derived. |
| `pclosure.certified-codepath-adoption` | P closure | `closed_guarded_codepath_adoption` | Switch live particle builders only after the source spectral measure payload emits R_Q(P), the interval certificate proves the full map, and the compressed trunk is promoted beyond guarded candidate metadata. |
| `charged.determinant-normalization-transport` | Charged leptons | `closed_corpus_limited_charged_end_to_end_no_go` | Keep charged masses suppressed. Reopen only for the theorem-grade sector-isolated normalized charged determinant-line attachment named by #546; the gate flips only when N_det has the certified singleton interval [0,0]. |
| `quark.selected-class-vs-global-classification` | Quarks | `selected_class_closed_global_classification_no_go` | Keep every public exact-quark claim explicitly selected-class. Reopen only for a new source-emitted global public-frame classifier artifact. |
| `neutrino.pmns-status-and-absolute-rows` | Neutrinos | `rejected_candidate_source_basis_and_kernel_open` | Derive a source-closed neutrino operator, stable physical charged-lepton left basis, and mass-eigenstate label/order rule without oscillation-target feedback; freeze that construction before evaluating a later likelihood. Keep all present weighted-cycle, bridge, and exact-adapter numbers rejected or compare-only. |
| `qcd.strong-cp-angle` | Strong CP | `open_theta_qcd_bar_theta_vanishing_gap` | Keep strong CP explicit as an open branch. Reopen only for a theorem-grade descent from exact quark/Yukawa phase data to the determinant-line phase contribution, together with a theorem fixing the topological-angle contribution and proving the physical strong-CP phase vanishes on the realized branch. |
| `calibration.direct-top-bridge` | D11/top codomain | `closed_corpus_limited_codomain_no_go` | Keep Q007TP compare-only while the theorem row remains anchored on Q007TP4. Reopen only for a concrete source-side extraction-response kernel. |
| `hadron.production-backend-systematics` | Hadrons | `source_backend_absent_empirical_policy_emitted` | Keep source-only hadron rows suppressed. Use empirical hadron closure rows only through the documented e+e- spectral payload. Promote source-only hadron rows only after a working OPH hadron backend emits production hadron output, Ward-projected spectral data, and systematics. |
| `hadron.empirical-ee-spectral-closure` | Hadrons | `payload_populated_endpoint_evaluated_gap_anchor_localized` | Emit the source-side electroweak scheme bridge for a0(P) that produces the certified anchor-gap interval; refine the payload with experiment-level tables when a finer compilation is ingested. |
| `d10.repair-tuple-selection` | D10/D11 electroweak masses | `conditional_selection_theorem_axioms_named` | Derive the descent axioms A2 and A3 from the realized carrier package (the color-singlet projection weight and the coherent neutral color sum), which closes the selection unconditionally and unblocks the D10 promotion review. |
| `d10.same-scheme-anchor-bridge` | P closure / D10 electromagnetic endpoint | `structure_resolved_reduces_to_source_hadron_backend` | A source-only anchor bridge reduces to the OPH hadronic spectral measure, blocked on the hadron backend (#425). No source-only scheme-bridge theorem exists on the current corpus; #545 stays open as that reduction. The empirical-class bridge is the built payload/endpoint pair. |

## Claim Policy

- The compressed P trunk is an audit/candidate artifact until the endpoint and certificate gates close.
- Claim gates are handled as coupled bundles rather than isolated one-off fixes.
- The particle pipeline must keep compare-only, continuation, selected-class, and theorem-grade rows mechanically distinct.
- Golden-ratio torus or resonance language is not a derivation input unless a separate representation-to-spectrum theorem is supplied.
