# Postdiction Ledger

Generated deterministically by `scripts/build_postdiction_ledger.py`; the JSON artifact is `runs/status/postdiction_ledger.json`.

Numeric values and measured references on this page are read live from the cited parent artifacts. Structural rows are derived from validated structured parents, and direct algebraic corollaries are identified. The ledger promotes nothing and changes no solve path. Interval rows report containment of the compare-only witness; conditional rows carry their declared premises; chart coordinates keep their NOT_EVALUABLE physical-comparison status.

## Principal results

- The anchor-gap value 0.6379 closes the charged-lepton lane exactly on the measured triple, inside the retrospective accounting interval [0.6199, 0.6506]; the distance +0.0070 to the standard on-shell reference deficit 0.6309 is the live scheme term of the open anchor bridge (issue 545). The lepton scale is localized only under that recorded accounting packet. A source-emitted bridge value is a falsification target: the closure value would satisfy the conditional lane, while a value outside the interval refutes the declared decomposition.
- The target-anchored measured charged-lepton triple lies inside every outward-rounded diagnostic interval; the payload-coherent logarithmic half-width is 1.732 percent, with one-sided multiplicative widths -1.72 and +1.75 percent. The conditional eight-register triple sits 84 ppm from measurement with the architecture declared.
- Under the balanced-circulant and mass-ordering premises the measured electron and muon masses fix the tau mass inside [1776.968991, 1776.969063] MeV, 0.4336 sigma from measurement; the window is three orders of magnitude narrower than the measurement uncertainty, so improving tau-mass averages test the premise directly. The premise ancestry is declared and the row stays conditional.
- The conditional Higgs envelope [125.183, 125.232] GeV sits 0.70 sigma from the measured 125.13 +- 0.11 GeV, and the top envelope [172.28, 172.35] GeV sits 0.36 sigma from 172.1 +- 0.6 GeV, compare-only, conditional on the declared selection premises.
- The gauge sector is pinned before any numeric lane runs: the twelve-port trichotomy forces su(3)+su(2)+u(1), the gluing-class quotient gives the Z6 global form, and the exhaustive scan inside the declared exterior-response algebra selects the charge-conjugate rank-15 chiral anomaly-free pair and its one-generation hypercharge multiset. The finite steps are machine checked in Lean/Screen, with the hypothesis boundaries recorded below.

## Forced structure

The icosahedral screen results pin the gauge sector before any numeric lane runs. The finite steps are machine checked in the Lean workspace; the recorded hypothesis boundaries are the exact classical inputs and open premises of The Standard Model gauge paper.

| Result | Observed counterpart | Match | Receipts |
| --- | --- | --- | --- |
| Compact-Lie trichotomy on the twelve-port screen: a compact connected group with a group-level A5 action equivalent to P12 has Lie algebra u(1)^12, su(2)^2+u(1)^6, or su(3)+su(2)+u(1); the noncentral quintet and the inner-action closure each select su(3)+su(2)+u(1) | Standard Model gauge Lie algebra su(3)+su(2)+u(1) | `exact` | `Lean/Screen/A5OPH.lean`, `Lean/Screen/A5CharacterField.lean`, `Lean/Screen/A5SixAxes.lean` |
| The screen gluing-class quotient Lambda_+/(Lambda_1 + Lambda_5) is Z/6 with proper-rotation invariance and antipodal sign reversal, matching the global form (SU(3) x SU(2) x U(1))/Z6 | Standard Model global gauge-group form and its charge quantization pattern | `exact` | `Lean/Screen/Z6Exact.lean`, `Lean/Screen/A5OPH.lean` |
| Inside the declared 10-component exterior-response algebra, an exhaustive scan of all 1024 subsets selects exactly one unordered charge-conjugate pair of nonempty chiral anomaly-free rank-15 projectors. Primitive determinant balance fixes the block charges up to conjugation, and the selected representative has multiset {Q: 1/6 x6, u_c: -2/3 x3, d_c: 1/3 x3, L: -1/2 x2, e_c: 1 x1} | Standard Model one-generation hypercharge assignment | `exact` | `Lean/Screen/ExteriorSelection.lean`, `code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json`, `code/a5_closure/manifests/matter_menu_spectral_ledger_reference.json` |
| A5-invariant readouts have port-independent group-averaged cap sums, so the per-cap ratio of any two averaged readouts is universal with zero spread | universality clause of the Einstein-branch coupling law | `structural` | `Lean/Screen/A5CouplingSymmetry.lean`, `Lean/Screen/A5PortAction.lean`, `Lean/Screen/PortFrameGram.lean` |
| On the declared unbroken Maxwell action and deconfined phase branch, the quadratic operator has zero hard mass parameter and two transverse classical modes with characteristic surface k^2=0 | massless classical electromagnetic propagation | `conditional structural` | `code/particles/runs/status/carrier_mode_acceptance.json` |
| On the declared pure Yang-Mills quadratic branch before nonperturbative confinement, every color generator has two transverse perturbative modes and zero hard quadratic mass parameter | perturbative color-gauge kernel before confinement | `conditional structural` | `code/particles/runs/status/carrier_mode_acceptance.json` |
| On the declared pure Einstein-Hilbert linearization about a suitable Ricci-flat background, the transverse-traceless quadratic operator has zero hard mass parameter and two classical modes with null characteristic | two massless classical gravitational-wave polarizations | `conditional structural` | `code/particles/runs/status/carrier_mode_acceptance.json` |
| The selected direct-sum current algebra has adjoint branch dimensions 8, 3, and 1. Its adjoint therefore contains no mixed (3,2,-5/6) (+) (bar3,2,+5/6) X/Y generator, so the ordinary minimal simple-GUT X/Y exchange channel is absent | no observed proton decay through the minimal simple-GUT X/Y channel | `structural channel exclusion` | `code/a5_closure/receipts/port_current_inner_reference.receipt.json` |

Lean declaration bindings:

- `gauge_lie_algebra`: `A5OPH`: `sum_eq_eleven`, `sum_eq_twelve`, `action_trivial_of_card_le_four`, `sum_not_mem_excluded`, `quintet_noncentral`; `A5CharacterField`: `multiplicities_equal_of_galoisStable`, `centreDim_mem_trichotomy_list`; `A5SixAxes`: `two_transitive`, `V5_irreducible`, `no_three_plus_three_split`
- `global_form_z6`: `Z6Exact`: `gauge_eq_kernel`, `residue_surjective`, `representative_formula`; `A5OPH`: `screen_gluing_class`, `q_perm_invariant`, `q_neg`
- `hypercharge_spectrum`: `ExteriorSelection`: `selection_unique`, `parity_sectors_survive`, `conj_exchanges_survivors`, `witten_automatic`
- `coupling_universality`: `A5CouplingSymmetry`: `groupAverage_port_independent`, `coupling_ratio_universal`; `A5PortAction`: `transitive_on_ports`; `PortFrameGram`: `degree_five`, `gram_sq`

Hypothesis boundaries:

- `gauge_lie_algebra`: the compact-simple classification, the torus/cocharacter step of the rationality lemma, and irreducibility of the five-dimensional summand stay declared classical inputs on paper; the physical inner current action is the open premise of issues 567 and 599
- `global_form_z6`: the quotient isomorphism and both invariance clauses are machine checked; the identification with the physical global form rides on the same inner current action premise as the Lie-algebra row
- `hypercharge_spectrum`: the selection is exhaustive inside the declared exterior algebra; completeness beyond that algebra, selection of one charge-conjugate representative, light-sector attachment, family multiplicity, scalar content, and laboratory identification remain separate
- `coupling_universality`: reduces the universality clause to A5-equivariance of the implemented source law; no coupling value is implied
- `maxwell_classical_massless_kernel`: the Maxwell action, positive kinetic coefficient, field content, and phase are supplied branch data; no photon Hilbert space, positive-residue pole, or universal zero-mass particle theorem is emitted
- `yang_mills_classical_massless_kernel`: this is not a free asymptotic-gluon claim and supplies neither a continuum Yang-Mills gap nor a hadron mass
- `einstein_classical_massless_kernel`: the action and background are supplied branch data; no graviton Hilbert space, quantum pole, or exclusion of additional massive modes is emitted
- `simple_gut_xy_channel_absent`: the corollary applies to the source-bound selected finite current algebra. General proton stability does not follow; higher-dimensional baryon violation, scalar mediators, and other ultraviolet gauge mechanisms are not excluded

## Fine-structure lane

- `alpha_em^-1` Thomson endpoint: `136.3827548` in `[136.3670481, 136.3984652]` against CODATA `137.0359992` (compare-only). Payload release `knt19_pinned_v1`.
- Recorded retrospective same-scheme accounting interval `[0.6199, 0.6506]` inverse-alpha units; the standard reference deficit sits inside that interval.
- Independent-class verdict: `MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE`; evaluated independent classes: `0`.
- Reading: one retrospective KNT19 accounting row is arithmetically compatible with the recorded same-scheme interval. The multi-class HVP test is not evaluable because no independent frozen class is present. Containment does not identify the physical source of the gap or close source-only transport
- Blocking issues: #425, #545

## Charged leptons

- Closure target (T1_empirical_closure): the anchor-gap value `0.6379` closes the lane exactly on the measured triple (inversion machine-checked); the distance `+0.0070` to the on-shell reference deficit `0.6309` is the live scheme term of the bridge. The certified width floor is the scheme-band ambiguity; no budget is shrunk without the source bridge.
- MCPR conditional triple (T2): electron `-84.1 ppm`, muon `-84 ppm`, tau `-84 ppm` against the PDG witness triple; the eight-register architecture is a declared model input.
- Kappa interval, rectangle (T1_empirical_closure): outward-rounded target-anchored diagnostic intervals with logarithmic half-width `6.554%` and one-sided multiplicative widths `-6.34%` / `+6.77%`; the witness triple lies inside every interval.
- Kappa interval, coherent closure (T1_empirical_closure): outward-rounded target-anchored diagnostic intervals with logarithmic half-width `1.732%` and one-sided multiplicative widths `-1.72%` / `+1.75%`; the witness triple lies inside every interval.
  - Width reduction over the rectangle: `3.78x`; premise: payload-coherent anchor-gap premise, declared.
- Koide conditional tau (T2_conditional): under the balanced-circulant and mass-ordering premises the measured electron and muon masses fix the tau mass inside `[1776.968991, 1776.969063]` MeV, `0.4336` sigma from the measured `1776.93 +- 0.09` MeV; the premise ancestry is declared and improving tau-mass averages test the premise directly.

## Electroweak sector

| Quantity | Conditional central | Envelope | Measured | Delta/sigma | Status |
| --- | ---: | --- | --- | ---: | --- |
| `mH_gev` | `125.20748` | `[125.18329, 125.23167]` | `125.13 +- 0.11` (PDG 2025) | `0.704` | compare-only |
| `mt_pole_gev` | `172.31492` | `[172.27749, 172.35236]` | `172.1 +- 0.6` (PDG 2025 direct-average context row) | `0.358` | compare-only |
| `MW_chart_gev` | `80.373315` | `[80.369217, 80.377413]` | chart coordinate | n/a | NOT_EVALUABLE |
| `MZ_chart_gev` | `91.193124` | `[91.187978, 91.198269]` | chart coordinate | n/a | NOT_EVALUABLE |

W/Z rows are running/tree chart coordinates; no physical comparison is defined until the chart-to-pole map is complete. The Higgs and top rows are conditional on the declared selection axioms.

## Quarks

- Absolute masses (source_only_nonidentifiability_obstruction_transport): No absolute quark mass is emitted: the two-modulus spread fiber survives every certified structure transport, so the six absolute masses are non-identifiable from the corpus, by theorem rather than by omission (issues #591).
- Down-type register-Clebsch route, rejected (T2_conditional_rejected_candidate): `ms/md = 22.97` against FLAG 2024 (Nf=2+1+1: 19.94, Nf=2+1: 20.36); all six generation assignments are rejected by the retrospective conservative gate. The diagnostic `sqrt(md/ms) = 0.2086` is not a derived Cabibbo angle. Premise: a cross-sector register relation, independent Yukawa coefficient identification, and a physical generation order; the pairing receipt supplies channel compatibility only. The target-free F1/F2 scan fixes only the unordered multiset. All six assignments fail the retrospective conservative FLAG gate. The displayed GST value is sqrt(md/ms) under an assumed texture, not a derived CKM angle; a simultaneous diagonal mass ansatz would instead give the identity CKM matrix.

## Hadrons

- Correction engine payload: `Delta alpha_had^(5)(M_Z^2) = 0.027609 +- 0.000112` from `knt19_pinned_piecewise_v1` (pin factor `1.03176`). The published-compilation payload is the correction engine of the fine-structure lane; source-only hadron rows stay suppressed pending the source spectral measure (issue 425).
- QCD solver: `SOLVER_COMPILED_AND_SMOKE_BLOCKED_INVOCATION_GATED_ON_SOURCE_PARAMETERS`; invocation is gated on the source-side parameter emissions recorded in the standby receipt.

## Neutrinos

- dimensionless PMNS and mass-splitting-ratio comparisons live on the results status surface; the absolute attachment stays compare-only (`code/particles/RESULTS_STATUS.md`).
