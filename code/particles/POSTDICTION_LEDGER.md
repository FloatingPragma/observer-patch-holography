# Postdiction Ledger

Generated: `2026-07-23T06:54:31Z` by `scripts/build_postdiction_ledger.py`; the JSON artifact is `runs/status/postdiction_ledger.json`.

Every value and every measured reference on this page is read live from the cited certified artifact. The ledger promotes nothing, changes no solve path, and introduces no number of its own. Interval rows report containment of the compare-only witness; conditional rows carry their declared premises; chart coordinates keep their NOT_EVALUABLE physical-comparison status.

## Principal results

- The anchor-gap value 0.6379 closes the charged-lepton lane exactly on the measured triple, inside the certified band [0.6199, 0.6506]; the distance +0.0070 to the standard on-shell reference deficit 0.6309 is the live scheme term of the open anchor bridge (issue 545). The lepton scale is localized to the width of the scheme band, and a source-emitted bridge value is a falsification target: the closure value confirms, a value outside the band refutes.
- The measured charged-lepton triple lies inside every certified interval; the payload-coherent half-width is 1.73 percent per lepton, and the conditional eight-register triple sits 84 ppm from measurement with the architecture declared.
- The conditional Higgs envelope [125.183, 125.232] GeV sits 0.70 sigma from the measured 125.13 +- 0.11 GeV, and the top envelope [172.28, 172.35] GeV sits 0.36 sigma from 172.1 +- 0.6 GeV, compare-only, conditional on the declared selection axioms.
- The gauge sector is pinned before any numeric lane runs: the twelve-port trichotomy forces su(3)+su(2)+u(1), the gluing-class quotient gives the Z6 global form, and the matter lift realizes the exact one-generation hypercharge multiset, with the finite steps machine checked in Lean/Screen and the hypothesis boundaries recorded below.

## Forced structure

The icosahedral screen results pin the gauge sector before any numeric lane runs. The finite steps are machine checked in the Lean workspace; the recorded hypothesis boundaries are the exact classical inputs and open premises of the compact paper.

| Result | Observed counterpart | Match | Receipts |
| --- | --- | --- | --- |
| Compact-Lie trichotomy on the twelve-port screen: a compact connected group with a group-level A5 action equivalent to P12 has Lie algebra u(1)^12, su(2)^2+u(1)^6, or su(3)+su(2)+u(1); the noncentral quintet and the inner-action closure each select su(3)+su(2)+u(1) | Standard Model gauge Lie algebra su(3)+su(2)+u(1) | `exact` | `Lean/Screen/A5OPH.lean`, `Lean/Screen/A5CharacterField.lean`, `Lean/Screen/A5SixAxes.lean` |
| The screen gluing-class quotient Lambda_+/(Lambda_1 + Lambda_5) is Z/6 with proper-rotation invariance and antipodal sign reversal, matching the global form (SU(3) x SU(2) x U(1))/Z6 | Standard Model global gauge-group form and its charge quantization pattern | `exact` | `Lean/Screen/Z6Exact.lean`, `Lean/Screen/A5OPH.lean` |
| The super-Tannakian matter lift realizes the one-generation left-chiral hypercharge multiset {Q: 1/6 x6, u_c: -2/3 x3, d_c: 1/3 x3, L: -1/2 x2, e_c: 1 x1} | Standard Model one-generation hypercharge assignment | `exact` | `code/a5_closure/receipts/super_tannakian_matter_reference.receipt.json` |
| A5-invariant readouts have port-independent group-averaged cap sums, so the per-cap ratio of any two averaged readouts is universal with zero spread | universality clause of the Einstein-branch coupling law | `structural` | `Lean/Screen/A5CouplingSymmetry.lean`, `Lean/Screen/A5PortAction.lean`, `Lean/Screen/PortFrameGram.lean` |

Hypothesis boundaries:

- `gauge_lie_algebra`: the compact-simple classification, the torus/cocharacter step of the rationality lemma, and irreducibility of the five-dimensional summand stay declared classical inputs on paper; the physical inner current action is the open premise of issues 567 and 599
- `global_form_z6`: the quotient isomorphism and both invariance clauses are machine checked; the identification with the physical global form rides on the same inner current action premise as the Lie-algebra row
- `hypercharge_spectrum`: conditional on the declared super-Tannakian lift premises recorded in the receipt claim boundary (issue 314); the generation count is an input, not a consequence
- `coupling_universality`: reduces the universality clause to A5-equivariance of the implemented source law; no coupling value is implied

## Fine-structure lane

- `alpha_em^-1` Thomson endpoint: `136.3827548` in `[136.3670481, 136.3984652]` against CODATA `137.0359992` (compare-only). Payload release `knt19_pinned_v1`.
- Certified same-scheme anchor gap `[0.6199, 0.6506]` inverse-alpha units; the standard reference deficit sits inside the certified interval.
- Reading: the endpoint carries the full measured deficit into the certified same-scheme anchor gap; the standard reference deficit sits inside the certified interval, so the residual is the one-loop anchor running deficit, not a payload or solve defect
- Blocking issues: #425, #545

## Charged leptons

- Closure target (T1_empirical_closure): the anchor-gap value `0.6379` closes the lane exactly on the measured triple (inversion machine-checked); the distance `+0.0070` to the on-shell reference deficit `0.6309` is the live scheme term of the bridge. The certified width floor is the scheme-band ambiguity; no budget is shrunk without the source bridge.
- MCPR conditional triple (T2): electron `-84.1 ppm`, muon `-84 ppm`, tau `-84 ppm` against the PDG witness triple; the eight-register architecture is a declared model input.
- Kappa interval, rectangle (T1_empirical_closure): certified relative half-widths electron `6.56%`, muon `6.56%`, tau `6.56%`; the witness triple lies inside every interval.
- Kappa interval, coherent closure (T1_empirical_closure): certified relative half-widths electron `1.73%`, muon `1.73%`, tau `1.73%`; the witness triple lies inside every interval.
  - Width reduction over the rectangle: `3.78x`; premise: payload-coherent anchor-gap premise, declared.

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
- Down-type texture, conditional (T2_conditional): Cabibbo `0.2086` against `0.225`; `ms/md = 22.97` against `19.9`. Premise: generation register order (issue 569); pairing and weight set selected by the Clebsch selection artifact. Conditional texture rows: the register order premise is open and the recorded tensions stay in the normalization_tension block of the parent.

## Hadrons

- Correction engine payload: `Delta alpha_had^(5)(M_Z^2) = 0.027609 +- 0.000112` from `knt19_pinned_piecewise_v1` (pin factor `1.03176`). The published-compilation payload is the correction engine of the fine-structure lane; source-only hadron rows stay suppressed pending the source spectral measure (issue 425).
- QCD solver: `SOLVER_COMPILED_AND_SMOKE_BLOCKED_INVOCATION_GATED_ON_SOURCE_PARAMETERS`; invocation is gated on the source-side parameter emissions recorded in the standby receipt.

## Neutrinos

- dimensionless PMNS and mass-splitting-ratio comparisons live on the results status surface; the absolute attachment stays compare-only (`code/particles/RESULTS_STATUS.md`).
