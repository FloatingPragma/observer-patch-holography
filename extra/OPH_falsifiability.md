# OPH Falsifiability Surface: Tests, Tensions, And Permanent Records

## Standing Of The Hypothesis

OPH is stated under the strange-loop principles (SLH), SL-0 through
SL-5, canonically recorded in [STRANGE_LOOP_PRINCIPLES.md](../STRANGE_LOOP_PRINCIPLES.md). The principles are declared
architecture postulates. Evidence bears on the principle set through the closure
ledger and the compression scorecard, never through restatement.

The two substrate constants are handled in three layers
([STRANGE_LOOP_PRINCIPLES.md](../STRANGE_LOOP_PRINCIPLES.md)). The theory layer
takes zero quantitative inputs: P is defined as the solution of
P = φ + √π/A_T(P), interval-certified to exist and be unique on its interval,
and N is defined as the solution of N = F(N), with F specified and its
construction open (CL-7). Measured values sit at the test layer, as the
comparisons the forced outputs face, and at the working layer, where open
lanes borrow located values:

- **P = 1.630968209403959…** is the working value of the substrate pixel,
  located from measured α: CODATA α⁻¹ = 137.035999177(21), basin half-width
  in P ≈ 2.0×10⁻¹². The working P is located from the measured fine-structure
  constant via principle SL-3 (self-simulation constant identification), and
  it is an input to every downstream row.
- **N = 3.31×10¹²²** is the working value of the record capacity, located
  under SL-4 from measured Λ with the scale bridge. It is likewise a counted
  input.

Every working borrow is tied to a named generator in
[CONSISTENCY_STACK.md](../CONSISTENCY_STACK.md) (G1 hadronic transport, G2
capacity readback, G3 electroweak repair); a generator that lands retires its
borrow.

Neither α nor Λ nor G appears on any OPH output list. A measured value used to
locate a basin is an input at every surface, permanently
([STRANGE_LOOP_PRINCIPLES.md](../STRANGE_LOOP_PRINCIPLES.md), rule 1). The current evidential standing is the
compression scorecard ([COMPRESSION_SCORECARD.md](../COMPRESSION_SCORECARD.md)):
M − K is deeply negative. The program's task is to move rows upward, and this
document lists every way it can lose while trying.

The exposure documented here is structural. Under the
uniqueness lemmas of the consistency stack
([CONSISTENCY_STACK.md](../CONSISTENCY_STACK.md), L1–L3), the SLH principle set
admits at most one (P, N): there is no vacuum landscape to relocate into, and
a blind closure computation landing outside its basin falsifies the
formulation permanently. The uniqueness result is global on the declared
physical domain: exactly one fixed point per readout map
(`code/P_derivation/runtime/p_global_uniqueness_certificate_2026-07-14.json`).
Maximal falsifiability is a corollary of that uniqueness.

This surface lists **33 hard OPH-killing outcomes**: 16 structural, 7
discriminating physical, and 10 shared with the Standard Model. Each executed
test stays here permanently with its verdict. Every numeric comparison shows
the relative deviation and the pull in measurement sigma; choosing one unit per
row is the sigma-washing pattern and is retired.

Useful background:

- [OPH Textbooks](https://learn.floatingpragma.io/)
- [Reverse Engineering Reality](https://oph-book.floatingpragma.io/)
- [Compact proof of OPH](compact_proof_of_oph.pdf)

## Quick Dictionary

| Term | Plain meaning |
| --- | --- |
| Falsifiable | A real test can make the claim lose. |
| Observer patch | The part of reality one observer, device, or local system can access. |
| Overlap | The shared part two patches can both check. |
| Repair | The settling step when two shared descriptions disagree. |
| Fixed point | The stable answer after checking stops changing the result. |
| Closure residual | The gap in an equation SLH requires to hold exactly. Under SL-0 it is a measurement of the principle set. |
| Basin | The interval a measurement locates for a substrate parameter. Locating a basin is stage 1 of the four-stage protocol; it is never a hit. |
| Contraction certificate | An interval-grade proof that a fixed-point map has one solution in a declared interval. Stage 2. A numerical witness is not a certificate. |
| Frozen prediction | A registry entry (quantity, formula, value with uncertainty, dependency-cone hash, external timestamp) committed before the comparison data exists or is looked at. |
| Benchmark check | A comparison with known data. Weaker than a frozen prediction. |
| Calibration row | A row that consumed its measured target. It is an input and can neither support nor falsify. |
| Standard Model | The particle-physics rulebook for known particles and forces, except gravity. |
| Gauge group | The symmetry bookkeeping behind forces: allowed relabeling inside the force description that leaves physics unchanged. |
| Hypercharge | A Standard Model charge label used to build electric charge. |
| Proton decay | A hypothetical proton breakup. OPH's product-group branch forbids the grand-unified gauge-boson route; so does the Standard Model. |
| \(P\) | The substrate pixel ratio. Working value located under SL-3 from measured α: 1.630968209403959… A counted input. The forward-model closure point \(P_{\mathrm{fwd}}\) = 1.630972095858897… (interval-certified) is a distinct numeral; their difference is the loop residual (CL-1/CL-2). |
| \(N\) | The record capacity. Working value located under SL-4 from measured Λ: about \(3.31\times10^{122}\). A counted input. |
| \(N_{\mathrm{CRC}}^{\mathrm{EW}}\) | The electroweak-bridge capacity, \(3.5323546226929906511\ldots\times10^{122}\). SL-4 admits one N, so the 6.6% gap to the Λ-located value is an open contradiction (CL-3) rather than a second capacity. |
| 24-tick repair lock | One global screen step equals 24 rounds of local repair, \(2(8+3+1)\). The ×2 write/check orientation (12→24) is a discrete structural selection with a recorded menu size, not a parameter-free deduction. |
| \(G_{\mathrm{OPH}}\) | The SI bridge display for Newton's constant, \(6.674299995910528\ldots\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}\). It consumes measured bridge quantities and carries the QCD/cesium-clock caveat. The structural claim is the dimensionless gravity theorem. |
| \(\chi_\nu\) | The coherent scalar susceptibility continuation. Currently a speculative continuation with no testable content; see the appendix. |
| Golden ratio | The number about 1.618 that appears in five-fold symmetry. SL-2 places the pixel off this balance point by one Gaussian-normalized observation width. |

## 1. Closure Tests

[CLOSURE_LEDGER.md](../CLOSURE_LEDGER.md) is the canonical record; this section
mirrors it. SLH requires these equations to hold exactly. Under SL-0 the
simulating structure and the simulated world are one closed system, so each
residual is a direct measurement of the principle set; there is no compare-only
branch for any row, and no row closes by relabeling.

| # | Closure | Required | Current status | Residual (relative) | Residual (measurement σ) |
|---|---|---|---|---|---|
| CL-1 | P-loop, source chain: \(P = \varphi + \sqrt{\pi}/A_T(P)\) with source-only transport | exact | forward map contracts to α⁻¹ = 136.994835177413 (interval-certified unique fixed point, width 7.2×10⁻²⁴), outside the SL-3 basin | 3.0×10⁻⁴ | ≈2×10⁶ |
| CL-2 | P-loop with gauge-width term, self-consistent map | exact | certified unique fixed point α⁻¹ = 137.035660136947 vs 137.035999177(21) | 2.5×10⁻⁶ | ≈1.6×10⁴ |
| CL-3 | One capacity (SL-4): \(N_{\mathrm{EW}} = \pi\exp(6\pi/(P\alpha_U))\) vs N from Λ | exact | 3.53×10¹²² vs 3.31×10¹²² | 6.6×10⁻² | n/a (model-dominated, not Λ error) |
| CL-4 | Hierarchy bridge at the one N: \(\alpha_U\log(N/\pi) = 6\pi/P\) | exact | 11.5546 vs 11.5573 at the Λ-located N | 2.3×10⁻⁴ (log scale) | n/a |
| CL-5 | Forward electroweak emission at the SL-3 pixel: \((M_W, M_Z)\) | exact | (80.330, 91.119) vs (80.3692(133), 91.1880(20)) | 4.9×10⁻⁴ / 7.6×10⁻⁴ | 2.9σ / ≈35σ |
| CL-6 | Printed pair identity: \(\alpha_{\mathrm{root}} = (P_{\mathrm{fwd}} - \varphi)/\sqrt{\pi}\) | exact | closed 2026-07-14: holds to ≥ 35 digits after converged precision-100 reruns; the converged α⁻¹ = 136.994835177413 supersedes the earlier printed tail 136.994835164622 beyond digit 9 | 3.1×10⁻³⁶ (report), 1.2×10⁻³⁸ (trunk) | n/a |
| CL-7 | Capacity readback map: \(F(N) = N\) | exact | F not yet constructed | n/a | n/a |

CL-1 and CL-2 are the same open term seen from two points on the chain: the
Ward-projected hadronic transport. CL-3 supersedes all "two capacities" and
"resonance at logarithmic depth" language: SL-4 admits one N. The previously
displayed 137.0359595 is a mixed-provenance display packet (the inner value
evaluated at P_fwd plus α_U evaluated at the SL-3 pixel); it is not a fixed
point of any single declared map and stays on record as a display packet only,
never a prediction or derivation of α, which is the input that located P. The
certified self-consistent gauge-width fixed point is α⁻¹ = 137.035660136947.

### The Decisive Test: Blind Hadronic Closure

The blind hadronic target is the decisive test of the program, with both
outcomes pre-declared before any payload construction:

- **Frozen target:** `falsification/frozen_targets/hadronic_closure_target_2026-07-14.json`,
  sha256 `7cedad0a7281c74ca0fb1105120c991aeab2f3c45bf86adbbfd560c6324fb985`.
  The target is externally timestamped: OpenTimestamps proofs for the target
  file and its commit record are committed alongside it in
  `falsification/frozen_targets/`, so the freeze date is verifiable outside
  the repository.
- **Required values:** the completed transport must supply
  Δ_source = 0.0414658… α⁻¹ units at the endpoint (equivalently
  S_required = 0.8954001…). The computation must not read the target file,
  CODATA α, the SL-3 pixel, or anything downstream of them; blindness is
  enforced by the dependency-cone audit.
- **Pass:** the blind payload lands inside the interval implied by CODATA α
  uncertainty (2.1×10⁻⁸ on α⁻¹). CL-1/CL-2 close; the α row moves to the
  scorecard outputs column; stage-4 sharpened digits are registered.
- **Fail:** SL-2/SL-3 as formulated are recorded as falsified in the closure
  ledger. Any principle revision is a new frozen registration.

The corpus's own tautological-realization theorem is the reason blindness is
non-negotiable: the payload class can represent any residual, so only a
target-blind construction carries information.

### Structural Kill Conditions (Theory And Computation)

Sixteen outcomes, each executable by an external party against the stated
assumptions:

| Test | Outcome that kills OPH | Why it kills OPH |
| --- | --- | --- |
| Repair confluence | A finite observer-patch network satisfies the OPH descent assumptions, and two accepted repair orders from the same initial quotient state settle into different observer-facing quotient normal forms with no declared holonomy or higher-gauge obstruction. | Descent only proves termination. OPH needs the local-diamond condition plus repair completeness for confluence. Same-boundary uniqueness also needs a preserved boundary/sector map with a unique consistent extension. If the final physical answer depends on repair order under those hypotheses, the proposed repair law fails as an OPH consensus mechanism. |
| Mismatch descent | A repair move accepted by the OPH rules makes the mismatch larger and no compensating potential decreases. | Repair is supposed to settle disagreements. If an allowed move can make the system less consistent, the core dynamics points the wrong way. |
| Relativity emergence | A model satisfies the controlled BW branch certificate, cap-pair extraction, regularized modular transport, support-readable modular covariance, round-cap rigidity, and KMS/BW normalization, then fails to produce the Lorentzian cap/conformal light-cone structure. | OPH claims Lorentz kinematics on the controlled geometric cap-pair branch. A countermodel satisfying that certificate while lacking the light-cone structure breaks the branch. |
| Einstein branch | The entropy and modular assumptions hold, the controlled Markov-collar and controlled modular remainders vanish in the stated limits, and the Einstein-gravity limit fails. | OPH claims ordinary gravity appears as a controlled scaling branch of the overlap rules. If the stated assumptions do not give that branch, the gravity reconstruction fails. |
| Gauge reconstruction | A cofinal tail satisfies the full compact-gauge refinement receipt, yet the zero-obstruction construction fails to produce a well-defined compact group from the stated tensor-generated transportable sector category and forgetful fiber. | This breaks the receipt-conditional classification/reconstruction stage before any Standard Model selection claim can be made. |
| MAR selection | A valid MAR-admissible one-Higgs low-energy sector package, satisfying the same anomaly, refinement, CP, and weak-sector clauses, is selected instead of the Standard Model package. | The selection claim is MAR-local; a different MAR-minimal package would mean the theory did not derive the realized Standard Model branch. |
| Charge lattice | The stated assumptions allow a different hypercharge lattice. | The hypercharge lattice (given the MAR clauses as inputs) is one of the two landed scorecard rows. If OPH allows another lattice, the observed charge pattern was not forced. |
| Generation count | A valid OPH construction has two or four light matter generations. | OPH claims exactly three generations on the realized branch. |
| Fermions and chirality | Chiral fermions have to be imported from the Standard Model by hand. | OPH says matter structure follows from observer consistency. If the key fermion structure has to be assumed, the reconstruction fails. |
| P-lane protocol | The stated fixed-point equations admit no interval-grade contraction certificate on the declared interval, admit multiple solutions there, or the source map reads the low-energy target. | Any "one and only P" claim requires the stage-2 certificate. The lane now holds one: `code/P_derivation/runtime/p_interval_contraction_certificate_2026-07-14.json` proves existence and uniqueness of the closure fixed point on its stated interval (direct Banach, mpmath.iv outward rounding, L ≤ 0.0724, edge-sum tails bounded), for both readout maps at the declared cutoffs and conventions. The global at-most-one statement is discharged on the declared physical domain (alpha_inv in [100, 200]): sup|g'| < 1 on all 256 certified pieces, both maps, zero exceptional set (`runtime/p_global_uniqueness_certificate_2026-07-14.json`). A proof of non-uniqueness on the declared domain or a target leak kills the lane rather than merely stalling it. |
| N-lane readback | The capacity readback map F is constructed and \(F(N)\neq N\) at the Λ-located N, or the same equations allow incompatible capacity branches. | The capacity closure is the second substrate constant. CL-7 is open; a constructed F that fails is a falsification, not an open item. |
| 24-tick repair lock | The observer-visible product-adjoint repair spectrum is not \(2(8+3+1)=24\), or the local repair contraction fails to line up with the global screen step. | The repair lock is the bridge between P and N. Its ×2 orientation is a recorded discrete selection; the spectrum content is not, and a different spectrum breaks the mechanism. |
| Hierarchy witness | The declared first-principles hierarchy map imports \(G\), \(\Lambda\), \(W\), \(Z\), Higgs, low-energy Thomson, cesium-clock, or hadronic inputs, or it misses \(v/E_\star\) on its declared branch: 2.0199803239725553e-17 at the SL-3 working pixel, 2.0198114078576331e-17 at the certified P_fwd. | This row avoids QCD and clock assumptions, which is why it is the compact proof's answer to the QCD-circularity criticism. It is conditional on the frozen MSSM-like one-loop coefficients \((33/5, 1, -3)\), a recorded structural selection; every α_U descendant, including \(\epsilon_H = 0\), inherits that condition. |
| Dimensionless gravity theorem | The dimensionless gravity relation fails from the stated assumptions. | OPH claims the gravity coupling is a branch readout. The SI display \(G_{\mathrm{OPH}}\) is a bridge display consuming measured quantities and cannot carry this test; the dimensionless theorem can. |
| Low-energy uniqueness | The OPH axioms plus MAR allow many inequivalent low-energy physics worlds with the same admissibility predicates. | If many incompatible worlds pass the same rules, the selection claim fails. |
| Quark source-spread theorem | A target-free OPH observable breaks the claimed free \((\mathbb R_{>0})^2\) action using only the stated source corpus. | That result would falsify the present non-identifiability theorem and reopen numeric quark reconstruction. |

### Discriminating Physical Kill Conditions

Seven outcomes on which OPH and the Standard Model (or its standard
extensions) part company:

| Test | Outcome that kills OPH | Why it kills OPH |
| --- | --- | --- |
| Blind hadronic closure | The blind Ward-projected payload lands outside the frozen-target interval. | SL-2/SL-3 as formulated are falsified. This is the decisive test above. |
| Maxwell carrier | With the Maxwell kinetic term, ordinary Lorentz vacuum, and no Higgs, Stückelberg, or medium mass all fixed, the reduced quadratic theory lacks its two transverse \(k^2=0\) modes. | This falsifies the scoped classical Maxwell carrier theorem. A photon-particle mass test applies only after a physical-Hilbert-space and positive-residue pole receipt is also supplied. |
| Color carrier | With the pure Yang–Mills kinetic term about the trivial flat connection and the perturbative/deconfined phase fixed, the reduced quadratic theory lacks its \(2\dim G\) transverse \(k^2=0\) modes. | This falsifies the scoped quadratic Yang–Mills theorem. Confinement can remove a free physical gluon pole and is not a counterexample to that theorem. |
| Einstein tensor carrier | With the pure two-derivative Einstein–Hilbert action, flat vacuum, and no extra fields or higher-curvature terms fixed, the reduced quadratic theory lacks its two transverse-traceless null modes. | This falsifies the scoped Einstein linearization theorem. A graviton-particle mass test applies only after quantization and a positive-residue physical pole have been established. |
| de Sitter capacity | Cosmological data break the claimed de Sitter capacity relation. | SL-4 reads N from Λ through this relation. If the capacity relation fails, the located working value and the screen bookkeeping fail with it. |
| Black-hole entropy | Black-hole entropy data contradict the OPH screen-counting rule. | OPH treats screens and records as physical bookkeeping. Wrong black-hole entropy breaks that bookkeeping. |
| Page curve | Quantum black-hole information data contradict the OPH record claim. | OPH says records and overlap consistency control what can be recovered. A wrong Page-curve pattern breaks that claim. |

## 2. Frozen Predictions

**The registry is currently near-empty.** Zero frozen, externally timestamped
prospective physical-prediction receipts exist. The one frozen artifact is the
hadronic closure target above, which is a target for a computation, not a
prediction for future data. Every claimed hit to date is a retrodiction or a
basin location; the four internal audits agree on this. Stating it plainly is
the precondition for fixing it.

Queued for the registry (each entry freezes quantity, exact formula, value with
uncertainty, dependency-cone hash, and external timestamp before the comparison
data exists or is looked at):

1. **Stage-4 sharpened P digits.** After an interval-grade contraction
   certificate exists for the completed P-loop map and its certified width
   beats CODATA's ±2.1×10⁻⁸ on α⁻¹, the sharpened digits are registered as a
   prediction for future metrology.
2. **Forward W/Z.** After the forward electroweak chain is repaired with
   preregistered revisions (declared suspects: the frozen MSSM one-loop
   coefficients, QT1–QT5, β_EW), the re-emitted \((M_W, M_Z)\) pair is frozen
   before comparison. The current forward pair is a live tension (CL-5), and
   the exact-match W/Z companion lane is calibration: it consumed its targets
   and can neither support nor falsify.
3. **DESI stance rows.** The fixed-N reading commits to
   \((w_0, w_a) = (-1, 0)\); a stance on the S8/growth row (currently
   disfavored, frozen anyway); any birefringence refinement carries its
   trials count (the existing 0.37501° expression was the 4th of 4 tried).
4. **Prospective neutrino branch.** A source-derived family kernel, cycle law,
   orientation, and physical charged-lepton basis, hash-locked without
   oscillation-target feedback, then evaluated once. OPH currently emits no
   neutrino prediction.
5. **Dark-sector joint run.** A source-eligible parent with frozen source map,
   applicability rule, nuisance policy, covariance treatment, and likelihood,
   evaluated against held-out CMB, BAO, SPARC, cluster, and weak-lensing data.

Retrospective comparisons that are not predictions and are not counted as
evidence: fixed-\(\mathbb{Z}_6\) RAR scatter 0.1328335 dex and error-aware BTFR
slope 3.8456543 ± 0.0858172 (slope 4 at 1.80σ), where the interpolation
function is the McGaugh–Lelli–Schombert empirical RAR fit and the SPARC
samples are not held out; analytic scalar tilt
\(n_s = 1 - P/48 = 0.9660214956\) at 0.267σ from the Planck summary with an
83-bin CAMB diagnostic of Δχ² = +0.83, where two formula candidates circulate
and the trials count is part of any future registration.

## 3. Live Tensions

Open quantitative conflicts, each attributed to the model under the
measurement-attribution rule (section 6):

| Tension | Values | Relative deviation | Pull (measurement σ) |
| --- | --- | --- | --- |
| One-N gap (CL-3) | \(N_{\mathrm{CRC}}^{\mathrm{EW}} = 3.53\times10^{122}\) vs Λ-located \(3.31\times10^{122}\) | 6.6×10⁻² | model-dominated; Λ error subdominant |
| Hierarchy bridge (CL-4) | 11.5546 vs 11.5573 at the Λ-located N | 2.3×10⁻⁴ (log scale) | same object as CL-3 |
| Forward W (CL-5) | 80.330 vs 80.3692(133) GeV | 4.9×10⁻⁴ | 2.9σ |
| Forward Z (CL-5) | 91.119 vs 91.1880(20) GeV | 7.6×10⁻⁴ | ≈35σ |
| S8 growth row | out-of-loop row of the dark-sector likelihood surface | n/a | +2.4σ against the model |
| Source-α gap (CL-1) | forward chain contracts to α⁻¹ = 136.994835177413 (certified) vs 137.035999177(21) | 3.0×10⁻⁴ | ≈2×10⁶σ |
| Gauge-width α gap (CL-2) | certified fixed point 137.035660136947 vs 137.035999177(21) | 2.5×10⁻⁶ | ≈1.6×10⁴σ |
| Bench bound on \(\chi_\nu\) | coupons at rest bound gravitational ΔS ≲ 10⁻¹⁰ | record-level contrasts are 10⁻²–1; the 1:1 record↔gravity identity is ruled out by 8–10 orders | n/a (null bound) |

Notes:

- SL-4 forbids two capacities. The one-N gap is carried as an open
  contradiction until the EW-bridge/capacity connection is derived; it is not
  resolved by naming two numbers.
- The forward W/Z pulls are the display of record for the electroweak lane. The
  exact-match pair is calibration and appears in no evidence count.
- The S8 row is the one out-of-loop comparison in the dark-sector surface; the
  in-loop rows (Ω_m closure against Planck inputs) are arithmetic identities
  and carry no evidential weight.
- The bench bound derives from the corpus's own collaborator analysis
  (`hoverboard-experimental-precursor/reply_to_bmu.md`): with
  \(F \approx 5.14\times10^{8}\cdot\Delta S\) N, a healthy-device record
  contrast of 10⁻²–1 would produce 10⁷–10⁸ N on an 80×60 mm coupon at rest.
  No such force exists, so the gravitational ΔS is bounded at ≲10⁻¹⁰. This is
  the strongest existing \(\chi_\nu\) constraint and it binds every future
  version of the force law before new hardware is built.

## 4. Executed Tests, Permanent Record

Every test ever run stays on this surface with its verdict. A failed row is
revised only by a new frozen registration, never by reclassification
([STRANGE_LOOP_PRINCIPLES.md](../STRANGE_LOOP_PRINCIPLES.md), rule 7).

| Test | Target and data | Result (relative deviation; pull) | Verdict |
| --- | --- | --- | --- |
| Neutrino weighted-cycle point | \(\sin^2\theta_{23} = 0.5820560367\), \(\delta_{CP} = -54.419^\circ\) (≡ 305.581°) vs NuFIT 6.1 correlated profiles | Δχ² = 20.12 (with tabulated atmospheric likelihood) and 18.44 (without); both exceed the two-parameter 3σ value 11.83 | **REJECTED.** The candidate was target-informed (exponent law ranked against the measured splitting ratio). The displayed mass tuple (0.01745, 0.01948, 0.05308 eV), sum 0.0900 eV, and Majorana readouts are coordinates on the rejected candidate. |
| Universal static dark law | fixed-\(\mathbb{Z}_6\) equation applied unchanged to the Milky-Way–Sun system vs the Cassini quadrupole summary | predicted \(Q_2 = 3.6201778\times10^{-26}\,\mathrm{s^{-2}}\), ≈19× the summary 1σ uncertainty above it; pull 19.2232σ (unit endpoint \(3.4021876\times10^{-26}\): 18.0122σ) | **EXCLUDED** as a universal law. The scoped old-settled-galaxy branch survives only with a source-derived applicability or screening rule, which does not yet exist. |
| IBM Z₅ high-precision run | five-state golden-ratio target \(\varphi^2 = 2.618033988749895\) | measured 2.5498, CI [2.4919, 2.6023]; relative deviation 2.6%; the target sits above the upper bound (≈2.4σ at the interval-implied width) | **The run excludes the target.** This is the one reconciled verdict; it replaces the pair of documents that called the same run confirmation and disproof. Standard QM predicts the ideal programmed target, so the run also carries no OPH-support weight; the exclusion stands as executed. |
| IBM S₃ six-state ratio test | preregistered ratio 2 | initial runs 1.70–1.93; relative deviation −3.5% to −15% | **FAILED** as preregistered. The reversed-layout reading 2.0348 (n = 1 configuration, found after the failures) is a new hypothesis requiring a fresh frozen run, not a rescue. |
| Bottom-quark normalization | \(m_b(m_b)\) forward emission vs PDG 4.186(6) GeV | emitted 6.03 GeV; relative deviation +44%; pull ≈3×10²σ | **OVERSHOOT stands.** Carried openly as the one-loop b–τ register-selection and threshold-packet gap. |
| Quark denominator-tuple lanes | six quark masses from the RSCC tuple grammar | tuple was the minimum of a 219,615-member post-hoc grammar; an ablation fit better | **SELF-REJECTED** by the quarks ledger (status: post-hoc, target-informed). The "all six within 11.5%" positives are void; no surface may cite them. |
| Charged-lepton phase rigidity | \(m_\mu/m_e\) from the \(\delta_e = 2/9\) phase | +9.8336 ppm; pull +442σ | **EXCLUDED as a prediction.** The phase is the Brannen/Koide fitted value; "Koide follows exactly" is an identity of the ansatz. The lane displays (deviation, pull) pairs wherever it appears. |
| Ringdown comb | α = 4 comb prediction | −5.9σ | **FAILED.** The α_eff = 16 variant is a new candidate for the frozen registry, scored separately; it does not repair this row. |

The remaining IBM circuit-family archive (recovery/CMI fingerprints, three-state
sanity, decoy and relabeling controls, bundle readings 2.5974/2.7383/2.5786 for
Z₅ and 1.872/2.030 for the six-state layouts) is retained as an engineering
record. Standard quantum mechanics predicts every programmed target in that
family, so those runs are hardware benchmarks on both sides of the ledger: they
provide no OPH support and, apart from the executed exclusions recorded above,
no OPH falsification. Any future quantum-hardware test must state, before
running, an observable where OPH and standard QM predict different values.

## 5. Null Checks Shared With The Standard Model

Ten outcomes that would kill OPH's realized branch but are equally forbidden in
the Standard Model (or standard GR). They are real constraints and
**non-discriminating**: their continued absence supports OPH exactly as much as
it supports the Standard Model, which is to say not at all.

| Null check | Outcome that kills the OPH branch |
| --- | --- |
| Gauge-mediated proton decay | A proton decay channel mediated by simple-grand-unified gauge bosons is observed. |
| Grand-unified leptoquark bosons | Simple-GUT leptoquark gauge bosons are discovered. |
| Fourth light generation | A fourth light chiral matter generation is discovered. |
| Fourth color | A real fourth QCD color degree of freedom is found. |
| Hypercharge outlier | A real elementary particle has charges outside the hypercharge lattice. |
| Fractional color singlet | A stable color-singlet particle has fractional electric charge. |
| Extra light Higgs sector | Additional light Higgs multiplets are discovered in the minimal regime. |
| Low-energy supersymmetry | A low-energy superpartner spectrum is discovered. |
| Lorentz violation | High-confidence low-energy Lorentz violation is measured where OPH requires Lorentz symmetry. |
| Equivalence principle | High-confidence EP violation is measured in the Einstein regime. |

Calibration rows do not appear anywhere on this surface as tests: the
exact-match W/Z pair and the Higgs/top D11 lane consumed their measured targets
(the D11 selectors are target interpolation, admitted in code) and are inputs.

## 6. Measurement-Attribution Rule

Closure residuals are attributed to the model. The position "any deviation must
come from additional influences on the measured value" is quantitatively
excluded: independent α determinations with different physics (electron g−2
with QED, and Rb/Cs photon-recoil interferometry) agree at the few-×10⁻¹⁰
relative level, while the loop gaps are 2.5×10⁻⁶ (gauge-width chain) and 3×10⁻⁴
(source chain): 10³–10⁶ times larger than the room available on the measurement
side. The corpus's own hadron-scheme bridge places the missing piece in the
model's hadronic transport. A measurement-side mechanism claim is admissible
only as a frozen registered prediction ([STRANGE_LOOP_PRINCIPLES.md](../STRANGE_LOOP_PRINCIPLES.md), rule 5).

## Appendix: Speculative Continuations

Rows here have no current testable content. They are not on the falsification
surface and no outcome concerning them bears on the core until the stated gaps
close.

### \(\chi_\nu\) Lift Response And The Hoverboard Lane

Recorded values: \(\chi_\nu^{\mathrm{can}} = \exp(-P/24) =
0.9343006394893864\ldots\); engineering chart \(\chi_\nu^{\mathrm{eng}} =
\chi_\nu^{\mathrm{can}}/N_{\mathrm{coh}}\), which for \(N_{\mathrm{coh}} =
10^{22}\) gives \(9.34\times10^{-23}\) to \(1.00\times10^{-22}\);
hoverboard-class areal load Σ = 200–600 kg m⁻² with
\(\Delta S_{\mathrm{coh}}^{\mathrm{can}} \simeq 2\times10^{-8}\) to
\(6\times10^{-8}\) for full response; null bound
\(|\chi_\nu^{\mathrm{eng}}| \le 4\pi G F_{\min}/(g^2 A_\perp|\Delta
S_{\mathrm{coh}}^{\mathrm{eng}}|)\).

**The lane has no testable content until two artifacts exist:** (1) an
operational, instrument-independent recipe for producing and measuring the
coherence contrast \(\Delta S_{\mathrm{coh}}\) on a real substrate; without
it, every null is attributable to not reaching the contrast, and
\(N_{\mathrm{coh}}\) free over 10²⁰–10²⁴ makes the engineering coefficient
anything; and (2) a reconciliation with equivalence-principle and fifth-force
bounds, since the assignment \(c_U = 0\) to ordinary matter currently exempts the
lane from every torsion-balance, Eötvös, and MICROSCOPE constraint by
construction, and that exemption requires a derivation, not a definition.

The decision rules of book chapter 14 (null without receipt → engineering
failure; receipt without lift → bound for that build surface only) are
superseded by one rule: **a null with receipts binds the continuation.** A
declared number of independent builds reaching their receipts with null force
retires the continuation. The bench bound of section 3 already binds: any
future force law must be consistent with gravitational ΔS ≲ 10⁻¹⁰ for coupons
at rest, derived before new hardware is built.

### Bench P-Body

The twelve-port icosahedral coupling target (edge : diagonal : opposite =
1 : 0.381966… : 0.145898…) is an analogue-engineering claim about a built
body, not a substrate measurement. Evidence requires public bundles with
controls; no bundle currently changes any physics row.

![OPH P hardware confirmation schematic](assets/oph_p_hardware_confirmation_animation.svg)

### Hardware Proof-Of-Work

The compute claim is a measured distributional lift \(B = p_Q/p_U\), not a
complexity-class theorem. The one existing scorebook is 512 synthetic
candidates with control ≈ candidate; the enrichment premise has no supporting
data. Engineering claim only.

### Coding-Theory Extensions

A bare OPH overlap net is a finite constraint code. QECC distance, min-cut
resilience, and Knill–Laflamme correction claims require an explicit code
subspace, logical operators, error family, and recovery certificate; only
certified topological-code branches may state them.

## Bottom Line

OPH is a principle set (SL-0…SL-5) with zero dials, zero quantitative inputs at
the theory layer, two working substrate values located from measurement and
counted in the scorecard, and a closure ledger that either goes to zero or
falsifies the principles as formulated.
The decisive experiment is the blind hadronic closure against the frozen
target; the live tensions are the one-N contradiction, the forward
electroweak pulls, the S8 row, the α-loop gaps, and the bench bound; the
executed failures stay on this page permanently; the SM-shared nulls
discriminate nothing; and the continuations in the appendix earn their way
back onto the falsification surface only by producing operational recipes and
frozen registrations. The theory can lose, and this document is the standing
list of exactly how.
