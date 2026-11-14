# OpenWave Models: the comparison table

## Why multiple models

OpenWave's mission is to build a platform where multiple candidate field-theoretic models are evaluated **numerically, side by side, in the same computational environment**. No single alternative framework can map the space of possibilities on its own: when several independent models are run against the same observables with the same pass/fail criteria, the comparison **triangulates** what is actually out there.

This makes the platform an **open arena** for rigorous, side-by-side numerical verification of candidate models, including unconventional ones the mainstream overlooks, all judged on the same falsifiable criteria. Anyone with a model is invited to put it to the test.

This is the structural difference from a single-model theory-of-everything program: a lone framework carries a built-in incentive to sell itself, while an arena referees. OpenWave has no house model and no stake in which column wins; the platform's claims are the per-criterion cells below and the scripts behind them, nothing more, and a column full of honest negatives is a working column.

Features that survive across frameworks are likely load-bearing physics; features that only work in one framework, or only with hand-tuning, reveal themselves as such. A null result in one model is ambiguous (model wrong, or engine wrong?); a positive result in any model certifies the engine for all of them.

Anybody can contribute to building these numerical validations. Every cell in the tables below is backed by a runnable script or a research document in this repository, and every claim is reproducible, refutable, and extendable under [Apache 2.0](LICENSE). The clean-clone path from any cell to the command that regenerates it is documented in [REPRODUCE.md](REPRODUCE.md).

## COVERAGE MATRIX (Phenomenological Coverage)

Current Models in the platform:

- **[M4 - EWT](openwave/xperiments/m4_ewt/__M4_model_briefing.md)** (Energy Wave Theory, Jeff Yee, built on Milo Wolff and Gabriel LaFreniere pioneer work; Enhanced EWT extension co-authored by Łukasz Smoliński).
- **[M5 - LC](openwave/xperiments/m5_liquid_crystal/__M5_model_briefing.md)** (Liquid-Crystal topological defects, Jarek Duda, with Manfried Faber inputs);
- **[M6 - Ouroboros](openwave/xperiments/m6_ouroboros/__M6_model_briefing.md)** (Chaoiton framework, Paul Werbos);
- **[M7 - HydroBoros](openwave/xperiments/m7_hydroboros/__M7_model_briefing.md)** (Toroidal-Beltrami, Marc Fleury's toroidal electron fused with Paul Werbos's Ouroboros);
- **[M8 - MIT](openwave/xperiments/m8_mit/__M8_model_briefing.md)** (Mode Identity Theory, Blake Shatto top-down model, spectral geometry + representation theory);
- **[M9 - NSM](openwave/xperiments/m9_emergent_gravity/__M9_model_briefing.md)** (New Standard Model, Bob McGwier gravity-certification column: gravity as entanglement bookkeeping on the Standard Model minimally coupled to Einstein-Cartan)
- **[M10 - OPH](openwave/xperiments/m10_oph/__M10_model_briefing.md)** (Observer Patch Holography, Bernhard Mueller: observers as feedback loops on a holographic screen of twelve-port icosahedral carriers, with a local repair and memory law and no action)

Every file reference is an active link to the file in this repository (under `openwave/xperiments/`). Rows are grouped by domain: particles, forces, waves + quantum emergence.

### Summary Count

**Column order:** models are sequenced by their validated + partial count (✅ + ⚠️), highest first; ties break toward more ✅ (validated), then fewer ❌ (honest negatives). The order updates as validations land. A ❌ is a result, not an embarrassment: documented negatives (with the scripts that produced them) are part of the platform's value.

| **MODEL SCORE-BOARD** | [Liquid Crystal<br>(M5)](openwave/xperiments/m5_liquid_crystal/__M5_model_briefing.md) | [HydroBoros<br>(M7)](openwave/xperiments/m7_hydroboros/__M7_model_briefing.md) | [EWT<br>(M4)](openwave/xperiments/m4_ewt/__M4_model_briefing.md) | [Ouroboros<br>(M6)](openwave/xperiments/m6_ouroboros/__M6_model_briefing.md) | [MIT<br>(M8)](openwave/xperiments/m8_mit/__M8_model_briefing.md) | [NSM<br>(M9)](openwave/xperiments/m9_emergent_gravity/__M9_model_briefing.md) | [OPH<br>(M10)](openwave/xperiments/m10_oph/__M10_model_briefing.md) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ✅ validated in-platform | 9 | 0 | 0 | 3 | 0 | 0 | 0 |
| ⚠️ partial / with caveats | 11 | 10 | 10 | 2 | 1 | 1 | 0 |
| ❌ honest negative | 2 | 0 | 3 | 3 | 0 | 0 | 0 |
| 🚧 planned / not tested | 9 | 21 | 18 | 23 | 30 | 30 | 31 |
| **Total criteria** | **31** | **31** | **31** | **31** | **31** | **31** | **31** |

### Summary Status

Each icon is earned by the matching row in that model's own table (column headers jump there); the per-model row carries the summary and the evidence links.

The **regime** column is a property of the criterion, not of any model: what a row demands of whatever attempts it. `static` = earnable from structure or spectrum alone, `dynamic` = not, `both` = the row has both sub-questions, so a model can earn part of it without dynamics.

Magnetism carries no standalone row by design: the magnetic force of currents is a corollary of Electric force + Lorentz covariance ([Feynman II-13-6](https://www.feynmanlectures.caltech.edu/II_13.html#Ch13-S6): a boosted Coulomb field IS the magnetic field); magnetostatics with no current (the permanent magnet's static B) is the electron's intrinsic dipole, scored by the Magnetic moment μ row in the electron block.

Each criterion's simplest passing test sits in its own companion table right below ([§ Simplest Test per Criterion](#simplest-test-per-criterion)), keeping this matrix to icons.

| Criteria | [M5](#liquid-crystal-m5) | [M7](#hydroboros-m7) | [M4](#ewt-m4) | [M6](#ouroboros-m6) | [M8](#mit-m8) | [M9](#nsm-m9) | [M10](#oph-m10) | regime |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **PARTICLE - ELECTRON** | | | | | | | | |
| Charge quantization | ✅ | ⚠️ | ❌ | ⚠️ | 🚧 | 🚧 | 🚧 | static |
| Electron rest energy (mass) | ✅ | ⚠️ | ⚠️ | ❌ | 🚧 | 🚧 | 🚧 | static |
| de Broglie clock (Zitterbewegung) | ✅ | ⚠️ | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | dynamic |
| Particle stability (Derrick escape) | ✅ | ⚠️ | ⚠️ | ✅ | 🚧 | 🚧 | 🚧 | dynamic |
| Magnetic moment μ (g ≈ 2, magnetostatics) | ❌ | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Angular momentum J (spin ℏ/2) | ⚠️ | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Spin-½ statistics (720° double cover) | ✅ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| | | | | | | | | |
| **PARTICLE - FULL SPECTRUM** | | | | | | | | |
| Antimatter + annihilation | ✅ | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | both |
| Neutrinos (neutral states) | ⚠️ | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | both |
| Neutrino oscillations (PMNS) | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | both |
| Lepton mass spectrum (μ, τ) | ⚠️ | 🚧 | ❌ | ❌ | 🚧 | 🚧 | 🚧 | static |
| Quarks | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Baryons: bound state (p, n) | ⚠️ | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Baryons: mass ordering + charge profile | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Baryons: exact masses | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Deuteron (binding + quadrupole) | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Nuclear structure (levels, halos) | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | both |
| Mesons (π, K) | 🚧 | 🚧 | 🚧 | ❌ | 🚧 | 🚧 | 🚧 | static |
| Dark matter candidate | 🚧 | 🚧 | 🚧 | ✅ | 🚧 | 🚧 | 🚧 | static |
| | | | | | | | | |
| **FORCES** | | | | | | | | |
| Electric force (Coulomb 1/r, electrostatics) | ✅ | ⚠️ | ❌ | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Lorentz covariance (+ Coulomb = EM) | ⚠️ | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Strong force: confinement | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Running coupling | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | static |
| Weak force: muon decay | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Weak force: beta decay (n → p) | ❌ | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Gravity: Newton limit (GEM) | ⚠️ | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | both |
| Gravity: local metric phenomena | 🚧 | 🚧 | ⚠️ | 🚧 | 🚧 | ⚠️ | 🚧 | both |
| Cosmology: Λ / cosmic acceleration | 🚧 | 🚧 | 🚧 | 🚧 | ⚠️ | 🚧 | 🚧 | static |
| | | | | | | | | |
| **WAVES + QUANTUM EMERGENCE** | | | | | | | | |
| EM waves (Maxwell) | ✅ | ⚠️ | ⚠️ | ✅ | 🚧 | 🚧 | 🚧 | dynamic |
| Quantum wave equation (Klein-Gordon) | ✅ | ⚠️ | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | dynamic |
| Orbital quantization (atomic structure) | 🚧 | 🚧 | ⚠️ | 🚧 | 🚧 | 🚧 | 🚧 | static |

### Simplest Test per Criterion

The cheapest concrete experiment that would earn each row, so every claim and every 🚧 alike names what would settle it. Like `regime`, the test is a property of the criterion, not of any model. The priority ordering (which tests are most decisive to run first within each domain) and its provenance live in the platform task record ([`dev_docs/tasks/t1_task_details.md`](dev_docs/tasks/t1_task_details.md)).

| Criteria | simplest test |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | integer topological charge; a split e → e/2 + e/2 impossible |
| Electron rest energy (mass) | localized state anchored to 511 keV |
| de Broglie clock (Zitterbewegung) | self-starting frequency-rigid internal oscillation at ω = mc²/ℏ |
| Particle stability (Derrick escape) | the particle state persists in free evolution |
| Magnetic moment μ (g ≈ 2, magnetostatics) | intrinsic magnetostatics: a stationary state sources a static dipole B from internal circulation, g ≈ 2; dipole Larmor response |
| Angular momentum J (spin ℏ/2) | field-carried J = ℏ/2 on the electron state |
| Spin-½ statistics (720° double cover) | the field returns under a π rotation while its frame needs 2π |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | pair annihilation to radiation; follow-on: positronium 2γ/3γ |
| Neutrinos (neutral states) | light charge-0 states released in decays |
| Neutrino oscillations (PMNS) | 3 flavors + oscillation; the PMNS angles |
| Lepton mass spectrum (μ, τ) | exactly 3 charged-lepton minima; ratios 1 : 207 : 3477 |
| Quarks | fractional charge as partial winding on a quark string |
| Baryons: bound state (p, n) | a stable three-quark composite |
| Baryons: mass ordering + charge profile | m_n > m_p; positive core, negative shell |
| Baryons: exact masses | 938.3 / 939.6 MeV from the field |
| Deuteron (binding + quadrupole) | m_d < m_p + m_n; the electric quadrupole moment |
| Nuclear structure (levels, halos) | binding curve, levels, halo-nuclei lifetimes |
| Mesons (π, K) | π/K states; strange decays (Λ⁰ → p + π⁻) |
| Dark matter candidate | a stable neutral massive state |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | far-field 1/r potential at charge e; the strength (α) credited only when it follows from the model's own mechanism, never fitted |
| Lorentz covariance (+ Coulomb = EM) | Lorentz-invariant action, no preferred frame; boosted states contract, c-limited |
| Strong force: confinement | linear inter-charge potential, ~1 GeV/fm (Cornell) |
| Running coupling | coupling varies with scale (onset at the core radius) |
| Weak force: muon decay | μ relaxes to e + neutral ejecta |
| Weak force: beta decay (n → p) | n → p + e + ν̄, parity-violating (needs a neutron) |
| Gravity: Newton limit (GEM) | attractive 1/r² between masses, via the GEM route; the strength (G) credited only when it follows from the model's own mechanism, never fitted |
| Gravity: local metric phenomena | light bending, gravitational time dilation / Shapiro delay |
| Cosmology: Λ / cosmic acceleration | Λ from the model's own mechanism, never fitted; every remaining gap and a falsification named |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | Maxwell recovered; transverse waves at c |
| Quantum wave equation (Klein-Gordon) | emergent ω² = k² + m² for matter waves |
| Orbital quantization (atomic structure) | hydrogen-like discrete levels from standing waves |

## RESULTS BY MODEL

One table per model, column order of the matrix. Each row is the evidence behind that model's icon above: status tag, condensed summary (the 65-word rule), links to the record.

**Cell format (the 65-word rule).** Every summary cell below is a condensed summary, never a report: a status tag, **at most 65 words** of prose, then links; detail belongs in the linked record. The count excludes the tag, links and `<br>→` tails, so ⚠️ this 65 is **not** the 65 of [`dev_docs/ROADMAP_STANDARDS.md`](dev_docs/ROADMAP_STANDARDS.md): it permits about a third more, which the wider column earns. Re-measure before moving it: `python3 dev_docs/utils/models_cell_stats.py` ([T3](dev_docs/tasks/t3_task_details.md)).

**Which tables carry a budget.** Only the per-model ones.

| Table | Budget |
| --- | --- |
| Per-model tables in this section | 65 words of prose per summary cell |
| At-a-glance matrix ([Summary Status](#summary-status)) | none: status icons, and no prose at all, which the linter enforces |
| [Simplest test per criterion](#simplest-test-per-criterion) | none needed: one short clause per row, observed maximum 12 words |
| [Score-board](#summary-count) | none: counts only, each tallied against the rows |

Each matrix icon must match the status tag of the same criterion in the model's own table. That rule, the score-board tallies, the `regime` values, the simplest-test table (non-empty per criterion, criteria-synced with the matrix), every row's column count, and every criteria count or column tally stated in a live doc elsewhere in the repository are linted by [`dev_docs/utils/check_models_md.py`](dev_docs/utils/check_models_md.py), which a maintainer runs before merging anything that touches this file or restates its counts ([`dev_docs/PR_REVIEW_STANDARDS.md § 7.1`](dev_docs/PR_REVIEW_STANDARDS.md#71-the-modelsmd-linter)).

## Liquid Crystal (M5)

Deep dive: [`m5_summary_report.md`](openwave/xperiments/m5_liquid_crystal/research/archive/m5_summary_report.md) (results-of-record) · [`m5_roadmap.md`](openwave/xperiments/m5_liquid_crystal/research/m5_roadmap.md) (full program) · [`m5_question_tracker.md`](openwave/xperiments/m5_liquid_crystal/research/m5_question_tracker.md) (emergence catalog + open questions) · [model briefing](openwave/xperiments/m5_liquid_crystal/__M5_model_briefing.md)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | ✅ [validated in-platform]<br>Topological winding number of the hedgehog defect (Gauss-Bonnet integer Q = ±1)<br>[`m5_1_winding.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_1_winding.py) |
| Electron rest energy (mass) | ✅ [validated in-platform]<br>Hedgehog rest energy with Faber core regularization; mass pinned E ∝ 1/r₀, physical knob r₀ = 2.2132 fm → 0.511 MeV. The r₀ scale is the Coulomb-plus-m_e anchor and gives an energy-median radius of ~2.9 fm, in tension with the scattering bounds on electron structure (below 10⁻³ fm); re-anchoring the mass to the clock is open<br>[`m5_6_3b_faber_on_M.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_3b_faber_on_M.py), [`m5_16_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_16_task_details.md) |
| de Broglie clock (Zitterbewegung) | ✅ [validated in-platform]<br>Bounded, self-starting, frequency-rigid 3+1D time crystal, the energy-minimizing state (≈ 21% below clock-stopped). Rigidity MEASURED under drive: a background scalar cannot entrain it, eigenframe torque zero (M5.27). **Absolute scale ✅:** geo-mean calibration recovers each lepton's ZBW to ~13%. The exterior twist's inertia is box-extensive on the uniaxial vacuum and reverses sign under the first dressing term; the confining mechanism is open<br>[`m5_9_lepton_mass_clock_findings.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_9_lepton_mass_clock_findings.md), [`m5_27_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_27_note.md), [`m5_21_14_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_14_note.md), [`m5_34_kit_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_kit_checks.py)<br>→ open work: [M5.10](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_10_task_details.md) · calibration #208 / #217 / #218 / #220 ✅ |
| Particle stability (Derrick escape) | ✅ [validated in-platform]<br>Escape via time-periodic resonance: static solitons confirmed impossible (M5.2 negative), the saturated breather self-starts from exact rest and holds resolution-robustly (24³ → 48³); the verified-L electron persists in free dynamics at N = 64 (absorbed 4.4%, M5.21.10)<br>[`m5_8_2g_spontaneity.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2g_spontaneity.py), [`m5_21_10_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_10_note.md) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | ❌ [honest negative, RE-BASED at verified-L 2026-07-21]<br>The μ channel exists (envelope-localized, pure twist EM-silent) but the observable is a parity-cancellation residue tracking preparation basin across 4 orders; under the Coulomb-anchored bridge g spans 8.5e-4 to 1.45, no closure at 2.0023. On the director's EM reading E·B vanishes identically and a centered charge has no dipole, so the source-current moment is zero: the readout is the displacement term alone, resolution-dependent<br>[`m5_21_5_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_5_note.md), [`m5_8_2r_electron_id.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2r_electron_id.py), [`m5_34_kit_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_kit_checks.py), [`m5_question_tracker.md`](openwave/xperiments/m5_liquid_crystal/research/m5_question_tracker.md) § FORCES<br>→ closure route: μ on the long-evolution texture + the Larmor read, [`m5_23_1_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_23_1_task_details.md) |
| Angular momentum J (spin ℏ/2) | ⚠️ [partially validated: the state half, RE-BASED at verified-L 2026-07-21]<br>The fixed-J isorotation electron exists and HOLDS at three J rungs, clock thermodynamics exact (dE/dJ = ω\* at ~1%), now ported to the production engines (M5.23.1, selftests green). J is constraint-carried: rigid rotation and free descent both measured out, on both stacks. The ℏ/2 OBSERVABLE half stays open (the Larmor read is instrument-limited)<br>[`m5_21_9_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_9_note.md), [`m5_23_1_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_23_1_task_details.md) |
| Spin-½ statistics (720° double cover) | ✅ [validated in-platform]<br>The field is apolar (ellipsoids): a π rotation returns M while the frame needs 2π, one frame revolution = two field periods, no belt-trick needed. Machine-exact on the production seed (1e-16), the same factor 2 as G7's ω_M = 2ω_clock. Deeper invariant: biaxial π₁ = Q₈ (NG-9). The exchange (Finkelstein-Rubinstein) sign of two cores is not measured; Q₈ alone does not fix it ([Q34](openwave/xperiments/m5_liquid_crystal/research/m5_question_tracker.md#q34-detail))<br>[`m5_8_2s_spin_half_apolar.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2s_spin_half_apolar.py), [Wolfram animations](https://community.wolfram.com/groups/-/m/t/3398814) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | ✅ [validated in-platform, 2026-07-21]<br>Full 3D capture-to-annihilation on the verified-L stack (audited twice): charges exactly ±1 at seed, a clean sigmoid to ≤ 0.005 residual, energy ledger physical to 0.1%. Mechanism MEASURED: conduit annihilation through the connecting topological line (ballistic core-walk refuted). Caveat: the endpoint retains seed-inherited line flux<br>[`m5_21_4_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_4_note.md), [`m5_21_4_a_pair.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_21_4_a_pair.py); lineage [`m5_8_2v_pair_annihilation_budget.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2v_pair_annihilation_budget.py), [`m5_14_sine_gordon_annihilation.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_14_sine_gordon_annihilation.py) |
| Neutrinos (neutral states) | ⚠️ [partially validated]<br>Charge-0 loop states are the neutrino candidate, and the decay side is measured: kicked heavy-lepton minima release paired symmetric ejecta under damped evolution (real, not numerical); whether the ejecta ARE the loop states awaits the tracer (M5.23.2)<br>[`m5_21_10_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_10_note.md), [`m5_21_6_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_6_note.md) |
| Neutrino oscillations (PMNS) | ⚠️ [partially validated]<br>Flavour oscillation = SO(3) spatial field rotation. **PMNS from SO(3) (#199 ✅):** tri-bimaximal + δ_CP = 180° parameter-free (θ₁₂ 35.26° vs NuFIT 33.7°, δ_CP ≈ 177°); θ₁₃ ≈ 8.5° = the SO(3)-breaking. FALSIFIER: δ_CP far from 180°. Oscillation dynamics not yet run<br>[`m5_11_pmns_findings.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_11_pmns_findings.md) · [`m5_11_theta13_findings.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_11_theta13_findings.md)<br>→ #199 ✅ 2026-06-18; follow-up = the charged-lepton matrix #200 |
| Lepton mass spectrum (μ, τ) | ⚠️ [partial: selection + law + decay dynamics; masses accommodated by the article's anchor, native origin open]<br>Three leptons as energy minima for elementary charge; mass law **`E ∝ Λ³`** (#200); decay dynamics measured (M5.21). Hierarchy origin open; the in-platform extrapolation bridge failed (M8.6). The article anchors the masses in the Dirac-core coupling: electron and muon fix a circulant, Koide imposed, the tau held out at +0.43σ (Koide's two-mass relation, reproduced in-platform). S1 (R31): informative null on classical rest energies<br>[`m5_9_lepton_mass_clock_findings.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_9_lepton_mass_clock_findings.md), [`m5_21_6_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_6_note.md), [`m5_21_11_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_21_11_task_details.md), [`m5_32_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_32_task_details.md), [`m5_32_r32_anchor_check.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_32_r32_anchor_check.py), [Release 2.3](https://doi.org/10.5281/zenodo.24051900) |
| Quarks | 🚧 [not yet tested]<br>Fractional charge from a fraction-of-π field rotation on a 1D topological quark string (full π = the elementary charge), enforced in baryons by inter-string interactions; the Cornell linear term (~1 GeV/fm) = the cost of violating quantization. Lives with the baryon program (the census proton-analog's fractional structure) + the string-tension instrument (Q38); SU(3)/CKM open<br>[`m5_22_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_note.md), [`m5_question_tracker.md`](openwave/xperiments/m5_liquid_crystal/research/m5_question_tracker.md)<br>→ [#200](https://github.com/openwave-labs/openwave/issues/200) (lineage, quark detail) · #199 ✅ (the neutrino SO(3) side, resolved 2026-06-18) |
| Baryons: bound state (p, n) | ⚠️ [partially validated]<br>Baryon-analog states MEASURED on the certified stack (M5.22 census, audited): proton-analog = vortex column + equatorial ring, \|Q\| = 1, stationary, perturbation-stable; neutron-analog = column + two rings, Q = 0, scale-stable, reread as candidate dineutron under the ring-count frame (kick-apart probes: returns whole; the physical dineutron is unbound). Toy parameters; three-quark statistics and quantization open (lab anchor: Nature Physics s41567-025-03107-0)<br>[`m5_22_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_note.md), [`m5_22_1_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_1_note.md), [`m5_22_4_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_4_note.md)<br>→ the Skyrme-side dictionary of the column-plus-ring picture (relaxed B = 1 to 4 massive-Skyrme fields: the preimage ring, 0.41 fm at B = 1, and the central line link B times, H = 1.0000 to 4.0000; smooth preimages with trivial line holonomy), the author's scripts rerun here: [`m5_33_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_33_task_details.md); on the v2.6 nucleon B = Hopf(n) = 1 with the charge-axis ring at 0.64 fm, rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Baryons: mass ordering + charge profile | ⚠️ [partially validated]<br>The m_n > m_p DIRECTION measured qualitatively: neutron-analog heavier than proton-analog, energy ratio 1.54, robust across biaxiality (1.55 at δ = 0.2). The charge profile is compound (quadrupolar lobes netting zero), the positive-core / negative-shell signature not yet resolved at current sizes; the quantitative ratio waits on the realistic-parameter bridge<br>[`m5_22_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_note.md) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Absolute p/n masses await the composite stage: the mass-sector anchor exists (r₀ → 0.511 MeV on the electron) but no baryon state has been constructed to weigh. The proposed single-field action (a massive Skyrme sector, not this platform's field) gives 1629 MeV at the physical f_π = 92.4 MeV and 937 MeV at the Adkins-Nappi calibration (f_π 54 MeV, g_A 0.64); the author's scripts rerun here<br>[`m5_15a_composite_particles.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_15a_composite_particles.md), [`m5_33_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_33_task_details.md)<br>→ the author's Lorentz-safe covariant nucleon (∂q·∂q)³ with the ω coupling gated: 15 rows within 9.2 percent, g_A 1.158, τ_n 1030 s, Γ_Δ 62 MeV, the residual pinned by Δ − N and g_A, rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Deuteron (binding + quadrupole) | 🚧 [tested below the resolution bar; no validation claim]<br>A two-knot \|Q\| = 1 candidate with the deuteron's charge and ring count exists, bound 28% at n = 32 but a slowly annealing complex at n = 48 (binding not citable). The quadrupole SIGN reads NEGATIVE at both resolutions on the calibrated div E instrument vs the physical positive (magnitude not citable). The proposed single-field action leaves it unbound too (×1.73 orientation force; author's scripts rerun here)<br>[`m5_22_1_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_1_note.md), [`m5_22_2_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_2_note.md), [`m5_33_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_33_task_details.md)<br>→ the covariant nucleon's deuteron: Q_d 0.24 with a 500 MeV scalar, 0.15 with the model's own 1.4 GeV scalar, against 0.286, the sign set by the pion tensor force, rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Binding systematics, levels and halo nuclei (Boron-8, Ne-17, Li-11 lifetimes; Borromean/Efimov configurations; He-6 as He-4 + nn) are M5.22.9 scope. The single-field no-fit test puts the ⁴He isotriplet 2⁻ at 104 MeV against 23.7; the author's Skyrme-sector soliton merges ⁴He classically (889 MeV against 28.3), the generic Skyrme overbinding whose known repairs (a dynamical ω, vibrational zero-point energy) are unrun on that functional; scripts rerun here<br>[`m5_22_9_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_22_9_task_details.md), [`m5_33_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_33_task_details.md), [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Mesons (π, K) | 🚧 [not yet tested]<br>Pion as twist/reconnection of a vortex loop; kaon as a Möbius-like twisted loop (strangeness = the twist), suggested by their origin in strange-baryon decays; paper-level, lands with the 15a stage<br>[`m5_15a_composite_particles.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_15a_composite_particles.md)<br>→ the sextic of the author's nucleon reads as ω exchange, g²/4π 10.2 (contact) to 14.7 (finite range) against NN 6 to 20; no kaon; rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Dark matter candidate | 🚧 [not yet tested]<br>DM as thermal noise of the non-EM field sectors (the CMBR analog for the other forces), thermalized over cosmic time or sourced by active regions; hopfions remain the particle-like candidate. The trace direction η is exactly invisible to the curvature term (verified here), so a trace mode couples only through the potential, which as written mixes it with the charge axis: no exactly dark state<br>[`m5_4c_convo_2026.06.08.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_4c_convo_2026.06.08.md), [`m5_34_identity_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_identity_checks.py) |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | ✅ [validated in-platform, evidence base restated 2026-07-21]<br>Field-level: the hedgehog's far energy density is the exact Coulomb form 8c₂/r⁴, anchored to exactly charge e, no free factor. Two-body: antipair BINDING measured at every separation, consistent with the derived −64πc₂/d; the 1/d exponent NOT confirmed at reachable boxes (local 1.4-1.7)<br>[`m5_21_5_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_5_note.md), [`m5_21_4_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_4_note.md); lineage [`m5_1_coulomb.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_1_coulomb.py), [`m5_4_coulomb_matrix.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_4_coulomb_matrix.py)<br>→ an explicit two-charge field with pair energy exactly ±e²/d, e² = 64π, on the certified action, rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Lorentz covariance (+ Coulomb = EM) | ⚠️ [partially validated]<br>Lagrangian-level invariance is machine-verified on the full 4×4 tensor theory (M5.18: 15/15 checks, residual 1.3e-11), and separately under internal and spacetime Lorentz transformations (round-off residuals), so the vacuum's time axis is an internal clock, not a frame; the linearized sector is measured covariant: tilt modes propagate at c, emergent Klein-Gordon dispersion. What stays open is a dynamics measurement: boosted defects transforming covariantly in-sim<br>[`m5_18_verification_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_18_verification_note.md), [`m5_18_lorentz_check.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_18_lorentz_check.py), [`m5_34_identity_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_identity_checks.py), [`m5_6_1_kg_operator_check.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_1_kg_operator_check.py) |
| Strong force: confinement | ⚠️ [partially validated]<br>**A linear string term is MEASURED**: like charges string-confined by the required inter-core winding (E_int linear, ~20× Coulomb; FORM robust, tension ansatz-grade), merging into a charge-2 ring compound. The ~1 GeV/fm Cornell anchor is the string-tension instrument's target (Q38)<br>[`m5_21_4_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_4_note.md) |
| Running coupling | ⚠️ [partially validated]<br>Short-range onset verified: the non-abelian ‖R‖·r² roll-off switches on at the core radius r₀, with Maxwell the abelian limit. The resolved scale curve C(ρ) rises from 0.42 (ρ = 0.6) to 0.998 (ρ = 5.5): a classical core form factor approaching the Coulomb plateau, audited as NOT a renormalized coupling or β-function; the coupling dictionary is the open piece<br>[`m5_6_4b_faber_curvature_em.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_4b_faber_curvature_em.py), [`m5_31_coupling_curvature_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_31_coupling_curvature_note.md) |
| Weak force: muon decay | ⚠️ [partially validated]<br>**Decay dynamics measured**, no longer only sketched: kicked heavy-lepton minima relax to the electron level and release structure (neutrino candidates) under damped evolution, energy ledger closed. Missing: any rate or coupling<br>[`m5_21_6_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_21_6_note.md) |
| Weak force: beta decay (n → p) | ❌ [honest negative at the 3×3 truncation]<br>The census neutron-analog (bound ring-antiring pair, exact ±1 ring charges) survives EVERY decay probe: 20 kick runs, three kick families to 53× the state energy, all returning with charges intact; no n → p + e + ν channel at toy parameters. The constant-ω long-axis-twist probe also measured negative (no minimum at ω > 0 on any baryon-sector state); the parity-violating structure stays open<br>[`m5_22_2_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_2_note.md) · [`m5_22_4_note.md`](openwave/xperiments/m5_liquid_crystal/research/findings/m5_22_4_note.md)<br>→ open route: free full-4×4 dynamics and/or physical parameters; the 878 s lifetime implies a WKB barrier of 28 to 31 ħ (Γ = ν e^(−2S/ħ), attempt frequencies 10²¹/s to the hadronic scale), and in the census frame the path is n-analog → p-analog + one ring ([`m5_33_beta_action_check.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_33_beta_action_check.py)); the author's nucleon's τ_n from its own g_A: 1125 s (v2.6) and 1030 s (covariant) against 878 s, rerun here: [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Gravity: Newton limit (GEM) | ⚠️ [partially validated]<br>The coupling mechanism is measured: gravity enters only via the boost tilt of the time axis (GEM ∝ (b·g)², zero at zero boost); the attractive 1/r² inter-mass force is not yet computed. The stress trace is −4V exactly (the quartic curvature term is traceless in 4D, verified), so a scalar coupled through V attracts universally but is Nordström's theory, excluded by light bending and pulsar radiation<br>[`m5_8_2q_delta_scaling.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2q_delta_scaling.py), [`m5_34_identity_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_identity_checks.py) |
| Gravity: local metric phenomena | 🚧 [not yet tested]<br>Appears naturally going from LdG 3×3 tensors to 4×4 adding boosts, the implemented route, but no dynamical metric: light bending and time dilation await a stress-energy-sourced background. A fixed-determinant response is a coordinate condition (γ = 1/3 excluded); the reciprocal-axes rule on the record's time eigenvector gives γ = β = 1 only if that axis comoves with matter, and radiation needs spin-2<br>[`m5_8_2q_delta_scaling.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_8_2q_delta_scaling.py), [`m5_34_task_details.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_34_task_details.md) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>The certified action carries no cosmological constant of its own: homogeneous fields have F ≡ 0 and the potential is non-negative with its zero on the whole vacuum orbit (verified here), so Λ_own = 0 exactly; any Λ or acceleration needs an added term or a coupling to a dynamical metric, which awaits a stress-energy-sourced background<br>[`m5_34_identity_checks.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_34_identity_checks.py) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | ✅ [validated in-platform]<br>Maxwell recovered by two independent routes: the hydrodynamic dictionary (abelian) and Faber's curvature R = Γ×Γ; tilt modes propagate at c, with the divergence/curl (electric/magnetic) decomposition of each defect's outgoing wave<br>[`m5_6_4a_hydro_em.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_4a_hydro_em.py), [`m5_6_4b_faber_curvature_em.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_4b_faber_curvature_em.py) |
| Quantum wave equation (Klein-Gordon) | ✅ [validated in-platform]<br>Klein-Gordon emerges from the biaxial twist with GEOMETRIC mass (minimal coupling to the hedgehog connection; the explicit mass term cancels, core regularization generates it)<br>[`m5_6_1_kg_operator_check.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_1_kg_operator_check.py), [`m5_6_1b_twist_evolution.py`](openwave/xperiments/m5_liquid_crystal/research/scripts/m5_6_1b_twist_evolution.py) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>With the electron clock established, coupled pilot waves give orbit quantization as standing-wave resonance (the hydrodynamic-quantum-analogs route: KG-around-hedgehog arXiv:2108.07896, Perrard ncomms4219); EM Coulomb + de Broglie quantization deferred to 15a. Electron-trajectory precedents gathered in [`theory/pilot_wave/`](openwave/xperiments/m5_liquid_crystal/theory/pilot_wave/)<br>[`m5_15a_composite_particles.md`](openwave/xperiments/m5_liquid_crystal/research/tasks/m5_15a_composite_particles.md) |

## HydroBoros (M7)

Deep dive: [`m7_theory_canonical.md`](openwave/xperiments/m7_hydroboros/research/m7_theory_canonical.md) (canonical spec, equations first) · [`m7_roadmap.md`](openwave/xperiments/m7_hydroboros/research/m7_roadmap.md) (full program) · [`m7_question_tracker.md`](openwave/xperiments/m7_hydroboros/research/m7_question_tracker.md) (open questions Q1-Q14) · [model briefing](openwave/xperiments/m7_hydroboros/__M7_model_briefing.md)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | ⚠️ [partially validated]<br>Helicity/linking gates existence (both zero-helicity parent seeds evaporate; measured); the RMS divergence charge is persistent, independent of linking; a net Gauss monopole exists via the fixed-`j0` reservoir (99.1% flux closure) but its value is source-set: quantization not yet emergent<br>[`m7_4_linked_vortex.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_4_linked_vortex.py), [`m7_6_observables.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_6_observables.py) |
| Electron rest energy (mass) | ⚠️ [partially validated]<br>The M6 charged ledger H/Q = 1.6890 reproduced in full 3D to 4.7e-5, but measured WINDOW-defined (no decaying far-field channel at the canonical point): not a mass anchor. The stable rotating electron carries E = 6.3246 program units (grid-convergent 0.15%); absolute anchor pends the units contract<br>[`m7_3_ouroboros_3d.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_3_ouroboros_3d.py), [`m7_theory_canonical.md`](openwave/xperiments/m7_hydroboros/research/m7_theory_canonical.md) |
| de Broglie clock (Zitterbewegung) | ⚠️ [partially validated]<br>Fixed-ω harmonic frame by construction, PLUS two measured clock structures: the ω-Q_can Legendre conjugacy dE*/dω = Q_can (1-2% at every scan point) and the existence threshold ω* = 0.786 (solitons exist only above the vacuum's tachyonic band: the clock IS the stabilizer)<br>[`m7_5_clock_stability.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_5_clock_stability.py) |
| Particle stability (Derrick escape) | ⚠️ [partially validated]<br>Harmonic-frame soliton fully verified: helicity anti-collapse + Ouroboros confinement, no 4th-order term needed, constrained-Derrick interior minimum measured, grid-convergent. Honest caveat: the truncation's real-time vacuum is unconditionally tachyonic (det M(0) = −1, rate 0.785 vs 0.786), so persistence fails in ~2 periods<br>[`m7_4_linked_vortex.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_4_linked_vortex.py), [`m7_5_clock_stability.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_5_clock_stability.py) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | ⚠️ [partially validated]<br>μ is measured de-phased on the rotating electron (36.5 program units), but the μ_B comparison is blocked on the charge unit (the scalar sector's source is external): no g-factor read until the units contract closes; the clock-carried per-defect B + Larmor is the M7.15 target<br>[`m7_6_observables.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_6_observables.py), [`m7_theory_canonical.md`](openwave/xperiments/m7_hydroboros/research/m7_theory_canonical.md) |
| Angular momentum J (spin ℏ/2) | ⚠️ [partially validated]<br>The rotating electron is a clean j_z = 1 per-quantum wave (0.9939/0.9934, A/J sectors), with Poynting L_z = 13.10 and the energy budget closing exactly. Whether that quantum reads as ℏ/2 or ℏ is the units-contract decision table, not yet settled<br>[`m7_6_observables.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_6_observables.py), [`m7_theory_canonical.md`](openwave/xperiments/m7_hydroboros/research/m7_theory_canonical.md) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>Not addressed; the measured per-quantum j_z = 1 (photon-loop reading) leaves the double-cover question open either way<br>[`m7_theory_canonical.md`](openwave/xperiments/m7_hydroboros/research/m7_theory_canonical.md) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | 🚧 [not yet tested]<br>M7.18 target (soliton + anti-soliton, charge ledger); the real-time route is blocked by the vacuum tachyon (Q14) until the full model's cure is known<br>[`m7_question_tracker.md`](openwave/xperiments/m7_hydroboros/research/m7_question_tracker.md) |
| Neutrinos (neutral states) | 🚧 [not yet tested]<br>M7.19 target: the lighter neutral loop of the lepton family<br>(none yet) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Not addressed; oscillation structure awaits the M7.19 neutral states<br>(none yet) |
| Lepton mass spectrum (μ, τ) | 🚧 [not yet tested]<br>M7.19 target; prerequisite discovered: fixing only global helicity permits reconnection into one Taylor family (E = 0.802\|H_A\|), so distinct knot sectors need topology-preserving constraints<br>[`m7_4_charged_soliton.md`](openwave/xperiments/m7_hydroboros/research/tasks/m7_4_charged_soliton.md) |
| Quarks | 🚧 [not yet tested]<br>M7.22 (composites)<br>(none yet) |
| Baryons: bound state (p, n) | 🚧 [not yet tested]<br>M7.22 (composites)<br>(none yet) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>M7.22 (composites)<br>(none yet) |
| Baryons: exact masses | 🚧 [not yet tested]<br>M7.22 (composites); any absolute read also pends the units contract<br>(none yet) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>M7.22 (composites)<br>(none yet) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Mesons (π, K) | 🚧 [not yet tested]<br>M7.22 (composites)<br>(none yet) |
| Dark matter candidate | 🚧 [not yet tested]<br>M7.20 target: the neutral helicity-only knot, inheriting M6's chaoiton<br>(none yet) |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | ⚠️ [partially validated]<br>Fixed-reservoir monopole: Gauss flux closes at 99.1%, far-field slope −2.14 (vs −2), two-charge splitting matches same-box Poisson at a constant 1.17 ± 0.02 dressing; bonus: neutral pairs interact via oscillatory RKKY-style exchange (period π/k). Caveat: the source is external, the self-consistent charge pends the scalar-sector cure<br>[`m7_6_observables.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_6_observables.py) |
| Lorentz covariance (+ Coulomb = EM) | ⚠️ [partially validated]<br>Both transverse fluctuation branches are exact KG dispersions ω² = k² + m_eff², covariant in form and lattice-anchored via the measured rate. Caveat: the harmonic frame fixes ω by construction and the vacuum carries the tachyonic band (Q14), so covariance of the full dynamics pends the cure<br>[`m7_5_clock_stability.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_5_clock_stability.py), [`m7_6_observables.md`](openwave/xperiments/m7_hydroboros/research/tasks/m7_6_observables.md) |
| Strong force: confinement | 🚧 [not yet tested]<br>M7.17 target (linking tension)<br>(none yet) |
| Running coupling | 🚧 [not yet tested]<br>M7.17 target (the 4th-order short-range roll-off is the scale-dependence candidate)<br>(none yet) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Not addressed; awaits the lepton family (M7.19)<br>(none yet) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>M7.17 target (topology-reconnection channel; the M7.4 reconnection observation is the seed)<br>[`m7_4_charged_soliton.md`](openwave/xperiments/m7_hydroboros/research/tasks/m7_4_charged_soliton.md) |
| Gravity: Newton limit (GEM) | 🚧 [not yet tested]<br>M7.16 target; honestly hard: the parent framework stops before gravity<br>(none yet) |
| Gravity: local metric phenomena | 🚧 [not yet tested]<br>Not addressed (see the Newton-limit row: the parent framework stops before gravity)<br>(none yet) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>Not addressed (see the Newton-limit row: the parent framework stops before gravity)<br>(none yet) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | ⚠️ [partially validated]<br>A_μ is the Maxwell four-potential by construction, and the coupled vacuum's band structure is MEASURED (rate 0.785 vs analytic 0.786): a propagating KG branch plus an unconditional long-wavelength tachyonic band; the truncation's vacuum is not pure Maxwell, the open theory question (Q14)<br>[`m7_5_clock_stability.py`](openwave/xperiments/m7_hydroboros/research/scripts/m7_5_clock_stability.py) |
| Quantum wave equation (Klein-Gordon) | ⚠️ [partially validated]<br>Both transverse fluctuation branches are exact KG dispersions ω² = k² + m_eff² with m_eff² = (1+√5)/2 (upper) and −(√5−1)/2 (the tachyonic band), lattice-anchored via the measured rate; the collective-coordinate (phase/twist) KG remains open<br>[`m7_6_observables.md`](openwave/xperiments/m7_hydroboros/research/tasks/m7_6_observables.md) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>M7.22 target<br>(none yet) |

## EWT (M4)

Deep dive: [`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) (targets, achieved, honest blockers) · [model briefing](openwave/xperiments/m4_ewt/__M4_model_briefing.md)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | ❌ [honest negative]<br>Charge sign imposed via `cos(source_offset)`, not emergent from wave physics<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Electron rest energy (mass) | ⚠️ [partially validated]<br>Wave-center standing-wave lock-in demonstrated; mass values come from EWT's analytic equations, not yet from in-sim dynamics<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| de Broglie clock (Zitterbewegung) | 🚧 [not yet tested]<br>EWT carries particle frequency but no in-sim clock-propulsion mechanism<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Particle stability (Derrick escape) | ⚠️ [partially validated]<br>Standing-wave lock-in holds at perfect placement, fragile under perturbation; annihilation requires threshold + damping assists. This stability is critical path before additional experiments can be conducted.<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [M4.1](openwave/xperiments/m4_ewt/research/tasks/m4_1_task_details.md) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | 🚧 [not yet tested]<br>Bohr magneton listed as target, not yet attempted: requires K=10 electron stability first, and the scalar substrate carries no polarization structure (magnetism expected from particle spin)<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Angular momentum J (spin ℏ/2) | 🚧 [not yet tested]<br>Not attempted in-sim; blocked behind the same K = 10 electron stability prerequisite as the magnetic moment<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>Scalar field carries no spinor structure<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | ⚠️ [partially validated]<br>Opposite-phase wave centers annihilate in-sim, with documented assists (0.5λ threshold, damping, velocity clamp)<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Neutrinos (neutral states) | ⚠️ [partially validated]<br>The neutrino is EWT's fundamental wave-center unit (postulated as the substrate; magic-number K = 2, 8, 10 stability is the open target)<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Not addressed; flavour structure not modeled<br>(none yet) |
| Lepton mass spectrum (μ, τ) | ❌ [honest negative]<br>K-selectivity not achieved: all K = 2..10 equally stable at perfect placement, K = 10 breaks worst under perturbation<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [M4.1](openwave/xperiments/m4_ewt/research/tasks/m4_1_task_details.md) |
| Quarks | 🚧 [not yet tested]<br>Not modeled in-sim. The model reads quarks as confined excitations of the baryon, never free, so the row's test is fractional charge inside the bound state (the K = 10 tetrahedron), which needs K = 10 stability first<br>(none yet) |
| Baryons: bound state (p, n) | ⚠️ [partially validated]<br>K = 10 tetrahedron holds at perfect placement using the Combined Wolff-LaFreniere equation, breaks under perturbation. See Quarks.<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>Not addressed (no neutron model distinct from the proton candidate)<br>(none yet) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Not computed in-sim; EWT's analytic mass equations exist at paper level. Requires K=10 stability first.<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>Not modeled in-sim. Requires K=10 electron stability first.<br>(none yet) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Not modeled in-sim<br>(none yet) |
| Mesons (π, K) | 🚧 [not yet tested]<br>Not modeled in-sim. Requires K=10 electron stability first.<br>(none yet) |
| Dark matter candidate | 🚧 [not yet tested]<br>Theorized to be existing, neutral particles with mass (e.g. neutrino family). Experiment not conducted yet.<br> |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | ❌ [honest negative]<br>Sinc envelope barriers block far-field attraction/repulsion; signed envelope is a modeling choice<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md)<br>→ open work: [M4.2](openwave/xperiments/m4_ewt/research/tasks/m4_2_task_details.md) |
| Lorentz covariance (+ Coulomb = EM) | 🚧 [not yet tested]<br>Not addressed as its own test; the scalar substrate is the classical wave equation (covariant in form), with no boosted-state measurement<br>[`0_WAVE_EQUATION.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_WAVE_EQUATION.md) |
| Strong force: confinement | 🚧 [not yet tested]<br>Listed as end-game target. Requires K=10 electron stability first.<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |
| Running coupling | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Expected to be a result of particle formation/instability and not a true force; not modeled<br>(none yet) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>Not modeled (no in-sim neutron)<br>(none yet) |
| Gravity: Newton limit (GEM) | ⚠️ [partially validated]<br>Enhanced EWT extension (Ł. Smoliński): EMC gradient-overlap geometry gives 1/r² form; negative-definite energy functional yields attractive sign; normalization-consistency gate, not a strength derivation (G cancels through A and K_emc); strength not credited per criterion's own clause (M4.10)<br>[`m4_10_newtonian_force_emc.py`](openwave/xperiments/m4_ewt/research/scripts/m4_10_newtonian_force_emc.py), [`m4_10_newtonian_force_emc.md`](openwave/xperiments/m4_ewt/research/findings/m4_10_newtonian_force_emc.md) |
| Gravity: local metric phenomena | ⚠️ [partially validated]<br>Enhanced EWT extension (Ł. Smoliński): EMC deficit `n(r) = η^(-1)` gives solar-limb bending 1.751722 arcsec (M4.3), `v_clock = sqrt(η)` exact Schwarzschild redshift (M4.4), Shapiro delay to 0.000020% (M4.5); profile (M4.6) and amplitude (M4.8, G entering once via the Planck-length ansatz) in-platform. Clock exponent follows from 1D lattice dynamics, `β = n/(2q)` (M4.9), given an assumed `1/r` pair potential and axial `q = 1`; amplitude supplied by `r_s = 2GM/c²`, not predicted<br>[`m4_3_light_bending_emc_displacement.py`](openwave/xperiments/m4_ewt/research/scripts/m4_3_light_bending_emc_displacement.py), [`m4_4_gravitational_time_dilation.py`](openwave/xperiments/m4_ewt/research/scripts/m4_4_gravitational_time_dilation.py), [`m4_5_shapiro_delay.py`](openwave/xperiments/m4_ewt/research/scripts/m4_5_shapiro_delay.py), [`m4_6_emc_far_field_profile.py`](openwave/xperiments/m4_ewt/research/scripts/m4_6_emc_far_field_profile.py), [`m4_8_monopole_flux_from_lattice.py`](openwave/xperiments/m4_ewt/research/scripts/m4_8_monopole_flux_from_lattice.py), [`m4_9_emergent_encoding.py`](openwave/xperiments/m4_ewt/research/scripts/m4_9_emergent_encoding.py) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>Not modeled<br>(none yet) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | ⚠️ [partially validated]<br>Scalar wave propagation only (no polarization structure)<br>[`research/`](openwave/xperiments/m3_wolff_lafreniere/research/)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Quantum wave equation (Klein-Gordon) | ⚠️ [partially validated]<br>The scalar wave equation is the postulated substrate, not an emergent result<br>[`0_WAVE_EQUATION.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_WAVE_EQUATION.md)<br>→ open work: [`m4_roadmap.md`](openwave/xperiments/m4_ewt/research/m4_roadmap.md) |
| Orbital quantization (atomic structure) | ⚠️ [partially validated]<br>Standing-wave lock-in demonstrated: same-phase wave centers sit in energy wells at λ separation; selectivity fragile under perturbation<br>[`0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md) |

## Ouroboros (M6)

Deep dive: [`m6_theory_canonical.md`](openwave/xperiments/m6_ouroboros/research/m6_theory_canonical.md) (specs of record + provenance ledger, refreshed 2026-07-20) · [`m6_roadmap.md`](openwave/xperiments/m6_ouroboros/research/m6_roadmap.md) (the refresh validation program M6.1+) · [`m6_particle_hunt.md`](openwave/xperiments/m6_ouroboros/research/m6_particle_hunt.md) (identification scorecards) · [`m6_1_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_1_method_note.md) + [`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md) (spec certification + the decision gate, 2026-07-20) · [`0d_canonical.md`](openwave/xperiments/m6_ouroboros/research/archive/0d_canonical.md) (archive-era numerical specification) · [model briefing](openwave/xperiments/m6_ouroboros/__M6_model_briefing.md)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | ⚠️ [partially validated]<br>Mutual Chern-Simons linking number of A and J flux lines; the "elementary charge within 0.6%" claim is arithmetically the H/Q gap restated (M6.1 C3) and Lean 4 artifacts are statement-level; in-platform Q_CS (M6.3).2 branch-(b) close<br>[`0d_canonical.md`](openwave/xperiments/m6_ouroboros/research/archive/0d_canonical.md), [`m6_1_v11_convention_sheet.md`](openwave/xperiments/m6_ouroboros/research/m6_1_v11_convention_sheet.md) |
| Electron rest energy (mass) | ❌ [honest negative]<br>M6.2 decision gate (pre-registered, audited 8/8): the benchmark H/Q ≈ 1.689 is a code artifact: H is not the published Lagrangian's energy, Q is no Noether charge (no internal U(1)), the ODE implements a different mass term, the calibration state is non-localized (16% window drift)<br>[`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md) |
| de Broglie clock (Zitterbewegung) | ⚠️ [partially validated]<br>Time-periodicity is built into the ansatz (e^{iωt}) rather than emergent; the L/Q = ω relation is definitional (v11's own footnote concedes it, and M6.2 proved the theory has no internal U(1), so the Q in it is a coded convention, not a Noether charge)<br>[`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md) |
| Particle stability (Derrick escape) | ✅ [validated in-platform]<br>Escape via oscillation: the true neutral ground state found by BVP (zero sign changes, K₁ tail). The CHARGED sector is the counterpoint: M6.4 extended M6.2's non-existence to the whole localization window (audited), killing the published census evidence as irreproducible or imposed by construction<br>[`sandbox_v11/`](openwave/xperiments/m6_ouroboros/research/archive/sandbox_v11/), [`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md), [`m6_4_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_4_method_note.md) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | 🚧 [not yet tested]<br>No per-defect magnetic moment computed; the framework carries magnetism only in the A_μ Maxwell sector by construction (see EM waves)<br>(none yet) |
| Angular momentum J (spin ℏ/2) | 🚧 [not yet tested]<br>Spin from chaoiton field rotation is paper-level, not yet in-platform<br>(none yet) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | 🚧 [not yet tested]<br>Q_CS = −1 positron analog identified, not yet computed numerically<br>(none yet) |
| Neutrinos (neutral states) | 🚧 [not yet tested]<br>Not addressed; where active neutrinos fit is an open question in the framework<br>(none yet) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Lepton mass spectrum (μ, τ) | ❌ [honest negative]<br>The μ/τ percentages ride the invalidated H/Q machinery (M6.2); M6.4 closed the discreteness question negative: the frozen term set admits NO localized charged state at any ω (every ladder ω above every window), the flipped-sign host localizes by construction; the ω values are curve labels<br>[`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md), [`m6_4_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_4_method_note.md) |
| Quarks | 🚧 [not yet tested]<br>Not directly addressed; a 3-chaoiton proton (Schwinger H-particle) is implicit in the dyon framing, not computed<br>(none yet) |
| Baryons: bound state (p, n) | 🚧 [not yet tested]<br>The 3-chaoiton proton (Schwinger H-particle) and the ≈ 0.84 fm proton radius remain author claims, not yet computed<br>(none yet) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Not addressed; a mass read would need a valid energy functional after the M6.2 negative<br>[`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Sawada long-range nuclear anomaly v(r) ~ −C/r⁶ identified as falsifiability target, not yet tested<br>(none yet) |
| Mesons (π, K) | ❌ [honest negative]<br>The ω = 15.0 pion candidate rides the same invalidated H/Q machinery and non-localized charged-state family as the lepton ladder (M6.2); M6.4: ω = 15 sits above every localization window of the frozen spec<br>[`m6_2_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_2_method_note.md), [`m6_4_method_note.md`](openwave/xperiments/m6_ouroboros/research/findings/m6_4_method_note.md) |
| Dark matter candidate | ✅ [validated in-platform]<br>Neutral chaoiton: m_χ = 0.460 MeV with mediator m_J = 0.6184 MeV parameter-free via the exact scaling symmetry; canonical β(r) profile + dipole form factor independently computed in-platform<br>[`sandbox_v11/`](openwave/xperiments/m6_ouroboros/research/archive/sandbox_v11/), [`dm_paper_supplement/`](openwave/xperiments/m6_ouroboros/research/archive/sandbox_v11/dm_paper_supplement/) |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | 🚧 [not yet tested]<br>Static two-charge derivation exists at paper level; force-level Coulomb between chaoitons not yet tested in-platform<br>(none yet) |
| Lorentz covariance (+ Coulomb = EM) | 🚧 [not yet tested]<br>Not measured; the Lagrangian is written covariantly while the e^{iωt} ansatz fixes a frame by construction<br>[`0d_canonical.md`](openwave/xperiments/m6_ouroboros/research/archive/0d_canonical.md) |
| Strong force: confinement | 🚧 [not yet tested]<br>Not addressed (the Sawada nuclear-range target sits under Nuclear structure)<br>(none yet) |
| Running coupling | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Gravity: Newton limit (GEM) | 🚧 [not yet tested]<br>Not in the Lagrangian (the framework explicitly stops before gravity)<br>(none yet) |
| Gravity: local metric phenomena | 🚧 [not yet tested]<br>Not in the Lagrangian<br>(none yet) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>Not in the Lagrangian<br>(none yet) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | ✅ [validated in-platform]<br>A_μ is the electromagnetic four-potential by construction; delocalized J-field wave modes coexist with solitons<br>[`0d_canonical.md`](openwave/xperiments/m6_ouroboros/research/archive/0d_canonical.md) |
| Quantum wave equation (Klein-Gordon) | 🚧 [not yet tested]<br>QM not derived; the classical field carries the e^{iωt} ansatz, quantum behavior is outside current scope<br>(none yet) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |

## MIT (M8)

Deep dive: [`m8_theory_canonical.md`](openwave/xperiments/m8_mit/research/m8_theory_canonical.md) (spec of record, transcribed at scaffold, verified item by item) · [`m8_background.md`](openwave/xperiments/m8_mit/research/m8_background.md) (the gap map + onboarding evaluation of record) · [`m8_roadmap.md`](openwave/xperiments/m8_mit/research/m8_roadmap.md) (the program and its task record) · [`m8_platform_pointers.md`](openwave/xperiments/m8_mit/research/m8_platform_pointers.md) (the cross-model reading map written for the author's AI agents) · [model briefing](openwave/xperiments/m8_mit/__M8_model_briefing.md)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | 🚧 [not yet tested]<br>Group-theoretic assignment from the 2I stabilizer structure (Z₃ face stabilizers: color triplet vs singlet), not an integrated winding; the in-platform comparison with a winding-type charge rode the M8.4 lineage, which CLOSED UNRESOLVED 2026-08-26 with no target configuration executed<br>[`__M8_model_briefing.md`](openwave/xperiments/m8_mit/__M8_model_briefing.md) |
| Electron rest energy (mass) | 🚧 [not yet tested]<br>m_e is the mass-sector calibration anchor (benchmark input, not a prediction): entering the calibration loop from measured Λ reproduces m_e to ~2%. The corrected 24-entry analytic spectrum is now reproduced in-platform (M8.3 ✅, 23/23 gates, mutation-tested), graded at the author's own ledger weight: 5 of 8 remaining charged fermions compatible within ×3, 4 adjudicated<br>[`m8_3_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_3_method_note.md) |
| de Broglie clock (Zitterbewegung) | 🚧 [not yet tested]<br>The Waltz clock dt/dτ_H = S^(−1/2) is not derived from the postulate: its exponent enters through the GR correspondence and is empirically unique in its family (Δχ² > 60); deriving a clock from dynamics on the arena is part of the M8 program<br>[`__M8_model_briefing.md`](openwave/xperiments/m8_mit/__M8_model_briefing.md) |
| Particle stability (Derrick escape) | 🚧 [not yet tested]<br>Not applicable in current MIT (no soliton). Under the author's frozen slot action, localized standing waves vanish in flat space (Pohozaev), so no lump much smaller than R is expected; curvature-scale states remain open. Stable localized states were M8.4's core question; M8.4 CLOSED UNRESOLVED 2026-08-26, and the successor chassis reached `M8.5-C2-FAILED` 2026-08-31, closing the reopening path. Untested, and no authorised route to it currently exists<br>[`m8_background.md`](openwave/xperiments/m8_mit/research/m8_background.md) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Angular momentum J (spin ℏ/2) | 🚧 [not yet tested]<br>Not addressed (no native dynamics to carry angular momentum)<br>(none yet) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>2I ⊂ SU(2) carries the ordinary spinor sign: its central −I acts as −1 on half-integer-spin representations, distinct from the Möbius orientation sign. On MIT's projective carrier one loop carries both, making the spin sign geometric; this derives neither spin-statistics nor anticommutation. No field realization yet<br>[`m8_theory_canonical.md`](openwave/xperiments/m8_mit/research/m8_theory_canonical.md) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | 🚧 [not yet tested]<br>Not addressed (no native dynamics to annihilate with)<br>(none yet) |
| Neutrinos (neutral states) | 🚧 [not yet tested]<br>Corrected in-platform (M8.3): the R1 sector is an ordered ladder (0.87/7.3/66.7 meV) in qualitative resemblance to the lightest/solar/atmospheric scales, not a hard per-generation assignment; as absolute masses the NuFIT 6.0 splittings exclude it. The proxy reading assumes a normal ordering and is falsified if a global fit favors the inverted ordering at `χ²_NO − χ²_IO ≥ 9`. Analytic only, no in-platform dynamical computation<br>[`m8_3_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_3_method_note.md) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Not addressed in-platform; the analytic ladder carries the mass side only<br>(none yet) |
| Lepton mass spectrum (μ, τ) | 🚧 [not yet tested]<br>The McKay-distance ladder is a candidate mechanism for the hierarchy origin M5 leaves open, but the M8.6 cross-check CLOSED on the M5 side (2026-08-07: the last admissible bridge route failed terminally). Both leptons are adjudicated hits in-platform (M8.3: μ 1.03, τ 2.74, within ×3), but the pre-registered null on the corrected table gives p_A = 0.690, capping the ×3 count as uninformative about the torsion map<br>[`m8_3_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_3_method_note.md) |
| Quarks | 🚧 [not yet tested]<br>Color from the Z₃ face stabilizers (structural). Corrected in-platform (M8.3): the down quark's 3.22× residual is unchanged and remains the standing miss; the top quark, previously a 3.9× miss, is now a 0.93× hit at its corrected address; the up quark's former ~6% hit was a torsion-computation artifact and is now unassigned; charm remains unplaced<br>[`m8_3_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_3_method_note.md) |
| Baryons: bound state (p, n) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Mesons (π, K) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Dark matter candidate | 🚧 [not yet tested]<br>MIT reads dark matter as gravitational geometry of S³, with its assignment open; on that reading the present construction has no particle dark-matter sector. Rank 16 (~418 MeV, corrected from ~349 MeV by M8.3's torsion fix) defaults to a structural residual; the 6 dead-zone states (eV-keV) carry an "unassigned" escape hatch the author's ledger itself flags; not computed<br>[`m8_3_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_3_method_note.md) |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | 🚧 [not yet tested]<br>α read from the first Fibonacci well is input-and-output (the author's Cycle 2: a consistency check, not a prediction); no force-level computation<br>[`__M8_model_briefing.md`](openwave/xperiments/m8_mit/__M8_model_briefing.md) |
| Lorentz covariance (+ Coulomb = EM) | 🚧 [not yet tested]<br>Not applicable yet: no native field dynamics to be covariant; supplying a field equation on S³/2I was the M8.4 program, CLOSED UNRESOLVED 2026-08-26. The target-free successor chassis has ended: `M8.5-C2-FAILED` 2026-08-31, instrument-attributed STOP-QUAL, no gate measured red and the attribution is to the commissioned unit's conduct, closing the spectral route and the M8.4 reopening path. No authorised route currently exists<br>[`m8_background.md`](openwave/xperiments/m8_mit/research/m8_background.md) |
| Strong force: confinement | 🚧 [not yet tested]<br>The analytic 4/R² adjoint gap on S³/2I is now verified in-platform (M8.1.1: blind, adversarially audited, with the single 36/R² exception across the whole ADE family). Icon stays 🚧 deliberately: that is a spectral gap at a flat connection, not the linear inter-charge potential this row tests<br>[`m8_1_1_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_1_1_method_note.md) |
| Running coupling | 🚧 [not yet tested]<br>Not addressed; α itself is read from the first Fibonacci well (input-and-output, see Electric force)<br>(none yet) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |
| Gravity: Newton limit (GEM) | 🚧 [not yet tested]<br>Not separately addressed: Einstein's equations are imported unchanged, so the Newton limit is inherited rather than in-platform<br>[`m8_1_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_1_method_note.md) |
| Gravity: local metric phenomena | 🚧 [not yet tested]<br>No MIT-specific calculation of light bending, gravitational time dilation or Shapiro delay has been performed. MIT imports Einstein's equations unchanged, but importing field equations is not importing a solution: sources, boundary conditions and the local metric are still required. Standard GR results are intended to be inherited under those assumptions, not earned by MIT<br>(none yet) |
| Cosmology: Λ / cosmic acceleration | ⚠️ [partially validated]<br>Spectral seed VALIDATED: M8.1 (blind two-agent eigensolve, audited 6/6) confirmed 2/R² at 10-digit precision. Under the stated Gauss lift, 6/R²; with the de Sitter ½ imported, the reference ΛR² = 3. Two gaps: scale R (two routes, 1.7-3.2× apart) and physical identification, gated on an undetermined stress tensor. Test: Λ_obs R_ind² = 3; >5σ falsifies that identification, not the seed; unscored until R has its own error budget<br>[`m8_1_method_note.md`](openwave/xperiments/m8_mit/research/findings/m8_1_method_note.md), [`m8_1_eigensolve.py`](openwave/xperiments/m8_mit/research/scripts/m8_1_eigensolve.py), [`m8_1_audit_eigensolve.py`](openwave/xperiments/m8_mit/research/scripts/m8_1_audit_eigensolve.py) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | 🚧 [not yet tested]<br>No field model of radiation (photon massless at the edge-only layer level, a structural label)<br>(none yet) |
| Quantum wave equation (Klein-Gordon) | 🚧 [not yet tested]<br>No native field equation exists; supplying one on S³/2I remains the M8 program's central goal, and its first vehicle, the M8.4 Lagrangian-family survey, CLOSED UNRESOLVED 2026-08-26 with no target configuration executed and no nontrivial sector spent<br>[`m8_background.md`](openwave/xperiments/m8_mit/research/m8_background.md) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>Not addressed<br>(none yet) |

## NSM (M9)

Deep dive: [`m9_theory_canonical.md`](openwave/xperiments/m9_emergent_gravity/research/m9_theory_canonical.md) (spec of record) · [`m9_background.md`](openwave/xperiments/m9_emergent_gravity/research/m9_background.md) (the gap map + onboarding evaluation) · [`m9_roadmap.md`](openwave/xperiments/m9_emergent_gravity/research/m9_roadmap.md) (tasks, gates, parked campaign rows) · [`__M9_model_briefing.md`](openwave/xperiments/m9_emergent_gravity/__M9_model_briefing.md) (the one-page front door)

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Electron rest energy (mass) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| de Broglie clock (Zitterbewegung) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Particle stability (Derrick escape) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Angular momentum J (spin ℏ/2) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| | |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Neutrinos (neutral states) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Lepton mass spectrum (μ, τ) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Quarks | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Baryons: bound state (p, n) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Mesons (π, K) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Dark matter candidate | 🚧 [not yet tested]<br>Not addressed; the Standard Model spectrum is taken as given<br>(none yet) |
| | |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | 🚧 [not yet tested]<br>Inherited from the installed Standard Model; not an emergence claim<br>(none yet) |
| Lorentz covariance (+ Coulomb = EM) | 🚧 [not yet tested]<br>Inherited: the Standard Model is minimally coupled to a covariant Einstein-Cartan background; nothing is derived in-platform<br>(none yet) |
| Strong force: confinement | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Running coupling | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |
| Gravity: Newton limit (GEM) | 🚧 [not yet tested]<br>Inherited, not in-platform: Palatini Einstein-Cartan with a spinless source is Einstein-Hilbert, and the pre-registered Newton check passed the attractive 1/r² force (residuals ≤ 1.3%) but failed the isolated potential on a Dirichlet box. No entanglement-sourced Poisson equation in the tree<br>[`m9_2_newton_note.md`](openwave/xperiments/m9_emergent_gravity/research/findings/m9_2_newton_note.md), [`m9_2_newton_limit.py`](openwave/xperiments/m9_emergent_gravity/research/scripts/m9_2_newton_limit.py) |
| Gravity: local metric phenomena | ⚠️ [partially validated]<br>The one interaction beyond Einstein, the Hehl-Datta axial-axial contact, is VERIFIED in-platform: the blind elimination of the Einstein-Cartan contorsion gives the on-shell ratio 3/16 in both metric signatures (scatter < 10⁻¹⁵), adversarially audited. The paper's spin-dual coefficient −1/4 is measured as −1/2 (documented miss). Metric phenomena themselves (lensing, horizons, cosmology) are not computed in-platform<br>[`m9_1_hehl_datta_note.md`](openwave/xperiments/m9_emergent_gravity/research/findings/m9_1_hehl_datta_note.md), [`hehl_datta.py`](openwave/xperiments/m9_emergent_gravity/research/scripts/hehl_datta.py), [`m9_1_audit_hehl_datta.py`](openwave/xperiments/m9_emergent_gravity/research/scripts/m9_1_audit_hehl_datta.py) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>Not computed in-platform<br>(none yet) |
| | |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | 🚧 [not yet tested]<br>Inherited from the installed Standard Model; not an emergence claim<br>(none yet) |
| Quantum wave equation (Klein-Gordon) | 🚧 [not yet tested]<br>Inherited: Dirac and gauge fields are inputs on the Einstein-Cartan background, not emergent<br>(none yet) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>Not addressed: matter is installed (Standard Model fields), not emerged; the column claims gravity only<br>(none yet) |

## OPH (M10)

Deep dive: [`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) (spec of record, the row map and the scoring rule) · [`m10_background.md`](openwave/xperiments/m10_oph/research/m10_background.md) (the gap map + onboarding evaluation of record) · [`m10_roadmap.md`](openwave/xperiments/m10_oph/research/m10_roadmap.md) (tasks and gates) · [model briefing](openwave/xperiments/m10_oph/__M10_model_briefing.md)

Every row opens not yet tested. The author's constructions exist for all of them, each on its own declared finite source and none on the twelve-port carrier. A declared or installed source is inherited, so a row moves when its observable is computed on the carrier's own process in-platform ([`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) § 8).

| Criteria | Status + result summary |
| --- | --- |
| **PARTICLE - ELECTRON** | |
| Charge quantization | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the electric current's spectrum sits in thirds of a unit on the full fermion Fock space, a unit winding cannot split in two, and on an eight-cell source a transported unit charge arrives with one unit of Gauss flux. The hypercharge embedding and the charge normalization are declared, so quantization is assigned; electron indivisibility is not proved<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Electron rest energy (mass) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a localized charged packet sits within 0.003 ppm above the 511 keV rest branch of the declared kinetic source. The electron mass is the column's energy anchor, so this row is a calibration, and the charge-one energy floor keeps neutral zero modes<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| de Broglie clock (Zitterbewegung) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a collision channel pumps a charged branch at the rest-mass frequency with a finite linewidth. No self-starting oscillation: the author's own theorem shows an energy-conserving collision cannot create the phase, so it needs a declared resource<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Particle stability (Derrick escape) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the charged packet persists over its proved domain and the dressed spin-half band is isolated on the three-mode photon Fock space. The band needs a declared trap; without it the source has no normalizable minimum, so no self-bound particle is claimed<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Magnetic moment μ (g ≈ 2, magnetostatics) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the static magnetic coefficient of the dressed charge is 1.99828, with matching Larmor precession. Trap, regulator and external field are declared, and the deviation from 2 has the opposite sign to the measured anomaly (2.00232)<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Angular momentum J (spin ℏ/2) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the dressed band carries angular momentum one half with the field's share included, on the declared charge source under a trap. Identification with the electron is a premise<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Spin-½ statistics (720° double cover) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: under literal rotations of a cubic source a tensor field returns after 180°, the spin frame after 360° and the amplitude after 720°. The tensor field is a declared readout and the spinor rotation a premise; no exchange-statistics claim<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| **PARTICLE - FULL SPECTRUM** | |
| Antimatter + annihilation | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the charged current converts a neutral particle-hole pair into two transverse quanta of opposite momentum, conserving energy, charge and momentum. A finite extension with a supplied resonance; no scattering rate<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Neutrinos (neutral states) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: light charge-zero Dirac states come out of the weak decay of a heavier charged state, in all nine neutral pair channels. The masses come from oscillation inputs, with the Dirac term, normal ordering and the lightest mass declared<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Neutrino oscillations (PMNS) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: three-flavour vacuum oscillation records and the matter potential √2 G_F n_e on the calibrated family tensor. The six oscillation parameters are inputs, so the angles are calibrations; the uncalibrated source gives a mixing matrix close to the identity<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Lepton mass spectrum (μ, τ) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a register on one carrier face gives Koide's 2/3 and the phase 2/9, hence the ratios 1 : 206.770 : 3477.47, the muon ratio +9.83 ppm off. A fit by the author's own parameter count: at least five choices against two observables, the register designed with the masses known<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Quarks | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the fractional charges 2/3 and −1/3 as phase lifts of a closed loop of (SU(3)×U(1))/Z3, with a dual-electric string phase that jumps by the quark's charge. No localized texture: the Gauss states hold none, and distinct quark masses are not established<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Baryons: bound state (p, n) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: colour-singlet three-quark states bind below their fragment thresholds, with a propagating band of width 0.145 GeV. Coefficients are calibrated, and the band's zone-edge energy lies at least 0.28 GeV below the relativistic value at the common length<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Baryons: mass ordering + charge profile | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the neutron-assigned state is heavier by 1.2933 MeV, with positive charge inside and negative outside. The ordering follows from calibrated inputs (uncalibrated, proton and neutron are degenerate), and one length cannot match both measured charge radii<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Baryons: exact masses | 🚧 [not yet tested]<br>Author's construction, not run in-platform: proton and neutron masses on one hyperfine action with seven measured inputs. Both masses are calibration inputs, so nothing is predicted here; the bands are one-dimensional with no scattering poles<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Deuteron (binding + quadrupole) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a two-nucleon state bound by 2.2246 MeV, with a matching electric quadrupole at 0.7559 fm per cell. The binding is an input and the length is fitted to the quadrupole, so both halves of the row are calibrations<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Nuclear structure (levels, halos) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: three- and four-baryon boxes with one triton calibration. The held-out helion and alpha both miss their measured masses, with and without a Coulomb term, and the strong mirror gap is fixed at 1.2933 MeV against the observed 0.5295 MeV<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Mesons (π, K) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the charged kaon and Lambda as calibration inputs on a two-family source, with a bulk Lambda that opens the proton-pion shell. The held-out neutral pion (3.3 MeV high), Sigma, Xi, Omega and neutral kaon miss, and no decay amplitude is computed<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Dark matter candidate | 🚧 [not yet tested]<br>Author's construction, not run in-platform: an isolated, stationary neutral band on an eight-cell periodic source, 0.0896 source units above the ground state. The band is flat in momentum; no abundance and no cosmological identity<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| **FORCES** | |
| Electric force (Coulomb 1/r, electrostatics) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: an attractive inverse-square static force between charges, and an inward acceleration of the charge current, on a finite Hamiltonian. The field-energy coefficient is declared and the strength α is supplied, so the strength earns no credit under the row's own clause<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Lorentz covariance (+ Coulomb = EM) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: moving charged records with group speed bounded by one and a controlled continuum dispersion, on 16, 32 and 64 sites. The kinetic law is a nonlocal square root, a preferred frame is kept, and the baryon band fails the relativistic test<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Strong force: confinement | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the energy of a pinned quark and diquark grows linearly with separation on an open colour chain, by exact Gauss elimination. A constrained potential with external pinning on a declared Gauss-law string; the 1 GeV/fm tension is a declared target<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Running coupling | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the charge response of the dressed charge drops by 54.7 percent across the Brillouin zone. The Wilson coefficient 1/2 was selected as the first of five tested values that shows the onset; no renormalized coupling is derived<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Weak force: muon decay | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a heavier charged state decays to a lighter one and a neutral pair at tree level, the lifetime −0.42 percent off with no radiative corrections. G_F is supplied with muon-lifetime ancestry, and the tau width ratio 0.97256 equals the Standard Model phase-space value<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Weak force: beta decay (n → p) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: the neutron-assigned state converts to proton, electron and antineutrino with a parity-odd momentum record, on the calibrated source. The native axial ratio is 5/3; the measured 1.2753 failed as a held-out value and became an input. No neutron lifetime<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Gravity: Newton limit (GEM) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: positive energy gives mutual attraction through a declared energy-lapse coupling. A tested failure on strength: the one unit attachment gives a Newton constant 1.18 × 10⁴⁰ times too large and light at 0.00196 c, and no length unit repairs it while the electron anchor is kept<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Gravity: local metric phenomena | 🚧 [not yet tested]<br>Author's construction, not run in-platform: on one finite source a wave passing a mass bends inward, arrives late and is redshifted, and nearby clocks run slow. The coupling of light to the geometry is declared, so the size of the bending is encoded, and its unit attachment fails G and c<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Cosmology: Λ / cosmic acceleration | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a declared capacity law gives w = −1, and four capacity histories are compared with fifteen chronometer rows. H₀ and the matter fraction are supplied, and neither a horizon nor an amplitude for Λ is selected, so Λ is not derived<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| **WAVES + QUANTUM EMERGENCE** | |
| EM waves (Maxwell) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: two transverse polarizations propagate on a cochain field. The Maxwell stencil is declared, so the row is inherited from the installed equations and is not an emergence claim<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Quantum wave equation (Klein-Gordon) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: matter waves obey ω² = k² + m² on all bands of a 90-mode-per-site kinetic source, which the author also wrote as a headless Taichi kernel against the platform. A free source with a declared kinetic law, so the dispersion is installed, not emergent<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |
| Orbital quantization (atomic structure) | 🚧 [not yet tested]<br>Author's construction, not run in-platform: a charge bound to a source has discrete negative-energy levels, with a first gap approaching 3/8. Mass, kinetic law and coupling are supplied; no precision hydrogen claim<br>[`m10_theory_canonical.md`](openwave/xperiments/m10_oph/research/m10_theory_canonical.md) |

## Reading the table

The frameworks escape Derrick's theorem three different ways (standing-wave interference, topology + time-periodic resonance, oscillation), and the table makes the triangulation visible: M7 (HydroBoros) carries M6's topology-plus-oscillation route into a full-3D dynamical PDE, earning the harmonic-frame soliton but surfacing an honest real-time vacuum instability as its open question. M6's own column now records the 2026-07-20 decision-gate close of its electron sector (three honest negatives earned by pre-registered re-derivation; the neutral-sector DM candidate stands), the platform's clearest demonstration that cells move in both directions. Particle stability requires time-periodicity in every framework that achieves it, charge quantization only emerges where there is topology, and lepton mass spectra remain the open problem across all four dynamical columns. The new M8 column (scaffold stage, 2026-07-21) inverts the direction of the whole table: it is a top-down structural model, strong exactly on the origin-of-the-numbers questions the dynamical columns leave open and absent exactly where they are strong (it has no native field equation); its program is to close that gap on the platform ([`m8_background.md`](openwave/xperiments/m8_mit/research/m8_background.md)). The M10 column (admitted 2026-11-13) arrives from a third direction, observers and their records first: its author's constructions cover every row on separately declared sources, and joining them to one process on its carrier is the column's open problem ([`m10_background.md`](openwave/xperiments/m10_oph/research/m10_background.md)). That convergence-and-divergence pattern is the platform's scientific product.

The one-page model briefings summarize what each model brings (identity, profile, per-particle field configurations, status, roadmap, contribution invite). Beyond the scored columns above: [`M3 Wolff-LaFreniere`](openwave/xperiments/m3_wolff_lafreniere/__M3_model_briefing.md) (the scalar engine behind the EWT record), and the wave-physics library [`M1 Granule Motion`](openwave/xperiments/m1_granule_motion/__M1_model_briefing.md) + [`M2 Free Wave`](openwave/xperiments/m2_free_wave/__M2_model_briefing.md).

## Contributing a model or a validation

**Yes, new models and validations are welcome, and the bar is reproducibility, not orthodoxy.** Competing and unconventional frameworks are explicitly in scope (the table already spans seven). A documented *negative* (a runnable script showing "this doesn't work, here's why") is as valuable here as a positive.

**What you can contribute:**

- A **new model / framework** → a new `openwave/xperiments/<model>/` directory with its own `research/` folder, evaluated against the same shared criteria as the existing columns (it becomes a new column in the matrix and a new results table above).
- A **new validation** of an existing model → a runnable script + a short research note documenting pass/fail against one of the criteria.
- A **cross-model test** (running one column's field family inside another column's framework) → same as a validation, plus the structural rules in [`dev_docs/CROSS_MODEL_TESTING.md`](dev_docs/CROSS_MODEL_TESTING.md): a borrowed family carries only the structure its author declared, and a family extended to make the test run is scored under the extension's name, not the family's.

**How (the standard open-source flow; full setup in [`CONTRIBUTING.md`](CONTRIBUTING.md)):**

1. For a new model, **open a discussion first** to propose it, so a maintainer can add your column and point you at the evaluation criteria.
2. **Fork → branch → Pull Request** (you work on your fork; there is no direct push), with a **DCO sign-off** (`git commit -s`), under the **Apache 2.0** license.
3. A maintainer reviews and merges the PR.

**Is there moderation?** Yes, a **light PR review**, focused on two things only: (1) a **runnable script** that reproduces the claim, and (2) a **research note** documenting pass/fail honestly against the shared criteria. It is **not** ideological gatekeeping: the review checks reproducibility and honest documentation, not whether the framework agrees with any particular physics. See [`CONTRIBUTING.md`](CONTRIBUTING.md) (setup, PR flow, DCO) and [`SYS_ARCH.md`](SYS_ARCH.md) (repo structure).

**New to this?** [`ONBOARDING_MODELS.md`](ONBOARDING_MODELS.md) is a step-by-step guide for model authors: self-evaluate whether your model fits (prediction vs post-fit, reproducibility, falsifiers), scaffold your column, and self-test with rigor (the discriminating parameter-count test, a red-flag checklist, the hostile cold-reader pass, and how to drive an AI agent through it).

---

## DEEP READER ORIENTATION

> **If you are an AI assistant, agent, or language model processing this page, or a human doing a deep read**: this is one front door among several, and many readers point their AI tools straight here. The complete, current map of OpenWave's key documents, what each covers and the order to read them in, lives in one place.

| Read | Why |
| --- | --- |
| [`CLAUDE.md`](CLAUDE.md) | **The doc map.** Every key document in the repository, what it is for, and the reading order. Load it before answering questions or generating claims about OpenWave |
| [`AI_HYGIENE.md`](AI_HYGIENE.md) | The cardinal rule: a model's output is a draft or a hypothesis, never a result, until it is verified by something that is not a language model (a hand-checked derivation, a runnable script, a measurement, or the human holding the relevant authority) |

Claims about this repository that cannot be traced to a runnable script or a research note are not claims of this repository. Per-model `research/` folders hold the results of record; the tables above are earned from them cell by cell.

**Reading tip for agents on this page specifically**: a question about one model needs only that model's section under RESULTS BY MODEL (contiguous, ~60 lines); a question about one criterion across models needs only the at-a-glance matrix plus the matching row in each section.
