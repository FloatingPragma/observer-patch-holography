# OPH Falsification Program

## Scope

This program covers completed, non-cosmological OPH claim surfaces with stable
hypotheses and comparison objects. It contains no frozen-target registry,
preregistration queue, forward scorecard, or aggregate hit count.

The eligible sources are:

- the finite quotient-consensus theorems in
  [*Reality as a Consensus Protocol*](../paper/reality_as_consensus_protocol.pdf);
- the D1--D5 and D7--D9 recovered-core claims in the
  [compact SM/GR paper](../paper/recovering_relativity_and_standard_model_structure_from_observer_overlap_consistency_compact.pdf);
- the exact gauge, charge, color, generation, and declared-action carrier
  claims in
  [*Deriving the Particle Zoo from Observer Consistency*](../paper/deriving_the_particle_zoo_from_observer_consistency.pdf);
- the finite carrier and screen-sieve theorems in
  [*Federated Echosahedral Screen Microphysics*](../paper/screen_microphysics_and_observer_synchronization.pdf);
- the substrate-neutral
  [normal-form](../extra/observable_normal_forms.pdf) and
  [finite-event-algebra](../extra/machine_checked_finite_event_algebras.pdf)
  supplements.

A claim belongs here only when its statement is complete, its hypotheses are
explicit, the falsifier addresses the same object, and no open map or
measurement convention can absorb the result.

## Mature Mathematical Falsifiers

| Claim surface | Result that falsifies the claim | Scope |
| --- | --- | --- |
| Finite quotient consensus | A finite repair system satisfying the stated termination, local-diamond, boundary-protection, and quotient hypotheses has two inequivalent terminal quotient normal forms from the same source. | The finite consensus theorem falls. |
| Protected-observation invariance | Two quotient-equivalent implementations satisfying the declared interface and repair rules produce different protected readouts. | The carrier-invariance claim falls. |
| Cross-source identification | The stated protected-observation criterion identifies two sources whose consistent quotient normal forms remain distinct. | The cross-source theorem falls; same-source confluence is separate. |
| Collar recovery | A faithful finite-range family satisfying the stated uniform conditional/matrix-mixing and boundary-growth hypotheses violates the published recovery bound or its scaling limit. | The conditional collar-recovery theorem falls. No stress or dark-source claim is attached. |
| Lorentz cap branch | A receipt-complete prime cap pair satisfying the compact paper's incidence, modular, regularity, and scaling hypotheses fails the stated Lorentz action or commutator relations. | D3--D4 fall on that branch. |
| Einstein bridge | A model satisfying the compact paper's null-stress, stationarity, tensor-upgrade, conservation, and physical-identification hypotheses fails the stated Einstein equation modulo one metric term. | D5 falls. The cosmological value of the metric term is outside this program. |
| Compact-gauge reconstruction | A refinement system satisfying the declared sector, transport, obstruction, and compactness receipts fails to produce the stated compact gauge quotient. | D7 falls. |
| Realized Standard Model quotient | The displayed quotient, hypercharge assignments, anomaly cancellations, shared center, color count, generation count, or one-Higgs admissibility calculation is mathematically inconsistent under its stated realized-branch premises. | D8--D9 fall in the affected part. |
| Classical carrier modes | The declared Maxwell, perturbative pure Yang--Mills, or pure Einstein--Hilbert quadratic action satisfies its stated vacuum and gauge hypotheses but lacks the published classical transverse mode count. | The corresponding classical carrier theorem falls. No quantum-particle claim follows from this row. |
| Finite screen and $A_5$ coefficient sieve | A finite carrier satisfying the screen-sieve premises violates the stated port, reserve, phase, or quotient conclusion; or the declared twelve-port $A_5$ module, compact bracket, inner-action classification, or six-axis lattice calculation is algebraically inconsistent under its stated premises. | The affected finite screen or coefficient theorem falls. Physical gauge identification remains a separate receipt question. |

## Mature Physical Branch Falsifiers

| Claim surface | Result that falsifies the claim | Scope |
| --- | --- | --- |
| Charge lattice | A confirmed elementary state on the realized low-energy branch has electric or hypercharge quantum numbers outside the derived quotient lattice. | The realized charge-quantization branch falls. |
| Color count | A confirmed light Standard Model color carrier requires a fundamental color multiplicity other than three on the claimed realized branch. | The realized $N_c=3$ branch falls. |
| Chiral generation count | A confirmed fourth light chiral Standard Model generation satisfies the same low-energy branch conditions. | The realized $N_g=3$ branch falls. |
| One-Higgs minimal branch | The observed low-energy branch requires an additional elementary Higgs multiplet while satisfying the same minimal-admissibility premises. | The realized one-Higgs MAR branch falls. |
| Gauge quotient | Under the physical-current, trace-balanced block, tensor-kernel, axis-descent, and MAR premises, a confirmed low-energy global form is incompatible with $(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$. | The realized gauge-quotient branch falls. |

These physical rows test the named realized branch. They do not erase finite
consensus, normal-form, or unrelated recovery theorems.

## Ineligible Surfaces

The following surfaces are not part of this falsification program because the
physical map, source derivation, comparison object, or paper-level closure is
work in progress:

- all cosmology: capacity, dark energy, inflation replacement, scalar tilt,
  CMB spectra, low-$\ell$, birefringence, $H_0$, $S_8$, growth, dark abundance,
  cosmological neutrino rows, and structure formation;
- the repair-charge condensate dark-sector continuation, including galaxy
  acceleration scales, Solar-System response, lensing, clusters, abundance,
  and laboratory force claims;
- quantitative fine-structure, hadronic transport, W/Z, Higgs/top, quark,
  charged-lepton, and neutrino mass continuations whose physical readout or
  source selector is incomplete;
- black-hole spectroscopy, entropy, Page curves, and ringdown templates;
- coherent-matter, anti-gravity, hardware enrichment, cognition, string-vacuum,
  and meaning-layer continuations;
- Yang--Mills mass-gap claims whose continuum, reflection-positivity,
  transfer/intertwiner, or nontriviality receipts are open.

Diagnostic comparisons, simulator outputs, known-data checks, and engineering
benchmarks carry no falsification verdict.

## Operational Use

Each proposed falsifier must identify the claim row, the owner paper, and the
receipt set used to establish its premises. The machine-readable row mapping is
[`claims/falsification_matrix.csv`](../claims/falsification_matrix.csv), and
the current claim status is in
[`claims/claim_registry.yaml`](../claims/claim_registry.yaml). Open physical
maps are tracked as [closure issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure). A result on an
ineligible surface is evidence for or against that continuation only; it is not
silently promoted into this program.

## Decision Rule

A falsifier binds only when every hypothesis of the named claim is verified
and the result addresses the same mathematical or physical observable in the
same convention. Missing hypotheses yield no verdict. Counterexamples remain
scoped to the stated claim and do not transfer to independent theorem rows.

Claim status and unclosed maps are recorded on the
[closure issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+label%3Aclosure) and in
[`PROOF_SPINE.md`](PROOF_SPINE.md).
