# Observer Patch Holography

**Physics from the records observers can share.**

[Main paper](https://philpapers.org/rec/MUEFOC) · [Book](https://oph-book.floatingpragma.ai/) · [Notebook](https://notebook.google.com/notebook/d5249760-6ce8-44a0-927b-ccf90402711a) · [Lean library](Lean/) · [Reproduce](REPRODUCE.md) · [Pragma Research](https://floatingpragma.ai/) · [Français](README_FR.md)

[![Mandatory Suite](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/mandatory-suite.yml)
[![Lean CI](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/lean-ci.yml)
[![Runnable observer model](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml/badge.svg)](https://github.com/FloatingPragma/observer-patch-holography/actions/workflows/observer-model.yml)

## Mission

The mission of Observer Patch Holography (OPH) is to create a consistent and
complete, formal and rigorous theory of reality. That theory has no starting
point: it is a closed loop, and it has to meet four conditions.

1. **Consistent.** The loop closes. Its fixed point exists, constructed
   explicitly and machine-checked, and a theory with a model contains no
   contradiction.
2. **Rigid.** Every dimensionless quantity is forced. One measured value sets
   the units and nothing else.
3. **Complete**, in the sense
   [Einstein, Podolsky and Rosen](https://doi.org/10.1103/PhysRev.47.777) gave
   the word in 1935. Every element of physical reality has a counterpart in
   the theory, and every measured quantity is derived.
4. **Correct.** Derived values agree with measurement, and new predictions are
   frozen before they are tested.

OPH does not claim to be the only theory with these properties.

## The idea

OPH proposes that physics is what bounded
observers can establish together. Finite observers compare records on shared
boundaries, repair their disagreements, and treat as real whatever survives
every comparison. Three problems that resist a classical account motivate the
programme: quantum measurement, correlations that no local classical model
reproduces, and the origin of an apparent spacetime. OPH seeks one explanation
for all three in the consistency of observers' records and operations, and it
pursues that explanation as an executable observer, machine-checked
mathematics and predictions fixed in advance.

**The Universe reads itself into existence.** OPH has one hypothesis, the
closure hypothesis: the Universe is a timeless, self-referential structure
whose observers recover the specification of their own world from their
records and build hardware that realizes that same history, including its
builders. Nothing lies outside this loop, so the recovered specification owes
no further justification. It has to run and reproduce what observers measure.
Space, time and matter are effective descriptions of the records, and "our
Universe" is the public record, what survives every comparison between
observers. The [basis](docs/BASIS.md) states the hypothesis, the no-outside
rule, the agreement condition and the executable specification; the
[working terms](docs/INTERPRETATION_AND_TERMS.md) explain the timeless view
and the clocks inside it.

## The executable hypothesis

The concrete bet is a program. One simple local repair and memory law,
running on a network of twelve-port carriers, should produce observers whose
shared records look like our world: interacting quantum fields in 3+1
dimensions, and gravity on large scales. Particles are the knots in that
network: protected defects that local repairs cannot undo.

The [hypothesis](docs/OPH_EXECUTABLE_HYPOTHESIS.md) says
what would count as success, and the
[completion program](docs/COMPLETION_PROGRAM.md) splits the work into gates
with stopping rules. Either the pieces join up on one running process, or a
proof shows that the whole class of laws cannot work.

## Highlights

- **The particle families are identified.** Every protected defect carries a
  six-valued flux class, and it is exactly the class the Standard Model gives
  its particles: left- and right-handed quarks, leptons with the Higgs, and
  neutral fields. Quarks and leptons come out as fermions, and no colour
  charge ever sits alone.
  [Particle paper](paper/deriving_the_particle_zoo_from_observer_consistency.pdf).
- **The charged-lepton mass ratios.** A generation register in the same repair
  law gives muon/electron 206.770 and tau/electron 3477.47 (measured 206.768
  and 3477.37), with no measured mass as input.
- **An early universe like ours.** The runs create every particle class in
  pairs, quark pairs annihilate, colour stays confined, and matter and
  antimatter come out equal, as before baryogenesis.
- **The gauge group with equal couplings.** Phases on the carrier's ports give
  (SU(3)×SU(2)×U(1))/Z6 with equal strengths at the cell scale, and neutrinos
  come out as Dirac particles.
- **An observer that runs.** A twelve-port observer keeps its records and
  uses them to choose its next repair. Overlapping observers read a common
  reference and its clock, a finite quantum process joins electric flux,
  charged matter and a scalar field, and refining screens give Gaussian
  field limits and waves of finite speed.
- **A Lorentzian scalar field from the observers' records.** On a declared
  radial extension of the refining screen, one action fixes correlations,
  response and a modular clock, and a smooth (1+3)-dimensional scalar field
  theory reproduces any finite set of joint readouts. The radial direction
  and the physical units are inputs.
  [Integrated verdict](model/CONSTRUCTION.md#integrated-effective-description).
- **Standard physics from finite agreement.** Under stated premises, finite
  theorems give public records that do not depend on the order of repairs,
  the Born rule and the Tsirelson bound, a carrier that is three-dimensional
  in its slow response, the Standard Model gauge Lie algebra with an
  anomaly-free matter representation, a (1+3)-dimensional causal order and the Einstein equation from null
  balance.
- **Energy readouts.** Complete interacting quantum sources give
  gauge-constrained energy records, including preparation failures, and a
  finite electric/magnetic field has a certified charged ground.
  [Charged-field attachment](model/native_abelian_field_join/PROOF.md).
- **Continuum fields and their measurements.** The ordered-spin construction has
  three continuous spatial dimensions, time, positive energy and finite propagation.
  Its position observable is conjugate to the same spatial momentum, with complete
  finite detectors and ordered source limits. Spectral measurements retain all outcomes.
  Proving Lorentz symmetry and surviving interactions on that same state is the
  central field requirement. A separate scaling gives free massless Weyl fields
  and the matter stress response used by the conditional gravity derivation.
  [Field constructions](docs/research/OPH_FIELD_RECONSTRUCTION.md#the-three-field-constructions-and-their-source-boundaries),
  [spectral readout](model/spin_shell_readout/README.md),
  [ten common-theory requirements](model/particle_joint_reconstruction/DEPENDENCIES.md).
- **Interacting fields and actual position records.** A declared
  [gauge-field construction](model/particle_observer_outputs/PROOF.md)
  has controlled local detector records and a causal bound uniform in
  physical refinement, including quantum links and all fermion occupations. A
  [joint field/reference instrument](model/particle_joint_translation_readout/PROOF.md)
  supplies conserved crystal momentum, position backaction and finite
  uncertainty. Physical chiral refinement and particle mass identification
  require their own same-state derivations. The
  [illustrated atlas](docs/research/OPH_PARTICLE_ATLAS_V2.html) explains
  these finite measurements and the remaining physical connections.
- **Thermal records, work and pressure.** One interacting quantum-flux
  control has a [certified thermal preparation](model/particle_field_state_response/PROOF.md)
  and complete two-time records, with every occupation and reference retained.
  Its [spatial strain response](model/particle_geometric_response/PROOF.md)
  calculates work and pressure from a declared geometry coupling on the
  same state and clock. The proper geometric interpretation and regional
  entropy and area premises required for Einstein dynamics are separate
  physical requirements.
- **One source, state, reference and clock.** The selected Wilson construction
  has a [common observable algebra](model/particle_primary_field_algebra/PROOF.md),
  a stationary limiting state, positive energy and convergent finite detector
  histories on one subsequence. Its energy currents obey exact finite
  conservation. A [controlled energy pulse](model/particle_primary_energy_response/PROOF.md)
  produces a nonzero original current response and positive cyclic work;
  its Ward identity includes the work supplied by the pulse.
  [Selected energy moments](model/particle_primary_weighted_spectral_tail/PROOF.md)
  obey a bound of 4000000/a, independent of link degree, where a is the source spacing.
  The positive finite-window probability deteriorates with capacity. An
  [equal-work control](model/particle_primary_equal_work_control/PROOF.md)
  on the same source has an elastic limit, so mean work alone cannot prove
  particle survival.
  The [native chirality audit](model/particle_primary_odd_sum_rules/PROOF.md)
  identifies a paired response that a physical chiral assignment must resolve.
  Lorentz covariance, surviving interacting particles and proper geometric
  stress remain requirements on this same source and state. The
  [ten-step dependency map](model/particle_joint_reconstruction/DEPENDENCIES.md)
  records the precise remaining physical requirements.
- **Mathematics you can check.** The [Lean library](Lean/) contains more than
  12900 public theorems and lemmas with no admitted proofs. Computational
  receipts bind their source inputs and include replay checks.
- **Tests fixed in advance.** Predictions and their rejection rules enter the
  [frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) before the
  eligible comparison data, with OpenTimestamps proofs, most of them anchored in
  Bitcoin; the
  dispersion prediction for the primitive port orbit is one of them. The
  [postdiction ledger](docs/POSTDICTION_LEDGER.md), the
  [falsification programme](docs/OPH_FALSIFICATION_PROGRAM.md) and the
  [cross-scale audit](docs/EMPIRICAL_STATUS.md) keep comparisons with
  measurements apart from tests of the mathematics.

**Agreement with the full range of observed measurement statistics has not
been established.** The theorems hold under stated premises, and the
executable observer is a candidate judged by declared acceptance criteria.

## Run the observer

From the repository root, with Python 3.11 or newer:

```sh
python3 -m model.observer run --trials 10000 --seed 42
```

This runs the minimal record-and-feedback experiment and prints forecast
scores and a result hash. The [model guide](model/README.md) covers the rest.

## Mathematical results

- **Agreement and stable records.** Under explicit termination and consistency
  conditions, local repairs reach one shared public result whatever their order.
  [Consensus paper](paper/reality_as_consensus_protocol.pdf).
- **Quantum observations.** Within the stated finite quantum algebra, verified
  results cover repeatable measurements, classical records that later
  operations cannot disturb, uncertainty and the Bell correlation bound, with
  a state that attains it. Native repair histories reconstruct a positive
  functional on the full 66-state matrix algebra, where the Born trace rule
  holds; native evolution does not prepare that state.
  [Finite quantum event algebras](extra/machine_checked_finite_event_algebras.pdf),
  [native reconstruction](model/native_quantum_reconstruction/PROOF.md).
- **Three spatial directions.** Optimal readout of a supplied three-dimensional
  response selects the twelve-port incidence, and the carrier's mean repair rule
  gives an exact three-dimensional slow response with Euclidean geometry, which
  an observer's current records recover.
  [Carrier construction](code/a5_closure/README.md).
- **Gauge symmetry.** Under stated premises on complete reversible response
  and observer transport, the twelve-port response algebra has the Standard
  Model gauge Lie type, and a specified matter representation gives an
  anomaly-free generation. Port phases realize the global group, whose centre
  is the flux alphabet of the particle defects.
  [Gauge derivation](paper/deriving_standard_model_gauge_structure_from_observer_overlap_consistency.pdf).
- **Thermodynamics.** Finite repair theorems give a Clausius inequality,
  fluctuation relations and an information Landauer inequality under a stated
  reference and energy calibration. Fair-coin seam repair is exact conditional
  resampling, so entropy never falls on binary populations while the complete
  tape keeps every bit. A native memory refresh realizes the machine-checked
  finite four-law core in internal cost units.
  [Thermodynamics account](paper/observers_are_all_you_need.pdf).
- **Information capacity.** Exact results count which records stay
  distinguishable after the allowed operations, bound that count by the
  available carrier, and follow it through refinement.
  [Record capacity](code/capacity_readback/README.md).
- **Entangled-pair coupling candidate.**
  [Alex Wissner-Gross](https://alexwg.org/publications/ChainTheoryAlpha.pdf)
  derives a finite spectral expression with inverse `137.035999189566`.
  The [independent OPH reconstruction](model/observer_particles/README.md)
  credits that construction and tests its complete-outcome and normalization
  assumptions. Its identification with the physical electromagnetic coupling
  and its proposed attachment to the screen-grain equation are conditional.
- **Matter and fields.** A supplied charged-scalar and Maxwell action with its
  geometry and quantization give an interacting Hamiltonian after gauge
  reduction, and controlled scalar refinement connects classical detector
  readouts to a nonlinear field equation. Gauge invariance, face locality and
  carrier rotations select the Maxwell action on the seams.
  A common metric dictionary is constructed.
  [Matter continuum](paper/observers_are_all_you_need.pdf).

These are conditional mathematical results. Making their premises hold for
the observations of one common model is the central construction problem.
The [premise map](model/PREMISES.md) records which premises the executable
specification supplies and which it takes as inputs.

## Papers

- **[Observers Who Read Their Universe into Existence](flagship/from_observer_consensus_to_standard_physics.pdf)**
  is the main synthesis: the closure hypothesis, the executable-model hypothesis and the tower of results built on it. Cite its released version, the
  [PhilPapers record](https://philpapers.org/rec/MUEFOC).
- **[An Executable Observer for Observer Patch Holography](paper/executable_observer_construction.pdf)**
  states the specification, the executed components and the ordered joins
  that remain.
- **The theorem library** holds the long proofs, from
  [consensus](paper/reality_as_consensus_protocol.pdf) and
  [observers and records](paper/observers_are_all_you_need.pdf) to
  [spacetime and Einstein dynamics](paper/recovering_observer_spacetime_and_einstein_dynamics_from_overlap_consistency.pdf)
  and [particles](paper/deriving_the_particle_zoo_from_observer_consistency.pdf).
  The [paper index](paper/README.md) and the [paper set](docs/PAPER_SET.md)
  list every paper and its role.
- **Companion papers** carry focused results. A lepton-mass relation, the
  Koide relation, holds exactly under a stated balance premise. The capacity
  programme asks whether the public capacity assigned to the Universe agrees
  with the one reconstructed from inside. These constant calculations,
  including the fine-structure root, have diagnostic status.

## Where things are

| Path | What it holds |
| --- | --- |
| [model/](model/README.md) | The executable observer: implementation, tests, retained receipts and protocol records |
| [flagship/](flagship/) | The synthesis paper |
| [paper/](paper/README.md) | The construction paper and the theorem library, with archived papers under [paper/archive/](paper/archive/README.md) |
| [extra/](extra/README.md), [cosmology/](cosmology/) | Companion papers |
| [Lean/](Lean/) | The machine-checked library |
| [code/](code/README.md) | Exact calculations, lane verifiers and receipts behind the papers |
| [evidence/](evidence/README.md) | Retained evidence packets and custody records |
| [claims/](claims/README.md) | Established claims and their dependency views |
| [tracking/](tracking/) | Premise, observation, prediction and instrument registers |
| [docs/](docs/README.md) | The basis, the working terms, policies and registers |
| [book/](book/) | The book |

## Reproduce and verify

The [reproduction guide](REPRODUCE.md) runs the verification suite from a
clean clone. The model verifier replays the receipt inventory declared in `model/verify.py`;
the regression suite also checks adversarial cases:

```sh
python -m pip install -r model/requirements.txt
python -m model.verify
python -m pytest -q -W error model
```

The guide also lists the standalone particle and source-field checks used by
the [observer workflow](.github/workflows/observer-model.yml).

The [Lean library](Lean/) builds with Lake, the [documentation index](docs/README.md)
lists the registers and policies, and the [book](book/) gives a longer road
in. The textbooks at [learn.floatingpragma.ai](https://learn.floatingpragma.ai/)
and the interactive simulation at
[simulation.floatingpragma.io](https://simulation.floatingpragma.io/) are
the gentlest introductions.

## Contributing

The task is to make one observation-and-repair process supply the effective
physics: the observer's clock, spatial interpretation and field measurements
on the same histories, compared with independent observations.
Contributions should resolve a named connection or expose a failure in the
shared model; the [construction guide](model/CONSTRUCTION.md), the
[repair-law guide](model/CANDIDATES.md) and the
[research issues](https://github.com/FloatingPragma/observer-patch-holography/issues?q=is%3Aissue+is%3Aopen+label%3Aprogram%3Aphysics)
set out that work. Paper releases follow the
[release policy](model/CONSTRUCTION.md#research-programme-and-paper-release-hold);
main-branch PDFs are source previews, and deposited releases retain their
original bytes.

## License

Software, including Lean and the model code, uses [Apache-2.0](code/LICENSE).
Papers, documentation, figures, data and ledgers use [CC BY-NC-SA 4.0](LICENSE),
with upstream Particle Data Group terms for the packaged particle values.
The [license map](LICENSE) gives the directory-specific terms.

<details>
<summary>Related projects and funding</summary>

[Cadence](https://github.com/muellerberndt/cadence) is the related learning-machine
project. Its capability evidence is separate from OPH's physics evidence.
[Pragma Research](https://floatingpragma.ai/) brings the projects together.

OPH's official TEO (Theory Of Everything) memecoin was deployed by the Base
community and lives on Base at
`0xd2a7055317D0c7b316319cDadaE592D0644A0b07`
([chart](https://dexscreener.com/base/0xba55cdae026b9dd63fcb50e82094d61472948e44527ca4791033497469228049)).
Its trading fees are routed to Bernhard Mueller to fund the project.
The token has no role in the theory. Other tokens using the OPH or TEO name
are unrelated. Nothing here is financial advice.

</details>

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->
