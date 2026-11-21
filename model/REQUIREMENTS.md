# Requirements for the native repair law

**Status:** open research problem\
**Main candidate:** [contrast_memory_spherical_v1](candidate_laws/contrast_memory_spherical_v1/README.md)\
**Candidates and their evaluation:** [CANDIDATE.md](candidate_laws/CANDIDATE.md)

## The criterion

The correct law is one executable repair and memory law that runs on the fixed spherical echosahedral federation
and fulfils the premises of the finite OPH theorems together, on one process. Nothing more is asked of it. The
finite theorems say how three directions of space, the rules of quantum measurement, the gauge structure,
thermodynamics and Einstein gravity come about once their premises hold, so a law does not have to produce them
separately.

| | |
| --- | --- |
| Fixed | the carrier: the spherical echosahedral federation |
| Fixed | the premises of the finite theorems |
| Free | every setting a designer can implement in the simulating universe |
| Forbidden | hardcoding what the theorems derive or what observers measure |

## 1. The fixed carrier

- **One echosahedron is one native observer.** Its twelve ports are its interfaces.
- Each echosahedron keeps its icosahedral incidence: 12 ports, 30 internal seams, 20 faces and the sixty rotations.
- The carrier is the spherical federation of 20, 80 and 320 cells, with its inter-cell seams and its refinement
  structure. [executable_model/](executable_model/README.md) implements it for every candidate.
- Larger observers are finite, possibly overlapping patches of echosahedra. Overlapping observers share registers
  and event identities, and shared records agree without duplication.
- The carrier is computational architecture and no physical space. No lattice, number of space dimensions, global
  clock, radial direction or spacetime metric enters the law.

## 2. The fixed premises

The premises of the finite theorems are the rows of the premise register
([PREMISE_REGISTER_V3.md](../docs/registers/PREMISE_REGISTER_V3.md), machine form
[premise_register.json](../tracking/premise_register.json)). For each premise the register says what a
construction must supply to discharge it. A law fulfils a premise when the same equipped federation satisfies it,
with the information or action at the interface the theorem names: when a conclusion concerns a simulated
observer's history, clock, causal parent or current, that observer can read it.

| Consequence | What the law and its design supply | What the theorem then fixes |
| --- | --- | --- |
| A common public result | a physical quotient, an accepted repair relation that descends and has joinable local forks, and equivalence of normal forms and consistency | one public normal form from each initial state |
| Three directions of space | the oriented native cell with its balanced binary sector and the equal-rate repair mean; the one-load law realizes the same eigenspaces | the slowest nonconstant response has rank three; a continuum description of patch records further uses a population of records, their authenticated order, local clocks and a compatible refinement family |
| The gauge Lie algebra | a faithful complete response that spans the twelve port directions, closes under ordered-response commutators, and implements proper rotations by accepted overlap histories | the compact Lie type $\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$; the global group and its action on matter use further premises |
| Quantum probabilities | a finite effect space with normalized positive record weights, one weight for every implementation of the same effect, and additivity for compatible effects | a unique density matrix and the Born rule |
| Einstein response | on one source, state and refinement family: an operational Lorentzian geometry, modular charges that measure the stress used in propagation, its Ward identity, the entropy-area normalization and fixed-volume stationarity, controlled small-region errors, enough null directions, and universal coupling with a common vacuum and a positive Newton constant | null balance with the Ward and Bianchi identities gives the semiclassical Einstein equation for that stress |
| Matter classes | a protected face-flux record on the five-level code; for the fermion grading also a realized spin doublet, a Gauss sector and antisymmetric occupation | the centre classes of the Standard Model; several species can share a class |
| The least-informative assignment | a feasible set, a faithful reference, positive weights and a cover that distinguishes feasible states, all with native provenance | an information projection within that family |

Existing theorems may be generalized or revised. Every such change is explicit, justified and formally verified.

### Together

The premises are fulfilled together, by one unchanged law on one process. Where two theorems refer to the same
object, the law supplies one record, one clock and one response for both. In the Einstein theorem the stress that
an observer infers from its modular flow is the stress carried by the excitations that propagate and interact. The
gauge currents act on those same excitations through the law's operations, and the detector, generator, clock and
currents that identify an excitation refer to one object. A concrete target follows: finite patches use their own
clocks and records to measure both the current response and the stress of one bound excitation.

## 3. Free: calibration

A designer may calibrate anything that can be implemented in the simulating universe, so that the observers'
records match ours: the repair and memory rule, couplings, preparations, feedback programs, service rates,
references, schedules, readouts and the depth of refinement. The test is implementability with finite local
resources on the federation. Every setting is declared with its value and provenance. A setting chosen from an
observation names that observation, and agreement with it is no independent prediction.

## 4. Forbidden: hardcoding outcomes

Nothing that the theorems derive or the observers measure is put in by hand: no particle masses or mass ratios,
no number of space dimensions, no spatial lattice or propagation graph, no spacetime metric, no global clock, no
gauge group and no named current matrices. Clock and wave behaviour obtained on a supplied spatial architecture
counts as hardcoding. It is at most a phenomenological benchmark ([benchmarks/](benchmarks/README.md)).

## 5. The prediction

Once a law fulfils the premises together, it is frozen and its observers' records are read: which bound objects
form, with which spins, charges and mass ratios. The matter theorem fixes the centre classes. The individual
spectrum, the representations and the dimensionless mass ratios are predictions of the frozen design wherever no
further theorem implies them. If the frozen law yields the observed spins, charges and mass ratios without using
them as inputs, that is independent evidence for identifying it with our Universe. A candidate is not discarded
for a particle its observers have not identified while its premises are open: the census of the main candidate
finds a bound pair and identifies no electron, and that decides nothing about the candidate.

## 6. Evaluation

`python -m model evaluate` scores every candidate in [candidate_laws/](candidate_laws/CANDIDATE.md):

- carrier checks on the candidate's own generator: the native cell, locality, the sixty rotations, the rank-three
  slowest response and the federations;
- the parameter check: every parameter declares its provenance, and none hardcodes an outcome;
- the premise matrix: every row of the premise register with the candidate's status (fulfilled, partial, open or
  fails) and its evidence; a failed premise names its counterexample;
- the joint requirements of the section "Together".

The [scorecard](evaluation/SCORECARD.md) collects the results, and [the evaluation](evaluation/README.md) says how a
candidate takes part. Each candidate delivers an executable specification with its permitted inputs, state space,
instruments and record semantics; the premise matrix with its evidence; formal proofs or reproducible finite
certificates; minimal counterexamples for failed premises; runs on single cells, overlapping patches and the
federations; and a list of its calibrations and open assumptions.

## 7. Search method

The search runs over native local operations built from the incidence of ports, seams and faces, the
contrast-memory and seam-repair operations, finite seam-current and record-memory registers, symmetry-equivariant
coherent interactions, record-conditioned bounded-support instruments, and the slow and fast spectral parts of the
cell. Promising hypotheses include conditional-expectation repair that preserves a dynamical sector, gauge brackets
synthesized from the incidence, wave generators derived from the incidence, and record populations selected jointly
by the slow and fast triplets. These are candidate mechanisms, never established features of the final law.

The search is counterexample-guided: generate a law, evaluate it against every applicable premise, extract the
smallest failing example, and exclude the failed construction from later searches.

## Success criterion

The decisive milestone is one law on the fixed federation that fulfils every premise the register assigns to the
law, jointly and with every calibration declared. Its frozen observers then recover the physical history and laws
of our Universe from their own records, without spacetime, forces or particles supplied in advance.

**Core principle: derive one universe from one native observer law, never a collection of individually successful
constructions joined by unproved assumptions.**
