# Observer Patch Holography

**Physics from the records observers can share.**

[Physics preprint](https://philpapers.org/rec/MUEFOC) · [Pragma Research](https://floatingpragma.io/research/) · [Cadence](https://floatingpragma.io/cadence/) · [Français](README_FR.md)

Observer Patch Holography (OPH) investigates whether familiar physics can be
reconstructed from bounded observers that reach agreement. An observer patch
has local state, a boundary with ports, readback, records, and feedback for
repairing disagreement. Within the model, facts become public when records
survive comparison across overlapping patches.

This repository contains the scientific work behind that idea: papers,
machine-checked proofs and executable models. It connects finite agreement
with quantum probability, spatial geometry and the symmetry pattern of the
known forces, with a conditional route to Einstein dynamics. Each construction
states the assumptions it consumes. Establishing a common physical
realization requires measurements that distinguish the proposed models.

## Start here

| To explore… | Start with… |
| --- | --- |
| The main physics argument | [*Finite Observer Consensus as a Reconstruction Principle*](https://philpapers.org/rec/MUEFOC) |
| The proofs and reproducible calculations | [Lean library](Lean/) and [reproduction guide](REPRODUCE.md) |
| The related work on learning machines | [Cadence](https://floatingpragma.io/cadence/) and [interactive demos](https://floatingpragma.io/demos/) |
| Explanations for a general audience | [Pragma Research blog](https://blog.floatingpragma.io/) |

## One Architecture, All Of Physics

The research program connects these parts of physics through one observer
architecture, with separate assumptions stated for each result.

- **Agreement between observers.** Under the stated termination and
  consistency conditions, local repairs reach a protected public record
  independent of their order. The finite theorems describe agreement,
  stability and refinement.
- **Quantum probability.** Finite event algebras support the Born probability
  rule, Lüders conditioning and exact Tsirelson results, including a Bell-state
  construction that attains 2√2. These are mathematical results on the stated
  event-algebra branch.
- **Geometry and the known forces.** A twelve-port response construction
  supplies a three-dimensional spatial readback. Complete reversible response
  and internal observer transport force the Standard Model gauge Lie type
  under the stated hypotheses. A supplied matter representation gives an
  anomaly-free fifteen-state generation; its physical realization and the
  global gauge group require additional structure.
- **Spacetime and gravity.** Declared source and reading laws support
  controlled causal-geometry constructions. The Einstein-equation implication
  assumes a common physical realization with the stated stress, entropy,
  continuum and scale data.
- **Classical and quantum matter.** A supplied charged-scalar/Maxwell action
  supports controlled nonlinear continuum motion in its real sector and an
  interacting quantum state space on a fixed mesh. Selecting the action from
  observer histories and identifying it with physical matter are separate
  requirements.
- **Constants as fixed-point problems.** Koide's relation holds exactly under
  a stated balance premise. The fine-structure calculation certifies a root
  of a declared closure map and carries diagnostic status. The capacity
  program asks whether the public capacity assigned to the universe agrees
  with the capacity reconstructed from within it. Connecting these
  constructions to measured constants requires physical identification.

The [main preprint](https://philpapers.org/rec/MUEFOC) develops these connections
and their precise scope. The [paper index](paper/) provides the specialist
accounts, including thermodynamics, field geometry and fixed-point
constructions for the constants.

<!-- PUBLIC-QUANTITATIVE-CLAIMS:BEGIN -->
<!-- Quantitative table suppressed while physical_establishment count is zero. -->
<!-- PUBLIC-QUANTITATIVE-CLAIMS:END -->

## Evidence you can inspect

The [Lean library](Lean/) contains more than 12300 public theorems and lemmas
with no admitted proofs. Exact certificates and reproducible simulations
accompany the mathematical arguments. The companion
[physics simulator](https://github.com/muellerberndt/oph-physics-sim) provides
executable observer dynamics and retained evidence.

The [axiom reference](docs/AXIOM_REFERENCE.md) states the three core axioms;
the [premise register](docs/PREMISE_REGISTER_V3.md) and
[claim registry](claims/claim_registry.yaml) record additional premises and
result-specific assumptions, including physical identifications and empirical
inputs.
The [postdiction ledger](docs/POSTDICTION_LEDGER.md) records comparisons with
measured values and their input ancestry; the
[frozen-prediction ladder](docs/FROZEN_PREDICTION_LADDER.md) records tests whose
conditions must be fixed before the comparison data are examined.
The [falsification program](docs/OPH_FALSIFICATION_PROGRAM.md) gives the
observations that would refute particular claims.

### Reproduce the finite core

After setting up the dependencies in [REPRODUCE.md](REPRODUCE.md), check the
claim graph and selected finite algebra, record-capacity and consensus results:

```bash
python3 tools/check_claim_registry.py
python3 -m pytest -q \
  code/a5_closure/test_audit.py \
  code/capacity_readback/test_correctable_public_record_capacity.py \
  code/capacity_readback/test_reversible_public_checkpoint_packet.py \
  code/consensus/test_reference_architecture_benchmark_suite.py \
  code/consensus/test_verified_tree_packet_net.py
```

The reproduction guide describes the wider checks and the scope of each
evidence family.

## Reconstruction map

<p align="center">
  <a href="assets/prediction-chain.svg" target="_blank" rel="noopener noreferrer">
    <img src="assets/prediction-chain.svg" alt="OPH reconstruction chain" width="92%">
  </a>
</p>

<p align="center"><sub>The OPH reconstruction map connects observer records, three-dimensional source geometry, causal order, clocks, fields and quantum states. Each arrow names a mathematical connection; the papers state the assumptions that join them into an effective physical description.</sub></p>

## Repository guide

| Path | Contents |
| --- | --- |
| [`flagship/`](flagship/) | Main standalone physics paper, TeX source and PDF |
| [`paper/`](paper/) | Core papers and publication index |
| [`Lean/`](Lean/) | Machine-checked mathematical development |
| [`code/`](code/) and [`evidence/`](evidence/) | Executable models, certificates and reproduction evidence |
| [`extra/`](extra/) and [`cosmology/`](cosmology/) | Focused mathematical and physical research |
| [`book/`](book/) | *Reverse Engineering Reality*, source and downloadable book |
| [`docs/`](docs/) | Reader policies and scientific ledgers |

## Contribute

OPH welcomes proofs, counterexamples, simulations, independent reviews and
readable explanations. Start with the [reproduction guide](REPRODUCE.md) and
the [selection ledger](docs/SELECTION_LEDGER.md), which states the premises
and scientific boundaries relevant to contributions.

## From observer patches to learning machines

[Cadence](https://github.com/muellerberndt/cadence) makes the settling idea
executable as a learning architecture. Its bounded software patches carry
local state and readback, with feedback that repairs prediction errors.
Its default brain, System 1, learns from experience as animal brains do and
not by backpropagation. Optional System 2 observers add recursive feedback
inside the same settlement.

The [Cadence preprint](https://philpapers.org/rec/MUECAP-2) describes the
architecture and experiments. Its computational case rests on learning
behavior and measured resource use. [Pragma Research](https://floatingpragma.io/)
connects this work to embodied AI.

## License

The repository uses split licensing. All software, including the Lean library, [`code/`](code), [`tools/`](tools), and the schemas in [`schemas/`](schemas), is licensed under [Apache-2.0](code/LICENSE). Papers, the book, documentation, figures, data, the generated ledgers in [`tracking/`](tracking), and the packaged particle data in [`pdg_data/`](pdg_data) are licensed under [CC BY-NC-SA 4.0](LICENSE); the tabulated values in `pdg_data/` carry their upstream Particle Data Group terms. The [LICENSE](LICENSE) file gives the per-directory map.
