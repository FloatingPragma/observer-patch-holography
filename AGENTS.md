# Repository guidance for agents

Research results belong in the papers, the Lean library and the code. The
[writing style guide](docs/policies/STYLE_GUIDE.md) binds all prose.

## Construction before uniqueness

The mission is: **Build an executable realization of the OPH observer algebra
whose compatible inter-observer measurements reproduce the observed Universe.**
Fields, particles and spacetime are effective descriptions of those records.
OPH does not claim to prove its own uniqueness or force its architectural
choices. Verify self-reference and the observable consequences on the same
declared construction.

Describe physics as **simulation-generated observables**. The executable
produces finite events, apparatus records and observer comparisons of their
interpretations. Agreed records support effective field and spacetime
descriptions; the program does not output an additional objective continuum
universe. Smooth continuum limits are mathematical descriptions with stated
domains and error bounds. Agreement alone neither fixes Born probabilities
nor proves a relativistic field limit. A record-copying operation copies an
orthogonal classical record, not an arbitrary quantum state.

Use the [construction contract](model/CONSTRUCTION.md#what-is-being-built)
and [executable-model hypothesis](docs/OPH_EXECUTABLE_HYPOTHESIS.md).
The carrier, code, repair law, scheduler, preparation and observer apparatus
may be declared design choices. Do not make proof that these choices are
unique or forced by OPH a prerequisite for building the model or calculating
its consequences. A component may be a primitive of the declared architecture;
claiming that a narrower source implements it requires an implementation proof.
Conditional finite theorems retain their stated hypotheses and domains.

Regulators are construction choices too. Reproducing a particular earlier
regulator, such as the discontinuous Wilson sign, is required only when
claiming equivalence to that program. An explicitly declared smooth or
finite-range alternative may instead prove its own complete records and
physical refinement. Keep its action and parameter ancestry explicit;
changing a regulator does not transfer a chirality or continuum theorem.

Freeze the declared construction and observation domain, then derive its
complete native observer outputs and consensus statistics, including retained
failures, no-operation records and quantum references. Establish the mappings
between actual preparations, interventions and shared records. Deriving a
consequence of a specified construction does not require deriving that
construction's uniqueness.

Field claims use the [primary source/state/reference/clock contract](model/particle_common_contract/PROOF.md)
and its typed consumers. Regenerate the actual member and result; matching
labels or hashes alone do not admit a different preparation or clock.
Comparison programs retain their stated domains until an explicit transfer
proof connects them. Fixing this construction does not prove its physical
continuum or observed particle dictionary.

Use the [joint field and particle dependency map](model/particle_joint_reconstruction/DEPENDENCIES.md)
and its machine register as the common acceptance list. Component controls
retain their source domains; their separate checklists do not add duplicate
physical requirements. Pure Yang-Mills gap claims require a pure-gauge
instance or a proved transfer. A photon-containing theory has no positive
global vacuum gap under its physical massless-photon interpretation.

Keep architecture choices distinct from empirical calibration. Settings
chosen using target masses or other comparison data are calibration inputs;
agreement with those inputs is not an independent prediction. Changes to the
construction or interpretation must preserve that input ancestry.

## The `docs/` directory

The full rules are in the [docs index](docs/README.md#rules-for-this-directory).

- The `docs/` root holds only the selected high-level programs, references
  and ledgers, such as `COMPLETION_PROGRAM.md`, `OPH_FALSIFICATION_PROGRAM.md`
  and `POSTDICTION_LEDGER.md`. Do not add a root document unless it is a new
  program or ledger of that kind.
- Never add a daily note, progress report, audit result, integration record,
  PR handoff or status summary anywhere under `docs/`, and never put a date in
  a file name or heading there. Record integration and review history in the
  GitHub pull request or issue.
- A note in `docs/research/` is allowed only when it states a current result
  that its code package verifies, without history, dated findings or progress
  narrative. Prefer extending the owning paper or package README, because
  separate notes go stale.
- Release papers do not cite files in `docs/`. Cite papers, Lean modules, code
  packages, machine registers and public preprints. The one exception is
  `docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md`.
- Code, receipts, reviews and registers bind code, Lean and receipts as
  evidence, never a note in `docs/`.
