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

The executable model in [model/](model/README.md) has a fixed carrier,
[model/executable_model](model/executable_model/README.md), and candidate laws in
[model/candidate_laws](model/candidate_laws/CANDIDATE.md). The main candidate is
[A-contrast-memory-spherical-v1](model/candidate_laws/contrast_memory_spherical_v1/README.md),
specified in [SPEC.md](model/candidate_laws/contrast_memory_spherical_v1/SPEC.md); the scalar channel is a
second. [REQUIREMENTS.md](model/REQUIREMENTS.md) states the criterion for the correct law: it runs on the
fixed federation and fulfils the premises of the finite theorems together, with free calibration and no
hardcoded outcomes. A new candidate gets a folder with a `candidate.py`, and `python -m model evaluate --write`
updates its evaluation. Use a candidate's law, carrier, emitter and reader for model claims. Independent
component proofs retain their own hypotheses and require an explicit source map.

Use the [specification](model/candidate_laws/contrast_memory_spherical_v1/SPEC.md)
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

Field claims use the [specification](model/candidate_laws/contrast_memory_spherical_v1/SPEC.md)
and its typed consumers. Regenerate the actual member and result; matching
labels or hashes alone do not admit a different preparation or clock.
Comparison programs retain their stated domains until an explicit transfer
proof connects them. Fixing this construction does not prove its physical
continuum or observed particle dictionary.

Component controls retain their source domains; their separate checklists do not add duplicate
physical requirements. Pure Yang-Mills gap claims require a pure-gauge
instance or a proved transfer. A photon-containing theory has no positive
global vacuum gap under its physical massless-photon interpretation.

Keep architecture choices distinct from empirical calibration. Settings
chosen using comparison data are calibration inputs;
agreement with those inputs is not an independent prediction. Changes to the
construction or interpretation must preserve that input ancestry.

## The `docs/` directory

The full rules are in the [docs index](docs/README.md#rules-for-this-directory).

- The `docs/` root holds only the selected high-level programs, references
  and ledgers, such as `OPH_FALSIFICATION_PROGRAM.md` and
  `FROZEN_PREDICTION_LADDER.md`. Do not add a root document unless it is a new
  program or ledger of that kind.
- Never add a daily note, progress report, audit result, integration record,
  PR handoff or status summary anywhere under `docs/`, and never put a date in
  a file name or heading there. Record integration and review history in the
  GitHub pull request or issue.
- No research note repeats what the code, the Lean library or the papers
  contain. Extend the owning paper or package README instead.
- Release papers do not cite files in `docs/`. Cite papers, Lean modules, code
  packages, machine registers and public preprints. The one exception is
  `docs/research/PUBLIC_DATA_PREDICTION_AUDIT.md`.
- Code, receipts, reviews and registers bind code, Lean and receipts as
  evidence, never a note in `docs/`.
