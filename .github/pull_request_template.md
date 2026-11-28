## Row

Row ID from docs/registers/OBSERVABLES.json or model/evaluation/SCORECARD.md, or the issue this pull request answers (#1033 audit, #1078 obstruction search):

## Evidence

What the change shows and which files hold the evidence (code, Lean, paper text):

For a theorem-backed observable, identify the shared record, the exact premise it closes, the theorem and its remaining hypotheses. State the consequence established on this construction, with the measurement bounds. [Premise closure](../AGENTS.md#closing-a-physical-premise)

## Replay

Exact commands that reproduce the result from a clean clone:

```sh

```

## Falsifier

The run outcome, fixed before the run, that would have rejected the claim:

## What failed

Approaches that did not work:

## Checks

- [ ] `python -m model evaluate --check` passes
- [ ] `python tools/check_observables.py` passes
- [ ] No outcome is put in by hand, and every calibration has declared provenance
- [ ] Every status change follows from the evidence above
- [ ] Coding agent and model used, if any (accepted models: the newest Claude Fable, GPT-6.2 Astra, Gemini Carbon):
