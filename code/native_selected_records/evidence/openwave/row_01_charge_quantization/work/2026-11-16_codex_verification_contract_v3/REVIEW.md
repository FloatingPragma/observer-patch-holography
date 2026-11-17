# Executable-verification audit and repairs

Reviewer/producer: Codex `/root/submission_review`, same group. Recorded
2026-11-16 Asia/Bangkok, with exact UTC registrations retained. This is known-
result reproduction and verification tooling, not prospective physical research.
The original source pin is `b6df2961f4ec8fc9a9f70a3c65e39da8dafa7e60`.

The current criteria were copied before process tests. Frozen CAUTIONS SHA-256
is `3f2ef5433513d8d6ed3796c64c2dfdb52cb2a6ca301c2f5bef1cf4a5df32c6f7`;
SUBMISSION_CHECKLIST is
`f8a08b120f8f298fccee025398341c210f67f58e1ca1c5198a3bf94262c8a501`.
They still matched the mutable workspace on completion. The binding audit
registration is [AUDIT_PREREG.json](AUDIT_PREREG.json). Static observations
preceded this registration; mutations and descendant runs followed their code
freezes. No source, criteria, outgoing V2, prior review or historical checklist
was edited.

## Findings on the original executables

1. **A real failed scientific gate could exit successfully.** In a disposable
   copy of the reservoir producer, the actual conditional transport rejection
   predicate was changed to reject the observed discrepancy. Its output had
   `conditional_transport_prediction=false` and `passed=false`, but process
   exit was 0. This is not a substituted report bit: the scientific predicate
   executed on the full calculated observations. The retained receipt is
   [reservoir_actual_failed_gate.json](original_processes/reservoir_actual_failed_gate.json).
   Its evidence digest refers to the unchanged full underlying scientific
   records; the mutation altered only the rejection condition.
2. **Optimization disabled scientific checks.** Each of the five outgoing
   payloads rejected a real scientific-predicate mutant under normal Python,
   yet accepted the same mutant with exit 0 under `-O`. The old source harness
   inherited `PYTHONOPTIMIZE` and therefore did not reliably preserve the
   original assertion checks in its subprocesses. This is an execution-contract
   defect, not evidence that any stored normal-mode scientific run failed.
3. **The old harness accepted insufficient verdicts.** Its parser looked only
   at the final top-level `passed` value. Separate harness probes accepted both
   failure-then-success JSON output and a false current gate map with no final
   `passed=false`. These probes expose the parser; they are not scientific
   competing laws. Original receipts are under `original_processes/harness_*`.
4. **Verification and generation needed distinct interfaces.** The five
   standalone scripts calculated stdout; the surrounding V2 integration did
   compare their output with saved summaries, but the scripts themselves did
   not offer a saved-result verification mode. Historical direct producers
   could overwrite records. In particular, the pair producer wrote its default
   named record files without a distinct first-record command. No observed
   historical replay silently rewrote the saved science: the primary wrappers
   generally already used temporary outputs and comparisons. A changed saved
   reservoir observation was rejected by its independent reader, with saved
   bytes and mtime unchanged. It would be incorrect to label the entire b6df
   replay as regeneration.

## Repair and validation

All five complete descendant scripts replace each assertion with an
unconditional `_require` call using the same predicate. Their output is held
until scientific checks pass, then compared with the complete saved typed JSON.
No numeric tolerance or scientific field omission is introduced. Default
verification performs no saved-data write; explicit `--record-to` creates a new
file exclusively. `adjudicate.py` independently normalizes the contract additions
and proves equality of each calculation AST with the original. The five files
have 197/338/200/438/114 lines for charge/energy/spin/common/pair, respectively;
the extra lines implement the inspectable contract, not an external archive
loader.

The reservoir direct-producer descendant checks the final verdict and exits
nonzero on a failed predicate, including after retaining an explicitly requested
failed record. Pair and reservoir dependencies also use unconditional checks.
Their default modes compare complete saved evidence. Static self-review found
that the first pair descendant did not additionally compare its saved receipt;
that narrower gap was repaired in a preserved, separately frozen amendment.
Its final mode validates receipt metadata, evidence digest and joint-column
count as well as the full source record. No scientific calculation changed.

The historical-source harness preserves the original packets and executes them
with `PYTHONOPTIMIZE=0`, inherited by child subprocesses. It validates every
returned JSON record, 34 explicit current-verdict paths, and full saved output
for all 28 jobs. Runtime exclusions are exact listed JSON paths. Legitimate
historical failed raw-positivity records and `external_G6=false` are retained;
they are not mistakenly treated as new scientific failures. Unexpected or
changed fields, missing records, current failed gates and nonzero process exits
reject. It checks custody plus before/after file bytes, sizes, mtimes and inodes.
It is an archive replay, not the small primary execution payload.

The retained results are:

| Check | Observed result |
| --- | --- |
| Five descendant scientific calculations | All five exactly match their saved V2 summaries; longest 24.25 s |
| Actual predicate mutants in those calculations | All ten reject: five normal and five optimized |
| Changed saved scientific values | All five reject; no bytes/mtime change |
| Explicit record / existing record / missing reference | New record succeeds; overwrite and missing reference reject |
| Direct reservoir producer | Full byte replay passes in 9.08 s; actual failed gate exits 1 under optimization and retains explicit failure receipt |
| Direct pair producer | Joint-composition mutant rejects under optimization; full source/receipt replay, changed receipt, first-record replay and overwrite tests pass |
| Historical source replay | All 28 jobs match saved scientific output exactly; 91.29 s internal / 92.56 s wall, maximum job 28.70 s |
| Historical source custody | 35 providers and 793 distinct manifest leaves pass; entire saved source tree unchanged |
| Distinct historical replay entry predicates | All 21 reject with actual `AssertionError`; repeated k/receipt arguments share the tested entry and all execute unmutated in the 28-job replay |
| Harness rejection behavior | Equal-length failure-then-success, false gate, changed observable and real job with changed saved baseline all reject; explicit historical failure is allowed |

`FINAL_VERDICT.json` checks 41 retained validation-process receipts and the 21
distinct historical entry mutations. It checks the rejection cause, not merely
whether any error occurred. Original defects remain in `ORIGINAL_AUDIT.json`.
The wrapper audit rebounded only the disposable reservoir headline's own
registration hash so the deliberately inverted real predicate could execute;
the original registration was untouched. All other mutated entry checks
executed against the unchanged original saved evidence. Mutation-generation
source and exact predicate fragments are retained; no corruption test receives
scientific emergence credit. No timeout was counted as rejection.

The environment was Python 3.13.7 with NumPy and one BLAS/OpenMP thread. Shell
startup emitted a pyenv rehash notice because its shim directory was read-only;
the actual verification subprocess receipts retain their own stderr and show
no such scientific-process warning. The suite can take over 60 seconds in
aggregate; each advertised scientific headline and each source job is bounded
and was measured below 60 seconds.

## Item and row dispositions

| Item | Rows | Historical V2 finding | Descendant disposition |
| --- | --- | --- | --- |
| E1 / failure exit and explicit verification | 1, 2, 7 | FAIL under newly explicit contract: demonstrated exit/optimization/parser gaps; historical positive results not erased | PASS for the exact repaired local interfaces and retained tests, subject to independent integration review and final outgoing binding |
| E1 / small payload and bounded run | 1, 2, 7 | Five original calculations readable; full archive never a substitute | PASS for the five complete descendants; source harness explicitly scoped as the longer archive suite |
| E2 / saved custody and replay | 1, 2, 7 | Existing b6df custody/replay remains valid at its scope | PASS for this local repaired replay of b6df; OPEN for a not-yet-created final V3 source export and new science jobs |
| D3 / registration | 1, 2, 7 | Historical registrations unchanged | This is disclosed reproduction/tool repair; no prospective discovery or upstream task-registration approval is inferred |
| E3/E4 / independent acceptance | 1, 2, 7 | Earlier scoped scientific reviews retained | Same-group execution audit only; external full G6 remains OPEN |
| F1/F2/F3/F4 / outgoing and release | 1, 2, 7 | Old sealed versions remain historical | New final integration, accessible source, completed checklist and owner release remain separate obligations |

The previous 22-PASS/6-OPEN historical checklist must not be reused unchanged
for the new current E1 or a new outgoing version. This packet repairs the
demonstrated local executable gaps; it does not by itself close other current
scientific findings, choose an upstream pin, satisfy external G6, publish a
source bundle or authorize release. No canonical RER or physical claim was
changed. New scientific descendants must be independently bound and checked;
their expected results must never replace the b6df contracts silently.
