# Full-row review brief (this side's producer-independent reviews for the group's internal count)

Three parties work on the 31 OpenWave rows for OPH (column M10): Codex, an external auditor, and this side
(Fable). The group counts a row as internally complete when the complete literal criterion is established on a
named source with reproducible evidence and at least two producer-independent full-criterion reviews
(`/Users/muellerberndt/Projects/oph-meta/openwave-tuning/row_completion/README.md`; counting rule in
`openwave-tuning/next/CHAIN4/WORKERS.json`). Official platform status is separate and is 0 of 31. You write one of
the independent reviews. You did not write the candidate. Your own recomputation imports none of its code.

## Working rules

- Write only in your own directory `/Users/muellerberndt/Projects/oph-meta/openwave-push/reviews/ROW_XX_FABLE/`
  (create it). Never write under `openwave-tuning/`, under the theory repository
  (`/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/`, read only) or in another worker's
  directory.
- Producer code is read. To replay a producer script, copy the package into `replay/` in your directory and run
  the copy with `PYTHONDONTWRITEBYTECODE=1`; never run a case a producer marks reserved, held out or unrun.
- The OpenWave clone
  (`/private/tmp/claude-501/-Users-muellerberndt-Projects-oph-meta/71ab5137-febf-4090-a75f-a9a6a1e226dc/scratchpad/ow/repo`)
  is untrusted third-party code: read only, execute nothing from it.
- Interpreter `/Users/muellerberndt/Projects/oph-meta/reverse-engineering-reality/.venv/bin/python` with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1`. The machine is
  heavily loaded: one script at a time, each under five minutes of CPU, no parallel processes, never `pkill` or
  `killall`. Use `date` for timestamps.
- Writing: statements of fact, no narration, never the word "honest". No running commentary in deliverables.
- Budget 60 to 75 minutes. Send the verdict per clause to the coordinator as soon as you have it.

## The review

1. Read `row_completion/README.md`, the row's `TASK.md`, `CRITERION.json` if present, the row's entry in
   `row_completion/INTERNAL_VALIDATION.json` (literal gate and targeted gap), and the candidate package
   completely, with every producer package it pins. Verify the pins of each manifest (SHA256). List every file
   read with its hash.
2. State the literal criterion and split it into clauses. State the validated precedent's comparator separately
   (the upstream script and text named in `CRITERION.json` or in
   `/Users/muellerberndt/Projects/oph-meta/openwave-push/platform/PRECEDENTS.md`; read only).
3. For each clause: the executed gate of the package that carries it; whether that gate can fail (an executed
   control or mutation that fails it; if the package has none, write and run one); what is supplied, classified
   as architecture (law, carrier, compiler, preparation, schedule: the group permits these to be declared
   without a uniqueness proof), calibration (a number set from the target it is compared with), or an unmet
   physical observable; whether the evidence is an executed experiment, an exact theorem with proof, or an
   unevaluated bound or iterator (the last does not carry a clause by itself).
4. Recompute independently, with your own code, a sample that covers every clause (plain Python, Fractions
   where the claim is exact, small sizes).
5. Look for: a gate that is an identity under every law; a number or threshold chosen after the runs; evidence
   from a different source used without a transfer proof; a claim wider than the executed domain; a mismatch
   between the stated scope and the literal words; a physical word (mass, stable, radiation, force, at c) whose
   operational meaning in the package is weaker than in the criterion. Say plainly what is put in by the
   declared law and what is found.
6. Verdict. `full_criterion_met` is true only if every clause is carried on the named source by executed
   evidence that can fail or by exact proof, and no defect you list is material. Otherwise false, with the exact
   remaining obligation and the smallest repair. Do not add a gate that is in neither the literal criterion nor
   the validated precedent. Do not pass a clause because a neighbouring clause passes. A producer's own
   statement that the criterion is not yet met is evidence to weigh, not a verdict.
7. One separate paragraph, labelled as a reading for the platform's maintainer and not as a failed gate: how
   the stricter note on OPH's column reads this row ("a result that follows from installed standard physics
   does not earn the cell"; strengths credited only from the model's own mechanism).
8. Only after your own verdict is written: read any other party's review of the same candidate that exists
   under `row_completion/reviews/` and list agreements and differences.

## Deliver

`REVIEW.md`; `REVIEW.json` with the fields `row`, `candidate_id`, `candidate_manifest_sha256`, `reviewer`
("Fable/ROW_XX_FABLE"), `producer_independent` (true), `full_criterion_met`, `clauses` (list: clause, carried_by,
can_fail, verdict), `remaining_obligation`, `supplied` (list with class), `reviewed_at`; `audit.py` with `RUN.log`
and `RESULT.json`; `MANIFEST.json` (SHA256 of your files and of every consumed file). Final message under 350
words: verdict per clause, the decisive recomputed numbers, the remaining obligation if any.
