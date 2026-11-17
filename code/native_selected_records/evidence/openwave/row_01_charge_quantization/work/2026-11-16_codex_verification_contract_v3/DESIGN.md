# Verification contract V3 audit and repair

This is an executable-verification audit and reproduction of already known scientific results, not prospective physical discovery. It addresses the frozen current E1 and CAUTIONS §11. The original b6df source, V2 outgoing package and all sealed records remain untouched. Only this work directory and disposable temporary copies may be changed.

## Fixed audit obligations

1. Inspect all five outgoing scripts and the 28 current source replay jobs, including primary wrappers and their execution dependencies. Distinguish actual scientific gate failure from altered metadata and integrity-only corruption.
2. Mutate an actual scientific condition or input in each outgoing calculation, and in each distinct primary wrapper/checker family, then run the real process and capture exit status. Use normal Python and explicitly test optimized Python where asserts protect gates. The expected result is nonzero. A timeout is a failure of the audit, not evidence that a gate rejected the mutant.
3. Mutate saved scientific values in disposable copies and execute verification. Require rejection by a scientific comparison (separate from custody-only rejection). Compare file bytes, sizes, mtimes and inodes before/after replay; temporary fresh output must not replace the saved expectation.
4. Inspect/test the orchestration parser against multiple JSON records and nested explicit failed gates, including an earlier failure followed by success. Historical negative records and external_G6=false are not scientific failures by themselves.
5. Repair descendants only. Actual calculation scripts must enforce their gates with nonzero process exit; a wrapper cannot mask broken direct entry points. Verification must compare to an existing saved baseline without writing it. First-record creation must be explicit and refuse overwrite. Retain failed original tests and the changes.
6. Reproduce unaffected outputs and apply the same gate/data mutants to corrected descendants. Limit a subprocess to 60 seconds. Report dependency/environment/runtime. Preserve every result and exact test code.

## Interpretation

These mutations test execution and evidence integrity, not source emergence or alternative physical laws. No new scientific award, source-family theorem, external G6 or release approval follows. Corrected code requires a new source pin and renewed affected checklist entries. Before-run freeze binds this design, current criteria and original entry points; subsequent implementation amendments retain their predecessors and are separately hashed before validation.
