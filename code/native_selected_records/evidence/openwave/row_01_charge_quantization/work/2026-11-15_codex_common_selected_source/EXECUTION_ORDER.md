# Execution and verification order

1. `common.py`, DESIGN, TASK and the original local LEDGER were hashed in
   PREREG.json at2026-11-15T17:01:54.506993Z, before any common experiment.
   V1 first run passed all11 gates in3.619seconds and is retained in first_run/.
2. Independent verification was written after the result. Its first attempted
   preregistration helper failed because of an incorrect working directory;
   the first audit nevertheless ran and passed. It is retrospective verification,
   not prospective discovery. That exact audit/output is in verification_v1/.
3. AUDIT_PREREG.json at17:05:18.991723Z froze the strengthened work-word/source-
   graph/higher-mode verifier before its next run. It passed in1.855seconds.
   Integrity mutations, geometry and replay controls followed, openly as checks.
4. V2 native-selector design/code was hashed in V2_PREREG.json at
   17:10:49.327560Z, before its first execution. Prior conditional component
   numbers were known. The new selector/controller run passed in6.800seconds.
5. selector_audit.py independently checked its words, coherent labels, failures
   and reach factors. The first strict replay wrapper compared JSON-sorted
   integer keys to loaded string keys and raised an assertion despite matching
   scientific values. Its source/error remain in verification_v1/. Canonicalizing
   JSON keys repaired the wrapper; scientific source/evidence stayed unchanged.
   The successful strict replay compares every value except measured runtimes.

Times are UTC. The task folder retains its start-date naming although local
Asia/Bangkok calendar time crossed midnight. File timestamps alone are not used
as preregistration evidence; the ledger/hash records and ordered tool execution
supply the local provenance. No platform task was registered or pin advanced.
