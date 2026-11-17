# Execution history

Attempt01 is retained under `attempts/01_invalid/`, including its original
frozen design/task/code, full compressed records and terminal traceback. It is
invalid evidence. The new Observer decoder, not the sealed provider, assigned
its mutable reconstruction state to the saved initial reference. Later writes
therefore changed the retained reference. This was identified during execution
before any success claim. All native stages finished and the raw record file
was saved; receipt serialization then failed on a NumPy int64 defect counter.

The correction gives the decoder independent working copies of the initial
load/slot lists and explicitly converts the matrix defect counter to a Python
int. Neither source law, input schedule, checkpoint denominator, coin tape,
gate nor threshold changes. The scope text also clarifies that the existing
outgoing-shell experiment is a local proposal, not an extra universal platform
requirement. The new freeze must be registered before the corrected rerun.
The sealed R-work and memory provider files remain unchanged.

The corrected preregistration hash is
`e62d9daf6e80535e3130d0fa73d3b563232350800f498feb61181bd19a95ee89`.
The parent registered that hash before the corrected run. The corrected
execution exits zero and passes all nine gates. The two forced timeout
controls are correctly reported failures, not execution defects. Counts and
complete records are in RECEIPT.json and RECORDS.json.gz; RUN.log retains the
terminal report. Headline replay is recorded separately in REPLAY.log.
