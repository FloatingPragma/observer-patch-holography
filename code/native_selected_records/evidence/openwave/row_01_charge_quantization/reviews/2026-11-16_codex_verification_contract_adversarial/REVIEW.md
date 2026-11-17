# Independent same-group verification-contract review

The repaired local interfaces pass this review. No scientific calculation or
claimed physical result was added by the repair. This review binds producer
manifest `b25c9c07328cc62fe7852261af837ddf5d33a03d3031dd74b54d56c74ac1dd42`
and its 146 leaves. Reviewer: Codex `/root/charge_review`; the producer is
Codex `/root/submission_review`. This is same-group review, not external G6.

The reviewer read the builder, descendants, direct-producer interfaces,
historical harness, mutation generators, retained rejection traces and final
adjudicator. The [independent audit](audit.py) transforms the **original** ASTs
forward, rather than importing the producer's reverse-normalization routine.
All five complete primary calculations match exactly after converting each
original assertion to the same unconditional predicate and replacing its final
print with the verification payload. All six direct-producer/dependency
calculation bodies also match; their explicitly changed command-line main
blocks are excluded from that second AST comparison and were reviewed directly.
No assertion remains in any of the five repaired primary scripts.

All 28 expected historical outputs equal the original clean b6df replay,
before runtime normalization. The retained new full replay agrees on every
scientific JSON field. Its 22 exclusions are exact runtime paths; its 34
required current-verdict paths pass. The full replay reports 91.294 seconds
internally, maximum job 28.695 seconds, and unchanged source bytes, sizes,
mtimes and inodes. This reviewer checked that complete saved result and
independently ran one actual historical job under an optimized parent with
`PYTHONOPTIMIZE=2`; the harness resets child optimization to zero and the job
passes without changing the archive. This review does not claim a second
independent execution of all 28 jobs.

All ten retained compact scientific mutants fail at their intended changed
predicate, in normal and optimized Python. All 21 historical-entry mutation
traces identify the intended scientific comparison or its named assertion;
none is credited for a custody error, import error or timeout. The initially
unequal-length earlier-failure parser test is supplemented by a retained
equal-length test that reaches the actual failed-verdict check. The audit also
independently reruns all five compact predicate mutants under optimization.

Ten fresh subprocess checks pass: five optimized predicate rejections, a clean
optimized charge calculation, exclusive first recording, overwrite rejection,
the real historical job and read-only final adjudication. Three additional
actual-parser checks reject an earlier current failure and a changed nested
scientific check while accepting an explicitly historical negative. The
first-record file and saved reference directory retain their bytes and file
metadata after verification/refused overwrite. The entire sealed producer
packet is unchanged. The independent audit took 4.944 seconds on Python 3.13.7
with one BLAS/OpenMP thread; the [saved audit](AUDIT.json) contains the exact
bindings and dispositions.

Direct reservoir recording correctly preserves an explicitly requested failed
scientific run and exits nonzero. Its default mode compares complete evidence
and receipt. The pair producer's separately frozen receipt amendment checks
receipt metadata, evidence digest and column count in addition to the full
source result, excluding only its runtime. The producer's final adjudicator is
now read-only by default and compares its saved verdict; explicit recording is
exclusive. Earlier versions and demonstrated original failures remain filed.

Two scope notes must remain visible in the next integration:

- The historical b6df direct entry points were not edited. The new harness
  safely replays that archive; the repaired descendants are separate files.
  New scientific packets require their own saved expectations and a new clean
  export. This review cannot close that future integration obligation.
- The charge descendant retains an original docstring saying the program does
  not read files. Its unchanged calculation is file-independent, but its new
  default verification interface reads saved JSON. That inherited sentence
  should be clarified when preparing the next outgoing version; it is not an
  unreported scientific dependency or an execution blocker.

No additional scientific or executable blocker was found for the exact repair
scope. E1 is supported for these repaired local interfaces. A final source pin,
clean replay including new science, exact outgoing binding, external review
and release authorization remain separate gates. No prior checklist or
official OpenWave mark is promoted here.

## Replay

From this directory, `python3 -B audit.py` performs the same focused checks and
compares [AUDIT.json](AUDIT.json), without replacing it. `--record-to NEW_FILE`
is the explicit exclusive first-record mode. The producer packet and original
local integration paths must be present. NumPy is required by the inspected
energy/common calculations. There is no network access.
