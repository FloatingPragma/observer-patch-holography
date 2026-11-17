# Frozen execution task

Read TASK.json for all numeric inputs. Run the 36 preparations generated in
experiment.py (B=2,3; six families; rotations0,3,7) under each of native and
swap/hold. Each trace has8192 settling attempts,256 calibration attempts,
2048 held-out attempts, then a native cap64 spanning-tree observation. All
source records, initial audit states, reader outputs and failures are retained.
The first fixed native stack trace alone supplies the calibration magnitude.
The control uses that same magnitude. No seed or window selection after results.

The exact census has325 states/46800 branches at N=24 and703 states/160284
branches at N=36. Its class closure, reachability and uniform stationary column
sums must all hold; rate equals k(N−k)/(N(N−1)). The local endpoint audit has
148 actual native calls. Main body calls are377856 native and377856 control;
probe counts depend on actual observed coins and are reported literally.

Every gate is specified in DESIGN.md and implemented in experiment.py.
Missing terminal readout or calibration is inconclusive and makes its positive
gate false. Finite-window misses remain in the artifact. A source or decoder
defect abort is retained as a failed run; amendments preserve original bytes.

Prior knowledge: the earlier adjacent-level theorem, census and decoder were
already executed and known. The reservoir application and transport predictions
were proved before this run; no numerical pilot of this new protocol has run.
The source/task/code are hashed in PREREG.json and the local LEDGER before the
first execution. This is prospective execution of predicted behavior, not a
claim that known mathematics was discovered blind. No values of nature enter.
