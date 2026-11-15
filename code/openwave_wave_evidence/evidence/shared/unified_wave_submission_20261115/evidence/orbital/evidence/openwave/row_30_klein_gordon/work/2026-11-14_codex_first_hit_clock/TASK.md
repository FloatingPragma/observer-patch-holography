# First-hit source clock task

Run `python3 headline.py`. Candidate MFR-first-hit-clock is the registered
auxiliary N+ source defined in `DESIGN.md`; it is unpinned. Exact pre-run
terms, machine time and hashes are in `PREREG.txt` and `PREREG.json`.

The experiment executes all48,600 endpoint/spectator branches and4096 actual
trials on each of two schedules. It keeps every attempted current and every
censor, computes the branch-derived capped endpoint law, and compares actual
Fourier readouts at three fixed momenta using the registered simultaneous99%
Hoeffding bound. The source probe and control use the same pinned endpoint
function. This task supplies a record-defined Euclidean transfer clock for
the free-field reconstruction; it is not a platform submission.

Values of nature: none. Fixed seed2026111430; no selected seed or exploratory
fit. The two schedules and all thresholds were registered before execution.
