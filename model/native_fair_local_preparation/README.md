# Preparing a local fair source

A recorded native repair word prepares each twelve-port cell close to a
uniform state within its actual occupation sectors. The exact error is at
most `2^-32` after 1931 attempts, independent of the federation size.
All sectors and raw records are retained. The error comparison allows
arbitrary initial correlations with outside systems and traces only the
mixing-record values. Local occupation weights
and the empty/full coherent corner are preserved. The
[proof](PROOF.md) specifies the admitted program, reference boundary and
source coupling required for coherent field preparation. The explicitly
supplied [occupation-wire coupling](COUPLING.md) also implements a complete
305-pulse conditional phase experiment on every nonconstant tree branch,
retaining the constant outcomes and all apparatus failures.

```sh
python -m model.native_fair_local_preparation.packet --check
python -m pytest -q -W error model/native_fair_local_preparation
```
