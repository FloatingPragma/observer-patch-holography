# Position measurements that retain the evolving source response

The same ordered-spin field admits a translation-covariant position
measurement on the closure of all its time-evolved field responses. The
detection sector and its dark complement are invariant under the actual
time generator. Momentum is the existing spatial translation generator,
and complete position outcomes obey the uncertainty inequality when their
second moments are finite. The [proof](PROOF.md) constructs the measurement
from the source history Gram. It supplies no sharp particle mass shell,
Lorentz covariance or native instrument realization.

The [finite source certificate](FINITE_PROOF.md) extends the source cube's
coordinate detector to its invariant128-dimensional odd-parity response
enclosure. Eight rank16 position effects retain every evolving response,
with an explicit even dark complement. Exact intervals give the complete
position probabilities under the same Hamiltonian. A nonzero matrix
entry proves that this finite detector requires all eight sites.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.spin_dynamical_localization --check model/spin_dynamical_localization/results.json
python -W error -m pytest -q model/spin_dynamical_localization
```
