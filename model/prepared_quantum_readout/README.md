# Preparing the native input for the complete phase/work readout

Nineteen fixed sweeps of the native thirty-seam repair prepare the input of
the [clocked phase/work experiment](../clocked_quantum_readout/README.md) to
an error below 3.72e-8. The 570-attempt prefix accepts any density on the
78-dimensional total-two carrier, including collision support and entanglement
with an untouched reference. It replaces exact diagonal I/66 preparation by
a finite source operation with an exact rational error bound.

The prefix ignores its records without selecting trials. The continuation
retains every failed herald, both herald signs, work record, phase bit and
unmeasured timeout. Its full cap-64 output contains 8,195 records. The
preparation error adds once to the complete output bound, independently of
the number of records or the useful herald probability. The separate
source-compiled mixer error is retained on its full 512-address carrier.

The [proof](PROOF.md) composes the existing exact deterministic-sweep theorem
with the complete finite Kraus instrument. The [receipt](results.json) binds
both sources, gives exact phase-event intervals at caps 4, 16 and 64 for C
and E, and adds preparation to the complete phase/work output certificate.
The executor additionally propagates arbitrary native/reference inputs
through the prefix and complete native suffix, without normalization or reset.

The total-two sector, native fair-coin quantum instrument semantics, coherent
service and work action/preparation/clock are premises. This preparation
theorem neither derives the Born rule from scalar consensus nor identifies
physical temperature or an equation of state. Continuations must ignore the
warm-up tape; selected herald laws require a separate conditional error bound.

```
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.prepared_quantum_readout.certificate --check
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/prepared_quantum_readout
```
