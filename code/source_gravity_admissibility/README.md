# Admissible source processes with identical clocks and different transport

The [proof](../../docs/research/SOURCE_GRAVITY_ADMISSIBILITY.md) strengthens the
two-field comparison: these are actual words in the
coherent twelve-port source, whose A1--A3 membership is constructed
separately. One stationary source and the entire resting-clock channel
are shared. Two fixed local coin policies give different finite
receiver probabilities. The missing spatial response is a function of
source strength; clock calibration and the typed A3 state objective do
not select it.

This is a non-entailment result, not a derivation of physical gravity,
Maxwell photons or selection by the integer-repair simulator. RG and CCG
are proposed services in the named mathematical source. The
[contract](CONTRACT.md) distinguishes the exact finite exit from the
long-wavelength interpretation and physical parent obligations.

From the repository root:

```sh
python -m source_gravity_admissibility.verify
python -m pytest -q code/source_gravity_admissibility
```

Set `PYTHONPATH=code` first. The checked-in `receipt.json` contains the
complete compact experiment census. To reproduce it:

```sh
python -m source_gravity_admissibility.build
```

The producer solves the grounded source equations exactly and executes
six-level native pulses reconstructed from the source response. The
checker imports neither that solver nor its circuit constructor: it
checks local source equations, independently reconstructs the response,
and uses Clifford projectors to build every reflected flight. It also
replays the retained native basis words on every input column. The
finite click gap and error intervals use exact rational arithmetic;
floating-point matrices are corroborating executions, not interval
proofs or experimental data. No external measurements are consumed.

The primary experiment uses five sites and all three Clifford flight
bases, four source strengths and both policies. Eight additional 27-site
three-dimensional executions retain massive and massless sectors.
Complete clock matrices and small-spacing Fourier matrices test the
channel and principal-symbol calculations. The exact one-axis massive
dispersion is also checked at 36 parameter settings with all four bands;
the velocity freedom is present before a continuum approximation.
An inhomogeneous native-word regression also checks the required
Hermitian gradient term and the quadratic local expansion error in all
three Clifford directions; removing that term fails the control.
Only the named receiver
read is a local physical report; full-matrix comparisons are software
checks, not an observer pooling spatial quantum states instantaneously.
