# Operational readouts and finite sparse evolution

`ObserverProtocol` consumes an actual `JointRecords` state space and its
native target map. `readout` reports all position probabilities, their port
and transit sums, return to the prepared port, applied-current probabilities,
and native/boundary/transit outcome probabilities for a uniform monitored
attempt. These outcome probabilities sum to one. They are conditional on an
attempt and do not describe a coherent hard-wall loss measurement.

`sparse_hamiltonian` assembles the half-hop matrix as CSR. `evolve` uses
`expm_multiply` on a normalized pure state, avoiding a dense density matrix.
`trajectory_step` first evolves this state, samples the declared commit choice,
and, on commit, samples a source basis label from its Born probabilities,
then a uniform seam, fair coin and C/E law. It uses `native_target` and retains
no-commit, native, boundary-rejected and transit-noop outcomes. Ensemble averaging
samples the same finite density channel as `JointRecords.step`. Individual
trajectories are not the density matrix. The tests compare an ensemble with the
exact one-step channel with a stated numerical tolerance.

`RecordWindow` projects onto the radius-r ball, measured in half-hops, from the
prepared position and zero records, within the specified integer cutoff. This
is a supplied finite projection distinct from the complete cutoff sector.
Transitions across the window boundary are omitted symmetrically in H;
native attempts targeting outside it are retained as rejections. No equivalence
of the finite window and full cutoff sector is assumed.

For the prepared basis state, powers of the complete and projected H agree
through order r. Each lifted incidence column has two unit entries and each
row at most d_max entries. Therefore `||D|| <= sqrt(2 d_max)` and
`||H|| <= sqrt(m²+2 d_max)` in both spaces. With
`x=|t| sqrt(m²+2 d_max)`, two exponential Taylor tails give

```text
||exp(-it H_full)|initial> - embed(exp(-it H_window)|initial>)||
  <= min(2, 2 exp(x) x^(r+1)/(r+1)!).
```

`initial_state_error_bound` implements this bound in logarithmic arithmetic.
It applies to one coherent interval from the specified preparation. It supplies
no bound for repeated monitored steps or cutoff removal. Tests compare it with
the complete triangle sector for radii 0 through 3.

The radius-4 carrier window at cutoff 1 and m=0.23 has 51 states and 151 CSR
nonzeros. At t=0.37 its prepared-state norm error bound against the complete
cutoff sector is 0.119693660. Its source-port probability is 0.464401425 and
transit probability 0.517052970. These are window observables, with the stated
bound and preparation, without a physical-particle interpretation.

The exact density monitoring control uses the oriented triangle, cutoff 1,
m=0.23 and eight intervals t=0.37. Commit probabilities 0, 0.29 and 1 give
final purities 1, 0.102994192 and 0.103104175. Source-return peak-to-peak ranges
across the eight post-step samples are 0.719932216, 0.547201907 and 0.318383960.
These finite oscillations and support leakage do not establish a rest clock
or persistent localized particle. No frequency is inferred from this short
sample sequence. Monitoring choices are supplied independently of a physical
measurement-rate calibration.

Implementation: [observer_protocol.py](../oph_openwave/observer_protocol.py).
Checks and controls: [test_observer_protocol.py](../tests/test_observer_protocol.py).
Run: `python -m oph_openwave.runner observer_protocol` when registered by the
candidate runner. Graphs, preparation, integer cutoff, window radius, hop clock,
monitoring, and trajectory seed are declared inputs. No OpenWave criterion is
fully established by these calculations.
