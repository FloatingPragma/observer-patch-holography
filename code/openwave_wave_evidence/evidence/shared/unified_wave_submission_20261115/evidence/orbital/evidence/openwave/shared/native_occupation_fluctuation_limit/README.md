# Gaussian native-record limit

Native occupation blocks, completed by independent native zero-current coin
records, converge in every fixed moment to a Gaussian with covariance L^+.
Rooted differences give the Dirichlet Green covariance. The complete finite
source proof, exact covariance error and scalar reflection-positivity domain
are in PROOF.md. The theorem fixes the graph before the block limit; no
uniform joint refinement rate or physical clock is asserted.

The exact finite execution uses744 native calls, retains every branch, and
sums later IID words through an explicitly declared transfer recurrence.
For the port0-minus-port1 contrast, the limiting variance is11/30. At2048
events the exact variance is0.36456705729166666 and fourth moment is
0.41950248585806954, compared with Gaussian0.4033333333333333. The unchanged
12.5% fourth-moment tolerance passes. At512 it failed with14.76% error; that
entire attempt and its preregistration remain in attempts/01_window_512.
The longer resource choice was preregistered before evaluation.

Replay is byte-identical in4.78 seconds:

```
python3 openwave/shared/native_occupation_fluctuation_limit/headline.py
```

This is an N+ source component. No full OpenWave row is passed. The readout,
preparation, auxiliary coin address, scheduler, normalization and block limit
are declared. Four-dimensional routing, physical time and interventions are
separate statements. Same-group theorem review is bound at
../../row_29_em_waves/reviews/2026-11-14_codex_occupation_samegroup.
