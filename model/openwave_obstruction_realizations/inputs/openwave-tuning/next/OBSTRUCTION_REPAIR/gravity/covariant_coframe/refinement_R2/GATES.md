# Representation refinement diagnostic, before results

Recorded2026-11-14 10:15:13 UTC; R2 results absent at recording.

The run doubles R/h from1 to2, coordinate ball radius from3 to6, and torus
side from40 to80. Thus proper R, nominal proper ball radius, qR and Gaussian
preparation width in units of R match. Digital ball volume is reported
separately; the ball-to-regulator ratio remains3 and is not a mesoscopic
R/ell convergence test. The source excitation is not re-fitted.

The run is a component diagnostic, not acceptance of the Einstein tower.
The following gate is fixed before the R2 output has been observed:

* The complete finite operator/entropy control must pass independently.
* Uniform-sigma residual must decrease under refinement.
* Interpret the sign of the first nonzero geometric response only if the
  raw uniform-sigma residual is at most10 percent of the nonuniform part
  |K(q)−K(0) W_B(q)/V_B|. This subtraction is used ONLY to size the error;
  the preparation kernel remains the raw K, with no fitted subtraction.
* Report all derivative stencils. Do not choose a stencil based on agreement
  with Einstein or discard the zero mode.
* Each array remains below the5000-side dense ceiling and the run below
  five minutes CPU. Longer wall time due CPU sharing is not increased CPU
  allowance. If interrupted for resources, retain partial outputs.

Physical volume equivalence of the CAR factor, independent common-ball
stationarity, a temporal/lapse response and physical G are separate gates.
