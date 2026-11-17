# Exact support proof and bounded numerical roundoff

The initial frozen implementation and its diagnostic both fail the explicit
bitwise exterior-probability-zero assertion. They remain unchanged. The raw
outside probability is7.124627860079694e-32, norm residual6.661338147750939e-15,
state residual4.454769886308441e-15 and record residual1.2490009027033011e-15.
The latter three satisfy the original2e-11 numerical tolerance.

The amended run adds an exact integer check on every padded arc. With
U=2sqrt1260 times u, J=SU, it verifies

    1260 WJ=2508 J−1260 U.

The source equation already gives WU=J, and the cycle's integer zero sums
and native reversal sign give Wv=v. Thus the invariant span and its fixed
support are exact consequences of all native source columns. No threshold
enters that proof. The numerical vector is also compared with the analytic
span at the existing per-coordinate tolerance2e-11. On the full524288-domain,
this implies outside probability at most524288(2e-11)^2=2.097152e-16. The
amendment uses this derived bound for the numerical exterior check. It does
not call the old bitwise gate a pass, alter source operations, fit a frequency,
change norm or other control gates, or remove any outcome. Receipt labels
retain the failed original gate and the exact/numerical support statements
separately. This is a newly declared numerical gate, with a bound derived
from the existing state tolerance rather than chosen from the observed leak.

Use run_amended.py and headline_amended.py for this stage. The original
run.py, headline.py, preregistration and two failed logs remain historical
replays. All executions set OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1 and
VECLIB_MAXIMUM_THREADS=1. No evidence has yet been run under this amendment.
