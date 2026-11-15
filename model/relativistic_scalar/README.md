# A source-executed canonical scalar field regulator

This construction executes an interacting scalar field on all eight vertices
of a three-dimensional cube. It allows particle creation, retains two local
measurement records and their feedback, and compares the **complete history**
with eight unrestricted oscillators using exact rational certificates.
The untruncated field has a finite-volume spectral gap bounded below by its
supplied mass. Its free spatial refinement approaches Klein--Gordon.
A controlled **interacting relativistic** observation domain is open.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.relativistic_scalar --check model/relativistic_scalar/results.json
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/relativistic_scalar/tests
```

The [derivation](DERIVATION.md) supplies the operator domains, observation
dictionary, error arguments and theorem reuse. The [receipt](results.json)
retains all 1,056 native/quantum histories at each cutoff, the independent
sparse-Hamiltonian comparison, source costs and exact probability intervals.
The [manifest](manifest.json) binds implementation, tests and results.

## The field and its observer protocol

The target action is the standard canonical scalar regulator

    H = (1/2) sum_x p_x^2 + (m^2/2) sum_x q_x^2
        + (c^2/(2a^2)) sum_<xy> (q_x-q_y)^2
        + (lambda/(24a^3)) sum_x q_x^4,
    q_x = a^(3/2) phi_x,  [q_x,p_y] = i delta_xy.

The executed cell has m=c=a=lambda=hbar=1 and twelve bonds. The oscillator
basis frequency is Omega=2. Its product blank is supplied as a preparation;
the protocol does not prepare or measure the interacting vacuum.

After evolution for Hamiltonian **model time** 1/128, an origin-site pointer
measures the projector onto (|0>+s i sqrt(24)|4>)/5 and its complement.
The first native repair current fixes s before the field experiment.
A committed one-bit result controls a phase on the origin's fourth
occupation. The same disturbed field continues for another 1/128 and is
read again. Its final bit controls the next native E/F repair. Both pointer
outcomes and both native repair coins remain in the record law.

All local evolution, pointer and feedback matrices are obtained from the
existing twelve-port source compiler. Their targets are used for independent
checks. The numerical global reference instead exponentiates the complete
sparse Hamiltonian. Cached local contractions save host computation while
the resource ledger counts every requested physical source word.

Unlike the [earlier Bose construction](../spatial_interaction/README.md),
this action has number-changing terms. Local occupation cutoffs therefore
need error bounds. The implementation compresses the **full** q^2, p^2 and
q^4 operators; taking powers of an already truncated coordinate would
define a different Hamiltonian.

## Retained comparisons

The origin retains eight occupation levels so that the same bounded local
instrument is used at both cutoffs. The other seven sites retain either
two or four levels. At 32 symmetric steps per evolution interval:

| Local occupation dimensions | State dimension | Occupation history TV bound | Ideal source-to-untruncated history TV bound |
| --- | ---: | ---: | ---: |
| 8 at origin, 2 elsewhere | 1,024 | <0.001131 | <0.001163 |
| 8 at origin, 4 elsewhere | 131,072 | <0.000618 | <0.000928 |

The last column adds the temporal product error to the exact occupation
comparison. The untruncated target is the **same eight-site spatial box**;
neither a spatial boundary limit nor an interacting continuum is included.
The bounds apply to complete joint probabilities, with disturbed states and
feedback, rather than a single fitted expectation. They cover the ideal
compiled operations; floating contraction errors, analogue precision and
physical service errors are separately excluded.

The matched interaction-deletion control changes the executed record law
by approximately 5.35e-5 in TV; deleting field feedback changes it by
approximately 0.00287. These numerical contrasts are diagnostics. A separate
bounded one-read witness at model time 1/4096 has an **analytic positive
interaction margin** for both the unrestricted action and the executed
cutoffs. Its observed contrast is approximately 1.22e-6. After charging the
source splitting error, its certified lower bound exceeds 6.1e-7 at the
larger cutoff. No conditional branch is normalized or postselected.

The unrestricted certificate uses exact finite-support oscillator
polynomials, with all Taylor remainder norms calculated by rational
arithmetic. Its degree-two approximation to the complete history has
retained-vector error about 2.07e-5. That approximation's raw probabilities
are deliberately left unnormalized; exact intervals contain the true law.
The occupation comparison adds the polynomial-state difference and both
remainder budgets, rather than treating agreement between finite runs as
a proof of convergence.

## What is established about relativistic physics

This is a spatial regulator for a real scalar field with a local quartic
interaction. The unrestricted equations have the discrete nonlinear
Klein--Gordon form. Positive mass and repulsive quartic coupling give
Hess V >= m^2 I. The existing harmonic-oscillator comparison theorem then
proves **gap >= m for every finite spatial box**, with unrestricted local
oscillators. The mass is supplied. This gap statement does not automatically
apply to the finite occupation matrices or an infinite-volume vacuum
representation.

The free spatial family has

    omega_a(k)^2 = m^2 + (4c^2/a^2) sum_j sin^2(a k_j/2).

For |k|<=K and aK<=1, its frequency differs from the relativistic dispersion
sqrt(m^2+c^2|k|^2) by at most c^2 a^2 K^4/(24m). The receipt includes free
evolution and covariance bounds for a=1,1/2,1/4,1/8 on the declared common
momentum band. These are analytic free-field comparisons, separate from
the finite-cell executions.

The next physical obligation is to compare **interacting** observation
statistics over a common spatial and energy domain, with matched mass and
coupling and quantified Lorentz-breaking corrections. Free dispersion and
finite-cell interaction measurements do not by themselves close that step.
In particular, holding the bare parameters fixed while changing a need not
hold the observed mass fixed. A mass counterterm can also invalidate the
positive-Hessian hypothesis used by this gap proof; it must be rechecked.

Positive scalar phi4 is an effective-theory candidate, not a demonstrated
nontrivial ultraviolet completion. The four-dimensional triviality result
and its precise relevance are discussed in the [derivation](DERIVATION.md).
The current action, spacing, preparation, pointer, native attachment and
model-time schedule are apparatus choices. Their success does not select
them from native repair or constitute an independent prediction of nature.

Blank initialization, coherent code services, storage, routing and pointer
readout are specified apparatus. Source resource counts do not include
those physical preparation/readout costs or establish efficient hardware.
