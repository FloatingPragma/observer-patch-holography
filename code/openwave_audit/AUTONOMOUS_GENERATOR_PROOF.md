# A joint generator and its complete compression instrument

The declared carrier has eight position cells, four signed spin coordinates,
190 coexact oscillators and three noncompact harmonic coordinates. Integer
Gauss constraints fix the longitudinal field. The positive quadratic field
Hamiltonian includes the exact position-dependent Coulomb energy. Each
positive route operator T_a combines the position permutation with its actual
field Weyl translation. The reverse current and reversed path define its
adjoint. The executable independently checks all48 signed paths and their
Gauss equations from the archived geometry.

In the shared physical time t=h u, the declared generator is

    H = H_field/h + m beta + sum_a alpha_a(T_a^dagger-T_a)/(2 i h).

Positive translation acts on exp(ipx) amplitudes by exp(-iph); this
orientation gives the kinetic symbol sin(ph)/h. Each alpha commutes with
its route. The matter term is bounded self-adjoint with norm at most
M=|m|+3/h. Bounded perturbation defines the joint self-adjoint operator on
the field domain. Its conserved expectation gives

    0 <= <H_field>-E0 <= h(<H>-E0/h+M),

where E0 is the coexact vacuum energy. The coexact occupation and harmonic
momentum have corresponding bounds. The harmonic coordinate is unbounded;
this theorem supplies no stationary charged state or filled-particle vacuum.
The ordered routed macro and this continuous generator are distinct sources.
An exact spectral frame belongs to the theorem; floating coordinates require
separate certified comparison bounds.

For a unitary U=alpha T, define the block involution
S=[[0,U],[U^dagger,0]] and the isometry J psi=(psi,i psi)/sqrt(2). Its
compression is A=J^dagger S J=(U^dagger-U)/(2i), a Hermitian contraction.
A finite pulse followed by the complete ancilla measurement has successful
Kraus operator

    K_epsilon = cos(epsilon) I - i sin(epsilon) A,
    I-K_epsilon^dagger K_epsilon = sin(epsilon)^2(I-A^2).

For an eigenvalue lambda of A, take X in{+1,-1} with mean lambda.
Taylor's integral remainder and cancellation of the linear term give

    |E exp(-i epsilon X)-exp(-i epsilon lambda)|
       <= epsilon^2 Var(X)/2 <= epsilon^2/2.

Functional calculus and contraction/unitary telescoping therefore bound
the N-pulse amplitude error by N epsilon^2/2. Every first-failure outcome
is retained, with combined probability at most N epsilon^2. For an input
purified with any reference, the successful block contributes at most the
amplitude error to trace distance; the failure blocks contribute half their
probability. The complete flagged instrument is within N epsilon^2 of the
ideal unitary instrument. Three axes with epsilon=t/(Nh) contribute at
most3t^2/(Nh^2), before native implementation and field approximation errors.

The code checks exact block algebra with nontrivial holonomy and112 scalar
repeated-instrument cases. Its zero-eigenvalue control detects omission of
the failure records. These are route and instrument-algebra checks, with
no claim of a full native field execution.

The unbounded field/matter Lie product has strong convergence. A finite
native error certificate requires an additional field-domain estimate.
One admissible comparison method expands bounded hopping in cell time:
the norm remainder after K Dyson terms is at most(3u)^(K+1)/(K+1)!.
For two unitary sources, twice that tail plus a comparison of each retained
history suffices. The history comparison must be uniform over all ordered
continuous free-time partitions and at most K actual Weyl kicks. A fixed
per-step native error floor also requires scaling of precision and reset
caps with the number of pulses.

The audited source and exact-frame scope are bound by the receipt under
`openwave-tuning/common_dynamics/autonomous_dirac/`, with manifests
44c008892484f8d346fe5f8e29cf29a5a185159e7f748e295cce86b3457756d1 and
71efd702dd1254e33439dd2b048d8de0853b0101b2717a259339ccb9dd499a2f.
The executable uses that explicitly supplied sibling workspace. Run from
the repository root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  .venv/bin/python code/openwave_audit/autonomous_generator_review.py
```

The output is `autonomous_generator_results.json`. Its complete-row flag
is false: source-domain and compression results do not establish a charged
particle spectrum, the quantitative full native transfer or an interacting
continuum limit.
