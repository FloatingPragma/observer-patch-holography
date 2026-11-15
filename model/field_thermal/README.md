# Source preparation of an interacting field reference

This package executes the preparation of a faithful Gibbs marginal for the
same [complete-Fock Weyl/contact action](../full_fock_field/README.md), then
executes its Hamiltonian continuation and a complete occupation read. All
particle numbers, the purifier, the two read outputs and every first failed
proper-code export remain present. This supplies a bounded preparation
instance for CP-02/09. It does not select the action, temperature, physical
vacuum, clock calibration or geometric G1--G6 premises.

The executable domain is one or two two-component sites: at most four field
modes, four purifier bits and one read pointer. Its complete field dimension
is sixteen. This restriction comes from the existing complex-state compiler;
it is not a particle-number projection. It does not prepare the separate
8,128-dimensional fixed-number field reference or an infinite-volume state.

## A preparation with all outputs retained

For the declared Hamiltonian `H` and `beta>0`, set

```
rho = exp(-beta H)/Z,
R = exp(-beta H/2)/sqrt(Z),
|Omega> = sum_(i,j) R_ij |i>_field |j>_purifier.
```

The purification is normalized and its field marginal is `RR*=rho>0`.
It uses the positive matrix square root, so no arbitrary choice within
degenerate energy eigenspaces is part of the target. The existing source
compiler prepares the absolute amplitudes and then all their complex phases.
The word includes its scalar phase and zero-angle service requests.

Each actual six-dimensional primitive is followed by the existing
code-plus-failure export, which resets the sender. Every complement basis
value has a separate Kraus vector. They are combined incoherently within
the common first-failure event; the basis value is not an extra observed
record or a surviving coherent source register. A complement terminates
the program in that event's slot. Summing success and failure Gram matrices
preserves every untouched reference marginal in the exact unitary model;
this follows by telescoping the complete local export isometries. The
implementation retains all failure Kraus vectors, including numerical
near-zero vectors, and never renormalizes the success block. A deliberately
leaking-unitary test checks that failures and reference coherences really
survive. Valid-code failures are zero for the exact source generators.
Floating residual admission is not an outward or hardware certificate.
The preparation error includes its orthogonal abort blocks. The continuation
error separately compares its complete output, including its aborts, with
exact Hamiltonian evolution and readout of the same executed preparation
success block. For comparison with the ideal Gibbs experiment these errors
add; the earlier preparation aborts remain unchanged terminal outputs.

At `beta=1`, `c=1/8`, `U=1` on two sites, preparation requests **10,249
source pulses** and **20,498 load/export service events**, eight initialized
work bits and four retained purifier bits. The explicit model-time `1/16`,
four-step field continuation plus number pointer requests another **1,411
pulses**, **2,822 service events**, and one initialized pointer. Initial
blanks and physical preparation/storage costs are declared resources;
they are not assigned calibrated energy or duration here.

All three bond directions are tested against independent Jordan--Wigner
matrices and a matrix exponential. The contact interaction produces a
nonzero onsite Wick defect, so this is not a free Gaussian reference.

## The reference algebra matters

The field functional is faithful. The joint purification is pure and is
not faithful on the full field/purifier algebra. Reading the purifier and
conditioning on one result produces a rank-one field state, not the Gibbs
reference. Retaining all purifier records and taking the field restriction
recovers the original marginal. These are distinct contracts, and no
conditional outcome is used to promote faithfulness.

Likewise, `[H,rho]=0` makes the field marginal stationary under the declared
Hamiltonian model time. It does not make the joint purification static, or
make a number measurement preserve the reference. The read's backaction and
its complete outcomes remain part of the protocol.

For any subset of the field modes, `car_restriction` implements the full
signed factorization `F(A direct_sum B)=F(A) tensor F(B)`. Its matrix entries
contain the wedge permutation sign for ordering the kept modes first.
The Kraus family is complete, agrees with the actual embedded local CAR
algebra, and composes under nested restrictions. Dropping these signs fails
the odd-generator comparison for nonconsecutive kept modes.

Every restriction of a faithful full-Fock density to a nonzero local factor
is faithful: for nonzero `v`, sum the strictly positive quadratic forms
`<v,b|rho|v,b>` over an environment basis. Hence the existing
[reference-observation maps](../reference_observation/README.md) apply
without removing any local number sector. Their normalized square-root
factors cancel in a nested tower. The executable checks cover all matrix
units and complex untouched-reference columns. These observation/Petz maps
are mathematical maps here; compiling the reference preparation does not
automatically compile those maps or make them local retractions.

## Entropy and energy on that same field

The finite Gibbs identity of the
[same-action entropy analysis](../field_stress/README.md) has an
executed bounded full-Fock preparation here. For every normalized field density
`sigma` on this complete algebra,

```
S(sigma)-S(rho) = beta Tr[(sigma-rho)H] - D(sigma || rho).
```

This follows by substituting `-log rho=beta H+log Z`; it does not require
`sigma` and `rho` to commute. For
`sigma_epsilon=(1-epsilon)rho+epsilon tau`, `0<=epsilon<1`, the logarithm's
resolvent derivative gives

```
0 <= D(sigma_epsilon || rho)
   <= epsilon^2 ||tau-rho||_HS^2 / [2(1-epsilon) lambda_min(rho)].
```

Indeed the segment's minimum eigenvalue is at least
`(1-epsilon)lambda_min(rho)`, and its relative-entropy second derivative is
bounded by that denominator's inverse times `||tau-rho||_HS^2`. Integrating
twice gives the factor `1/2`. The noncommuting contact-pair control tests
this identity and its remainder on the actual interacting action. The
related Lean Gibbs ledger identities have their declared finite-energy
domain; this quantum matrix execution is not a kernel-checked theorem.

Finite Gibbs free-energy stationarity is not fixed-proper-volume
generalized-entropy stationarity, an area law, a Newton normalization, or
a physical rest mass. The current preparation closes none of those joins.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.field_thermal.thermal
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/field_thermal
```
