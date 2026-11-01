# Compatible reference observations

One global faithful reference determines compatible local observation maps
without assuming that its regional densities factor. This uses the standard
finite Petz transpose construction. It provides a UCP observation family
with an exact nested composition law. It does not turn those observations
into conditional expectations or implement their reference square roots on
the native source.

## Finite construction and proof

Let `T:B(H_S)->B(H_B)` be a completely positive trace-preserving state
restriction, with Kraus matrices `K_a`. Fix faithful normalized `rho_S` and
faithful `rho_B=T(rho_S)` on their declared supports. Define

```
E_BS(X) = rho_B^(-1/2) T(rho_S^(1/2) X rho_S^(1/2)) rho_B^(-1/2).
```

Its Kraus matrices are `L_a=rho_B^(-1/2) K_a rho_S^(1/2)`. Thus it is
completely positive with arbitrary untouched references. Their `L_a L_a*`
sum is `I_B`, so it is unital. Cyclicity and trace preservation give
`Tr[rho_B E_BS(X)]=Tr[rho_S X]` for every matrix, not just commuting ones.
The adjoint `R_SB(Y)=sum_a L_a* Y L_a` is a trace-preserving recovery channel
and maps `rho_B` to `rho_S`. It need not recover any other state; preparing
the correlated global reference is part of this declared channel's cost.

For a second restriction `T_AB` with `rho_A=T_AB(rho_B)`, use the same
reference reductions and set `T_AS=T_AB T_BS`. The inner positive square
roots cancel exactly:

```
E_AB(E_BS(X)) = E_AS(X),     R_SB(R_BA(Y)) = R_SA(Y).
```

This establishes the entire nested tower from the state restrictions and
one reference. No product, Markov or modular-covariance hypothesis is needed
for this composition identity. Incompatible intermediate references break
it. Complete POVMs remain complete under `E`; all outcomes and arbitrary
reference systems survive without success normalization.

## What this changes in the regional interface

The [regional retraction criterion](../regional_modular/README.md) excludes
an exact Gibbs-preserving conditional expectation for the stated interacting
field region. The observation above has a different contract: in general
`E_BS(A tensor I)` is not `A`. Its effects can depend on the global reference
and need not be native local detectors. It is not the ambient tracial
projection of `CPRegionalNet`, and is not supplied as an inhabitant of that
unchanged Lean interface.

Likewise, reference preservation does not imply
`E_BS(rho_S^(is) X rho_S^(-is)) = rho_B^(is) E_BS(X) rho_B^(-is)`.
The tests retain failures of both retraction and modular covariance on a
correlated faithful state. For a product reference and ordinary tensor
partial trace, the construction reduces to `id tensor phi_rho_exterior`,
and these stronger properties hold. Thus the compatible observation tower
can be built before solving the separate common-clock and local-access joins.

The [same-source clock comparison](../modular_clock_comparison/README.md)
now implements relative modular corrections for the actual interacting
two-site Gibbs field. Exact spectral intervals also prove that any
uncorrected same-parameter modular intertwiner on this reference has
stationary regional range. The active corrections provide transitive clock
comparisons with explicit disturbance and complete source outputs; they
do not turn this observation map into a covariant retraction or calibrate
physical time.

## Actual CAR domain and controls

The implementation reuses the signed CAR restriction with every permitted
local number sector. It checks the construction on the same Weyl/contact
action restricted to a four-site three-axis star, on its complete 28-state
two-fermion sector. A nested four-mode to two-mode observation retains
interacting Gibbs correlations and agrees with direct restriction. The
subregion and Gibbs preparation are supplied; this is a bounded control of
the general finite proof, not a continuum limit or physical clock.

Controls also cover arbitrary complex input/reference operators, complete
POVMs, recovery of the declared reference, a product positive control,
inconsistent references and singular support. No small eigenvalues are
clipped and no missing support is silently discarded. Floating powers and
tolerances check the implementation; the algebraic argument proves the
identities. The numerical constructor rejects whitening when either unitality
or reference preservation has operator-norm residual above `1e-10`.
An ill-conditioned faithful reference can trigger this execution limit;
the exact finite construction remains valid. A small eigenvalue alone does
not trigger rejection when the computed maps meet the residual budget.
Host matrix execution is bounded to dimension 512.

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/reference_observation
```
