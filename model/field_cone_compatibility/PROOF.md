# A common principal cone, clock and frame

## What the source supplies

The [native null decoder](../common_source/null_tomography.py) supplies twelve
spatial vectors `u_p` from the actual source Gram. Their common squared length
is `R²=phi+2`. For the **declared** lifts `k_p=(1,u_p)`, the map
`K -> (k_p^T K k_p)_p` on real symmetric four-matrices has rank nine and kernel
spanned by `G=diag(-R²,1,1,1)`. Its exact left inverse on the slice `K00=0`
proves this for arbitrary real matrix entries. The
[transition criterion](../native_lorentz_transition/README.md) reuses that
decoder for changes of chart. Neither construction declares source events
to lie on these null directions.

## Conditional field matching

At one point, let `a` be a positive definite symmetric spatial principal
tensor of the existing static wave equation. Its principal symbol on
frequency covectors is `-omega²+xi^T a xi`. The dual form on tangent
displacements, up to a positive conformal factor, is

```math
-dt^2+dx^T a^{-1}dx.
```

Supply an invertible spatial readback Jacobian `dx=J dz` and one positive
clock conversion `dt=alpha d tau`. Put `s=alpha²>0`. In observer coordinates
the field form is

```math
K=\operatorname{diag}(-s,J^T a^{-1}J).
```

The following conditions are equivalent:

1. All twelve tests `k_p^T K k_p=0` hold.
2. `K=(s/R²)G`.
3. `J^T a^-1 J=(s/R²)I`.
4. `a=(R²/s) J J^T`.

The decoder gives `K=lambda G` from the first condition. The time-time
entry fixes `lambda=s/R²>0`. The block identity proves the third condition;
inverting it proves the fourth. Every step reverses. Thus the finite tests
identify the entire principal cone, without checking infinitely many
directions. Positive definiteness and invertibility are necessary input
conditions, not consequences of nearly vanishing directional residuals.

This is a local statement about principal characteristics. Mass terms,
lower-order couplings, boundary conditions, state preparation and the full
field law are not determined by the cone. Even exact cone equality is
weaker than comparison of complete experiments. Conversely a conservative
upper bound on propagation speed does not supply equality of cones.

For example, with constant `a`, the free equations with any `m²>=0` have
the same principal cone and dispersion `omega²=xi^T a xi+m²`. Directional
matching alone cannot select even this free mass parameter, let alone the
particle content or interacting spectrum. This comparison does not change
the source action; it identifies information absent from the cone tests.

## Constructive compatibility and frame covariance

For every positive `a` and `s`, `J=(sqrt(s)/R) a^(1/2)` is an invertible
solution at a point. This proves nonvacuity. Such a chosen frame is not an
independently measured spatial readback, and a pointwise Jacobian field
need not integrate to one global coordinate chart.

The executable uses a nonsingular upper-triangular rational `J`, sets
`a=JJ^T` and `s=R²`, and obtains exactly `K=G`. It also applies the exact
algebraic boost control from `native_lorentz_transition`; the source does
not select that boost. For any invertible real four-matrix `L`,

```math
K=\lambda G\quad\Longrightarrow\quad
L^T K L=\lambda L^T G L.
```

Both field and observer forms must transform. If additionally
`L^T G L=mu G`, then the field factor becomes `lambda mu` in the fixed
native chart. The transition certificate proves `mu>0` for an invertible
null-preserving map. Future orientation and properness are separate tests;
null-cone equality alone does not choose them. A general frame change may
introduce mixed time-space entries, so it must not be forced back into a
static block-diagonal ansatz.

## The fixed polar dictionary for a radial product

The existing [radial product](../physical_space_time/OCCUPATION_FIELD.md#independent-native-records-give-a-radial-product-field)
uses `B=-d_r²+A` and has a valid conditional curved static `1+3` wave
description. Test a particular proposed identification with flat space:
`x=r*n`, where `n` is a unit angular direction and `r>0`, together with
one global positive clock scaling. At a fixed angle the angular tensor of
a product does not depend on `r`.

For the explicit round control `A=kappa*(-Delta_S2)`, `kappa>0`, its tangent
form is `-dt²+dr²+kappa^-1 dOmega²`. Since
`dr=n^T dx` and `|dn|²=r^-2 dx^T(I-nn^T)dx`, the spatial form in `x` is

```math
H_r=nn^T+\frac1{\kappa r^2}(I-nn^T).
```

Matching `diag(-s,H_r)` to `G` requires the radial eigenvalue `1=s/R²`,
and then the tangential eigenvalues require `kappa*r²=1`. Two distinct
positive radii cannot satisfy this with one `kappa` and one `s`. The exact
control uses `kappa=1`, `s=R²`, radii `1,2`, and `n=(0,0,1)`. The first
radius passes every native null test; the second has tangential entries
`1/4` and fails. A nonconstant positive angular tensor has the same
two-radius obstruction at a fixed angle: its tangential eigenvalues are
multiplied by `r^-2` while the radial eigenvalue stays one.

This conclusion is restricted to this fixed polar map and global clock
conversion. It does not obstruct local orthonormal frames, another
coordinate dictionary, or the product's curved Lorentzian interpretation.
A warped `r^-2 A` angular term would change the tested action. No such
replacement is made or inferred from the source here.

## The actual causal-field join

The [native controller dataflow](../native_causal_dataflow/README.md) has
authenticated endpoint, controller and storage dependencies. Its shared
FIFO makes its complete order a chain; it supplies no spatial placement or
proper clock. The golden-population cone limit has a different specified
population, communication rule and layer clock. The field construction's
active coherent clock and radial readback are further specified operations.

To use the criterion above, a common source construction must produce the
event placements, readback Jacobian and clock conversion and identify its
propagation with that same field principal form. Equality must hold across
the required charts or be replaced by a controlled approximation, with
compatible source histories and all relevant controller paths included.
The exact twelve-direction test can then discharge the algebraic matching
step. It cannot manufacture its input maps, prove order reflection, or
replace the state and complete-instrument comparison. No common
source-selected causal field, interacting relativistic continuum or
particle spectrum is established by this conditional certificate.
