# Maxwell evolution as a weak-field limit

The committed gradient `D` and curl `C` have sizes 30 by 12 and 20 by
30. Their exact identities are `CD=0`, `rank D=11`, `rank C=19` and
`1^T C=0`. Hence `ker C=range D` and `B0=range C={B:1^T B=0}`. The sixty proper carrier
rotations act by signed seam permutations `U` and ordinary face permutations
`V`, with `CU=VC`. These are the same tables and finite identities used in
[the local action classification](../maxwell_local_selection/PROOF.md).

Smooth actions and smooth conservative evolution need not be quadratic or
linear. Their leading response has the Maxwell form under the hypotheses
below, with a uniform error bound. The result concerns a declared finite
field model; it does not select that model from the native repair law.

## Smooth magnetic action

Let `S` be a real C3 function on a neighborhood of zero in seam-potential
space, invariant under local gauge translations `A -> A+D chi` and the
committed rotations. Suppose its Hessian **at zero** has strict face-sharing
support: two different seams can couple in that Hessian only when they
belong to one triangular face. This is a hypothesis about the vacuum
susceptibility; it places no locality restriction on higher derivatives.

Differentiating the invariances shows that `ell=dS(0)` is an invariant seam
vector and `H=d2S(0)` satisfies `HD=0`, `H=H^T` and `U^T H U=H`.
The exact finite classification gives

```
ell=0,                    H=beta C^T C.
```

No quadratic ansatz was used. A local minimum at zero gives `beta>=0`.
The additional requirement `H!=0` gives `beta>0`. A strict local minimum
on the gauge quotient alone does not imply this nonzero susceptibility:
`S(A)=sum_f (CA)_f^4/4` is a counterexample with zero Hessian.

Use the orthogonal representative `a` in `(ker C)^perp` of a gauge class.
On a quotient ball, if the operator norm of the third derivative is at most
`M3`, Taylor's integral remainder gives

```
|S(a)-S(0)-beta ||Ca||^2/2| <= M3 ||a||^3/6,
||grad S(a)-beta C^T C a|| <= M3 ||a||^2/2.
```

The gradient here is the quotient gradient, equivalently the transverse
seam gradient. If additionally `S(-a)=S(a)`, `S` is C4, and its fourth
derivative is bounded by `M4`, the two right sides improve to
`M4 ||a||^4/24` and `M4 ||a||^3/6`. The derivative bounds are explicit
assumptions about the neighborhood, rather than consequences of smoothness
with a selected numerical constant.

## Deriving the energy Hessian

Let `x=(E,B)` range over `X=R^30 x B0`. Let `W` be a real C3 energy on a
neighborhood of zero, invariant under `(E,B)->(UE,VB)` and the reversal
`J(E,B)=(E,-B)`. Assume that zero is a local minimum, and impose two
vacuum-susceptibility conditions:

1. The electric-electric Hessian is diagonal in the committed seam basis
   and is nonzero.
2. The pulled-back magnetic action `S(A)=W(0,CA)` has a nonzero Hessian
   with the strict face-sharing support stated above.

These conditions imply, on `X`,

```
dW(0)=0,       d2W(0)=K=diag(epsilon I_30,beta I_B0),
epsilon>0,     beta>0.
```

For the electric block, a signed permutation conjugates a diagonal matrix
by permuting its diagonal entries: the two signs square to one. The seam
action is transitive, so all diagonal entries equal `epsilon`. The minimum
and nonzero condition make it positive. This electric locality premise is
separate from the magnetic one; rotation covariance by itself permits
other electric quadratic forms.

Reversal removes the mixed electric-magnetic Hessian block. For the magnetic
block, the action classification gives `C^T H_B C=beta C^T C`. Since `C`
is onto `B0`, substituting preimages of any `B1,B2 in B0` gives
`B1^T H_B B2=beta B1^T B2`. This is equality of the restricted bilinear
form, so its representing self-adjoint operator on `B0` is `beta I_B0`.
It makes no claim about an ambient constant-flux extension. The executable
exhibits nineteen independent columns of `C` and their exact preimages;
no inverse of the rectangular curl is assumed.

The zero derivative also follows from the local minimum. Independently,
the averaged signed seam action is zero, while a rotation-invariant face
vector is constant and has zero projection to `B0`. Thus covariance
forbids either linear field term on `X`.

## Deriving linear evolution without a linear ansatz

Let `F:X->X` be a real C2 autonomous vector field on the energy neighborhood.
Require the following exact identities for every field there:

```
F_B(E,B)=-CE,                  dW(x)[F(x)]=0,
F(Jx)=-J F(x).
```

They state Faraday kinematics, conservation of the specified smooth energy,
and the specified electric-even/magnetic-odd time reversal. Differentiating
the conservation identity once at zero gives `K F(0)=0`. Since `K` is
positive definite on `X`, `F(0)=0`; equilibrium need not be assumed
separately. Differentiating again gives

```
K G+G^T K=0,                 G=dF(0).
```

Faraday fixes the lower row of `G`. Reversal gives `JGJ=-G`, eliminating
the electric-electric block. The cross term in the energy identity then
gives the uniquely determined derivative on closed flux:

```
G(E,B)=((beta/epsilon) C^T B,-CE).
```

Thus `F(x)=Gx+O(||x||^2)` without assuming a linear evolution law. The
two supplied positive susceptibilities and the evolution parameter remain
free. This result classifies the derivative of an autonomous field law;
it does not derive the quantum, Hamiltonian, or native implementation of
that law.

Faraday preserves closed magnetic flux exactly. The leading Maxwell law
also preserves `D^T E`. The hypotheses do **not** force exact nonlinear
Gauss conservation. If `D^T F_E=0` is imposed, it is preserved exactly by
the complete flow. A constructive nonlinear model below derives that
identity from its curl form, and a separate counterexample shows why it
cannot be inferred for every energy-conserving reversible flow.

## Uniform finite-time comparison

All norms in this section use `||x||_K=sqrt(x^T K x)` and the induced
multilinear operator norms. Suppose `W,F` are defined on an open neighborhood
of the closed radius-`R` ball, with uniform bounds

```
||d3W||_K <= M,               ||d2F||_K <= L.
```

Set `a=M R/3` and require `a<1`. Taylor's theorem on the ball gives

```
(1-a)||x||_K^2 <= 2(W(x)-W(0)) <= (1+a)||x||_K^2,
||F(x)-Gx||_K <= L ||x||_K^2/2.
```

If `||x0||_K<=r` and `J0 r^2<R^2`, where `J0=(1+a)/(1-a)`, conservation
and the first-exit argument keep the exact trajectory in the smaller ball
`||x(t)||_K^2<=J0 r^2`. Its invariant compact sublevel is inside the domain,
so the solution extends to every finite positive or negative time.
The Maxwell propagator `exp(tG)` is an isometry in the `K` norm. Duhamel's
formula therefore gives the uniform comparison

```
||x(t)-exp(tG)x0||_K <= (L/2) J0 r^2 |t|.
```

There is no exponential amplification factor. In normalized fields
`x0=eta z0`, the difference divided by `eta` tends to zero at least
linearly in `eta` on every fixed time interval. A caller must supply valid
derivative bounds for the entire ball; the rational bound API verifies
their arithmetic and confinement conditions, not those external estimates.
This is a finite-amplitude error estimate in the supplied clock, with no
uniform-in-infinite-time or continuum assertion.

## A nonlinear local model preserving all constraints

For supplied `epsilon,beta>0` and `lambda>=0`, define

```
W=epsilon ||E||^2/2+beta ||B||^2/2+lambda sum_f B_f^4/4,
dot E=epsilon^(-1) C^T(beta B+lambda B^3),
dot B=-CE.
```

The cube is componentwise. This smooth local model has exact energy
conservation, reversal, rotation covariance, closed-flux conservation and
`D^T dot E=0`. Its magnetic action is gauge invariant and its Hessian
`C^T diag(beta+3 lambda (CA)^2) C` is strictly face local at every field.
It has the selected Maxwell derivative and a nonzero cubic correction for
`lambda>0`. The finite assumptions permit that correction.

Write `r=||x0||_K`. Since `sum B_f^4<=||B||^4`, energy conservation gives
the global bound

```
||x(t)||_K^2 <= A2 = r^2+lambda r^4/(2 beta^2).
```

The carrier has `||C||_1=2` and `||C||_infinity=3`, hence `||C||_2<=sqrt(6)`.
The nonlinear remainder consequently satisfies

```
||F(x)-Gx||_K <= lambda sqrt(6/(epsilon beta^3)) ||x||_K^3.
```

The squared finite-time comparison bound is the exact rational expression

```
||x(t)-exp(tG)x0||_K^2 <=
    6 lambda^2 t^2 A2^3/(epsilon beta^3).
```

The API returns that exact bound and an outward rational square-root
enclosure. No numerical ODE trajectory is used as a proof. Uniform bounds
valid on a weighted ball are `M=6 lambda R/beta^2` for `d3W` and
`L=6 lambda sqrt(6/(epsilon beta^3)) R` for `d2F`; these also instantiate
the general theorem. At fixed time the quartic model has cubic absolute
error, or quadratic relative error after amplitude rescaling.

## A nonlinear Gauss counterexample

Let `L0=D D^T`, `z(B)=|C|^T B`, and
`A(B)=L0 diag(z(B))-diag(z(B))L0`. Then `A(B)` is skew, odd in `B`, and
covariant under the signed seam rotations: absolute incidence transports
`z` by the unsigned seam permutation, while diagonal conjugation removes
orientation signs. With the quadratic energy and unit susceptibilities,

```
F_E=C^T B+A(B)E,              F_B=-CE
```

preserves energy, reversal, proper rotations and closed flux. Its derivative
at zero is Maxwell. The executable supplies a closed `B` and a source-free
electric `E` with `D^T A(B)E!=0`. Thus the nonlinear flow can generate
Gauss charge even from initially source-free data. Exact nonlinear Gauss
conservation is an independent condition unless a compatible curl-form
construction is supplied. This counterexample is not the quartic model.

## Scope

Quadratic action and linear evolution are conclusions about the leading
response. Smoothness, a vacuum with nonzero susceptibilities, the two
specific locality conditions, carrier covariance, the field interpretation,
positive energy conservation, Faraday kinematics and time reversal are
explicit premises. They are not selected by deterministic native repair.
Physical energy coefficients, a physical clock, currents, continuum
refinement, photon poles and laboratory record maps require separate
attachments. The exact executable certificates and the analytic estimates
are not presented as a Lean-kernel proof.
