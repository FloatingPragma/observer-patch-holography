# A positive-energy continuum field net with all eight Weyl valleys

The original Weyl/contact lattice action has a subsequential continuum
field construction under `c(a)=v/a`, for fixed `v>0` and uniformly bounded
`0<=U(a)<=Umax`. All eight microscopic Weyl valleys are retained. The
construction supplies bounded smeared CAR fields, their complete
multi-time correlation hierarchy, continuous space and time translations,
positive energy, and exact graded locality outside a conservative cone.
It does not prove Lorentz covariance, surviving continuum interactions,
physical particle identification, or any completed Einstein premise.

This is an existence theorem on a common subsequence of the **spacing**
parameter, in addition to the previously constructed infinite-volume
ground at each spacing. It is not a claimed convergence rate for states,
a unique limit, or a numerical construction of its correlations. The
norm estimates below are uniform and explicit. No finite-volume gap or
ground-state clustering assumption is used.

## Inputs and the cell-average map

At every spacing choose the translation-invariant, gauge-invariant full
ground from `thermodynamic_ground`, with original Hamiltonian

```
H_a=dGamma(h_a)+U(a) sum_x n_xup n_xdown,
h_a(k)=-(v/a) sum_j sigma_j sin(k_j).
```

Keep its full quasi-local CAR algebra, original automorphism group
`tau_a,t`, positive ground GNS generator and full-CAR locality estimate.
No particle-number sector or momentum projection replaces that dynamics.
We use the conventional antilinear annihilator `c(f)=sum_xsigma
conjugate(f_xsigma)c_xsigma`, so the free field evolves as `c(exp(it h)f)`.
Let `C_a,n=a(n+[-1/2,1/2)^3)` and define the contraction

```
(J_a f)(n)=a^(-3/2) integral_(C_a,n) f(x) dx,
J_a^* u = a^(-3/2) u(n) on C_a,n,
P_a=J_a^*J_a.
```

Cellwise Cauchy--Schwarz proves `||J_a f||<=||f||`; `P_a` is the
orthogonal projection onto cellwise constants. For `f in H^1(R^3;C^2)`,

```
||P_a f-f|| <= a sum_j ||partial_j f||.
```

Indeed one-dimensional averaging on an interval of length `a` has error
at most `a||f'||`, by the fundamental theorem of calculus and
Cauchy--Schwarz. The three coordinate projections commute and are
contractions, so telescoping proves the displayed bound. Density extends
strong convergence `P_a -> I` to all `L^2` functions.

The cell map intertwines lattice shifts and translations by `a e_j`.
Writing `D_a,j f=[f(.+ae_j)-f(.-ae_j)]/(2a)`, the integral formula
`D_a,j f=(1/2a) integral_(-a)^a partial_j f(.+s e_j) ds` gives
`||D_a,j f||<=||partial_j f||` and `D_a,j f->partial_j f` in `L^2`.
For `f in H^3`, Taylor's integral remainder also gives
`||D_a,j f-partial_j f||<=a^2||partial_j^3 f||/6`.

For each `b in {0,1}^3`, let `M_b u(n)=(-1)^(b dot n)u(n)` and define
the physical smeared lattice field `C_a,b(f)=c(M_b J_a f)`. The convention
for linearity of `c(f)` is inherited from the source CAR algebra. These
fields still evolve under the full original Hamiltonian. Their kinetic
terms obey

```
h_a M_b J_a f = M_b J_a [i v sum_j (-1)^b_j sigma_j D_a,j f],
||h_a M_b J_a f|| <= v sum_j ||partial_j f||.
```

Thus every valley has the appropriate signed Weyl differential expression
`h_b=i v sum_j (-1)^b_j sigma_j partial_j`. This identifies the kinetic
part of the source, not a full interacting dispersion or a species.

## Orthogonal valleys without discarding the other modes

For equal valleys the CAR pairing tends to `<f,g>` since `P_a->I`.
For `b!=d`, choose a coordinate `j` where they differ. A one-step shift
changes the sign of `M_(b+d)`, so

```
|<M_b J_a f,M_d J_a g>|
 <= (a/2) [||partial_j f|| ||g||+||f|| ||partial_j g||].
```

To see this, subtract the simultaneously shifted pairing from the
original one; it equals twice the original pairing. Bound the two
differences using `||J_a(f(.+ae_j)-f)||<=a||partial_j f||`.
Consequently the limiting eight valley fields have the exact CAR over
`L^2(R^3;C^2) tensor C^8`, including zero cross-valley anticommutators.
All eight are constructed jointly. None of the unused lattice momenta
is removed from the finite-spacing Hilbert space or dynamics.

## A uniform time bound in the full CAR operator norm

For a lattice spinor `u`, define
`B(u)=sum_xsigma conjugate(u_xsigma) n_x,opposite c_xsigma` with the same linearity
convention as `c(u)`. On one site, `B(u_x)` maps the doubly occupied
state into a one-particle vector of norm `||u_x||` and vanishes on the
other number sectors. Therefore

```
{B(u_x),B(u_x)^*} <= ||u_x||^2 I.
```

Different-site odd operators anticommute. Summing the local inequalities
proves `||B(u)||<=||u||` for finite support and hence for all `ell^2`
spinors by norm completion. This is an operator bound before taking any
state expectation. The original commutator is a kinetic annihilator plus
`U B(u)` and gives

```
||[H_a,C_a,b(f)]|| <= v sum_j ||partial_j f||+Umax||f|| = L(f),
||tau_a,t(C_a,b(f))-tau_a,s(C_a,b(f))|| <= |t-s| L(f).
```

Local finite-support vectors are in the derivation domain. For general
cell-averaged `H^1` vectors, truncate in lattice `ell^2` at fixed `a`:
the kinetic operator is bounded there, and both annihilator and `B`
maps are norm-continuous. Closedness of the generator extends the formula.
Adjoints satisfy the same bounds. Products of finitely many fields have
the usual telescoping bound, with the product of the other field norms.

## One spacing subsequence for the complete field hierarchy

Choose a countable rational linear set dense in `H^1` and in `L^2`, with
compactly supported smooth members, and rational times. Include every
finite word in these fields and their adjoints, in all eight valleys.
Every expectation in the chosen ground at spacing `a` is bounded by the
product of the field norms. From any sequence `a_n->0`, diagonal
subsequence extraction therefore makes all these word expectations
converge simultaneously. No word chooses its own subsequence.

The norm bound in `f` and uniform time Lipschitz bound extend those limits
uniquely to arbitrary `L^2` smearing functions and real times. For a
general `L^2` function, first approximate it by the dense `H^1` set;
the time-extension error is uniform because every field has norm at most
its smearing norm. Multilinearity, adjoints, equal-time CAR and positivity
of every polynomial `X^*X` pass to the limit. The field norm inequality
also passes in its form `omega(X^* C(f,t)^* C(f,t) X)<=||f||^2 omega(X^*X)`.

The resulting positive functional on the algebraic *-algebra of field
symbols has a GNS representation by bounded operators `C_b(f,t)`.
The polynomial norm inequalities justify boundedness without assigning
a pre-existing norm to the free algebra. Equal-time CAR identify the
ordinary continuum CAR algebra faithfully inside this representation
(a nonzero unital representation of the simple CAR algebra is faithful).
Equivalently, the nonzero CAR identity itself already prevents collapse
of a nonzero smearing function to the zero operator.

In particular, a norm estimate passes to an **operator** identity by
testing it inside every polynomial context. If `||Y_a||<=epsilon_a->0`,
then `omega_a(X_a^*Y_a^*Y_a X_a)<=epsilon_a^2 omega_a(X_a^*X_a)`.
The limit operator `Y` vanishes on every polynomial vector and hence
vanishes everywhere. This is stronger than a zero vacuum expectation.

Common time shifts preserve every finite-spacing expectation and every
limit. They define isometric *-automorphisms of the generated C*-algebra
and unitaries fixing the GNS cyclic vector. Strong continuity follows
first on smooth smeared fields from the Lipschitz bound, then on the
dense polynomial domain and the Hilbert completion. All finite-degree
correlations are continuous and retain their product norm bounds.

## Continuous spatial translations and positive energy

For a physical translation `y`, choose a nearest lattice translation
`a m_a`, so `|y_j-a m_a,j|<=a/2`. Cell averaging commutes with that lattice
translation, while

```
||f(.+y)-f(.+a m_a)|| <= (a/2) sum_j ||partial_j f||.
```

There is an extra valley phase `(-1)^(b dot m_a)`. Choose the nearest
**even** lattice vector instead, with `|y_j-a m_a,j|<=a`; then that phase
is one simultaneously for all eight valleys. The same estimate holds
with `a` in place of `a/2`. Translation invariance of each lattice ground,
the field norm bounds and product telescoping give invariance of the
limit under common continuous spatial translations. The corresponding
unitaries are strongly continuous by `L^2` translation continuity, commute
with time translations and fix the ground vector. Cubic symmetries may
also be retained, but continuous rotations and boosts are not inferred.

Positive energy holds for **every** polynomial vector, not just a linear
field. Let `X` be any finite multi-time field polynomial and `X_a` its
lattice realization. Since the finite-spacing ground generator is
nonnegative, the bounded function
`omega_a(X_a^* tau_a,t(X_a))` has a positive spectral measure supported on
`[0,infinity)`. Its limit is uniform on compact `t` intervals by the
word construction and the time estimates. Against the inverse Fourier
transform of any smooth test function supported in negative energy,
dominated convergence therefore gives zero. The limiting time generator
has no negative spectral projection on any polynomial vector. These
vectors are dense, so the entire ground GNS generator is nonnegative.
No subtraction of a fitted chemical potential or new energy origin is
introduced in this passage.

## Exact graded locality in a conservative physical cone

First let `f` be supported in `[-r,r]^3` and `g` outside `[-R,R]^3`,
with `R>r`. Modulation by a valley sign does not change these supports.
Take lattice radii `r_a=ceil(r/a)+1`, `R_a=floor(R/a)-1` and
`0<a<=min(1,(R-r)/8)`. The cell-averaged fields are respectively inside
the first cube and outside the second. A finite-volume approximation to
`tau_a,t(C_a,b(f))` inside the second cube remains odd, so its
anticommutator with either `C_a,d(g)` or its adjoint is zero. The existing
full-CAR bound, with `s<=(2r_a+1)^3`, gives the following explicit bound
for `|t|<=T` and
`delta=(R-r)/2-48vT>0`:

```
||{tau_a,t(C_a,b(f)),C_a,d(g)^#}||
 <= min(2||f||||g||,
        192 v T (2r+5)^3 (2R+1)^2 ||f||||g|| a^-6 exp(-delta/a))
 <= min(2||f||||g||,
        192*7! v T (2r+5)^3 (2R+1)^2 ||f||||g|| a/delta^7).
```

Here `#` denotes either adjoint or no adjoint. The factor16 in the
first bound uses `R_a-r_a >= (R-r)/a-4`; `log 2>=1/2` supplies the
displayed exponent. The final inequality follows from
`exp(x)>=x^7/7!`. At `T=0` the anticommutator is exactly zero. Thus the
limit anticommutes whenever `R-r>96v|t|`.

For arbitrary separated compact supports, cover the first support by
finitely many small cubes and use a smooth partition of unity. If their
sup-norm distance is `d>96v|t|`, the cubes can be chosen small enough that
each inner/outer cube pair has a gap still larger than `96v|t|` and the
outer cube misses the other support. Apply the preceding bound to each
piece and sum. Time translation invariance gives the statement for any
two times. Algebraic graded commutation then extends from field generators
to the local polynomial algebras and their norm closures. This constructs
a translation-covariant field net with finite propagation speed `96v` in
the sup norm. This conservative bound is not a measured light speed,
an isotropic Lorentz cone, or an Einstein metric.

## Uniform low-momentum inverse response

There is also a useful bound before taking the spacing limit. At each
valley use the exact Fourier isometry from physical `|p|<=P` into lattice
momenta `k=pi b+a p`, with `aP<pi/2`. The images of all eight windows are
disjoint (under the torus convention). Translation invariance makes each
window reduce the full retarded Green operator and its Schur memory.
Its kinetic norm is at most `vP`. The previously proved exact formula
`G_a^-1=z-h_a-U r_a-Sigma_a`, with
`||Sigma_a||<=U^2/(4 Im z)` and `1/8<=r_a<=1/2`, gives

```
||G_a(z)^-1||_(window) <= |z|+vP+Umax/2+Umax^2/(4 Im z),
-Im G_a(z) >= Im(z)/( |z|+vP+Umax )^2 I_(window).
```

The second inequality is the same spectral Cauchy--Schwarz argument as
in `retarded_spatial_locality`, now applied only to this reducing window.
It proves invertibility with a bound uniform in spacing on compact upper
half-plane frequency domains. It neither takes the real-frequency limit
nor proves convergence of inverse operators from weak convergence alone.
For compactly Fourier-supported `f`, cell averaging differs from this
isometry only by the multiplier `product_j sinc(a p_j/2)`; its error is
at most `a^2 P^2||f||/24`. Thus both smearing conventions have the same
limiting field correlations on these packets.

The density bound used here is the previously proved one for `U>0`.
At `U=0` the filled negative-energy free ground has onsite covariance
`I/2`, hence `r=1/2`, giving the same bound directly.

## A controlled free-continuum corollary when the interaction vanishes

There is a full-CAR comparison independent of any ground-state limit.
Let `tau_kin,t` denote the free lattice dynamics at the same spacing.
Differentiating `tau_H,t-s(tau_kin,s(c(u)))` and integrating gives

```
||tau_H,t(c(u))-tau_kin,t(c(u))|| <= |U(a)t| ||u||.
```

Only the interaction derivation contributes. At each intermediate time
the free field is linear and its one-particle norm is still `||u||`, so
the bound on `B(u)` proves this Duhamel estimate. The domain is justified
by the same closed derivation argument used above; it is not an estimate
of the norm of an extensive interaction Hamiltonian.

For `f in H^3`, the central-difference estimate and a second Duhamel
comparison give the explicit all-valley bound

```
||tau_H,t(C_a,b(f))-C_a,b(exp(it h_b)f)||
 <= |t| [U(a)||f||+(v a^2/6) sum_j ||partial_j^3 f||].
```

The continuum free group preserves each of the displayed third-derivative
norms. Consequently, if `U(a)->0`, the continuum field dynamics is exactly
the free eight-valley Weyl dynamics, with norm convergence on these
smeared fields. No subsequence is needed for this dynamical comparison.

Its positive-energy ground is also fixed. A free annihilator supported
strictly in the positive one-particle spectral subspace has negative
GNS energy when applied to the ground and must annihilate it; a creator
supported strictly in the negative subspace must also annihilate it.
This follows by applying compactly supported spectral filters to the
time evolution, not by assuming individual momentum eigenvectors.
The zero-energy one-particle kernel is zero in `L^2`, since it is supported
only at `p=0`, a measure-zero set. Approximation by spectral bands therefore
extends the two annihilation rules to their complete subspaces. Moving
these annihilators through a finite CAR word determines every moment by
the CAR recurrence: odd moments vanish and even moments obey Wick's rule
with covariance `1_(h_b<0)` in each valley. This proves uniqueness of the
filled-sea ground state on the generated field algebra. All subsequential
hierarchies therefore coincide in this vanishing-interaction case.

This corollary is a controlled **free** continuum. For fixed `U0>0`, the
Duhamel bound retains `U0|t|||f||` and gives no vanishing interaction
error. It does not settle the required interacting continuum.

## What is still missing

The source interaction is retained throughout every finite-spacing
evolution, but a strict connected correlation has not been shown to
survive the spacing limit. PASS10's one-site non-Gaussian witness probes
shrinking physical regions and times; it cannot establish that claim.
The construction supplies neither a unique continuum vacuum nor numerical
limiting correlations. It proves no real-frequency particle classification,
scattering, Lorentz covariance, common metric or physical unit map.
The existing Einstein theorem still requires its same-family G1--G6
receipts: the finite propagation cone and positive energy here do not
provide a symmetric Lorentz stress, geometric modular flow, renormalized
area stationarity, curvature control or a Newton scale.

The finite calculations in this package audit cell averaging, the exact
onsite CAR inequality, valley signs and the quantitative estimates. They
are not simulations of the continuum state and do not replace the
compactness and GNS proof above.
