# Ground constraints for the original bounded-contact continuum refinement

For the original zero-chemical-potential action, the ground occupation
approaches one fermion per lattice site when `U/c` tends to zero. Every
fixed microscopic local correlation approaches the free filled-sea
correlation, with an explicit error bound. Along `c=v/a` and bounded
`U(a)->U0`, this fixes the limiting first energy moment of a continuum
wave packet, including the scalar term `U0/2`. It does not prove that the
macroscopic continuum field is Gaussian or that its interaction vanishes.

The state is exactly the translation-averaged full-ground state from
[thermodynamic_ground](../thermodynamic_ground/README.md), with the
hidden-spin symmetry proved in
[retarded_spatial_locality](../retarded_spatial_locality/README.md).
No chemical-potential subtraction, new action or replacement reference is
made. The free sea below is a variational comparison, not the interacting
ground. All integrals use normalized Brillouin measure `d^3k/(2pi)^3`.

## Energy and occupation defect

Write `h(k)=-c sum_j sigma_j sin(k_j)`, with eigenvalues `+-c s(k)` and
`s(k)=sqrt(sum_j sin(k_j)^2)`. Its zero set has measure zero. Let `P_+-`
be the two rank-one spectral projectors away from that set. Translation
invariance makes the one-body covariance a measurable matrix multiplier
`rho(k)`, with `0<=rho(k)<=I`. The free filled sea has covariance `P_-`.
Oddness of the sine vector gives its onsite covariance `I/2`; the free
CAR Wick rule therefore gives double occupation `d_free=1/4`.

Restrict this free state to each finite open box. The variational
principle bounds the full ground energy above by its trial expectation;
missing bonds contribute only a vanishing boundary/volume error. The
existing energy-density limit then gives

```
0 <= e_kin-e_free+U*d = e_U-e_free <= U/4,
e_free=-c integral s(k),             d=omega(n_up n_down).
```

For completeness, the first inequality does not assume that the
interacting covariance commutes with `h`. Set

```
w(k)=tr(P_+ rho)+tr(P_- (I-rho)),       0<=w(k)<=2.
D=integral w(k),
e_kin-e_free=c integral s(k) w(k).
```

Off-diagonal covariance terms drop out of the kinetic trace. Both terms
in the energy excess are nonnegative, so `d<=1/4` and
`integral s w <= U/(4c)`.

For `0<delta<=1`, the set `s(k)<=delta` is contained in the Cartesian
product of the three sets `|sin k_j|<=delta`. In one dimension that set
has normalized measure `2 arcsin(delta)/pi<=delta`: convexity puts
`arcsin` below the chord from zero to one. Consequently its three-dimensional
volume is at most `delta^3`. Splitting the occupation defect gives

```
D <= min(2, 2 delta^3 + U/(4c delta)).
```

When `0<U<=24c`, choose `delta=(U/(24c))^(1/4)` to obtain
`D<=8(U/(24c))^(3/4)`. The calculator uses exact rational thresholds
instead of evaluating fractional powers. Its finite dyadic search is a
calculator choice; the theorem holds for every threshold in the stated
interval.

The particle density obeys
`|n-1|=|integral tr(rho-P_-)|<=D`. Combining this with the earlier
`1/4<=n<=1` and onsite symmetry `rho_x=rI` gives

```
max(1/4,1-D)<=n<=1,      |r-1/2|<=min(3/8,D/2).
```

Thus `n->1` and `r->1/2` whenever `U/c->0`. For example, at
`v=1/8,U=1,a=1/65536`, the rational choice `delta=1/16` proves
`1023/1024<=n<=1`. This is a bound, not a computed density.

## Every fixed microscopic local word converges

This statement does not require importing a uniqueness theorem for the
free ground. For a unit basis vector `f=delta_(x,sigma)`, the free-sea
annihilators are `c(P_+ f)` and `c^dagger(P_- f)`. Each has squared GNS
norm in the interacting state at most `D`, by the positive Fourier
integrands defining `D`. The projected vectors need not have finite
support: CAR norm continuity defines their operators on the same
quasi-local algebra.

Write any local unit basis field `A_j=c_xsigma` or its adjoint as
`A_j=b_j+d_j^dagger`, where both `b_j,d_j` annihilate the free sea.
Both have operator norm at most one, and their interacting GNS norms are
at most `sqrt(D)`. In the expectation of a word of even length `m`,
the leftmost `d_1^dagger` term has magnitude at most `sqrt(D)` by
Cauchy-Schwarz. Move `b_1` to the right using exact CAR. Its final term
has the same bound; each of the `m-1` intervening anticommutators is a
scalar of magnitude at most one. This is the free Wick recursion plus
two controlled errors. Induction yields

```
|omega(A_1...A_m)-omega_free(A_1...A_m)| <= C_m sqrt(D),
C_0=0,       C_m=2+(m-1)C_(m-2).
```

In particular `C_2=2` and `C_4=8`, so `|d-1/4|<=8 sqrt(D)`.
Odd words have zero expectation in both gauge-invariant states.
Every fixed finite local matrix is a finite linear combination of such
words; hence all its expectations converge. Onsite symmetry makes the
four probabilities `1-n+d,n/2-d,n/2-d,d`. Their trace-norm distance from
`I_4/4` is at most `min(2,2D+32 sqrt(D))`.

## What passes to continuum energy measures

For a physical smooth wave packet use cell averages
`(J_a f)_n=a^(-3/2) integral_(a(n+[-1/2,1/2)^3)) f(x) dx`. They are
contractions and satisfy `J_a^*J_a->I` strongly. Central differences
and translation continuity give
`h_a J_a f-J_a h_W f->0` for `f` in `H^1`, where
`h_W=i v sum_j sigma_j partial_j`. Multiplication by a lattice valley
phase replaces each derivative coefficient by `(-1)^b_j`.

The exact addition/conjugate-removal spectral measure has first moment

```
integral E dmu_(a,f)(E)
 = <J_a f,(h_a+U(a) r_a)J_a f>.
```

If `U(a)->U0<infinity`, the result above gives the limit
`<f,h_W f>+(U0/2)||f||^2`. The existing identity
`L V=V h+U J` and `||J||<=1` bounds its second moment by
`(||h_a J_a f||+U(a)||J_a f||)^2`, uniformly in `a`.
Therefore the first moment also passes to any weakly convergent spectral
subsequence: its omitted first-moment tail is at most `C/R`. This is a
full spectral sum rule, not a mass-shell equation or the low-energy
principal symbol. In particular the scalar term is retained, not
subtracted to create a preferred dispersion relation.

Two negative controls prevent invalid conclusions. Probability measures
`(1-n^-2) delta_0+n^-2 delta_n` converge weakly to `delta_0`, although
each has second moment one and its variance tends to one. A memory Gram
can therefore disappear from a weak spectral limit despite uniform second
moments. Separately, an occupation defect supported in a momentum cube
of volume `n^-3` has vanishing integrated defect, but a normalized packet
supported in that cube still has defect expectation one. Fixed microscopic
local convergence is not uniform convergence on the scaling wave packets.

These results fix ground inputs and energy sum rules for the continuum
construction. They neither identify physical particles nor settle
non-Gaussian continuum interactions, Lorentz covariance, real-frequency
thresholds or poles, or the joint Einstein premises.

## Replay

```
.venv/bin/python -W error -m pytest -q model/continuum_ground_constraints
.venv/bin/python -W error -m model.continuum_ground_constraints --write
.venv/bin/python -W error -m model.continuum_ground_constraints --check
```

The receipt binds the existing ground/response proofs and the exact
calculators. The infinite-state argument is the proof above; finite
examples are not asserted to compute the limiting state.
