# The centered bounded-contact continuum is free

For the source-selected centered ground and `c(a)=v/a`, with fixed `v>0`
and `0<=U(a)<=Umax<infinity`, **every continuum field net constructed by
the retained all-word limit is the free eight-valley Weyl vacuum net**.
This includes its entire field hierarchy, not only its two-point function.
Its particle and hole excitation energies are `v|p|` and its relativistic
rest masses are zero. These modes are not identified with measured
particle species. This branch therefore does not supply the required
interacting theory or its massive sectors. Other coupling trajectories,
field scalings and microscopic modifications are outside this theorem.

The result uses the new centered state and clock `K=H-(U/2)N`; it does
not transfer that state to the original uncentered Hamiltonian. No
external renormalization theorem is imported. The proof combines the
existing source inequalities with a uniform microscopic-time estimate,
spectral positivity, and an explicit vacuum cyclicity argument.

## Precisely retained inputs

The frozen [centered ground construction](../centered_source_ground/README.md)
supplies translation-invariant, gauge-invariant, particle-hole-invariant
ground states at every spacing, positive centered GNS generators `L_a`,
and exact onsite one-particle spin covariance `I_2/2`. Its kinetic occupation defect satisfies

```
Delta_a <= min(2,11u_a^2/8),       u_a=U(a)/c(a).
```

The [full continuum construction](../continuum_field_net/README.md)
supplies bounded all-time fields `C_b(f,t)` for all eight valleys
`b in {0,1}^3`, positive energy on their full cyclic GNS space, the
equal-time CAR, and graded locality outside a finite propagation cone.
The limiting fields come from `c(M_b J_a f)` with cell-average contraction
`J_a` and `M_b(n)=(-1)^(b dot n)`. Its exact band isometry differs from
`J_a` by `O(a^2 P^2)` on a fixed physical momentum window `|p|<=P`.

The [composite algebra](../continuum_composite_channel/README.md) gives
the canonical onsite field `d=(2n_opposite-1)c`, with `||d(f)||=||f||`.
For the centered action,

```
[K,c(f)]=-c(h_a f)-(U/2)d(f),
h_a(b pi+a p)=-c(a) sum_j (-1)^b_j sigma_j sin(a p_j).
```

The full one-particle symbol tends on bounded physical windows to
`h_b(p)=-v sum_j (-1)^b_j sigma_j p_j`. In coordinate space this is
`h_b=i v sum_j (-1)^b_j sigma_j partial_j`. The source hopping is
`i c sigma_j/2` plus its adjoint; there is no missing factor two.

The existing CAR occupation argument also applies to the centered
ground: its expectation of any word of `m` unit site fields differs
from the free filled-sea expectation by at most `C_m sqrt(Delta_a)`,
where `C_0=0` and `C_m=2+(m-1)C_(m-2)`. In particular `C_2=2` and
`C_4=8`. This follows by moving the positive-band annihilator or
negative-band creator through the word, bounding its ground-vector
norm by `sqrt(Delta_a)`, and iterating the CAR. It is used only on
microscopic site fields in the uniform estimate below.

## Uniform microscopic-time comparison of the composite response

Measure microscopic time in units `1/c(a)`, so the normalized action is
`T+u_a(D-N/2)`, with normalized kinetic speed one. Write

```
F_a(s;f)=omega_a({tau_(s/c(a))(d(f)),d(f)^*}).
```

This is the characteristic function of the complete positive signed
response measure of `d(f)`, after frequency rescaling `E -> E/c(a)`.
Its total mass is `||f||^2`. The positive and negative frequencies
represent addition and removal; each physical GNS excitation energy
is nonnegative. Let `F_free(s;f)` use the same lattice field and free
filled-sea state, with speed one.

On a CAR polynomial of degree at most `m`, the centered onsite
interaction group has frequencies in `[-m/2,m/2]`. Each single field
has frequencies `+1/2,-1/2`, and multiplication adds frequencies.
The Bernstein inequality for this bounded-frequency trigonometric
polynomial therefore bounds its derivation by `(m/2)||A||`. Free
quadratic evolution preserves degree. Duhamel on the degree-three
field `d(f)` gives operator error at most `3u_a|s| ||f||/2`, hence
anticommutator error at most `3u_a|s| ||f||^2`.

For completeness, the state comparison is uniform in the lattice
smearing `f`, including slowly varying physical packets. Write
`d=2b-c`, with `b=c_opposite^* c_opposite c`. A graded anticommutator
of two cubic words expands into at most nine single CAR contractions
leaving words of degree four. The cubic/linear terms give at most
three contractions leaving degree two. The linear/linear term is a
scalar and has no state discrepancy. With `C_4=8,C_2=2`, the resulting
coefficients are `288,12,12`. For the normalized free hopping,
the row and column absolute sums of the one-particle propagator are
at most `B(s)=exp(3|s|)`. Sum the two remaining propagator factors
in absolute value and use the Schur bound on the contracted factor.
This yields `(300 B^3+12 B)sqrt(Delta_a)||f||^2`. Thus

```
|F_a(s;f)-F_free(s;f)|
 <= [3u_a|s|+312 sqrt(Delta_a) exp(9|s|)] ||f||^2.       (1)
```

The estimate keeps all lattice modes and both charge branches. Finite
support smearings prove it first; exact CAR norm continuity extends it
to every square-summable `f`. It uses no pointwise GNS fiber.

## Free cubic phase space and escape from fixed energy windows

In the free filled sea the onsite one-particle spin covariance is `I_2/2`. Normal ordering
therefore cancels the linear term exactly:
`d=2 :c_opposite^* c_opposite c:`. Applied to the vacuum or its adjoint,
it creates three positive-energy excitations. Their energy is

```
e(k_1)+e(k_2)+e(k_3),        e(k)=sqrt(sum_j sin^2(k_j)),
```

with the appropriate total-momentum constraint. Choose measurable
normalized spinors for the two rank-one band projections. For each
charge branch the unsymmetrized three-creation amplitude is twice a
sum over the two external spin components of the smearing, multiplied
by three band-spinor components of absolute value at most one.
The squared coefficient is four; antisymmetrizing the identical pair
costs at most a factor two in squared norm, and Cauchy--Schwarz for
the two external spin terms costs another factor two. Thus the bound
is sixteen for each charge branch. For
`0<epsilon<=1`, the normalized volume of `{e(k)<=epsilon}` is at
most `epsilon^3`: it lies inside the product of the three sine
sublevel sets, each of normalized length at most `epsilon`.
Two independent internal momenta must lie in these sublevel sets;
discarding the third restriction only enlarges the integral. After the
Haar-measure-preserving change from three momenta to external momentum q
and two internal momenta, integration of `||f_hat(q)||^2` gives `||f||^2`
by Plancherel. Hence,
uniformly in the external smearing and its momentum,

```
mu_free,f([-epsilon,epsilon]) <=32 epsilon^6 ||f||^2.    (2)
```

For larger epsilon the same inequality follows from total mass.
In particular there is no zero-energy atom. The exact full-Fock dimer
control in `algebra.py` illustrates the cancellation: the linear field
has excitation energy one, whereas the composite has energy three,
with weight one half in each charge branch. This control is not used
as a substitute for the phase-space integral.

Gaussian smoothing of the complete positive response measure turns
(1)--(2) into an explicit uniform bound. With Gaussian width epsilon,
the free integral is at most `1536 epsilon^6 ||f||^2`; integration of
(1) contributes at most
`[3u_a/epsilon+624 sqrt(Delta_a) exp(81/(2epsilon^2))]||f||^2`.
For a physical energy window `[-R,R]`, multiply their sum by
`exp(R^2/(2 c(a)^2 epsilon^2))`. Thus the squared norm of either
positive-GNS charge-branch projection of `d(f)Omega_a` onto energies
`[0,R]` is bounded by

```
eta_a(R)^2 ||f||^2,
eta_a(R)^2 <= min(1, exp(R^2/(2c(a)^2 epsilon^2))
  [1536 epsilon^6+3u_a/epsilon
   +624 sqrt(Delta_a) exp(81/(2epsilon^2))]).             (3)
```

First let `a->0` with epsilon fixed, then let epsilon decrease to zero.
Because `u_a<=Umax a/v`, (3) proves `eta_a(R)->0` for every fixed `R`,
uniformly over all unit lattice smearings. The calculator also supplies
explicit simultaneous choices. This is escape of the composite's
energy to infinity, not disappearance of its exact equal-time CAR norm.

## The limiting initial state is the free filled sea

For a physical momentum window where the positive kinetic energy is
at least delta and at most P, define the annihilation and composite
maps into the centered GNS space by
`A_a f=c(f)Omega_a`, `D_a f=d(f)Omega_a`, using the conjugate input
space for these antilinear maps. In the following map equations and
their exponentials, h denotes the induced conjugate one-particle
operator; the final physical smearing formulas use the original h.
The exact source equation is

```
L_a A_a+A_a h_a=-(U(a)/2)D_a.
```

At each fixed spacing these bounded maps lie in the generator domain,
since `h_a` is bounded and the right side is bounded. Integrate
`exp(-tL_a) A_a exp(-t h_a)` on this window. Its boundary at infinity
vanishes because `h_a>=delta` and `L_a>=0`. Splitting the output
energy at R gives

```
||1_[0,R](L_a) A_a|| <= Umax eta_a(R)/(2 delta),
||1_(R,infinity)(L_a) A_a|| <= Umax/(2(R+delta)).         (4)
```

Taking first `a->0`, then `R->infinity`, kills the positive-band
annihilator on the limiting vacuum. The negative-band creator obeys
the identical estimate for holes. On bounded physical momentum windows
away from zero, the exact lattice band projectors converge uniformly
to those of `h_b`; the cell-average/band-isometry comparison transfers
(4) to the retained fields. Such windows exhaust the one-particle
space up to measure-zero nodes. The field norm bound extends the
annihilation statements by density to every positive/negative band.

The CAR then uniquely determine the pure free filled-sea state on the
entire initial-time CAR algebra: move its annihilating parts to the
vacuum and iterate the CAR, giving Wick's rule at every order. This is
now a statement about physical packets, not only microscopic site words.
Let `H_0` denote this initial algebra's cyclic subspace in the full
limiting GNS representation.

## Both charge branches have free vacuum time evolution

On an arbitrary kinetic window `||h_a||<=P`, Duhamel in the preceding
map equation compares `exp(itL_a) A_a` to `A_a exp(-it h_a)`. The
low-output-energy integral is bounded by `|t| eta_a(R)`. For `L_a>R`,
integrate by parts using `L_a^(-1)1_(R,infinity)(L_a)`; its norm is at
most `1/R`. The two boundary terms and the differentiated input group
give `(2+|t|P)/R`. Consequently

```
||(C_a(f,t)-C_a(exp(it h_a)f,0))Omega_a||
 <= (Umax/2)[|t| eta_a(R)+(2+|t|P)/R] ||f||.             (5)
```

Here the map notation absorbs the annihilator's conjugate linearity.
The creator equation gives the same estimate for the adjoint difference;
it is essential to keep both physical charge branches. Band-isometry
comparison, `a->0` followed by `R->infinity`, and density yield for
every valley and every `f` in the limiting one-particle space

```
A_b(t,f)Omega=0,       A_b(t,f)^*Omega=0,
A_b(t,f)=C_b(f,t)-C_b(exp(it h_b)f,0).                   (6)
```

These are vector identities, not an inference from a scalar retarded
function alone. A finite countercontrol in `algebra.py` shows why (6)
by itself would not yet prove the full hierarchy: an empty two-mode
vacuum can have free one-field propagation while a quartic interaction
changes its two-particle energy. The spatial cyclicity step is needed.

## Outside cyclicity without assuming Lorentz symmetry

Let O be any bounded spatial region and let `P_b,+,P_b,-` be the free
Weyl spectral projections. The ranges
`P_b,+ C_c^infinity(O^c)` and `P_b,- C_c^infinity(O^c)` are dense in
their respective one-particle bands. Here is a direct proof.

If a band vector is orthogonal to the indicated range, it vanishes
outside the closure of O. Compact spatial support makes its Fourier
transform F entire in all three complex variables. The band equation
on real momentum reads `M_b(p)F(p)=+/-|p|F(p)`, with
`M_b(p)=-sum_j(-1)^b_j sigma_j p_j`. Fix real `(p_2,p_3)!=0` and put
`r^2=p_2^2+p_3^2`. The identity in `z=p_1` continues from the real line
as

```
M_b(z,p_2,p_3) F(z,p_2,p_3)=+/-sqrt(z^2+r^2) F(z,p_2,p_3).
```

Continue once around the simple branch point `z=i r`, but not `-i r`.
The square root changes sign; M and the entire F remain single-valued.
Subtracting the two identities forces F to vanish. Real transverse
momenta with r>0 exhaust an open dense set, so F is zero everywhere.
This proves the density claim separately for every one of the eight
valleys, without deleting any of them.

The pure filled-sea representation is fermionic Fock space on
`ran P_+ direct-sum conjugate(ran P_-)`, with the projections taken
over all valleys. Initial outside fields and their adjoints applied to
Omega create the dense hole and particle one-particle sets just proved.
Products create their exterior products plus lower-order contractions.
Induction on word length subtracts those lower-order terms. Thus
outside initial CAR polynomials applied to Omega are dense in H_0.
This is the required outside cyclicity. No spacetime CPT, Poincare
covariance or imported Reeh--Schlieder theorem is assumed.

## Locality promotes the vacuum identities to the entire net

Take compact smooth f. The free Weyl equation has finite propagation
speed v: its conserved current obeys `|j dot n|<=v|u|^2`, and integration
outside an expanding ball preserves zero exterior initial data. The
already constructed interacting candidate has its retained finite
propagation cone. Therefore both terms of A_b(t,f) in (6) graded-commute
with every initial field supported outside one bounded enlarged region.
For every homogeneous polynomial B of those outside fields,

```
A_b(t,f) B Omega = (+/-) B A_b(t,f) Omega =0.
```

Outside cyclicity makes these vectors dense in H_0, so the bounded
operator A vanishes on H_0. The adjoint argument proves that H_0
reduces every all-time field and that its action there is the free one.
The full continuum GNS space was defined to be cyclic for these
all-time fields. Hence it equals H_0, and A is the zero operator on
the entire space. Norm density extends this equality beyond compact f.

Thus every field and every ordered correlation function is exactly the
free eight-valley Weyl field and its filled-sea vacuum. The centered
bounded-U branch has no surviving connected scattering interaction.
Its positive-energy particles and holes have massless dispersion
`E(p)=v|p|`; no positive relativistic rest mass is generated by U/2.
The composite's equal-time CAR in a larger auxiliary limit is compatible
with this result because its energy escapes and its time group was
never asserted to extend continuously to that enlarged representation.

## Replay and the remaining physical question

```
.venv/bin/python -W error -m pytest -q model/continuum_free_limit
.venv/bin/python -W error -m model.continuum_free_limit --write
.venv/bin/python -W error -m model.continuum_free_limit --check
```

The receipt authenticates the frozen source inputs, exact finite negative
controls, uniform bound calculators and this proof. It does not pretend
that finitely many numerical checks prove the infinite hierarchy.
The full proof above supplies that conclusion. The interacting massive
physical endpoint now requires a different justified coupling/refinement
trajectory or a changed microscopic construction. This theorem selects
neither such a trajectory nor a physical species dictionary, gauge sector,
detector calibration or common Einstein geometry.
