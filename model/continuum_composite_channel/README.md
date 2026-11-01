# The resolved composite force channel under physical refinement

The original full Weyl/contact action has an exact second, nonlinear CAR
family. Under the bounded-contact physical refinement its mixed equal-time
anticommutators with the original continuum fields vanish as operators in
a joint limiting GNS representation. This resolves the algebra of the
contact-force channel; it does not prove that its energy remains finite,
that it has continuous limiting time evolution, or that interaction
survives. Algebraic independence of two CAR families is not statistical
independence of their state and is not a particle-species identification.

All results use the original full ground at each spacing, the same
`c(a)=v/a`, `0<=U(a)<=Umax`, and all eight valleys from
[`continuum_field_net`](../continuum_field_net/README.md). Nothing changes
the Hamiltonian, chemical potential, ground selector, or clock. The
occupation defect `D_a` from
[`continuum_ground_constraints`](../continuum_ground_constraints/README.md)
satisfies `D_a->0` and `|r_a-1/2|<=D_a/2`. For `U=0`, use the filled free
sea and `D_a=0`. A stronger rate for the actual ground is not assumed.

## An exact canonical composite and the contact equation

At a site, with the conventional antilinear annihilator, put

```
b_sigma = n_opposite c_sigma,
d_sigma = (2 n_opposite-1)c_sigma = 2 b_sigma-c_sigma.
```

The unitary `Z=exp(i*pi*n_up*n_down)=I-2*n_up*n_down` obeys
`Z c_sigma Z = -d_sigma`. Thus `d` is an exact canonical CAR family,
on every finite set and in the full quasi-local algebra. It is an odd
local polynomial of degree three, not an additional microscopic mode.
For a lattice spinor `u`, write `B(u)=sum conjugate(u_xsigma)b_xsigma`
and `d(u)=2B(u)-c(u)`. In particular `||d(u)||=||u||`.
The exact contact equation is

```
[H,c(u)] = -c((h+U/2)u) - (U/2)d(u).
```

The scalar `U/2` is recorded in the original energy convention. This
identity is not a Hartree truncation: the last operator is retained.

Define the local quadratic operators, for any `u,v`,

```
Q(u,v) = {B(u),c(v)^*} = {B(u),B(v)^*}
       = sum_x dGamma(q_x),
q_x = (u_x^* v_x) I_2 - v_x u_x^*,
P(u,v) = {B(u),c(v)}
       = sum_x [conjugate(u_xdown v_xup-u_xup v_xdown)] c_xup c_xdown.
```

With `epsilon=[[0,1],[-1,0]]`,
`q_x=(epsilon conjugate(u_x))(epsilon conjugate(v_x))^*`, so it has rank
at most one and trace `u_x^*v_x`. Consequently

```
||q_x|| = ||q_x||_HS = ||u_x||||v_x||,
||Q(u,v)|| <= ||u||||v||,
||P(u,v)|| <= ||u||||v||.
```

The mixed CAR defects are exactly

```
{d(u),c(v)^*}=2Q(u,v)-<u,v>I,
{d(u),c(v)}=2P(u,v).
```

Their operator norms need not vanish at finite spacing. It is their
controlled action in the same-ground limiting representation that does.

## Uniform local Wick comparison controls averaged fluctuations

The previously proved local Wick estimate is uniform in the positions of
its site fields, including sites whose separations grow with refinement.
It also holds for any unit spinor at a site: its relative-free annihilator
has squared ground norm at most `D_a`, by the same positive covariance
trace bound. No sum over a macroscopic smearing was used in that estimate.
For every product of four such unit fields the error from the free filled
sea is at most `8 sqrt(D_a)`.

Let `L=||u||||v||` and `I4=sum_x ||u_x||^2||v_x||^2`. The rank-one
factorization expresses each summand of `Q` as one product of two local
fields, with total absolute coefficient at most `L`. Applying the four-word
bound separately to the two site sums gives

```
|omega_a(Q^*Q)-omega_free(Q^*Q)| <= 8 sqrt(D_a) L^2.
```

In the free gauge-invariant quasifree sea with covariance `Pi`,

```
omega_free(Q^*Q)-|omega_free(Q)|^2
 = Tr[Pi q^*(I-Pi)q] <= Tr(q^*q) = I4.
```

Here `q` is the block-diagonal one-particle operator with blocks `q_x`.
The exact means are `omega_a(Q)=r_a<u,v>` and
`omega_free(Q)=<u,v>/2`. Centering at the latter and using
`0<=1/2-r_a<=D_a/2` therefore proves

```
omega_a[(2Q-<u,v>I)^*(2Q-<u,v>I)]
 <= 4 [I4+(8 sqrt(D_a)+D_a/2)L^2].                  (1)
```

For `P`, gauge invariance makes both means zero. The free pair covariance
is the positive contraction `wedge^2 Pi`, so
`omega_free(P^*P)<=sum_x |u_xdown v_xup-u_xup v_xdown|^2<=I4`.
The same four-word comparison gives the bound (1) for `2P` as well.
Finite support proves these identities first; the `ell^1` coefficient
bounds extend them by operator norm to the stated smearings.

For the centered-cell maps in every pair of valleys,
`u=M_b J_a f`, `v=M_e J_a g`, with bounded smooth compactly supported
`f,g`,

```
L <= ||f||||g||,
I4 <= a^3 min(||f||_infinity^2 ||g||^2,
             ||g||_infinity^2 ||f||^2),
max_x ||q_x|| <= a^3 ||f||_infinity ||g||_infinity.
```

The same last bound holds for the pair coefficient. Thus (1) tends to
zero. This is an averaged fluctuation result, not a claim of a uniform
macroscopic Wick estimate for arbitrary smeared field words.

## Degree-controlled dynamics makes the fluctuation an operator result

An expectation tending to zero would not suffice on its own. We next
control commutators with every finite-time physical field, without a
clustering hypothesis or a factor exponential in `v|t|/a`.

A CAR polynomial of degree at most `m` is a linear combination of products
of at most `m` linear creation or annihilation operators, together with
the operator-norm closure of that class. Its coefficients may depend on
spacing and volume. Norm integrals in the Dyson construction belong to
this closure, and all degree bounds below pass to it. If `q=q^*` is a bounded one-particle
operator and `A` has degree at most `m`, then

```
||[dGamma(q),A]|| <= m ||q|| ||A||.                 (2)
```

In finite volume diagonalize `q`. The orbit
`F(t)=exp(it dGamma(q)) A exp(-it dGamma(q))` has frequencies in
`[-m||q||,m||q||]` and `||F(t)||=||A||`. The scalar Bernstein bound,
applied after any norm-one linear functional, proves (2). One direct
form of that bound, for type `R>0`, is the absolutely convergent Riesz
interpolation formula

```
F'(0) = (4R/pi^2) sum_(k in Z)
        (-1)^k F((2k+1)pi/(2R))/(2k+1)^2.
```

For an exponential it follows from the Fourier series of the triangular
wave on `[-R,R]`; linearity gives every finite-frequency orbit. The sum
of absolute coefficients is `R`, which proves (2) also for an
operator-valued orbit. Norm approximation removes finite-support
restrictions. For non-Hermitian `q`, splitting into Hermitian real and
imaginary parts gives the safe constant `2m||q||||A||`.
An onsite pair sum with coefficient supremum `s` has the same safe bound
`2ms||A||`, by diagonalizing the two Hermitian quadratic pairing
generators; their one-quasiparticle frequencies have magnitude at most
`s`. This argument allows Bogoliubov mixing and preserves polynomial
degree. Centering a quadratic operator does not change its commutator.

Write `D=sum_x n_xup n_xdown`. Every nonzero matrix element of a degree-m
monomial connects occupation configurations whose values of `D` differ
by at most `m`. The same interpolation proof gives

```
||[D,A]|| <= m||A||,
degree([D,A]) <= m+2.                              (3)
```

The degree assertion also follows directly from the Leibniz rule and
`[D,c_sigma]=-n_opposite c_sigma`. Free evolution preserves degree and
norm. Thus (3) applies to each free-interaction-picture contact generator,
with no norm estimate on the extensive operator `D` itself.

For a time interval `s`, a degree-m input `A`, and `N>=1`, truncate the
interaction-picture Dyson expansion after orders `0,...,N-1`. Its output
`A_N(s)` has degree at most `m+2(N-1)`. The integral remainder has its
last exact interacting propagator on the outside, where it is isometric;
the remaining nested contact derivations are bounded by (3). Hence

```
||tau_s(A)-A_N(s)||
 <= ||A|| (2 Umax |s|)^N (m/2)_N/N!.               (4)
```

This derivation is first finite-volume. The uniform bound, the previously
proved full-CAR thermodynamic dynamics, and norm approximation pass (4)
to the infinite system. No differentiability of a limiting continuum
field is asserted.

Here is a completely explicit finite-time approximation. Divide
`|t|<=T` into `K=max(1,ceil(8 Umax T))` steps, so
`2Umax |s|<=1/4`. Choose an integer `ell>=1` and at a step of input degree
`m` take `N=m+ell`. Since
`sum_n (m/2)_n (1/2)^n/n! = 2^(m/2)`, (4) is at most
`2^(-ell)||A||`. After the `K` steps, the output `A_app` has

```
M = 3^K (m_initial+ell-1)-(ell-1),
degree(A_app) <= M,
||A_app|| <= (1+2^(-ell))^K ||A||,
||tau_t(A)-A_app|| <= [(1+2^(-ell))^K-1] ||A||.      (5)
```

These estimates can be large, but for each fixed physical time they are
uniform in spacing and spatial volume. Each step's error is propagated
by a norm-preserving exact evolution; the approximant norm bound closes
the displayed recurrence. If `Umax=0`, simply use the exact free field.

For any of the averaged quadratic defects `Y_a` above, its norm is
bounded uniformly and its quadratic coefficient supremum tends to zero
as `a^3`. Equations (2) and (5) give

```
||[Y_a,tau_t(A)]||
 <= 2 M s_a (1+2^(-ell))^K ||A||
    + 2 ||Y_a|| [(1+2^(-ell))^K-1] ||A||.          (6)
```

For the mixed CAR defects absorb their factor two into `s_a`. First let
`a->0` at fixed `ell`, then `ell->infinity`. Thus the commutator tends
to zero for every finite-time smeared field, and also for every finite
product of such fields by the Leibniz rule. Taking
`2^(-ell)<=a^3/(2K)` gives a concrete `O(a^3 log(1/a))` bound at fixed
`T,m_initial,Umax`, with the coefficient in (6). The argument also covers
an initial composite `d(u)`, which has degree three and norm `||u||`.

Finally, for every such bounded polynomial `X_a`,

```
||Y_a X_a Omega_a||
 <= ||[Y_a,X_a]|| + ||X_a|| ||Y_a Omega_a|| -> 0.
```

Together with (1), this proves vanishing in every polynomial GNS
context. The conclusion is an operator identity in the joint limit,
not merely a vanishing vacuum expectation.

## A joint continuum algebra with a resolved composite channel

Adjoin the bounded equal-time composite fields
`d(M_b J_a f)` to the countable extraction used for the existing
continuum field net. The same diagonal argument extends the spacing
subsequence while preserving every already convergent original-field
correlation. The composite norm bound and the cell-average/valley
estimates extend their smearings to all `L^2` functions. Restriction to the
original field algebra preserves its correlations and its original cyclic
GNS subrepresentation. Its strongly continuous positive-energy time group
is not automatically extended to the entire enlarged joint GNS space. Include every
mixed polynomial with the original fields at arbitrary physical times.
The preceding context estimates imply, as operator identities,

```
{C_b(f),C_e(g)^*} = delta_be <f,g>I,
{D_b(f),D_e(g)^*} = delta_be <f,g>I,
{D_b(f),C_e(g)^*} = {D_b(f),C_e(g)} = 0,
{C_b(f),C_e(g)} = {D_b(f),D_e(g)} = 0.
```

Thus the original and composite families generate the CAR algebra over
two copies of the eight-valley one-particle space. The nonlinear source
composite has become a distinct canonical channel under refinement. No
microscopic degree of freedom was added, and the state on this enlarged
algebra is not asserted to be a tensor product or quasifree state.

The original exact equation couples the contact force through `(U/2)d`.
It does not follow that its time integral has a nonzero continuum limit.
The derivative bound on `d(M_b J_a f)` contains kinetic terms of size
`1/a`; the argument has not supplied uniform time continuity for the
composite channel or justified passing time derivatives through the
spacing limit. It can therefore carry spectral weight escaping to
infinite energy, as in the retained moment-escape controls. Conversely,
its nonzero canonical norm prevents simply deleting it in a claimed
finite-spacing operator identity. A proof that its bounded-frequency
spectral weight vanishes or survives is still required to select
Hartree/free or interacting dynamics.

## An actual positive-energy control for loss of the memory channel

The unresolved energy issue has an exact dynamical control. On two
canonical modes `a,b`, choose the one-particle Hamiltonian

```
h=[[E0,kappa],[kappa,Lambda]],
E0>0, Lambda>E0, kappa>0, E0*Lambda>kappa^2.
```

Its second quantization on all four Fock states has the empty positive-
energy ground. The exact equation is `[H,a]=-E0*a-kappa*b`, and the
compressed selfenergy is `Sigma(z)=kappa^2/(z-Lambda)` with constant
memory Gram `kappa^2`. Nevertheless, for `|z|<=Z<Lambda`,

```
|Sigma(z)| <= kappa^2/(Lambda-Z),
||tau_t(a)-exp(-it E0)a||
 <= (2 kappa+kappa^2|t|)/(Lambda-E0).
```

To prove the second bound, put `Delta=Lambda-E0`. The lower eigenvalue
is `E0-delta`, where `0<=delta<=kappa^2/Delta`, and the mixing angle
satisfies `sin(theta)<=kappa/Delta`. Separating the two eigenphases bounds
the one-particle vector error by `delta |t|+2 sin(theta)`. The CAR norm
identity transfers this directly to the annihilator norm. Thus the
finite-energy field becomes free as `Lambda->infinity` while the exact
independent canonical force channel and its Gram remain nonzero.
This control is a separate quadratic model, not a claim that it describes
the actual OPH continuum. It shows why the equal-time channel theorem
alone cannot settle the latter dynamics.

This result provides neither a Lorentz stress tensor nor geometric cap
modular flow. The existing `EinsteinRegisterData` remains the consumer:
its diamond-kernel identification, null matching, area stationarity and
physical normalization cannot be obtained by relabeling this new channel.
The fixed-spacing composite force algebra and its continuum canonical
resolution do not close any of those remaining geometric premises.
