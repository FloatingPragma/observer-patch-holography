# A spatially local inverse response on the interacting infinite lattice

The original zero-chemical-potential Weyl/contact action has a
translation-invariant ground with a **nonzero retarded interaction memory**.
Its full spatial Green operator and inverse are exponentially local on every
compact frequency domain separated from the real axis. Consequently its
inverse response has an analytic spatial symbol and a controlled derivative
expansion there. This supplies a spatial response interface on the actual
infinite interacting ground, rather than extrapolating the special cube.

The result is at fixed lattice spacing and `Im z>=eta>0`. It neither takes
`eta->0` nor computes the selected ground's response values. It establishes
no real mass shell, common Lorentz metric, scattering theorem or Einstein
balance. Its constants are deliberately conservative, not predictions of a
physical correlation length.

## Source, state and Hilbert spaces

Use the original action on `Z^3`, with two CAR modes per site,

```
H = dGamma(h) + U sum_x n_xup n_xdown,
h_(x,x+e_j) = i c sigma_j/2,       c>0, U>=0.
```

The [thermodynamic-ground theorem](../thermodynamic_ground/README.md) supplies
the state `omega` by spatially averaging normalized **full** finite-box
ground projectors before taking a weak-star subsequence. It proves
stationarity and a positive self-adjoint GNS generator `L_omega`, with
`L_omega Omega=0`. At `U>0` its density is between `1/4` and `1` per site.
Its [full-CAR locality estimate](../full_fock_field/README.md) is uniform in
volume, particle number and the onsite coupling. No new reference state,
chemical potential, regulator in the action, or geometry law is introduced.
The imaginary part of `z` regulates only the Laplace transform of a response.

Let `K=ell^2(Z^3;C^2)`. On the doubled GNS Hilbert space define

```
L = L_omega directsum -conjugate(L_omega),
V_(x,sigma) = pi(c_xsigma^dagger) Omega
               directsum conjugate(pi(c_xsigma) Omega).
```

The conjugate second summand retains the correct removal indices. The CAR
identity gives `V^dagger V=I`, so these columns extend to an isometry from
all of `K`. Define `b_xsigma=n_x,opposite c_xsigma` and define `J` from `b`
by the same doubled construction. With the one-site density convention
`rho_x[sigma,tau]=omega(c_xtau^dagger c_xsigma)`, the exact local identities are

```
{b_sigma,b_tau^dagger}
  ={c_sigma,b_tau^dagger}={b_sigma,c_tau^dagger},
R_x = tr(rho_x) I-rho_x,
V^dagger J = J^dagger V = J^dagger J = R = directsum_x R_x.
```

Different-site anticommutators vanish. `algebra.py` checks every spin pair
on the complete four-state onsite Fock space and all different-site pairs
on the complete sixteen-state two-site space. In particular, off-diagonal
spin coherence is retained. Since `0<=rho_x<=I`, the two eigenvalues of
`R_x` are the two eigenvalues of `rho_x` in reversed order. Thus `0<=R<=I`
and `||J||<=1`.

On finite-support vectors the original commutator gives

```
L V = V h + U J,              ||L V|| <= sqrt(3)c+U.
```

The kinetic bound follows from `h(k)=-c sum_j sigma_j sin(k_j)` and the
Pauli anticommutators. Closure of `L` extends this identity and its domain
to every vector in `K`; finite-support approximants converge both as
vectors and after applying `L V`.

## Exact memory, including a nonzero infinite-ground lower bound

Put `P=VV^dagger`, `Q=I-P`,

```
M1=V^dagger L V=h+UR,
W=Q L V=U(J-VR),
B=Q L Q restricted to Q Dom(L).
```

The off-diagonal operator `W V^dagger+V W^dagger` is bounded and
self-adjoint. Subtracting it from `L` gives a self-adjoint operator reduced
by `P,Q`; hence `B` is self-adjoint. This justifies the Schur decomposition
even on the infinite lattice; no finite-dimensional inverse is substituted.
For `Im z>0`,

```
G(z)=V^dagger(z-L)^(-1)V=[zI-M1-Sigma(z)]^(-1),
Sigma(z)=W^dagger(z-B)^(-1)W,
W^dagger W=U^2(R-R^2),
||Sigma(z)||<=U^2/(4 Im z),       -Im Sigma(z)>=0.
```

The same ground has a strictly nonzero memory when `U>0`. The existing
source frame is `R_x=(i sigma_1)^x1(i sigma_2)^x2(i sigma_3)^x3`.
Pairwise anticommutation of the three factors gives
`R_(x+t)=sign(x,t) R_t R_x`, with central sign `+1` or `-1`. The code checks
all 4096 residue-class pairs; exact period four makes this exhaustive for
all integer coordinates. Therefore spatial translations normalize the
global hidden `SU(2)` action `c_x -> R_x^dagger S R_x c_x`.

Each full-box ground projector commutes with that symmetry, and the
exterior vacuum is invariant. A translate of such a state is invariant as
well, because the symmetry is normalized by translations. Averages and
weak-star limits preserve it. Schur's lemma on each two-component site
then gives `rho_x=r I`; translation invariance makes `r` constant. The
earlier density theorem implies `1/8<=r<=1/2`, so

```
(7U^2/64) I <= W^dagger W=U^2 r(1-r) I <= (U^2/4) I.
```

By the spectral theorem, `i y Sigma(i y)` converges **strongly** to
`W^dagger W` as `y->infinity`. In particular, this response cannot have
identically zero memory, or equal the response of a number-conserving static
quadratic Hamiltonian on these same CAR modes. The comparison keeps the
source's number conservation; a pairing Hamiltonian's normal subblock alone
can have memory after its anomalous components are eliminated. This result
does not exclude auxiliary-mode representations, prove many-body scattering, identify hidden spin with
Lorentz spin, or imply that local ground density matrices are faithful.

## A nonzero connected four-point function on that same ground

The full multi-time field, unlike its two-point response alone, has a
quantitative non-Gaussian witness. Fix a site and spin, put
`a=c_xsigma`, `n=n_x,opposite`, and define

```
A(t)={tau_t(a),a^dagger},       D(t)=A(t)-omega(A(t)) I,
F(t)=omega(D(t)^dagger D(t)).
```

There is no onsite kinetic term. The exact local commutator consequently
gives `A(0)=I` and `D'(0)=-i U(n-r I)`. Since `n^2=n`,

```
lim_(t->0) F(t)/t^2 = U^2 r(1-r) >= 7 U^2/64.
```

This is also an explicit finite-time bound. Write `delta=[H,.]` on local
operators and `M=9c^2+12cU+U^2`. The six incident hopping terms give
`||delta(a)||<=3c+U`, `||delta(n)||<=6c`, and

```
||delta^2(a)|| <= 3c(3c+U)+U(6c+3c+U) = M.
```

These are bounds in the full CAR algebra, independent of volume. The
second-order integral remainder gives
`||tau_t(a)-a-it delta(a)||<=Mt^2/2`. Taking the anticommutator and
subtracting its expectation yields
`||D(t)+itU(n-r I)||<=2Mt^2`. Therefore

```
sqrt(F(t)) >= U |t| sqrt(7)/8 - 2Mt^2,
F(t) >= U^2 t^2/64 > 0,       0<|t|<=U/(16M), U>0.
```

To exhibit the connected correlation explicitly, for odd operators with
zero one-point functions define

```
kappa4(B1,B2,B3,B4) = omega(B1 B2 B3 B4)
 - omega(B1 B2) omega(B3 B4)
 + omega(B1 B3) omega(B2 B4)
 - omega(B1 B4) omega(B2 B3).
```

Set `(b1,b2,b3,b4)=(a,tau_t(a)^dagger,tau_t(a),a^dagger)`.
Then `F(t)` is the sum of `kappa4` for the four words `1234`, `1243`,
`2134`, `2143`. The first disconnected products sum to
`omega(A(t)^dagger)omega(A(t))`; every crossed contraction cancels
pairwise. In particular no anomalous pairing is omitted. One of these
four connected functions has magnitude at least `U^2 t^2/256` on the
displayed nonzero time interval.

A quasifree state with linear fermionic fields obeys Wick factorization
at every time. This strict witness excludes any such representation of
the complete multi-time field, even with pairing or auxiliary Gaussian
modes. It does not exclude a nonlinear encoding in a larger Gaussian
system. More strongly, any quadratic evolution on an enlarged CAR algebra
maps a linearly embedded observed mode to a linear combination of creation
and annihilation operators. Its anticommutator with the original mode is
therefore scalar, so `D(t)=0` in every reference state. The strict witness
excludes such a linear quadratic realization even with a non-Gaussian
reference. The witness establishes genuine non-Gaussian dynamics at fixed
spacing, not scattering, a continuum interaction, a mass shell or a
stress/geometry identification. `connected.py` independently checks the
local commutator, the four-word identity and its negative controls.

## The complete ordered Euclidean hierarchy at fixed spacing

The same selected ground supplies every ordered local Euclidean n-point
function, not only the retarded two-point response. Let `L_omega>=0` be
its already proved positive GNS generator, with `L_omega Omega=0`.
For bounded local CAR operators `A_1,...,A_n`, keep their displayed order
and take `s_1<=...<=s_n`. Define

```
S(A_1,s_1;...;A_n,s_n)
 = <Omega, A_1 exp[-(s_2-s_1)L_omega] A_2 ...
                  exp[-(s_n-s_(n-1))L_omega] A_n Omega>.
```

The represented operators `pi(A_j)` are written as `A_j` here. Spectral
calculus supplies every contraction, including the identity for a zero
gap, and gives `|S|<=product_j ||A_j||`. This defines all orders on one
state and one generator, with no separate subsequence for each n-point
function. It is invariant under common shifts of all Euclidean times.
It is continuous in nonnegative time gaps, and holomorphic when each
complex gap has strictly positive real part. These facts follow from
strong continuity and bounded holomorphic spectral calculus; they do
not require boundedness of the full infinite-volume generator.

The finite approximants must use **one shared ground and translation
label for the entire word**. In each finite-box summand, translate every
`A_j` by the same `x`, use the same full ground projector `rho_L`, and
insert the same positive generator `H_L-E_L` between every pair of
operators. The exterior Hamiltonian is zero. Average these whole
correlations over `x`, precisely as in `thermodynamic_ground`; multiplying
separately averaged factors would describe a different hierarchy.

To prove convergence of these finite correlations along the existing
subsequence, use the positive-generator Poisson identity, for `s>0`,

```
exp(-s L) = integral_R [s/(pi(s^2+t^2))] exp(it L) dt,   L>=0.
```

The integral is a strong integral of bounded operators; its scalar
spectral identity is `exp(-s|lambda|)`, which equals `exp(-s lambda)`
only on the nonnegative spectrum. Insert this formula independently
for every strictly positive gap; omit zero gaps. Fubini is justified
by the product probability density and the uniform bound
`product_j ||A_j||`. The resulting real-time integrand is

```
omega(A_1 tau_(t_1)(A_2) tau_(t_1+t_2)(A_3) ...
                         tau_(t_1+...+t_(n-1))(A_n)).
```

Its finite comparator uses the same shared translation label and finite
Hamiltonian throughout. On a compact box of the `t_j`, every cumulative
time is bounded. The earlier full-CAR local-dynamics estimate controls
each factor uniformly for good translates; telescoping the finite
product sums those errors. Bad translates have vanishing fraction and
each entire integrand is bounded by the same product of norms. Finally,
the infinite-dynamics products form a norm-compact set on that time box,
so weak-star state convergence is uniform there by a finite epsilon-net.
This proves uniform compact-time convergence of the whole integrand.
Dominated convergence against the fixed product Poisson measure proves
convergence of every ordered Euclidean n-point function to `S`.

There is an explicit integral-tail control. For each positive gap `s_j`
and cutoff `T_j>0`,

```
probability(|t_j|>T_j) = (2/pi) arctan(s_j/T_j)
                     <= 2s_j/(pi T_j).
```

The product measure outside the cutoff box has probability at most the
sum of these bounds, capped by one. Multiply by `product ||A_j||` for
the tail of one correlation, or by twice that product when comparing
finite and infinite correlations. The calculator returns the rational
upper bound obtained from `pi>3`. This is an integration-tail rate,
not a rate for the unknown ground-state subsequence.

There is exact time-reflection positivity for these ordered words. For
a positive-time word `F=(A_1,t_1;...;A_n,t_n)`, with
`0<=t_1<=...<=t_n`, put

```
Psi_F = exp(-t_1 L_omega) A_1 exp[-(t_2-t_1)L_omega] A_2 ... A_n Omega.
```

Reflecting the left word means reversing its order and taking adjoints:
its negative-time list is
`(A_n^*,-t_n;...;A_1^*,-t_1)`. Concatenating this reflected list with a
right positive-time word gives exactly `<Psi_F,Psi_G>`; the central
semigroup has duration `t_1+s_1`. Consequently every finite matrix of
reflected correlations is positive semidefinite, and
`sum_ij conjugate(z_i) z_j S(theta F_i,F_j)=||sum_j z_j Psi_Fj||^2>=0`.
This holds for arbitrary local operators, including odd CAR words and
complex coefficients. No CAR factors are sorted or commuted during
reflection. An extra bosonic permutation rule or an inserted fermionic
sign would generally destroy this Gram identity.

`euclidean.py` checks the reflection identity by two separate finite
constructions: ordered reflected products and half-word vectors. Its
rational two-mode witness checks all 31 principal minors of a five-word
Gram matrix, a negative odd-word cross entry `-1/64`, and the exact
polynomial norm `209/256`. These checks audit order and sign conventions;
the infinite-volume argument is the spectral and convergence proof
above. The hierarchy is not claimed to be Gaussian: the preceding
nonzero connected real-time witness belongs to this same state and
dynamics. No `O(4)` covariance, bosonic permutation symmetry, spatial
continuum, or complete relativistic Osterwalder--Schrader axioms are
inferred. The individual limiting values and a unique ground are still
not computed.

## A uniformly invertible spatial Green operator

Fix `Im z>=eta>0`, `|z|<=Z`, with `Z>=eta`. For a unit `f in K`, put
`v=Vf`. The spectral resolvent identity and Cauchy-Schwarz give

```
-Im <f,G(z)f> = Im(z) ||(z-L)^(-1)v||^2,
1 = |< (conjugate(z)-L)v, (z-L)^(-1)v >|,
||(conjugate(z)-L)v|| <= Z+sqrt(3)c+U <= Z+2c+U.
```

Consequently, with `B0=(Z+2c+U)^2/eta`,

```
-Im G(z) >= I/B0,       ||G(z)||<=1/eta,
||G(z)^(-1)|| <= B0.
```

For completeness, the imaginary-part inequality bounds both `G` and
`G^dagger` below by `1/B0`: apply it to their quadratic forms and take
absolute values. The range of `G` is closed, while `ker(G^dagger)=0`
makes that range dense. Thus it is surjective as well as injective.
The lower bound is not merely a finite-volume determinant argument.

## Direct spatial decay from the original dynamics

Use the site distance `d=||x-y||_infinity`. If `d>=2`, let `R0=d-1`.
The radius-`R0` approximation to `tau_t(c_xsigma)` is odd and excludes
site `y`, so its anticommutator with `c_ytau^dagger` is exactly zero.
The prior boundary bound therefore bounds a two-by-two response block by
four times

```
B_R0(t) <= 6c(2R0+1)^2 t 2^(-R0) [exp(48ct)-1].
```

Split the Laplace integral at `T=R0/(192c)`. On the first part,
`log(2)>=1/2` and `48ct<=R0/4` give the bound
`24c(2R0+1)^2 exp(-R0/4)/eta^2`. On the tail use the uniform block
bound four, giving `4 exp[-eta R0/(192c)]/eta`. Define

```
mu=min(1/8, eta/(384c)),       A=192c/eta^2+8/eta.
```

Since `exp(2mu)<2`, both parts give, for all `d` including zero and one,

```
||G_xy(z)|| <= A(1+d)^2 exp(-2mu d).
```

At zero and one use the previously proved operator bound `1/eta`.
This estimate uses no decay or clustering assumption on the ground state:
the distant anticommutator is controlled in operator norm before taking
its expectation.

## Exponential inverse locality with bounded weights

There are `24d^2+2` sites at distance `d>=1`. Let

```
S(mu)=sum_(d>=1)(24d^2+2)d(1+d)^2 exp(-mu d).
```

The calculator uses the explicit rational upper bound `S_hat` obtained by
replacing `exp(-mu)` by `q=1/(1+mu)`. In terms of
`T_n(q)=sum_(d>=1)d^n q^d`,

```
S_hat=24T_5+48T_4+26T_3+4T_2+2T_1,
alpha=min(mu, 1/(2 B0 A S_hat)) > 0.
```

Eulerian polynomial identities evaluate every `T_n` exactly. For any
**bounded**, real, one-Lipschitz function `f` on the lattice, conjugation
by `D_f=exp(alpha f)` changes a block of `G` by at most
`alpha A d(1+d)^2 exp(-mu d)`. The row and column Schur bounds imply

```
||D_f G D_f^(-1)-G|| <= alpha A S_hat <= 1/(2B0).
```

A Neumann series and `||G^(-1)||<=B0` now give
`||D_f G^(-1)D_f^(-1)||<=2B0`. For each desired pair `x,y`, choose the
bounded function `f(u)=min(||u-y||_infinity,R)` with `R>=d`. It follows that

```
||(G(z)^(-1))_xy|| <= 2B0 exp(-alpha d).
```

No unbounded similarity transformation has been assumed. Finally,
`Sigma=zI-h-UR-G^(-1)` differs from the inverse by onsite and nearest-
neighbor terms only, so

```
||Sigma_xy(z)|| <= (Z+U+2B0+2c) exp(-alpha d).
```

## Spatial and frequency derivative expansions

Translation invariance makes these convolution kernels. Their absolutely
convergent Fourier series define two-by-two symbols. The inverse symbol
`K(k,z)=G(k,z)^(-1)` is analytic for `||Im k||_1<alpha`, by normal
convergence; on real `k` it is the inverse of the Green symbol. For a real
increment `q`, truncating the Taylor expansion about any real `k0` at
total degree `p` gives

```
||remainder|| <= 2B0 ||q||_1^(p+1) M_(p+1)(alpha)/(p+1)!,
M_n(a)=sum_(r in Z^3)||r||_infinity^n exp(-a||r||_infinity).
```

The origin contributes one only when `n=0`. The code bounds `M_n` by
`24T_(n+2)(1/(1+a))+2T_n(1/(1+a))`, plus that origin term.
The remainder follows from the integral remainder of `exp(-i q.r)`;
it has no exponentially growing factor for real `q`. All derivatives
are absolutely convergent sums with corresponding moment bounds.

There is also a complex-frequency expansion. A circle of radius `eta/2`
around a point in the original domain lies in
`Im z>=eta/2, |z|<=Z+eta/2`. Let `B_star` be the inverse bound for that
larger domain. For real `k` and `|w|<eta/2`, Cauchy's integral formula gives

```
||frequency remainder after degree p||
 <= B_star (2|w|/eta)^(p+1)/(1-2|w|/eta).
```

Thus the local spatial coefficients and their frequency dependence exist
with explicit error bounds on the actual interacting infinite ground.
Their numerical values are not computed, because the underlying ground
functional remains selected only by a subsequence. The constants worsen
as `eta` tends to zero. These controlled complex-frequency symbols cannot
yet be inserted as real hyperbolic action coefficients into
`wave_density_metric` or `geometric_matter_join.source_tetrad`. That step
still requires real-frequency control, common matter coefficients and
the full action normalization, followed by the existing joint modular,
stress, area and Einstein premises.

## Boundary of the current refinement estimates

For a possible refinement `c=v/a` with fixed `v,U,Z,eta>0` and lattice
spacing `a->0`, these specific bounds give `mu=O(a)`, `A=O(a^-1)`,
`B0=O(a^-2)`, `S_hat=O(a^-6)` and `alpha=O(a^9)`. Thus the proved
physical-momentum analyticity width `alpha/a=O(a^8)` collapses. The
connected-witness time window `U/(16M)=O(a^2)` also shrinks. These are
limitations of the retained estimates, not a no-go theorem for a
refinement of the action. Fixed-spacing hierarchy existence and time
reflection positivity do not supply uniform continuum estimates.

## Replay

`python -W error -m pytest -q model/retarded_spatial_locality` checks the
exact source identities, complete frame residues, rational estimates and
domain controls. `python -W error -m model.retarded_spatial_locality --check`
replays the packet in a fresh process. The proof above uses standard
self-adjoint spectral theory, the block Schur identity, the operator Schur
test, norm-convergent Neumann series and Fourier/Cauchy analysis; those
steps are stated rather than inferred from finite tests.
