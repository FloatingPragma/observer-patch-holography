# Actual regional modular flow from the free lattice source

For every fixed bounded ball, the exact finite regional density matrices of
the `U=0` source have modular flows converging to the modular flow of the
same eight-Weyl vacuum. Convergence holds on every fixed one-particle probe,
uniformly on compact modular-time intervals, and hence on every finite CAR
word. No spectral gap away from covariance eigenvalues zero or one is
assumed. The unbounded modular generators converge on an explicit graph
core, with matched lattice inputs. A family of finite, nontrivial relative
entropies converges as well.

This strengthens FA94's convergence of selected stress responses to an
actual regional modular object. The region is fixed throughout this
theorem. A shrinking-ball/area/geometry limit, nonzero-U composite limit,
or a native implementation of these density-matrix operations is not proved.

## Source, region and the eight-valley comparison maps

Take precisely FA94's infinite free source ground, `U=0,c=v/a`, and physical
coordinate `x0=v*t`. Its one-particle filled-band projector is `Pi_a` on
`ell2(Z3;C2)`. Let `D=B_R` be a fixed ball with `R>0`. The finite site set
`D_a={n:a n in D}` defines the projection `E_a` and finite space
`H_a=E_a ell2`. Its exact regional covariance is

```
A_a=E_a Pi_a E_a | H_a.
```

The continuum space is `H=direct_sum_(b in {0,1}3) L2(D;C2)` with zero
extension outside D. In the actual source frames of FA87, write `Pi_b`
for each Weyl negative-band projector and
`C=direct_sum_b chi_D Pi_b chi_D | H`. Let `J_a` be the centered-cell
average contraction and `M_b u(n)=(-1)^(b.n)u(n)`. Define

```
T_a f=E_a sum_b M_b J_a f_b,             ||T_a||<=sqrt(8)<3.       (1)
```

One may instead conjugate each summand by FA87's constant spin frame;
all statements are unchanged. The eight images are not orthogonal at
finite a. FA94 and the oscillatory cell estimate in `continuum_field_net`
give, for fixed f,g in H,

```
<T_a f,T_a g> -> <f,g>,
R_a(f):=A_a T_a f-T_a C f -> 0 in H_a.                          (2)
```

Here is a direct proof of the new regional part. Begin with smooth probes
compactly contained in D; the input cell averages then lie in D_a once a
is small. FA94 controls each full-space projector intertwiner. Multiplying
its output by E_a is a contraction. The remaining difference is
`E_a J_a (1-chi_D) Pi_b f_b`, supported on cells crossing the boundary.
Their exterior portion lies in `R<=|x|<R+a`. If a<=R its volume is at
most `16 pi R^2 a`, so for a single such probe

```
||R_a(f)|| <= (a^2/epsilon)|||p|^3 f||
       +epsilon^(3/2)||f||_1+8R sqrt(a)||Pi_b f||_infinity.       (3)
```

The output supremum is finite by H2 Sobolev control. For example Fourier
Cauchy--Schwarz gives, for any length l>0,
`||Pi_b f||_infinity <= l^(-3/2)(||f||+l^2|||p|^2 f||)`;
the omitted exact factor `1/sqrt(8pi)` is less than one.
Choosing epsilon proportional to a gives an explicit O(sqrt(a)) bound
for each fixed collared smooth probe. Sum (3) over the eight components.
Density and the uniform bounds in (1) extend (2) to all H. The same
boundary-layer argument and the oscillatory cell estimate prove the
first line of (2), initially for collared probes and then by density.
No aliasing or band-limited sampling hypothesis is used.

## Both regional covariances have no endpoint eigenvectors

If `<f,Cf>=0`, zero extension of f belongs to the positive Weyl band;
if `<f,(1-C)f>=0`, it belongs to the negative band. Either would be a
compactly supported pure-band spinor. Such a nonzero spinor does not exist.
Its Fourier transform is entire. Fix real transverse momenta with
`r^2=p2^2+p3^2>0`. In the first complex momentum z the band equation is
`h(z)F(z)=+/-sqrt(z^2+r^2)F(z)`. Continue around exactly one simple branch
point. The entire left side and F are unchanged, while the square root
changes sign. Hence F vanishes, first on this complex line and then for
all momenta by analyticity. Apply this to every valley independently.
This is the compact-support argument already used for FA82's outside
cyclicity, now applied to the regional covariance.

The lattice argument is equally explicit, and includes all doubler zeroes.
A finite-support lattice spinor has a trigonometric-polynomial Fourier
transform. Fix real k2,k3 with `s=sin(k2)^2+sin(k3)^2>0`. The band equation
uses `sqrt(sin(z)^2+s)`. In w=exp(iz) its branch polynomial is

```
w^4-(2+4s)w^2+1,
discriminant=4096 s^2(s+1)^2>0.                                (4)
```

Its roots are nonzero and simple; equivalently `2sin(z)cos(z)` is nonzero
at each zero because `sin(z)^2=-s` and `cos(z)^2=1+s`. Continuation around
one such zero again forces the entire trigonometric polynomial to vanish.
The transverse momenta with s>0 are dense. Consequently

```
ker C=ker(1-C)=0,             0<A_a<1 on each nonempty H_a.       (5)
```

The strict finite inequalities give an actual faithful finite regional
density matrix. They give no lower eigenvalue bound uniform in a.

## Functional calculus with an explicit finite-residual bound

Since A_a and C are positive contractions, (2) implies
`p(A_a)T_a f-T_a p(C)f ->0` for every polynomial p, by telescoping powers.
Uniform polynomial approximation extends this to every continuous
function on [0,1]. This is convergence between H and the varying H_a;
it does not assert that `T_a^*T_a` tends strongly to the identity.

For example let
`k_delta(x)=log((1-x+delta)/(x+delta))`, delta>0, and
`r=(1+2delta)^(-1)`. The exact uniformly convergent series is

```
k_delta(x)=2 sum_(j>=0) r^(2j+1)(1-2x)^(2j+1)/(2j+1).
```

After N terms its uniform remainder is bounded by
`tau_N=2r^(2N+1)/[(2N+1)(1-r^2)]`. With M=3, the regularized-log
intertwiner obeys

```
||k_delta(A_a)T_a f-T_a k_delta(C)f||
 <=2M||f|| tau_N
   +[4r/(1-r^2)] max_(0<=j<2N-1)||R_a((1-2C)^j f)||.             (6)
```

Thus every fixed delta has controlled convergence; the right side uses
finitely many actual covariance residuals. Equation (3) supplies a direct
spacing rate on collared smooth probes. It is not asserted to apply with
the same smoothness constants to all iterated sharply cut-off probes
appearing in (6). For those probes (2), density and (6) give convergence,
without an invented uniform power law in a.

The embedded complement matters. If W_a is a one-valley cell isometry
into an ambient L2 space and P_a=W_a W_a^*, then for continuous F

```
F(W_a A_a W_a^*)=W_a F(A_a)W_a^*+F(0)(I-P_a).                   (7)
```

The second term is generally nonzero for k_delta. It is an artificial
zero-covariance sector, not a regional lattice modular mode. We apply
functional calculus on H_a and H themselves, as in (6); we never assign
the padded zero modes a physical logarithm.

## Removing the regulator: the actual modular flow and generator

Let `k(x)=log((1-x)/x)` on (0,1). By (5), k(C) is a selfadjoint operator
with a dense domain, and k(A_a) is an ordinary finite matrix. Put
`w(x)=x(1-x)`. The space `D_core=range(w(C))` is dense and is a graph
core for k(C): cut any vector in its graph domain to
`[1/n,1-1/n]`, where division by w is bounded, and then take n to infinity.

For f=w(C)g use the matched lattice vector `f_a=w(A_a)T_a g`.
It converges to T_a f, since

```
||f_a-T_a f|| <= ||R_a(g)||+||R_a(Cg)||.                         (8)
```

The functions `G(x)=w(x)k(x)` and
`F_t(x)=w(x)exp(-it k(x))`, assigned value zero at both endpoints,
are continuous on [0,1]. Hence

```
k(A_a) f_a - T_a k(C) f ->0,
exp(-it k(A_a)) f_a - T_a exp(-it k(C)) f ->0.                   (9)
```

This is genuine graph convergence on the stated core, not a claim that
`k(A_a)T_a f` converges for every f in the logarithm's domain.
Moreover `F_t'(x)=(1-2x+it)exp(-it k(x))`; therefore its Lipschitz
constant is at most `1+T` for |t|<=T. The degree-N Bernstein polynomial
has error at most `(1+T)/(2sqrt(N))`, and the sum of the absolute values
of its monomial coefficients is at most `3^N/4`. Power telescoping gives

```
sup_(|t|<=T)||exp(-it k(A_a)) f_a-T_a exp(-it k(C)) f||
 <=M(1+T)||g||/sqrt(N)
    +(N 3^N/4) max_(0<=j<N)||R_a(C^j g)||.                       (10)
```

Use (8) to replace f_a on the left by T_a f. Density of D_core and the
unitarity of both groups extend this convergence to **every fixed f in H**,
uniformly on every compact modular-time interval. Endpoint cutoffs provide
an equivalent proof: on the core the endpoint spectral norm is at most
`epsilon||g||`, and it tends to zero on every H vector by (5).

There is also a true joint delta/a limit. The elementary inequality
`log(1+y)<=y` gives

```
w(x)|k_delta(x)-k(x)|<=delta.
```

The extra core-generator error is at most `M delta||g||`; the extra
core-flow error on |t|<=T is at most `M T delta||g||`. Thus any positive
delta(a)->0 works with the preceding core comparison, and density extends
the regularized-flow result to all fixed probes. No relation between a
and delta and no uniform covariance gap is required. Operator-norm
convergence over all probes, or convergence of bare logarithms on all
vectors, is not asserted.

## Modular CAR automorphisms and finite relative entropy

At finite a the regional state is exactly
`rho_a=Z_a^(-1) exp(-dGamma(k(A_a)))`. With the repository's antilinear
annihilator convention its modular automorphism is

```
sigma_t^a(c_a(f))=rho_a^(it)c_a(f)rho_a^(-it)
                 =c_a(exp(-it k(A_a))f).                       (11)
```

For the continuum gauge-invariant quasifree state, the same relation is
the modular automorphism of its regional von Neumann algebra. The precise
standard input is Araki, *On Quasifree States of CAR and Bogoliubov
Automorphisms*, Publ. RIMS 6 (1970/71), 385--442,
[Corollary 4.10, Lemma 4.11 and Theorem 3, pp.393--404](https://ems.press/content/serial-article-files/41673),
[DOI](https://doi.org/10.2977/prims/1195193913).
In the doubled self-dual notation the covariance blocks are C and
`1-conjugate(C)` on H and its conjugate space. Condition (5) excludes
endpoint eigenvectors, so Corollary 4.10 gives a cyclic separating vector.
The 1/2 eigenspace has doubled, hence even or infinite, multiplicity;
Lemma 4.11 supplies the factor condition in Theorem 3. Equivalently, Wick
contraction verifies the KMS condition on spectrally compact analytic
CAR words, whose norm-dense continuation determines the group. This is
the same regional Weyl vacuum as FA87, not a new thermal state.

The CAR identity `||c(f)||=||f||` turns (9)--(10) into the precise comparison
`||c_a(exp(-it k(A_a))T_a f)-c_a(T_a exp(-it k(C))f)||->0`.
Both operators in this subtraction belong to the same source algebra;
operators in unrelated Hilbert-space representations are not subtracted.
Telescoping products gives convergence for every fixed finite
word, including adjoints and finitely many modular times. The asymptotic
isometry in (2), the covariance relation, and Wick's rule also give all
such state-correlation limits. No convergence in a single arbitrarily
identified many-body Hilbert space is presumed.

There is a concrete finite-entropy consequence. Choose nonzero
`f=w(C)g/||w(C)g||`, the corresponding normalized f_a, and
`U_theta=exp(i theta c^*(f)c(f))`. These are actual bounded local unitaries.
For all sufficiently small a the matching core vector is nonzero, by
asymptotic isometry. For `omega_theta(X)=omega(U_theta^* X U_theta)`,

```
D(omega_theta||omega)
 =2(1-cos(theta)) [<f,Cf><f,k(C)f>-<f,C k(C)f>] < infinity.       (12)
```

The finite regional matrix states obey the identical formula with A_a
and f_a. Their relative entropies converge to (12) by (2) and the graph
limit (9). The bracket is nonnegative because k is decreasing: it is
minus the covariance of x and k(x) in f's spectral probability measure.
There is no trace-class or absolute-continuum-entropy assumption.

For clarity, (12) can be derived without a thermodynamic-limit entropy
interchange. In the faithful regional GNS representation let Omega be
the reference vector and Psi=U_theta Omega. The relative Tomita operator
`S_(Omega,Psi)` equals `U_theta S_Omega` on its core, since it sends
`A U_theta Omega` to `A^* Omega`. Its relative modular operator is
therefore `Delta_Omega`. The defining relative entropy equals
`-<Psi,log(Delta_Omega)Psi>
 =i (d/dt) omega(U_theta^* sigma_t(U_theta))|_(t=0)`.
The derivative exists because f belongs to Dom k(C); (11) and Wick's
rule evaluate it as (12). Full finite-Fock controls independently verify
the modular sign and this entropy expression.

## Scope and verification

This supplies actual finite-density regional modular flows, a dense
generator graph core, regulator removal and convergent finite relative
entropy readouts for the explicitly free source. It strengthens the
matter/modular part of G5 while retaining the existing G1--G6 theorem
stack. It supplies no shrinking-region chart/area limit, boosted family
uniformity, Newton coefficient, generalized-entropy stationarity or
geometric null balance. Nonzero-U regional modular convergence and
interacting/massive particle identification are separate requirements.

Run `python -W error -m pytest model/lattice_modular_limit -q`, and
`python -m model.lattice_modular_limit --write` / `--check`.
The receipt binds only the new package and the narrow mathematical inputs;
inherited broad custody is recorded separately and is not silently repaired.
The exact historical FA94 receipt at checkpoint `ba312dbf0` is bundled under
`historical/`, checked against its fixed SHA-256, and source-pinned. Replay
therefore works in shallow checkouts without consulting Git history.
