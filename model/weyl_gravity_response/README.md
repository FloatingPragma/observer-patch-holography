# Explicit stress directions and finite modular tomography

The same free Weyl vacuum used in FA87 admits an explicit nine-dimensional
family of local normal perturbations. Their stress tangents span all real
symmetric traceless tensors at a point. Nine geometrically normalized,
boosted-ball modular charges distinguish these directions at finite radius,
with an explicit inverse bound. The vacuum-relative entropy first law holds
on this whole variation space, with a uniform quadratic remainder.

This closes a **matter-response and finite variation-space premise** beyond
the existence of a stress tensor and a modular group. It does not equate a
geometry charge to stress, supply an area term, establish generalized-entropy
stationarity, or instantiate the complete Einstein register.

## Existing target and source

`Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean` already proves null
and timelike tensor determination. We do not reprove that tensor theorem.
`RegisterSurface.lean` still requires an actual variation space, a common
stress for entropy and directional responses, the geometric null balance,
and the other named supplies. This package constructs the former matter
data; the geometric balance remains a distinct premise. `SmallBall.lean`
already fixes the coefficient `8*pi^2*R^4/15`; FA87 identifies its continuum
matter kernel. Here the kernel acts on an explicit, nondegenerate family.

The input is the complete centered free eight-Weyl continuum hierarchy,
with `c=v/a`, bounded nonnegative `U(a)`, and the ground of `H_a-U(a)N/2`.
We use the local spin-frame map in
[FA87](../weyl_geometric_modular/README.md). Excite one `h=+sigma.p` summand;
the other seven remain in the same product vacuum. No valley is removed
from the theory or the total stress. Coordinates are `x0=v*t`, source
`hbar=1`, and `eta=diag(1,-1,-1,-1)`. Lowering both indices agrees with
the geometry stack's opposite overall signature. No dimensional or Newton
calibration is added.

## A shell which interpolates the required Cauchy jets

Let `L>0`, and let `P_-=(I-sigma.p/|p|)/2` be the filled-band covariance.
Away from its light cone the negative-band time kernel is

```
K_-(t,z)=-i(t I+sigma.z)/(2*pi^2*(|z|^2-t^2)^2).             (1)
```

For completeness, Fourier transform `exp(i|p|t)/(2|p|)` with an exponential
regulator. Its spacelike limit is `1/[4*pi^2*(|z|^2-t^2)]`.
Applying `-i partial_t+i sigma.grad` gives (1), including its sign and
normalization. Direct differentiation checks the Weyl equation
`partial_t K+sigma.grad K=0`. We only integrate (1) at a positive spacelike
distance, so no light-cone distribution is evaluated pointwise.

Put `s=|y|/L`, `omega=y/|y|`, and set

```
phi(s)=s^2(s-1)^4(2-s)^4       for 1<=s<=2, zero otherwise,
g_n(y)=phi(s) sigma.omega [a+sum_j omega_j b_j].             (2)
```

This compact shell function is `C^3` across its two boundaries and is in
`H^3`. Choose any unit real `n` and unit spinor `u` with
`sigma.n u=-u`, and choose

```
a=-210*pi*i*u,       b_j=(945*pi/4)(n_j I-sigma_j)u.         (3)
```

The elementary angular integrals are
`integral omega_j=0` and `integral omega_i omega_j=(4*pi/3)delta_ij`.
The exact radial integrals are

```
A=integral_1^2 phi(s)/s ds=1/420,
B=integral_1^2 phi(s)/s^2 ds=1/630,
J=integral_1^2 s^2 phi(s)^2 ds=37949/669278610.              (4)
```

Write `F_n=P_-g_n` at time zero, and `S=sum sigma_j b_j`. Differentiating
(1) under the shell integral gives

```
F_n(0)=(2i/pi) A a=u,
partial_j F_n(0)=-(2i/(3*pi*L)) B(sigma_j S-4b_j)
               =i n_j u/L,
partial_0 F_n(0)=i u/L.                                  (5)
```

The last equality is the Weyl equation, not an independent chosen jet.
In fact the map `b -> (sigma_j S-4b_j)_j` is invertible: its inverse on
data `d_j` is `b_j=-(d_j+sigma_j sum_k sigma_k d_k)/4`.
Thus the shell construction supplies arbitrary spatial first jets as well
as the displayed null-wave jets.

The shell norm is independent of `n`:

```
N_L^2=||g_n||^2
 =4*pi*L^3*J [||a||^2+(1/3)sum ||b_j||^2]
 =(150088295/4249388)*pi^3*L^3.                           (6)
```

Here `sum ||(n_j I-sigma_j)u||^2=6`. All constants in (3)--(6), every
chosen helicity vector, and the Lorentz frame matrices are checked exactly.

## Positive local states and the sign of their stress

Let `f_n=g_n/N_L`, `r_n=<f_n,P_-f_n>`, and
`B_n=c^*(f_n)c(f_n)-r_n I`. The number operator is a projection, so
`B_n=B_n^*`, `||B_n||<=1`, and `omega(B_n)=0`. The number `r_n` is the
explicit covariance integral of the known shell and filled-band projector;
no guessed occupation or thermal state is used.

The shell at time zero is spacelike to the whole diamond
`D_{L/4}={|x0|+|x|<L/4}`. The even operator `B_n` commutes with that entire
field algebra. For real parameters `theta_j` with `a=||theta||_1<1`, put
`B(theta)=sum theta_j B_nj` and define on this common local algebra

```
rho_theta(X)=omega((I+B(theta)) X).                       (7)
```

Commutation, `I+B(theta)>=0`, and zero mean prove positivity and unit
normalization. These are normal states. A global vector implementation
uses `sqrt(I+B(theta)) Omega`, but the assertion that (7) is affine in
theta and has vacuum as its midpoint is local. A nontrivial affine two-sided
global state family cannot have a pure global vacuum as its midpoint; no
such affine global claim is made here. The displayed global vector
extension is generally nonlinear in theta.

Inside the diamond, free propagation of the shell itself is zero.
Therefore `P_+ exp(-ith)g=-P_- exp(-ith)g` there. Wick's rule gives the
exact local covariance tangent

```
delta_n rho(x,y)=-F_n(x) F_n(y)^*/N_L^2.                  (8)
```

This minus sign is essential: the shell removes a negative-band wave in
the local covariance. The plane-wave jet (5) has momentum
`p=(-1,n)/L` and Weyl current `j=(1,-n)`. Its classical contravariant
stress is `(p*j^T+j*p^T)/2=-(1/L)j*j^T`. Combining this with (8) and
lowering both indices gives the **positive** tangent

```
delta_n T_{mu nu}(0)=gamma_L ell_mu ell_nu,
ell=(1,n),       gamma_L=4249388/(150088295*pi^3*L^4)>0.    (9)
```

All other field summands contribute zero variation. This is the same
normal-ordered stress as FA87; there is no added classical matter source.

## Nine independent directions and the common first law

Choose the six axis directions and three rational directions

```
±e1, ±e2, ±e3, (3/5,4/5,0), (3/5,0,4/5), (0,3/5,4/5).   (10)
```

In coordinates `(T00,T01,T02,T03,T11,T22,T12,T13,T23)`, with
`T33=T00-T11-T22`, the nine columns `ell_j ell_j^T` have determinant
`13824/15625`. Thus (9) spans every real symmetric traceless stress
variation at the center. The metric-trace direction is absent, exactly
as required by the free massless stress; no full ten-component span is
claimed. In particular the existing null tomography can interrogate a
nondegenerate matter family, rather than a single vacuum with zero stress.

The relative-entropy argument is uniform over this entire variation space.
For `a>0`, the two normalized local states
`rho_±(X)=omega((I±B(theta)/a)X)` have vacuum as their equal-weight mean,
and `rho_theta=((1+a)rho_++(1-a)rho_-)/2`. Normal-channel monotonicity,
as imported and scoped in FA87, gives

```
0<=S(rho_theta||omega)<=[(1+a)log(1+a)+(1-a)log(1-a)]/2
                     <=a^2/[2(1-r0^2)],   a<=r0<1.        (11)
```

Restriction to any smaller ball algebra preserves this bound. With the
FA87 vacuum-relative definition `Delta S=Delta<K>-S(rho||omega)`, (11)
proves a Frechet first law on `V=R^9`, not just one directional derivative.
The energy/modular pairings and entropy differential all use (8)--(9).
The central/area contribution has not been identified; setting a matter
scalar offset to zero does not establish the geometric edge premise.

## Actual finite boosted-ball responses

For each `n_i` choose the unit future timelike normal
`u_i=(5/4,3 n_i/4)`. The boost has spatial block `I+n_i n_i^T/4` and
off-diagonal block `3n_i/4`. Its expansion of `|x0|+|x|` is at most two.
Thus all nine boosted diamonds of radius `R<=L/8` lie in the same
`D_{L/4}`, and the same states (7) apply to all of them.

Let `K_i(R)` denote the FA87 local modular pairing for that boosted ball.
Its stress is `T(u_i,u_i)` in the proper rest coordinates and its common
normalization is `2*pi`. Consequently

```
partial_theta_j Delta<K_i(R)>|0
 =(8*pi^2*R^4/15)*gamma_L*D_ij + remainder_ij,
D_ij=(5/4+(3/4)n_i.n_j)^2.                               (12)
```

Exact rational arithmetic gives
`det D=24407490807/2000000000` and `||D^(-1)||_infinity=68/5`.
The finite radius remainder below proves independence for genuine balls,
not merely a formal limit at their centers.

### Uniform quantitative bound

For `|t|+|x|<=L/4`, shell separation gives
`q=|x-y|^2-t^2>=L^2/2`. In the complex polydisc of radius `L/64` in each
of four coordinates, its perturbation is at most
`(19/128+1/1024)L^2<L^2/4`, and the numerator matrix of (1) has norm
less than `3L`. Hence the kernel norm is at most `24/(pi^2 L^3)`.
Operator-valued Cauchy estimates and
`||f_n||_1<=sqrt(28*pi/3)*L^(3/2)` imply

```
||partial^alpha(F_n/N_L)|| <= H_m L^(-3/2-m),
H_m=16*m!*64^m,      m=|alpha|<=3.                       (13)
```

The simple numerical majorant uses `pi>3` and `sqrt(28)<6`.
Every component of the Wick stress is a symmetric bilinear with one
derivative. Its coordinate Hessian is bounded by
`3H1 H2+H0 H3`, times `L^-6`. Both the timelike normal and each spatial
boost column have absolute component sum below three. Transforming the
two tensor indices, two derivatives, and the spatial Hessian norm gives

```
||Hess_spatial delta_theta T(u_i,u_i)||
 <=195689447424*||theta||_1/L^6.                          (14)
```

FA87's exact radial Taylor bound therefore gives
`|remainder_ij|<=4*pi^2*195689447424*R^6/(35*L^6)`.
After dividing (12) by `(8*pi^2 R^4/15)gamma_L`, each error entry is at
most `C_B(R/L)^2`, where the entirely rational majorant is

```
C_B=(3/14)*195689447424*64/(4249388/150088295).             (15)
```

Here only `pi<4` was used. For every `0<R/L<=2^-28`,
`9*(68/5)*C_B*(R/L)^2<1/4`. The Neumann series proves invertibility of
the actual normalized modular-response matrix, with inverse infinity norm
at most `272/15`. The small certified radius reflects conservative
uniform estimates; no claim of optimality or experimental feasibility is made.

There is also a finite entropy estimator. Use
`[Delta S_i(rho_(epsilon e_j))-Delta S_i(rho_(-epsilon e_j))]/(2epsilon)`
with `epsilon=(R/L)^6`. Taking `r0=1/2` in (11) bounds its error from the
modular derivative by `2epsilon/3`. Its normalized error is at most
`[5/(4249388/150088295)](R/L)^2`. Including this term still gives a
Neumann budget below `1/4` at `R/L=2^-28`. Thus nine finite, positive-state,
vacuum-relative entropy differences also have full rank. Along, for
example, `R_m=L/(2^28(m+1))`, both error channels vanish as `O(R_m^2/L^2)`.

## Scope and replay

The constructed finite variation space supplies nondegenerate matter
responses, a common entropy/stress source, and controlled continuum
small-ball and finite-perturbation limits relevant to G2--G5. The native
observer's ability to implement these continuum shell operations is not
derived. There is no joint lattice/collar refinement estimate, area or
Newton normalization, geometric null-stress balance, generalized-entropy
stationarity, universal coupling to missing matter, or Einstein evolution.
Those fields of `EinsteinRegisterData` are not filled by this calculation.

Run `python -m pytest -q -W error model/weyl_gravity_response`, followed by
`python -m model.weyl_gravity_response --write` and `--check`.
The new receipt recomputes exact shell, frame, positivity-control and error
certificates. It distinguishes its direct checks from the inherited
custody status of FA87; no old receipt or ledger is rewritten.
