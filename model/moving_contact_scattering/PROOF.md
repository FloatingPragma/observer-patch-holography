# Uniform moving contact scattering

Use the same two-body source, critical contact and global clock as the
band-edge and moving-bound-pair constructions:

```
c_a = sqrt(3)/a^2,
1/U_a = a^2 W0/2 - a^3 gamma/pi,   gamma>0,
J_a = 6/a^2 - H_a,
W0 = <1/sum_j cos(q_j)^2>,        2/3 <= W0 < 5.
```

The full torus average defines W0 without measured inputs. The original
action H_a has the unchanged Weyl hopping and repulsive onsite contact.
The scaling, gamma and reflected clock are specified family choices.
The scattering and bound-band statements use the two-particle sector. The
final sections address its additive extension to all particle numbers.

Let `K=aP`, `kappa=max_j |K_j|`, and let s denote energy below the exact
free upper edge `B_K=2c_a R`, where

```
C_j=cos(K_j/2), S_j=sin(K_j/2), R^2=sum_j C_j^2,
alpha=R/sqrt(3), h=a^2 s, t=R-h/(2sqrt(3)),
lambda=R^2-t^2=alpha h-h^2/12.
```

The uniform domain is `kappa<=1/4`, `0<=h<=1/128`. Thus
`1-kappa^2/4<=alpha<=1`, `t>=3/2`, and `0<=lambda<=h`.
For example `t>=(5/3)(1-1/64-1/768)>3/2` uses `sqrt(3)>5/3`.

## Exact source Green function and all edge charts

The contact vector is the constant normalized spin singlet. On each
relative-momentum fiber the exact four-spin resolvent gives

```
g_K(E) = < E(E^2-4c_a^2 e^2) /
 [E^4-4c_a^2(d^2+e^2)E^2+16c_a^4(d.e)^2] >,
d_j=C_j sin(q_j), e_j=S_j cos(q_j).
```

At real `E=2c_a t` near the upper edge this becomes

```
g_K(E+i0) = t/(2c_a) W_t(-lambda+i0),
W_t(z) = pi^(-3) integral_[-1,1]^3 J(u) N(u)/(F(u)+z) du,
J(u) = product_j (1-u_j^2)^(-1/2),
N(u) = 1-sum_j S_j^2 u_j^2/t^2,
F(u) = sum_j cos(K_j) u_j^2
       + [sum_j C_j S_j u_j sqrt(1-u_j^2)]^2/t^2.
```

This formula includes all eight edge extrema. Partition each q coordinate
into its two half-circles centered at `+pi/2` and `-pi/2`. In the chart
with sign sigma_j write `q_j=sigma_j*pi/2+r_j`, `|r_j|<=pi/2`.
Then `sin(q_j)=sigma_j cos(r_j)` and
`cos(q_j)=-sigma_j sin(r_j)`. Every component product d_j e_j is
`-C_j S_j cos(r_j)sin(r_j)`, independent of sigma_j. The complete Green
integrand is therefore identical on all eight charts. Their total
normalized measure is `pi^(-3) d^3r`. Setting `u_j=sin(r_j)` gives J above.
There is no selected valley or missing factor of eight.

The `+i0` orientation is also fixed by the source resolvent. For a small
positive imaginary part of t, the imaginary coefficient in its denominator
is `2t[1-(d.e)^2/t^4]`. It is positive: `d^2<=3`,
`e^2<=3kappa^2/4`, `t^4>=81/16` give
`(d.e)^2/t^4 <= (4/9)kappa^2 < 1`. Thus the boundary distribution is
the principal value minus `i*pi` times the positive coarea density.
Dependence of F on t does not reverse or remove this sign.

## Uniform radial change of variable

Put `Q=|u|^2`. On the entire cube,

```
(1-kappa^2/2)Q <= F(u) <= Q,       0<N(u)<=1.
```

Indeed `F=Q-2e^2+(d.e)^2/t^2`, `e^2<=kappa^2 Q/4`, and
`(d.e)^2/t^2<=4e^2/3`. The same inequalities give

```
|W_t(0)-W0|
 <= kappa^2 [W0/2+1/9]/(1-kappa^2/2) < 3kappa^2.       (1)
```

For the singular part write `u=r omega`, `|omega|=1`, `x=r^2` and
`F=r^2 f(x,omega)`. Define

```
w(x)=sum_j C_j S_j omega_j sqrt(1-x omega_j^2),
f(x)=sum_j cos(K_j) omega_j^2+w(x)^2/t^2.
```

On `0<=x<=1/4`, Cauchy--Schwarz and
`sum |omega_j|^3<=1`, `sum |omega_j|^5<=1` give
`|w|,|w'|,|w''|<=kappa`. Hence

```
31/32 <= f <= 1,  |f'|<=kappa^2,  |f''|<=2kappa^2.
```

Set `y=r sqrt(f)` and `z=y^2=x f(x)`. With `D=f+x f'`,

```
D>=61/64>7/8,      |D'|<=5kappa^2/2<=5/32<1/4.
```

Thus x maps monotonically to z. The common transformed ball `0<=y<=b=1/4`
has a unique inverse in `r<1/2` for every direction. The global lower bound
on F ensures that this ball contains every singular point of
`F-lambda` and that its complement has `F>=b^2`. The transformed radius b
is independent of K and s; its preimage in the u coordinates depends on them.

The radial measure becomes

```
J(u) N(u) r^2 dr = y^2 A_direction(y^2) dy,
A_direction(z) = J(sqrt(x)omega) N(sqrt(x)omega) / [sqrt(f(x)) D(x)].
```

Let A(z) be its normalized spherical average. Direct bounds on this
explicit expression suffice; no unproved smoothness constant is used.
On `x<=1/4`, `J<=2`, `|J'|<=2`, `N<=1`, `|N'|<=1`,
`f^(-1/2)<=2`, `|(f^(-1/2))'|<=1`, `D^(-1)<=2`, and
`|(D^(-1))'|<=1`. The product derivative is at most 24, and
`dx/dz=1/D<=2`. Therefore

```
|A'(z)| <= 48,   A(0)<=2,   |A(0)-1|<=2kappa^2.        (2)
```

For the last inequality, `A_direction(0)=f(0)^(-3/2)` and the derivative
of `v^(-3/2)` has magnitude below two for `31/32<=v<=1`.
All derivatives are uniform in K,t,omega. Smoothness follows from positive
square roots and the uniform positive inverse derivative D. Angular
integration preserves these derivative bounds.

## The interior principal value

The common-ball contribution is exactly

```
W_ball(-lambda+i0) = (4/pi^2) integral_0^b
                      y^2 A(y^2)/(y^2-lambda+i0) dy.
```

For `0<lambda<=1/128=b^2/8`, subtract its value at lambda=0.
The real part is

```
(4lambda/pi^2) PV integral_0^b A(y^2)/(y^2-lambda) dy.
```

The constant part uses the exact identity

```
PV integral_0^b 1/(y^2-lambda) dy
 = -atanh(sqrt(lambda)/b)/sqrt(lambda).
```

Since `atanh(v)/v<=1/(1-v^2)` and `A(0)<=2`, its magnitude after
multiplication by lambda is at most `64lambda/7`. For
`B(z)=A(z)-A(0)`, subtract B(lambda) at the singularity. The quotient
`[B(y^2)-B(lambda)]/(y^2-lambda)` is continuous and bounded by 48.
The additional constant term has magnitude at most `48lambda/[b(1-lambda/b^2)]`.
Together these bound the B integral by `96/7`. This subtraction, with
the uniform derivative bound (2), justifies the principal value directly.
Multiplying by `4lambda/pi^2` gives the real common-ball bound
`(640/63)lambda`, using `pi>3`.

Outside the common ball, both denominators are separated from zero:
`F>=b^2` and `F-lambda>=7F/8`. The original cube measure has mass one,
and `N<=1`. Consequently the real tail difference is at most

```
integral N lambda/[F(F-lambda)] dmu <= 512lambda.        (3)
```

The imaginary part of the common ball is exactly
`-(2/pi)sqrt(lambda)A(lambda)`. There is no imaginary tail.
Equation (2), `sqrt(lambda)<1/8`, and `pi>3` bound its difference from
`-(2/pi)sqrt(lambda)` by `kappa^2/6+4lambda`. Combining the estimates gives

```
|W_t(-lambda+i0)-W_t(0)+(2i/pi)sqrt(lambda)|
 <= 528lambda+kappa^2.                                 (4)
```

At lambda=0 the boundary value is the convergent threshold integral,
and the same bound follows directly. The preceding estimates also prove
continuous extension to this endpoint, uniformly over the declared K,t
domain. Neither (3) nor (4) is inferred from an exterior resolvent bound.

## Uniform source expansion and nontrivial scattering

Write `beta=t/sqrt(3)=alpha-h/6`. Then `0<beta<=1`,
`1-beta<=kappa^2/4+h/6`, and

```
0 <= 1-beta sqrt(alpha-h/12) <= kappa^2/2+h/4.
```

Combining (1) and (4), multiplying by `a^2 beta/2`, and using
`sqrt(h)<1/8` bounds the source Green function by

```
|(g_K(B_K-s+i0)-a^2 W0/2)/a^3+i sqrt(s)/pi|
 <= a (3 Pmax^2+265 S),                                (5)
```

uniformly for `|P_j|<=Pmax`, `0<=s<=S`, `aPmax<=1/4`, `a^2S<=1/128`.
The unrounded momentum and energy coefficients are respectively
`21/8+1/48=127/48<3` and `3173/12+1/96=25385/96<265`.
They include the K-dependent edge and the factor beta; neither is dropped.

Let `epsilon_G=a(3Pmax^2+265S)`. For `a gamma<=1` and
`44epsilon_G<=7gamma`, the same positive contact coupling obeys

```
D_a(P,s)/a^3 = (-gamma+i sqrt(s))/pi + error,
|error|<=epsilon_G,  |D_a/a^3|>=7gamma/44.
```

The existing rank-one stationary formula is `S_H=conj(D_a)/D_a` on the
contact cyclic channel. The ratio inequality therefore proves

```
sup_(P,s) |S_H-(gamma+i sqrt(s))/(gamma-i sqrt(s))|
 <= (88/7) epsilon_G/gamma.                            (6)
```

The comparison `J_a=6/a^2-H_a` reverses time, so `S_J=conj(S_H)` and its
limit is `(gamma-i sqrt(s))/(gamma+i sqrt(s))`. At s=0 both limits are
one, defined by the continuous threshold extension. Gamma is strictly
positive; no zero-gamma threshold assertion is included.

For s>0 the coarea density is strictly positive, since J,N,f,D are
positive. The contact cyclic spectral measure can consequently be flattened
to ds, and total momentum K=aP to dP, by their positive square-root
densities. Equation (6) is operator-norm convergence of the resulting
multiplication operators on `L2([-Pmax,Pmax]^3 x (0,S), dP ds)`. It controls
moving contact-channel wave packets, including arbitrary untouched
references. Noncontact reducing channels have identity scattering and are
retained. This comparison does not supply a full continuum momentum-space
wave-operator embedding or a many-particle asymptotic-completeness theorem.

## Bound states and scattering use one generator

On the same family, the moving-bound-pair theorem gives the bound-band
limit `|P|^2/4-gamma^2`. Its premise is
`256a(gamma^2+Pmax^2)<=gamma`. No momentum-dependent scalar is subtracted
from H. The free scattering energy under the global J_a is

```
6/a^2-B_a(aP)+s -> |P|^2/4+s,
|6/a^2-B_a(aP)-|P|^2/4| <= a^2 Pmax^4/4.
```

The last estimate uses the same [moving-edge rationalization](../moving_bound_pair/README.md):
`6/a^2-B_a(aP)=2 sum_j sin(aP_j/2)^2/[a^2(1+alpha)]`.
The sine remainder gives the error
`a^2[sum_j P_j^4+(sum_j P_j^2)^2]/48`, which is at most
`a^2 Pmax^4/4`. Its domain `aPmax<=1/2` includes the domain used here.

The error bounds (5)--(6) and the existing bound propagator therefore
compose on the same source, coupling and clock. For `gamma=Pmax=S=1`
and `a=1/16384`, the uniform scattering error is below 1/4, and the
bound-band unit-time propagator error is below 0.028. At `s=gamma^2`
the limiting S_J is `-i`; the source scattering is separated from identity
by more than one throughout the momentum window. This proves a surviving
moving two-body interaction and a bound branch in three spatial dimensions.
Its dispersion is nonrelativistic. It does not construct a relativistic
many-body quantum field theory, fix gamma, or identify physical particles.

## The additive Fock generator has no empty ground vacuum

Let h_a denote the original one-particle Weyl hopping on
`ell2(Z^3) tensor C^2`. Its Fourier eigenvalues are
`+-c_a sqrt(sum_j sin(k_j)^2)`, so

```
0 <= j_a := 3/a^2-h_a <= 6/a^2.
```

The additive extension of the same reflected clock considered on fermionic
Fock space is

```
J_(a,N) = dGamma(j_a) - U_a D,
D = sum_x n_(x,up) n_(x,down),
-U_a floor(N/2) <= J_(a,N) <= 6N/a^2.                 (7)
```

On an occupation configuration every doublon consumes two particles, hence
`0<=D<=floor(N/2)`. This proves (7), including the odd-N and empty sectors.
Each fixed-N operator is bounded and self-adjoint. Their direct sum is
self-adjoint on the domain
`{psi: sum_N ||J_(a,N) psi_N||^2<infinity}`. Finite particle-number vectors
with finite spatial support form a core: truncate the particle-number sum
in the graph norm and then approximate each bounded sector by finite spatial
support. Thus this construction supplies a mathematically defined all-N
generator. Its lower bound is linear in N, uniformly over finite spatial
volumes at each fixed a. The constant `U_a/2` depends on a and diverges as
`a->0`; (7) supplies no uniform continuum many-body stability estimate.

For `a gamma<=1/4096`, the full Fock operator is unbounded below. The proof
uses actual normalizable states of this same source. Set `Pmax=gamma/4`.
The moving-pole premise follows from `272a gamma<=1`; its pole offset is
at least `gamma^2/4`. The exact reflected bound energy on this momentum box
is therefore at most

```
3 gamma^2/64 + a^2 gamma^4/1024 - gamma^2/4
 = -13 gamma^2/64 + a^2 gamma^4/1024 < -gamma^2/8.      (8)
```

Take any normalized square-integrable total-momentum wavepacket supported
inside this box and apply the normalized antisymmetric bound-fiber isometry
from [the moving-pair construction](../moving_bound_pair/README.md).
The result is a normalized two-particle vector psi satisfying (8), rather
than a single nonnormalizable momentum eigenstate. The contact is positive
and obeys `U_a<=4/a^2`, since

```
1/U_a >= a^2(1-a gamma)/3 >= a^2/4.
```

Equation (7) then gives `||J_(a,2)||<=12/a^2`. Finite spatial supports are
dense in the antisymmetric two-particle Hilbert space. Choose a normalized
finite-support approximation phi with
`||phi-psi||<=a^2 gamma^2/768`. The expectation error is bounded by
`2||J_(a,2)|| ||phi-psi||<=gamma^2/32`. In particular,

```
<phi,J_(a,2) phi> < -gamma^2/16.                       (9)
```

The norm tolerance is explicit; this density argument does not calculate
a support radius. It suffices for the following variational obstruction.
Enclose both spin modes of phi's support in a finite set of sites Lambda.
For any integer m, translate Lambda to m disjoint site sets and form the
fermionic wedge product of the m translated phi vectors. Orthogonality of
their one-particle mode spaces makes this vector normalized. Each core has
exactly two particles. A hopping term between distinct cores, or between
a core and the empty exterior, changes those local particle numbers and
has zero expectation. Internal hopping expectations equal those of phi;
the onsite contact is additive on disjoint site sets. Therefore the full
2m-particle expectation is exactly m times (9):

```
inf spectrum(J_a on full Fock space) = -infinity,
<Phi_m,J_(a,2m) Phi_m> < -m gamma^2/16.                (10)
```

This is compatible with the fixed-N stability bound (7). It excludes the
empty Fock vector, whose energy is zero, as a ground vacuum of this
additive reflected generator. It is not an instability at fixed particle
number, a proof of continuum collapse, or an exclusion of different
finite-density representations.

For a modified generator `J_a+mu_a N`, an empty ground requires
`2mu_a>=delta_a(0)`, where delta_a(0) is the positive bound-pole offset at
zero total momentum. Indeed the bound-band energy there is exactly
`-delta_a(0)+2mu_a`; if it is negative, continuity gives a nonzero momentum
neighborhood supporting a negative normalizable wavepacket. Continuity
follows from the strict exterior-resolvent derivative and the isolated
simple bound pole. Since `delta_a(0)->gamma^2`, any finite limit of mu_a
must satisfy

```
mu >= gamma^2/2.                                      (11)
```

Equation (11) is necessary only. It supplies no lower bound for higher
particle sectors. A distinct sufficient finite-spacing bound is
`mu_a>=U_a/2`, directly from (7); this diverges as `a->0` and does not give
a finite-energy continuum vacuum on the identified pair branch.

The exact finite control retains all 256 Fock states on a four-site Weyl
ring, with two disjoint two-site cores and both inter-core hopping links.
The occupied core mode addresses interleave in the global ordering.
Each normalized fixed-two-particle core has energy `-4/3`, and their graded
fermion product has energy `-8/3`. Omitting the fermionic reordering signs
gives `-2`. Independent Jordan--Wigner matrices check the Hamiltonian,
every particle-number sector, and the exactly vanishing inter-core hopping
expectation. These finite parameters test the product identity; the
critical-family negative packet is certified separately by (8)--(9).

## The pair branch has no exact affine Lorentz-shell identification

The limiting pair branch, including a finite number shift, is
`epsilon(P)=|P|^2/4-gamma^2+2mu`. Consider any fixed affine clock map
`E=alpha epsilon+beta`, `alpha!=0`, and fixed invertible affine momentum
map `p=L P+p0`. An exact massive or massless Lorentz shell would require

```
[alpha(|P|^2/4-gamma^2+2mu)+beta]^2
  -v^2 |L P+p0|^2 -m^2 v^4 = 0                      (12)
```

on an open momentum set, with `v>0` and `m>=0`. This is a polynomial in P.
A polynomial vanishing on a nonempty open set vanishes identically. On
`P=r e1`, however, its coefficient of `r^4` is `alpha^2/16>0`;
the momentum and mass terms have degree at most two. This contradicts
(12). Thus this particular branch under the stated affine identification
is neither an exact massive nor an exact massless Lorentz shell. The
argument permits its ordinary nonrelativistic low-momentum interpretation
and makes no claim about arbitrary nonlinear dictionaries, other source
branches, or different continuum clock constructions.
