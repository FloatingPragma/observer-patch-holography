# An isolated zero pair in the bounded radial lattice source

The operator is the supplied bounded radial texture of
[bounded_defect_index](../bounded_defect_index/PROOF.md), with the same
original two-spin hopping and every lattice site retained:

```
h_(n,n+e_j)=iv sigma_j/(2a),
B_n=eta*m sum_j(-1)^n_j (n_j/|n|) sigma_j for n!=0,
B_0=0,          H=h+B,          eta=+1 or -1.           (1)
```

The theorem holds on `ell^2(Z^3,C^2)` and on each inversion-symmetric
open cube `{-R,...,R}^3`, R>=1. The finite operator is the compression
of (1), with its original links and mass values. Its two origin rows
and all six adjacent sites are present.

The source parameters a,v,m are positive. The sharp sufficient condition
is `m>sqrt(3)*v/a`. The executable interface admits a positive rational
kinetic bound k satisfying

```
k^2>=3(v/a)^2,       m>k,       d=m-k>0,
b^2=3v^2/(2a^2),    r=k/m<1.                           (2)
```

This is a strong mass condition on the full lattice band. It does not
hold under a->0 at fixed positive physical m and v. No such continuum
limit or classification of every nonzero bound level is inferred.

## Removing the origin leaves an invertible complement

The infinite free multiplier is
`-(v/a) sum_j sigma_j sin(k_j)`, so the exact norm of h is
`sqrt(3)*v/a`. A finite compression cannot increase that norm.
Let P be the rank-two origin projector and Q=I-P. In the decomposition
P plus Q,

```
H=[[0,V*],[V,D]],
V=QhP,          D=QHQ=B_Q+h_Q.                         (3)
```

On Q, `B_Q^2=m^2 I`, so `||B_Q psi||=m||psi||`.
The triangle inequality gives
`||D psi||>=(m-k)||psi||=d||psi||`. Since D is self-adjoint, its range
is closed and its orthogonal complement is ker(D)=0; D is onto.
Consequently

```
sigma(D) subset (-infinity,-d] union [d,infinity),
||D^(-1)||<=1/d.                                       (4)
```

The origin couples to six distinct neighbors. Each link has square
Gram matrix `v^2 I_2/(4a^2)`, which gives the exact full-source identity

```
V*V=b^2 I_2,               b^2=3v^2/(2a^2).             (5)
```

These bounds have no dependence on the finite cube radius.

## Exact kernel and the complete central gap

Spatial inversion `S psi(n)=psi(-n)` anticommutes with H and fixes
each origin row. It preserves Q, makes D odd, and sends V into the
odd subspace of Q. Its inverse D^(-1) reverses inversion parity.
Therefore

```
V*D^(-1)V=0.                                          (6)
```

Every zero vector of (3) has the form

```
Phi xi = (xi,-D^(-1)V xi),       xi in C^2.              (7)
```

Equation (6) proves that every such vector is a zero vector; (4) proves
that there are no others. The kernel has exactly dimension two, and
every zero vector is inversion-even. This is an exhaustion theorem on
the stated strong mass domain, in both finite and infinite volume.

For a real E with `0<|E|<d`, the complement resolvent exists. The origin
Schur complement is

```
F(E)=-E I_2-V*(D-E)^(-1)V
    =-E[I_2+V*(D-E)^(-1)D^(-1)V].                     (8)
```

The second equality uses (6) and the resolvent identity. On sigma(D),
the scalar function `1/[lambda(lambda-E)]` is positive for |E|<d.
Hence the bracket in (8) is at least I_2, and F(E) is invertible.
The bounded block inverse proves that E is in the resolvent set of
the full operator. Thus

```
sigma(H) intersect (-d,d)={0},          dim ker(H)=2.   (9)
```

This excludes all additional spectral values in the central interval,
not only additional eigenvalues detected in a finite calculation.
The condition is m>k; a separation assumption m>2k is unnecessary.

## Local zero-mode density and a contact-overlap lower bound

Define the staggered sign `D_s(n)=(-1)^(n1+n2+n3)`, and `J=i sigma_2`.
The source antiunitary `C=D_s J K`, where K is complex conjugation,
anticommutes with H. The hopping changes sign under D_s after spin
time reversal, while the spin mass changes sign under spin time reversal.
The antiunitary `Theta=S C` commutes with H and has square -1.
The corresponding canonical particle-hole map
`c_n -> D_s(n) J c_n^dagger` preserves the quadratic and centered-contact
source. These are exact source identities, with their full parity factors.

Let

```
A=-D^(-1)V,        G=Phi*Phi=I_2+A*A,
W=Phi G^(-1/2),    Z=WW*.                              (10)
```

The columns of W form an orthonormal basis of the exact kernel. The
origin rows of Phi are I_2. By the uniqueness in (7), applying Theta
to Phi gives `Theta Phi xi=Phi J conjugate(xi)`. Inversion-evenness
then gives the pointwise identity

```
D_s(n) J conjugate(Phi_n)=Phi_n J.                     (11)
```

If its first column is u, the second is `-D_s(n) J conjugate(u)`.
These two columns are orthogonal and have the same norm. Each local
Gram matrix `Phi_n* Phi_n` is a scalar multiple of I_2, so G is scalar
as well. After normalization there are nonnegative numbers rho_n with

```
W_n*W_n=W_n W_n*=rho_n I_2,
sum_n rho_n=1,      rho_(-n)=rho_n,
|det W_n|^2=rho_n^2.                                  (12)
```

Since `||A||<=b/d`, the origin weight and inverse participation obey

```
rho_0>=1/(1+b^2/d^2),
I=sum_n rho_n^2 >=rho_0^2>0,       I<=1.               (13)
```

The coefficient of the quartic zero-mode term in the original onsite
contact compression is U I, by the determinant identity in (12).
The bound (13) is independent of the cube volume. This statement does
not bound the interaction with all virtual Fock states; an interacting
ground and charged spectral gap require that separate estimate.

## The origin number-deviation identity

Let R_- and R_+ be the negative and positive one-body spectral projections.
Inversion interchanges them, and the origin is fixed. Since
`R_-+R_+=I-Z`, equation (12) gives

```
R_-(0,0)=R_+(0,0)=(1-rho_0)I_2/2.                    (14)
```

Fill every negative mode and choose an empty, singly occupied, or full
zero pair. These define the corresponding free finite Slater references
and, in infinite volume, their local quasifree CAR expectations. The
single-mode cases can be mixed without changing the expectations below.
For the actual local observable
`D0=(n_0-1)^2=1-n_(0,up)-n_(0,down)+2n_(0,up)n_(0,down)`, Wick contraction
of the two-row covariance in (14) gives

```
<D0>_empty=<D0>_full=(1+rho_0^2)/2,
<D0>_single=(1-rho_0^2)/2,
<D0>_charged-<D0>_single=rho_0^2.                      (15)
```

These are noninteracting reference identities. Projection-angle bounds
are needed to transfer them to actual interacting states.

## Local walk expansion and quantitative volume convergence

The complement inverse has a norm-convergent series. Put

```
T=-B_Q^(-1)h_Q,       A_0=-B_Q^(-1)V,
A=sum_(j>=0) T^j A_0,
||T||<=r=k/m<1,       ||A_0||<=b/m.                    (16)
```

The first N terms have operator-norm remainder at most

```
||A-sum_(j=0)^(N-1) T^j A_0|| <=(b/d)r^N.             (17)
```

The seed A_0 is supported at the six neighbors of the origin; each T
adds at most one nearest-neighbor step. Thus `T^j A_0` is supported
within infinity-norm radius j+1. For L>=1, terms with j<L-1 make no
contribution outside radius L. Since `||G^(-1/2)||<=1`, every normalized
exact zero vector satisfies

```
sum_(|n|_infinity>=L)|psi(n)|^2
 <=min(1,(b^2/d^2)r^(2(L-1))).                        (18)
```

Let A_R, Phi_R and Z_R be the same objects in the symmetric cube of
radius R, extended by zero into the infinite source carrier. Its first
R walk terms are identical to those in (16): none of their paths leaves
the cube. Both remaining tails have the bound in (17), giving

```
||A_R-A|| <=(2b/d)r^R.                                (19)
```

For any unit vector Phi xi in the full zero space, `||xi||<=1`, and
Phi_R xi belongs to the finite zero space. Their distance is bounded
by the right side of (19). The same argument works in reverse.
For two orthogonal projections, the norm of their difference is the
maximum of the two subspace-distance norms. This follows by squaring
their difference and using the two diagonal blocks relative to either
projection. Consequently

```
||Z_R-Z|| <=min(1,(2b/d)r^R).                          (20)
```

Equations (17)-(20) are quantitative operator estimates. They do not
assume convergence of an eigensolver or drop paths that return to the
cube after leaving it; those paths are retained in the bounded tails.

## Exact retained examples and physical domain

At a=v=1 and k=2, m=4 gives d=2, `rho_0>=8/11` and `I>=64/121`.
The m=64 example gives d=62, `rho_0>=7688/7691`,
`I>=(7688/7691)^2` and r=1/32. Both examples use the same supplied
bounded radial family and original source clock.

The theorem proves an exact isolated quadratic zero pair and its
finite-to-infinite-volume attachment in the strong mass regime. It
does not derive the radial field from the source interaction, establish
a volume-uniform interacting Fock gap, calibrate physical mass units,
or identify these modes with observed particles. Native repair and
Einstein dynamics are outside these spectral and reference-state bounds.
