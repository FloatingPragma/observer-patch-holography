# Exact zeros of the bounded source hedgehog

The operator is the bounded original-lattice profile of
[topological_defect_spectrum](../topological_defect_spectrum/PROOF.md),
without a change of mass profile or hopping convention. For a,v,m>0 and
orientation eta=+1 or -1, on `ell^2(Z^3,C^2)` set

```
(h0 psi)(n)=iv/(2a) sum_j sigma_j[psi(n+e_j)-psi(n-e_j)],
B_n=eta*m sum_j (-1)^n_j (n_j/|n|) sigma_j  for n!=0,
B_0=0,                       H=h0+B.                 (1)
```

This is the complete original two-spin carrier. The corresponding folded
description retains all eight valleys. The quadratic texture is supplied;
its selection by the interacting source and its autonomous formation
energy are not premises proved here.

## The inversion-graded Fredholm index

Let `R_inv psi(n)=psi(-n)` and `P_+/-=(I+/-R_inv)/2`. This is a
nonlocal spatial-inversion grading. The original central difference and
the stated mass both change sign under R_inv. Thus, for every t in [0,1],

```
H_t=t*h0+B,       R_inv H_t R_inv=-H_t,
A_t=P_- H_t P_+ : range(P_+) -> range(P_-).             (2)
```

Each H_t is bounded and self-adjoint. The path is norm-continuous at
every fixed a; there is no domain change at t=0.

The essential gap remains open on this entire path. To see this directly,
write u(n)=n/|n| for n!=0, u(0)=0, and epsilon_j(n)=(-1)^n_j.
The forward-edge coefficient of `{h0,B}` is `i*eta*m*v/(2a)` times

```
T_j(n)=-epsilon_j(n) Delta_j u_j I
       +sum_(k!=j) epsilon_k(n) Delta_j u_k sigma_j sigma_k,
Delta_j u=u(n+e_j)-u(n).                               (3)
```

The exact two-spin identity is `T_j* T_j=|Delta_j u|^2 I`.
Every coefficient therefore tends to zero at spatial infinity. Finite
truncation of the coefficients approximates this finite-hop operator in
norm by finite-rank operators. Also `B^2-m^2 I` has support only at the
origin. Hence

```
H_t^2=t^2*h0^2+m^2 I+compact,
sigma_ess(H_t^2)=[m^2,m^2+3*t^2*v^2/a^2].              (4)
```

Here `h0^2` has scalar Fourier multiplier
`(v/a)^2 sum_j sin(k_j)^2`. Essential spectral mapping and the exact
inversion anticommutation give

```
sigma_ess(H_t)=[-sqrt(m^2+3*t^2*v^2/a^2),-m]
               union [m,sqrt(m^2+3*t^2*v^2/a^2)].       (5)
```

In particular H_t and its off-diagonal block A_t are Fredholm. The
ordinary complex Fredholm index of A_t is constant along (2).

At t=0 this index is explicit. Partition every nonzero site into the
inversion pair `{n,-n}`. On its four rows, H_0 is `diag(B_n,-B_n)`.
The even/odd pair basis turns this into

```
[[0,B_n],[B_n,0]],       B_n^2=m^2 I_2.                 (6)
```

The direct sum of all nonzero-pair off-diagonal blocks is invertible
with inverse norm 1/m. The origin contributes two even rows, no odd rows,
and zero mass. Thus A_0 is onto, has a two-dimensional kernel and no
cokernel. Consequently

```
ind(A_1)=dim ker(A_1)-dim ker(A_1*)=2.                  (7)
```

There are at least two exact square-summable, inversion-even zero modes
of H at every positive a,v,m. Its kernel is finite-dimensional, because
zero is outside the essential spectrum. Equation (7) does not exclude
additional even and odd zero modes.

Both eta=+1 and eta=-1 give index +2 in (7). Spatial inversion is a
different grading from the continuum internal Clifford grading, whose
hedgehog and antihedgehog indices in the parent theorem are +2 and -2.
The reflection index is not identified with that internal index, the
hedgehog degree, an electromagnetic charge or Lorentz spin.

## The allowed perturbation class and its countercontrol

If V is compact, bounded and self-adjoint and
`R_inv V R_inv=-V`, then `H+sV` remains Fredholm and off diagonal for
every s in [0,1]. Its essential spectrum is unchanged. The same index
argument therefore preserves the lower bound of two exact zero modes,
without a smallness restriction on the compact perturbation.

The index also persists for a bounded self-adjoint, inversion-odd V
with `||V||<m`, whether or not V is compact. In the Calkin algebra the
inverse of H has norm 1/m by (5). The Neumann estimate for `H+sV`
keeps zero outside its essential spectrum throughout the path, and
the off-diagonal index is again constant.

The anticommutation condition is essential. At t=0, the compact scalar
term `(m/3) P_origin` changes the two origin energies from zero to m/3,
leaves every other site mass and the essential spectrum unchanged, and
has no zero eigenvector. This perturbation is even under inversion and
breaks (2). For any t, a uniform onsite scalar term lambda I moves each
original zero eigenvector to energy lambda. These are fixed-clock
Hamiltonian perturbations; an open essential gap alone does not pin an
energy level at zero. No claim is made that the actual native repair
channel realizes or preserves the admitted operator perturbation class.

## Explicit exponential localization

The following estimate applies to every normalized exact zero mode of
(1). It also applies after an arbitrary self-adjoint perturbation with
both matrix indices inside a finite cube K_S, once that cube is included
in the reference radius. Index preservation for that perturbation still
requires the anticommutation condition above. No exponential estimate is
transferred to arbitrary nonlocal compact perturbations.

Define

```
R=max(1,ceil(12v/(am)),S),     K_R={n:|n|_infinity<=R},
D=Z^3\K_R,                   q=1+am/(12v).             (8)
```

For the original operator S=0. Let h_D be the compression of h0 to D,
and B_D the unchanged diagonal mass there. The exterior
operator squares as

```
H_D^2=h_D^2+m^2 I+{h_D,B_D},        h_D^2>=0.           (9)
```

This is the square of the compressed operator, not the compression of
the full square. It introduces no omitted boundary cross-term. For every
retained edge in D, both endpoint lengths are at least R+1. The elementary
normalization bound `|x/|x|-y/|y||<=2|x-y|/min(|x|,|y|)` and (3) bound
each anticommutator edge by `mv/[a(R+1)]`. Summing the six edges, using
the block Schur bound, gives

```
||{h_D,B_D}||<=6mv/[a(R+1)],
H_D^2>=m^2-6mv/[a(R+1)]>=m^2/2,
||H_D^(-1)||<=2/m.                                     (10)
```

Put `d(n)=(|n|_infinity-R)_+` and
`chi(n)=min(1,d(n)/R)`. The scalar cutoff is zero on K_R, equals one
beyond K_(2R), and changes by at most 1/R on an edge. If H psi=0,
including a permitted perturbation supported in K_R, then

```
H_D(chi psi)=P_D[H,chi]psi=f,
||f||<=3v/(aR)<=m/4,             supp(f) subset {1<=d<=R}. (11)
```

The mass commutes with chi. The finite-core perturbation contributes
nothing, because chi is zero at both of its matrix indices. Formula
(11) retains the hopping across the inner cutoff boundary.

For an integer M>=1 let `w_M(n)=q^min(d(n),M)`. Adjacent weights have
ratio between 1/q and q. Only the hopping changes under conjugation;
its six block norms are v/(2a). Hence

```
||w_M H_D w_M^(-1)-H_D||<=3v(q-1)/a=m/4.
```

The inverse bound (10) and a Neumann series imply
`||(w_M H_D w_M^(-1))^(-1)||<=4/m`, uniformly in M.
Using the support and norm of f in (11) gives

```
||w_M chi psi||<= (4/m) q^R [3v/(aR)] <= q^R.
```

Monotone convergence as M tends to infinity proves the same bound with
weight q^d. On `|n|_infinity>=2R+L`, chi=1 and d>=R+L. Thus, for every
integer L>=0,

```
sum_(|n|_infinity>=2R+L) |psi(n)|^2 <= q^(-2L).         (12)
```

All constants are explicit at every positive spacing. On the stated
small-spacing domain `am<=12v`, the integral inequality
`log(1+x)>=x/2` for 0<=x<=1 yields

```
q^(-2L)<=exp[-m(aL)/(12v)].                            (13)
```

For the unperturbed profile, `2aR<=24v/m+2a`. Therefore the physical
reference radius and probability decay length stay bounded as a tends
to zero with v and m fixed. The estimate is an upper bound on localization
length; it is not an exact wavefunction or mode-count classification.

## Scope of the source conclusion

The result gives exact localized zero modes of the bounded original
lattice profile at each spacing, with explicit stability under the
stated spectral anticommutation. It uses the original two spin rows,
six hopping links and all folded valleys. Finite controls test the source
sign convention, full parity spaces, hopping homotopy, compact odd
perturbations, scalar countercontrols and exterior localization bounds.
Those finite controls are independent checks of the proof and are not
used as an infinite-volume eigenvalue extrapolation.

The texture, its site-centered inversion symmetry, its amplitude and
its clock remain specified inputs. Native repair protection, autonomous
texture formation, full interacting continuum particle states, physical
charges and a common Einstein stress tensor are not established here.
