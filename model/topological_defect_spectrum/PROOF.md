# Spectra attached to the native topological mass family

The matrices and onsite mass channels in this theorem are those of
[native_mass_channels](../native_mass_channels/README.md). All eight source
Weyl valleys, with two spin rows each, are retained. The added premise is
that the stated quadratic mass field is present. Its amplitude, radial
profile, dynamical formation and energy are not derived from the contact
interaction by this theorem. The conditional spectrum is determined once
that field and the source clock are specified.

## Bulk spectrum and its field content

Let gamma_j be the three native folded kinetic matrices, beta_j the three
spin-mass matrices, and beta_0 the scalar mass. They act on C^16. The
six matrices gamma_j,beta_j, j=1,2,3, anticommute pairwise and square to I.
Beta_0 anticommutes with every gamma and commutes with every spin beta.
For the continuum velocity v>0 and a uniform mass field,

```
H_bulk(p)=v sum_j p_j gamma_j+m0 beta_0+sum_j m_j beta_j,
r=|m|,
E(p)=+/-sqrt(v^2 |p|^2+(m0+r)^2),
     +/-sqrt(v^2 |p|^2+(m0-r)^2).                       (1)
```

Each generic signed level has multiplicity four. Thus the two gap
magnitudes |m0+r| and |m0-r| each occur in two four-component Dirac copies.
Equal levels have their multiplicities added. A pure spin mass has four
Dirac copies of mass magnitude r. With both masses zero there are eight
Weyl cones, four of each chirality. At m0=r>0, two Dirac copies are
massless and two have mass magnitude 2r. These are quadratic field species
and source energy scales; they are not a Standard Model assignment or
interaction-renormalized pole masses.

The tuned massive-and-massless case closes the uniform bulk gap. The
hedgehog Fredholm argument below instead uses m0=0 and an open essential
gap. These rows do not construct one regime with both a massless bulk
and the stated gapped-defect protection.

The source lattice counterpart replaces v p_j by (v/a)sin(a p_j).
Its full mass classification and uniform band multiplicities are exact in
the parent package. No node is removed in obtaining (1).

## Clifford algebra, degree and the comparison grading

Define

```
rho_j=i gamma_j beta_j,
b_j=(beta_j+i gamma_j)/2,
F=sum_j b_j* b_j,
Gamma=rho_1 rho_2 rho_3=i gamma_1 gamma_2 gamma_3 beta_1 beta_2 beta_3.
```

The b_j satisfy the three-mode canonical anticommutation relations;
`b_j* b_j=(I-rho_j)/2`. Their Fock representation has a two-dimensional
multiplicity space: the joint vacuum projector
`P0=product_j(I+rho_j)/2` has rank two. The form-number subspaces have
ranks 2,6,6,2. Gamma equals (-1)^F. It anticommutes with the massless
kinetic operator and each spin mass, so it grades every operator below.
This grading is a proved comparison symmetry, not an electromagnetic charge.

In the original valley basis, rho_j is exactly Pauli Y on valley bit j,
tensored with the identity on the other bits and original two spin rows.
The two vacuum vectors have valley amplitudes i^|b|/sqrt(8), with the two
independent spin columns. The exact algebra checks use the original folded
matrices rather than a replacement irreducible Clifford representation.

For a nonzero spin mass on a surrounding sphere, its direction n defines
an integer degree. A radial hedgehog n(x)=x/|x| has degree +1; its negative
has degree -1. These are the degree classes of the parent's central gapped
mass component. The explicit zero modes below connect those two classes
to this actual sixteen-row Dirac operator.

## A complete bounded hedgehog spectrum

Set m0=0, m>0, and consider the prescribed operator on R^3

```
H=-iv gamma.dot(grad)+m beta.dot(x/|x|).                  (2)
```

The value of the mass at the single point x=0 has no effect on the
multiplication operator. Its coefficient is bounded and Hermitian.
Consequently H is self-adjoint on H^1(R^3,C^16), as a bounded perturbation
of the free Dirac operator. There is no adjustable point boundary condition.
The mass direction has a singular core; smoothing that core changes the
excited spectrum and is treated separately below.

Put W(x)=m|x|, A_j=v partial_j+partial_j W and
`d_W=sum_j b_j* A_j`. Since the gradient is exact, d_W^2=0, and

```
H=d_W+d_W*,
H^2=sum_j A_j* A_j+2v sum_ij b_i* (partial_i partial_j W) b_j.
```

This is first a quadratic-form identity on smooth functions away from
zero and then on the full form domain by closure. In three dimensions
`Hess W=(m/r)(I-n n^T)` and `Delta W=2m/r`, without a delta at the origin:
the boundary flux through a radius-epsilon sphere is O(epsilon^2).
Hardy's inequality controls the 1/r terms. With
`b_n=sum_j n_j b_j`, `N_n=b_n* b_n`, the result is

```
H^2=-v^2 Delta+m^2+(2vm/r)(F-1-N_n).                    (3)
```

This scalar radial coefficient and matrix angular term retain every form
sector. In degree zero, (3) is the scalar Coulomb operator

```
L0=-v^2 Delta+m^2-2vm/r                                 (4)
```

with multiplicity two. Its complete sub-threshold eigenvalues are

```
lambda_n=m^2(1-1/n^2),  n=1,2,...,
orbital multiplicity n^2.
```

For completeness, separate into spherical harmonics of degree l. If
`lambda=m^2-v^2 k^2<m^2`, the regular radial solution has the form
`u(r)=r^(l+1) exp(-kr) L(2kr)`. Its equation is

```
z L''+(2l+2-z)L'+[m/(vk)-l-1]L=0.
```

Regularity at zero in the source form domain and square integrability at
infinity require polynomial termination. Hence
`m/(vk)-l-1` is a nonnegative integer and `k=m/(vn)`, with
`0<=l<=n-1`. Summing 2l+1 gives n^2. This also excludes an omitted
sub-threshold scalar radial solution.

On degree two, F-1-N_n=1-N_n is positive semidefinite; on degree three
it is I. Both sectors therefore have H^2>=m^2. For a degree-one
eigenvector with 0<lambda<m^2, d_W maps it into a degree-two eigenvector
at the same lambda, which must vanish. Thus
`psi=d_W d_W* psi/lambda`: every such vector comes from a degree-zero
eigenvector. Conversely d_W maps each positive degree-zero eigenvector
injectively into degree one, with norm factor sqrt(lambda). These statements
follow for the closed differential and its spectral subspaces, so they
are not based on a truncated form calculation.

At lambda=0 the positive identity preceding (3) first forces every
A_j psi=0, hence `psi=exp(-mr/v) chi` for a constant sixteen-row vector
chi. The Hessian term then requires the transverse b operators to annihilate
chi for almost every radial direction. Intersecting these transverse
conditions gives `b_j chi=0` for all j. Therefore the zero space has
exactly dimension two, with normalized vectors

```
psi_0(x)=(m^3/(pi v^3))^(1/2) exp(-mr/v) chi,
chi in range(P0),  ||chi||=1.                           (5)
```

If phi_nlm is a normalized scalar eigenfunction of (4), tensored with
either vacuum column, the two signed eigenvectors for n>=2 are

```
[phi_nlm +/- H phi_nlm/sqrt(lambda_n)]/sqrt(2).
```

Thus the complete spectrum inside the bulk gap is

| Level | Energy | Multiplicity for each displayed energy |
| --- | --- | --- |
| Topological core | 0 | 2 |
| n>=2 | +m sqrt(1-1/n^2) | 2n^2 |
| n>=2 | -m sqrt(1-1/n^2) | 2n^2 |

The matrix Coulomb perturbation in (3) is relatively compact with respect
to -Delta: local compactness uses the Hardy bound and Rellich compactness;
the coefficient tends uniformly to zero at infinity. Hence the essential
spectrum of H^2 is [m^2,infinity). The grading Gamma pairs positive and
negative energies, giving essential spectrum
`(-infinity,-m] union [m,infinity)` for H. This identifies the full spectral
set and every discrete level in its gap. It does not assert a multiplicity
formula for continuum scattering channels in the inhomogeneous background.

The first positive levels in units of m are sqrt(3)/2, 2sqrt(2)/3,
sqrt(15)/4, 2sqrt(6)/5, and sqrt(35)/6, with respective multiplicities
8,18,32,50,72. They accumulate at the bulk threshold m. These values depend
on the constant-magnitude radial profile in (2); topology alone does not
fix them.

## Smooth cores and the protected zero-mode count

For the smooth profile

```
W_l=m(sqrt(r^2+l^2)-l),
m(x)=grad W_l=m x/sqrt(r^2+l^2),  l>0,                  (6)
```

the Hessian is strictly positive at every finite x. Its radial eigenvalue
is `m l^2/(r^2+l^2)^(3/2)` and its two transverse eigenvalues are
`m/sqrt(r^2+l^2)`. The same positive factorization proves that the kernel
is exactly `exp(-W_l/v) range(P0)`, again of dimension two. The essential
threshold is m. The Coulomb excited-level formula is not transferred to
this different profile.

Changing every spin mass sign gives the degree -1 antihedgehog. Its two
normalizable zero modes lie in degree three, so Gamma=-1 on its kernel.
The chiral Fredholm index is therefore +2 for the stated hedgehog and -2
for its antihedgehog. These signs use the explicit Gamma convention above.
A continuous bounded chiral perturbation that keeps zero outside the
essential spectrum throughout its path preserves this index by Fredholm
deformation. The statement
does not assume every native repair branch preserves that grading, nor
that arbitrary Hamiltonian perturbations are chiral. No claim for arbitrary
higher-degree profiles is inferred from these two explicit solutions.

## Complete linear-core continuum spectrum

An independent solvable control uses the unbounded supplied texture
`m_j(x)=kappa x_j`, with kappa>0. Its potential is `W=kappa r^2/2` and

```
H_linear^2=-v^2 Delta+kappa^2 r^2-v kappa sum_j rho_j
          =2v kappa (N_oscillator+F).                  (7)
```

The Hermite functions and the exact finite Clifford representation give a
complete orthonormal basis for H^2. The closure of H on their finite span
is self-adjoint, with compact resolvent. There are two Gaussian zero modes.
For integer q>=1 the signed energies and multiplicities are

```
E_q=+/-sqrt(2v kappa q),
g_q=sum_(f=0)^min(3,q) binom(3,f) binom(q-f+2,2)=4q^2+2.
```

H exchanges the two grading sectors within each nonzero H^2 eigenspace;
its positive and negative multiplicities are equal. The total squared-level
multiplicity is 2g_q, including the twofold Clifford multiplicity. The
Gaussian width is sqrt(v/kappa). This confining texture has no finite
asymptotic bulk mass and is not substituted for (2).

## Exact zero modes and full spectral reduction on the original lattice

Use the original source hopping `h_(n,n+e_j)=iv sigma_j/(2a)` on Z^3,
with the supplied onsite native spin mass

```
B_n=kappa a sum_j n_j (-1)^n_j sigma_j.                 (8)
```

The common lattice clock and contact convention are unchanged. Define
`epsilon_j(n)=(-1)^n_j`, the central difference D_a, and the unit phase

```
p(n)=product_j [1+i epsilon_j(n)]/sqrt(8).
```

Since `p(n+/-e_j)=-i epsilon_j(n)p(n)`, conjugation by this phase gives

```
p* h p=sum_j sigma_j B_j,
B_j=epsilon_j A_j,   A_j=v D_(a,j)+kappa a n_j,
B_j*=B_j,           h^2 unitarily equals sum_j A_j* A_j. (9)
```

The B_j commute on different coordinates. No continuum approximation or
projection onto a subset of valleys is used in (9).

Let z=v/(kappa a^2), and let I_n denote the modified Bessel function.
The Fourier coefficients of exp(z cos k) satisfy
`I_(n+1)-I_(n-1)=-2n I_n/z`. Parseval gives
`sum_(n in Z) I_n(z)^2=I_0(2z)`. Consequently the two exact normalized
lattice zero vectors are

```
psi(n)=p(n) product_j I_(n_j)(z) / I_0(2z)^(3/2) xi,
xi in C^2, ||xi||=1.                                  (10)
```

The recurrence cancels every hopping/mass term individually. Fourier
transformation turns A_j phi=0 into a first-order periodic equation with
unique solution proportional to exp(z cos k_j), so (9) proves that (10)
exhausts the zero space. The three Bessel identities here follow directly
from the stated Fourier series and need no numerical small-residual premise.

In one dimension B is the self-adjoint Jacobi operator

```
B_(n,n)=kappa a n(-1)^n,
B_(n,n+1)=(-1)^n v/(2a).                               (11)
```

Its diagonal part has compact resolvent and its hopping is bounded, so B
has compact resolvent. Every eigenvalue is simple: the discrete Wronskian
of two square-summable eigenvectors is constant, and tends to zero at
infinity; all hopping coefficients are nonzero. Reflection n->-n
anticommutes with B. Its spectrum is therefore
`0, +/-b_1, +/-b_2,...`, with `0<b_1<b_2<...` and no missing multiplicity.

Diagonalize the three commuting copies in (9). For each ordered triple
of signed one-dimensional eigenvalues, the remaining two-spin matrix is
`sigma.dot(b)`. This gives the complete original-lattice spectrum

```
+/-sqrt(b_i^2+b_j^2+b_k^2),  i,j,k>=0, b_0=0,           (12)
```

retaining every independent sign of a nonzero b, every coordinate
permutation, and any coincident sums. The all-zero triple has multiplicity
two. The twofold spin row at a nonzero triple supplies one level of each
sign. Equation (12) is a full operator spectral reduction, not an omitted
valley count or a finite-box eigenvalue fit.

## Certified numerical lattice levels

For v=kappa=a=1, lattice_bounds.py certifies the first six b values using
exact rational arithmetic. The finite witness interval is [-32,32]. On
its complement, the diagonal has modulus at least 33 and the hopping has
norm at most one. Thus the tail resolvent for |E|<32 has norm at most
`1/(32-|E|)`. The two boundary couplings have joint norm 1/2. The
Feshbach self-energy is bounded by `1/[4(32-|E|)]`.

Exact tridiagonal LDL inertia proves a spectral gap around every coarse
interval endpoint larger than that bound. Multiplying the two boundary
couplings by s in [0,1] then produces no endpoint crossing. The isolated
finite-block spectral count is therefore exactly the infinite-source count.
The central interval contains only the exact zero; every subsequent listed
interval contains exactly the indicated positive eigenvalue.

Each retained rational vector is an untrusted eigensolver proposal. Its
residual is recomputed exactly after zero extension, including both omitted
boundary entries. The normalized residual bounds distance to the full
self-adjoint spectrum. Combining it with the preceding unique-eigenvalue
count gives the narrow enclosure. No floating eigenvalue or assumed
boundary convergence is used as the acceptance test.

These lattice intervals and the continuum analytic levels are separate
calculations at their stated spacings and profiles. Equality of their
finite-spacing excited spectra is not asserted.

## Original-lattice attachment of the bounded hedgehog levels

The bounded continuum hedgehog (2) also attaches directly to the original
lattice, independently of the linear confining control. At spacing a>0
retain the original hopping `iv sigma_j/(2a)` and set

```
B_n=m sum_j (n_j/|n|)(-1)^n_j sigma_j  for n!=0,
B_0=0,                      H_a=h_(0,a)+B.             (13)
```

This is a bounded self-adjoint operator on the complete original
`ell^2(Z^3,C^2)` carrier. The value at the single lattice origin is part
of this declared profile. It is not discarded as a measure-zero site.

First determine the essential gap at every fixed a. The onsite square
is `B_n^2=m^2 I` except at the origin. Write `u(n)=n/|n|`, with u(0)=0.
For a forward edge in direction j, the hopping-mass anticommutator has
coefficient `imv/(2a)` times

```
-epsilon_j(n)[u_j(n+e_j)-u_j(n)] I
+sum_(k!=j) epsilon_k(n)[u_k(n+e_j)-u_k(n)] sigma_j sigma_k.
```

Every coefficient tends to zero as |n| tends to infinity. A finite sum
of shifts times such multiplication operators is compact: cut each
coefficient off outside a finite cube, obtaining finite-rank operators
that converge in norm. The origin contribution to `B^2-m^2 I` is also
finite rank. Consequently

```
H_a^2=h_(0,a)^2+m^2 I+compact.
```

The Fourier multiplier of `h_(0,a)^2` is the scalar
`(v/a)^2 sum_j sin(k_j)^2`, whose essential range is `[0,3v^2/a^2]`.
Compact-perturbation invariance gives
`sigma_ess(H_a^2)=[m^2,m^2+3v^2/a^2]`. Spatial inversion n->-n
anticommutes with both terms of H_a, including their declared value at
the origin. Essential spectral mapping and this sign symmetry give

```
sigma_ess(H_a)=[-sqrt(m^2+3v^2/a^2),-m]
               union [m,sqrt(m^2+3v^2/a^2)].            (14)
```

Thus every spectrum point in (-m,m) is an isolated finite-multiplicity
eigenvalue, with possible accumulation only at the gap endpoints.

To compare continuum and lattice eigenspaces, write a continuum test
column as its eight two-spin rows f_b, b in {0,1}^3, and define

```
(J_a f)(n)=a^(3/2) sum_b (-1)^(b.n) f_b(a n).           (15)
```

There is no removed valley or chosen sublattice in this map. For smooth
compactly supported fields f,g, its norm/inner-product sum has diagonal
valley terms converging to the continuum Riemann integral. A cross term
has the form `a^3 sum_n (-1)^(d.n) F(a n)`, with nonzero binary d.
Choose j with d_j=1 and shift the sum by e_j. Twice this sum becomes

```
a^3 sum_n (-1)^(d.n)[F(a n)-F(a(n+e_j))].
```

The mean-value bound makes it O(a): the compact support involves at most
C a^(-3) sites and each difference is bounded by a times a fixed first
derivative. Therefore every finite test-field Gram matrix converges to
its full sixteen-row continuum Gram matrix.

For `f in C_c^infinity(R^3\{0},C^16)`, onsite multiplication intertwines
exactly. Multiplication by `(-1)^n_j` flips valley bit j, which is exactly
the original beta_j action. At n=0 both sampled mass actions vanish,
because f vanishes near the origin. The hopping intertwines with the
central-difference derivative: its valley sign is `(-1)^b_j`, and
`gamma_j=-Z_j tensor sigma_j` is the source matrix. Taylor's formula gives

```
||H_a J_a f-J_a H f||_(ell^2) <= C_f a^2.              (16)
```

Indeed the central-difference remainder is bounded by `a^2/6` times a
uniform third derivative on a fixed compact support enlargement. Each
sample carries `a^(3/2)`, there are O(a^(-3)) sites, and the finite sum
over all eight valleys preserves that order. Equation (16) is a test-field
consistency estimate; its constant depends on the chosen smooth field.

The test domain just used is a graph core for (2). Its graph norm is
equivalent to H^1 because the free Dirac square is `-v^2 Delta` and the
mass is bounded. To remove the origin from a compact smooth function,
multiply by a cutoff that is zero in radius delta and one beyond radius
2delta. The additional derivative has squared integral O(delta), while
the removed function and ordinary derivative have squared integrals
O(delta^3). Spatial truncation and ordinary smoothing then prove the
core statement. This is the zero H^1 capacity of a point in three
dimensions, here with an explicit cutoff estimate.

Let E be any continuum gap eigenvalue and d its finite multiplicity.
Choose an orthonormal basis of its eigenspace. Graph-core approximation
gives d compact smooth columns away from zero, as close to this basis
and with residual `||(H-E)f_j||` as small as desired. For these fixed
columns, Gram convergence makes their sampled span d-dimensional for
all sufficiently small a. Equation (16) and norm convergence for each
`(H-E)f_j` make the operator norm of `(H_a-E)` on that sampled span
arbitrarily small relative to its norm.

More explicitly, choose the core errors so that the test-field Gram
matrix has minimum eigenvalue at least 3/4 and each continuum residual
is below delta. For sufficiently small a the sampled Gram has minimum
eigenvalue at least 1/2 and each sampled residual is below 2delta.
Cauchy-Schwarz then bounds the relative residual on the entire sampled
span by `2sqrt(2d) delta`. For any epsilon>0 with
`(E-epsilon,E+epsilon)` inside (-m,m), choose this bound below epsilon.
If the lattice spectral projection onto that interval had rank less
than d, its restriction to the sampled span would have a nonzero kernel.
The spectral theorem would give residual at least epsilon on a vector
in that kernel, a contradiction. In view of (14), the projection counts
actual discrete lattice eigenvalues, including their multiplicities.

Thus, for every such fixed E and epsilon, every sufficiently small a
has at least d lattice eigenvalues in the stated interval. In particular
there are at least two near the continuum zero level and at least
`2n^2` near each signed Coulomb level for fixed n>=2. This is a direct
source-lattice spectral inclusion with preserved lower multiplicities.
It does not exclude additional lattice gap levels, prove equality of
multiplicities, establish an exact bounded-lattice zero mode, or give a
rate uniform in the continuum level or its core approximation. The
linear-lattice exact Bessel zeros belong to their separate profile.

## Interactions and physical interpretation

The [finite interacting core theorem](CORE_PROOF.md) retains the original
centered contact law in a finite box of (8). It fills every negative
one-body level, keeps both zero orbitals and every virtual Fock state,
and controls the four interacting core energies by a gap and a full-norm
Feshbach estimate. This supplies an actual interacting defect energy
splitting on its finite stated domain. It does not turn the quadratic
continuum results into a constructed interacting four-dimensional continuum.

A fermionic field mode trapped by a texture and a dynamical particle made
from the texture are different objects. The latter's rest energy includes
the field creating the texture and its response. This package has no
measured electromagnetic/color/weak charge dictionary, universal mass-unit
calibration, autonomous defect formation or repair protection, or common
Einstein stress tensor. Its conditional premises are explicit: the mass
field, its profile, and the clock. The spectra and the degree/zero-mode
association follow from those premises on the complete source carrier.
