# Winding degree and the full-valley defect index

The source kinetic and spin-mass matrices are the sixteen-row folded
matrices of [native_mass_channels](../native_mass_channels/README.md).
All eight Weyl valleys and the two original spin rows are retained.
The texture is a supplied field; its selection and formation energy are
separate from the spectral theorem.

## A degree family with an explicit configuration space

For a positive integer k use stereographic coordinates on the unit sphere
and the rational map `w -> w^k`. In Cartesian coordinates, write

```
a=x+i*y, b=r-z, A=a^k, B=b^k,
N_k(x)=(2 Re(A conjugate(B)), 2 Im(A conjugate(B)), |A|^2-|B|^2)
       /(|A|^2+|B|^2).
```

The north-pole value is `(0,0,1)`. This is a smooth sphere map including
both stereographic poles. On the radial extension, its first and second
derivatives are respectively `O(1/r)` and `O(1/r^2)`. It has unit length.
For orientation eta=+1 or -1, set `N=eta*N_k` and `d=eta*k`.

The sphere degree is d. For example, the normalized area integral in the
stereographic radial variable t is

```
degree(N_k)=2*k^2 integral_0^infinity t^(2k-1)/(1+t^(2k))^2 dt
           =k integral_0^infinity du/(1+u)^2=k.
```

Here `u=t^(2k)` and the last integrand has primitive `-1/(1+u)`.
Target reversal multiplies the degree by -1. Under spatial inversion,
the transverse components acquire `(-1)^k` and the third component changes
sign. Thus `N(-x)=-N(x)` exactly when k is odd. Degree is preserved under
continuous nonvanishing boundary-mass deformations, independently of any
particle interpretation.

## Exact comparison with the Callias operator

Let gamma_j and beta_j be the actual folded generators, with
`Gamma=i gamma1 gamma2 gamma3 beta1 beta2 beta3`. The calculation constructs
an explicit sixteen-row unitary U satisfying

```
U^dagger gamma_j U = sigma_j tensor I2 tensor sigma_x tensor I2,
U^dagger beta_j U  = I2 tensor sigma_j tensor sigma_y tensor I2,
U^dagger Gamma U   = I2 tensor I2 tensor sigma_z tensor I2.
```

The construction uses the three CAR annihilators `(beta_j+i gamma_j)/2`,
their two-dimensional joint vacuum, and all eight occupation words. It
checks every generator and both unitarity products exactly. The last
factor has dimension two; no valley projection is made.

For v,m>0 consider the bounded prescribed continuum operator

```
H=-iv gamma.dot(grad)+m beta.dot(N(x)),  domain H1(R3,C16).
```

In the comparison basis its positive-to-negative chiral block consists
of two copies of `i L^dagger`, where
`L=v sigma.dot(grad)+Phi` and `Phi=m N(x).sigma` acts on the distinct
two-dimensional internal factor. In particular the physical spatial
derivative Pauli matrices commute with Phi in L.

For a smooth Hermitian bounded Phi, invertible outside a compact set,
with first derivatives `O(1/r)` and second derivatives `O(1/r^2)`, use
the [Callias index formula of Gesztesy and Waurick](https://arxiv.org/abs/1506.05144),
equation (0.1) and Theorem 10.2. Positive v can be absorbed in Phi without
changing its normalized exterior direction. The exact Pauli trace is

```
tr[(n.sigma)(u.sigma)(v.sigma)]=2i*n.(u cross v).
```

Substituting this in the three-dimensional boundary formula gives
`index(L)=-degree(N)`. Taking the adjoint reverses the index; the constant
factor i does not change it. Both comparison copies therefore give

```
index(Pminus H Pplus)=2*d,
dim ker(H)>=2*abs(d).                                  (1)
```

This sign agrees with the explicit degree +1 hedgehog zero modes, which
have Gamma=+1. The positive-to-negative operator for d<0 has negative
index. Equation (1) does not exclude additional zero modes of opposite
chiralities and does not classify the nonzero bound levels.

The displayed constant-norm radial extension is singular at the origin.
To apply the smooth theorem, multiply it near zero by a smooth radial
function flat at zero and equal to one outside a ball. All derivatives
then exist and are bounded; the mass is invertible at infinity and has
the stated exterior decay of derivatives. The difference from the bounded
singular coefficient has compact support. Multiplication by this difference
is compact from H1 to L2 by local Rellich compactness. Along the bounded
chiral interpolation the off-diagonal operators remain Fredholm, so their
indices agree. Both full operators are self-adjoint on H1 by bounded
perturbation of the free Dirac operator. The value at the single origin
does not change either multiplication operator.

The imported analytic theorem is identified explicitly; the finite code
checks its matrix dictionary, oriented trace and degree calculation.
It does not numerically prove Fredholm theory or claim a Lean formalization
of that analytic result.

## A different exact index on the original lattice

On the original two-spin lattice use the unchanged links
`h_(n,n+e_j)=iv sigma_j/(2a)` and the prescribed onsite texture

```
B_n=m sum_j (-1)^n_j N_j(n) sigma_j,  n!=0;
B_0=0.
```

For odd k, spatial inversion S anticommutes with H and fixes the two
origin rows. The argument of [radial_defect_gap](../radial_defect_gap/PROOF.md)
extends to every such unit-length odd texture: it uses no radial differential
equation. For a rational kinetic bound K with `K^2>=3(v/a)^2` and m>K,
let P project onto the origin and Q=I-P. Then

```
H=[[0,V^dagger],[V,D]],  D=QHQ,
||D^(-1)||<=1/(m-K),  V^dagger V=3v^2/(2a^2)*I2,
V^dagger D^(-1)V=0.
```

The last equality follows from inversion parity. Thus every zero vector
is `(xi,-D^(-1)V xi)` and the kernel is exactly two-dimensional. For
`0<|E|<m-K`, its Schur complement is

```
-E*[I2+V^dagger(D-E)^(-1)D^(-1)V].
```

The bracket is positive definite since `lambda*(lambda-E)>0` on the
complement spectrum. Hence there is no additional spectrum in that interval.
The origin probability, contact-overlap and local Neumann bounds of the
parent theorem also transfer: their remaining hypotheses are unit onsite
mass, nearest-neighbor hopping and the same two-row origin map.

These statements hold on the full infinite lattice and every symmetric
finite open cube. For m=4 and a=v=1, K=2 gives a nonzero gap at least two
and exactly two zeros for degrees ±1, ±3, ±5 and every other odd degree.
They do not identify this nonlocal inversion grading with the continuum
internal Gamma. Fixed m and v with a approaching zero leave the strong-mass
domain. Thus a degree-three continuum index of six is compatible with a
strong-mass lattice kernel of two. There is no theorem here transporting
that kernel unchanged through the continuum limit.

Even k fails the inversion-odd premise. No inversion-protected kernel count
is assigned to those cases. The numerical search retains them as controls.

## Finite search and physical scope

The declared search uses d=1,2,3,4, cube radii 1,2,3, m=1 or 4 and a=v=1.
It diagonalizes every finite matrix, retaining all eigenvalues. The window
`|E|<m/4` is fixed before the calculation. Its rank, origin weight and outer
boundary weight are whole-projector traces, so arbitrary rotations among
degenerate eigenvectors do not select favorable modes. Numerical values
and residual admission are diagnostics, not exact spectral enclosures.

After the initial 24 cases, variation of the degree-three window rank
motivates four further volume controls: d=1 and 3, m=1 and radii 4 and 5.
These adaptive trials are labeled separately and all four are retained.
The total count is 28; no preselection or prospective physical test is
claimed. Degree three has window rank six at both larger volumes, with
two numerical zeros and doublets at ±0.1657820223 and ±0.1451903447,
respectively. Its boundary projector trace falls from 0.385777 to 0.154790.
This finite trend does not prove an infinite-volume bound-state limit.

Negative orientations satisfy `H[-d]=-Ds H[d] Ds`, where
`Ds(n)=(-1)^(n1+n2+n3)`. The source particle-hole symmetry pairs ±E, giving
the same complete spectral multiset for the two orientations. A vanishing
finite eigenvalue or a small window energy alone is not a continuum bound
particle; boundary weight and volume dependence are retained for that reason.

The index establishes a protected net zero-mode multiplicity for the
supplied continuum texture. Identifying that multiplicity with spin, family,
color or observed particles requires the same-current field intertwiner.
Degree three is one admitted input among infinitely many integer degrees.
Neither the native mass family nor this index calculation selects it, fixes
Yukawa couplings, supplies texture formation energy, or proves that native
repair preserves the protecting symmetry.
