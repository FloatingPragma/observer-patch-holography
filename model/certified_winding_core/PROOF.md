# Certified finite winding spectra and interacting poles

The action is the literal two-spin cubic-lattice hopping and prescribed
rational-map mass texture of [winding_sector_search](../winding_sector_search/PROOF.md).
The spacing and hopping velocity are one, the mass is one or four, and the
open cubes and texture degrees are exactly the 28 retained search cases.
Every finite site and spin row is included. These choices specify a source
family; the theorem does not select its parameters or texture from repair.

## Rational enclosure of the actual source

At each nonpolar integer site, compute `(x+i y)^k` with Gaussian integers.
The other homogeneous coordinate is `(sqrt(x²+y²+z²)-z)^k`. An integer
square root encloses the radical by rational endpoints with denominator
`2^64`. Interval multiplication and division enclose all three normalized
texture components, including numerator/denominator dependence by outward
interval arithmetic. The denominator is strictly positive. Polar values
and the origin are assigned exactly as in the source.

Round each onsite coefficient to denominator `S=2^26` and retain its
maximum endpoint error, rounded upward to denominator `2^48`. This last
rounding keeps a common denominator without weakening containment. The
nearest-neighbor hopping coefficients are exactly representable. The
resulting matrix Hq is exactly Hermitian. If the three coefficient errors
at site x are e_j(x), the exact Frobenius error obeys

```
||H-Hq||² <= ||H-Hq||_F² <= 2 sum_(x,j) e_j(x)².       (1)
```

Thus the reference matrix is the actual source, rather than a floating
matrix whose construction error is ignored.

## A complete eigenbasis certificate

A numerical eigensolver proposes a square complex Q and sorted real
diagonal Lambda. Round every component to denominator S. The proposed
columns need not be exactly orthonormal. Compute

```
G=Q†Q,       R=Hq Q-Q Lambda
```

with exact integer products and common dyadic denominators. Frobenius
norms bound the operator norms. All integer-product partial sums and
recombinations are bounded before execution; squared sums use arbitrary-
precision Python integers. To accelerate dense products, `exact_dot`
splits each signed integer into two base-8192 limbs. Each binary64 limb
dot has every product and partial sum strictly below 2^53 in magnitude,
so all are exactly representable integers. Int64 recombination has a
separate bound below 2^63. This relies on ordinary IEEE binary64 dot
arithmetic, not on a small floating-point residual or an eigensolver
accuracy guarantee. Independent controls compare these products with
arbitrary-precision integer multiplication.

Let `||G-I||<=g<1`, `||R||<=r`, `||Lambda||<=L` and `||H-Hq||<=s`.
The polar factor `U=Q G^(-1/2)` is square and unitary. Moreover

```
||G^(-1/2)|| <= 1/(1-g),
||G^(-1/2)-I|| <= g/(1-g),     ||Q|| <= 1+g.
```

Using `Hq U-U Lambda = R G^(-1/2)+Q[Lambda,G^(-1/2)]` gives

```
||H-U Lambda U†|| <= s+[r+2L(1+g)g]/(1-g).             (2)
```

Weyl's ordered eigenvalue inequality applies to these two Hermitian
matrices. Every ordered eigenvalue of H lies within (2) of the respective
diagonal entry of Lambda, including multiplicities and ties. This is an
exhaustion statement; a single small residual vector would not suffice.

For portable replay the certificate checks the exact inequalities

```
g <= 2n/S,       r <= 2(m+2)n/S,
s <= n/S,       L <= m+2,
```

where n is the full one-body dimension. The receipt uses twice the
resulting bound around dyadic centers proposed by the retained search.
Every fresh certified interval must lie inside its respective retained
target interval. Those numerical targets have no evidentiary role until
this independent admission succeeds. Raw basis vectors and platform-
dependent residuals do not define the receipt's exact claims. Every trial
is retained, including the four explicitly adaptive volume controls.

## An exact kernel outside the strong-mass domain

For odd degree the exact texture is inversion-odd, so spatial inversion S
anticommutes with H. On an odd cube, each non-origin site is paired with
its opposite. The two origin spin rows therefore give
`dim H_plus-dim H_minus=2`. The off-diagonal rectangular chiral block has
index two for every mass: rank-nullity gives

```
dim ker(H)|plus - dim ker(H)|minus = 2.                 (3)
```

If exactly two intervals from the complete spectral enclosure contain
zero, (3) forces exactly two zero eigenvalues, both inversion-even. Every
remaining interval is separated from zero, giving a strict one-body gap
d. This finite argument needs no assumption `m>sqrt3`. For even degree,
two zero-containing intervals do not certify a kernel. Such cases retain
an unresolved kernel count. If all intervals avoid zero, invertibility
is certified independently of parity.

Neither this finite inversion index nor its spectral gap is the continuum
Callias index. No volume-independent or zero-spacing limit follows.

## The same full-Fock contact

Suppose (3) and the exhaustion certificate give the two-dimensional kernel.
The source symmetries `S_stagger i sigma2 K` and its product with spatial
inversion give the same particle-hole and Kramers maps as in
[radial_interacting_core](../radial_interacting_core/PROOF.md). Choose a
symplectic basis of the two zero modes. Since inversion is identity on
this kernel, each site block satisfies `W_x†W_x=rho_x I2`, with
`sum rho_x=1`. No inverse of the origin complement is required here, and
the origin weight need not be positive. Consequently

```
I=sum rho_x² >= 1/V,                                  (4)
```

where V is the number of sites. For m>2 the independent strong-mass bound
also gives `I>=[1+3/(2(m-2)²)]^(-2)`; the larger lower bound is retained.
The paired filled-sea Wick calculation in the parent proof depends only
on these symmetries and the exhausted kernel. It therefore gives the
same four-state compression of the centered contact,

```
P W P = C+U I[n1 n2-(n1+n2)/2].                       (5)
```

There are V-1 occupied negative modes in the sea. The quartet has particle
numbers V-1,V,V,V+1, and its full-Fock complement has free energy at least d.
Keep the actual interaction `U sum_x(n_up-1/2)(n_down-1/2)` on every Fock
sector. Its norm is `w=UV/4`. If `w<d/4`, put

```
epsilon=w²/(d-2w),       t=w/(d-2w).
```

The full complementary Feshbach inverse, rather than a discarded virtual
sector, bounds each low eigenvalue correction by epsilon. If
`U I_lower/2-2epsilon>0` and `3t<1/2`, the exact ground is a Kramers pair
at particle number V. The two unique neighboring charged states have an
equal gap in `[U I_lower/2-2epsilon, U/2+2epsilon]`. Their positive
addition/removal pole matrices have rank two and error at most3t from
half the exact zero projector. These are the same finite full-Fock
consequences as the parent theorem with its one-body premises certified
on this source. No matrix of size `2^(2V)` is claimed to be executed.

The supplied contact grid is `1/65536, 2^-24, 2^-40`, independent of any
target particle masses. Every failed sufficient inequality is retained
as inconclusive, not as a negative physical mass or an exclusion. This
establishes finite interacting poles on admitted inputs. It supplies no
moving-defect rest mass, physical current identification, native vacuum
preparation or interacting relativistic continuum.
