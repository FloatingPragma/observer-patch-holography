# Infinite-photon confined particle theorem

## Declared source and its domain

Retain the complete matter oscillator shells n<=2, their twenty spin/spatial
coordinates, mass16, field coefficient1/8, radial form factor and trap of
`openwave_quantum_rotational`. Replace its finite photon cap by the full
three-mode bosonic Fock space. In Gaussian field coordinates,

    B(x)=B0+sum_i M_i sqrt(2)x_i/8,
    K(x)=sqrt(16^2+B(x)^2),
    H=K+N_matter+N_photon.

B is self-adjoint finite matrix multiplication on its maximal domain.
The closed nonnegative quadratic-form sum defines H. Finite Hermite
polynomials form a form core: K has at most linear coordinate growth and
is controlled by N_photon+1. Compactness follows from H>=16+N_photon
and finite matter dimension. The CAR lift uses the same pointwise positive
matrix K and retains the field blocks in their original order. Thus the
full source obeys H>=16 N_ferm+dGamma(N_matter)+N_photon. Vacuum energy is
zero and every higher fermion sector starts at32 or above.

Complete oscillator shells and the vector photon oscillator give exact
rotation covariance. The native Clifford carrier has eight identical
reduced spin/spatial copies after taking the complete absolute Dirac kernel.
All copies and both signed Dirac branches are retained.

## Exact isolated spin-half band

Let P0 be the two spin states with matter and photons in their vacua. The
scalar inequality sqrt(16^2+x)<=16+x/32 gives, on this entire doublet,

    P0 H P0 <= (16+998453/21233664) P0 < (16+1/16)P0.

The value follows from the exact matrices B0 and M_i and vacuum oscillator
contractions. The min-max principle gives two eigenvalues below17. Since
H>=16+N_matter+N_photon, its third eigenvalue is at least17. Denote the
entire two-dimensional low-energy space by G. The restriction P0:G->P0 is
injective: a nonzero vector in its kernel would have energy at least17.
This restriction is an equivariant bijection. Consequently G carries the
irreducible spin-half representation, its exact Casimir is3/4, and H is a
scalar E_infinity on G by Schur's lemma. In particular the two eigenvalues
are exactly degenerate. Temple's bound can therefore use the next distinct
energy17 for any trial vector, without choosing a numerical Jz eigenvector.
The full eight-copy carrier has sixteen ground states. No thermodynamic or
spatial continuum assertion follows from this fixed matter regulator.

## Exact moment shortcut and global square-root bound

For x>=0 Taylor's integral remainder gives

    0 <= sqrt(16^2+x)-p6(x) <= (33/2048) x^7/16^13,
    p6(x)=16+sum_(k=1)^6 binom(1/2,k)x^k/16^(2k-1).

The seventh derivative is positive and decreases on the positive axis.
Functional calculus consequently bounds the vector error for any polynomial
photon vector phi by

    ||(K-p6(B^2))phi|| <= 33 ||B^14 phi||/(2048*16^13).

This estimate needs no global bound on the unbounded B. Each B application
changes total photon number by at most one. For a trial supported through
cap C, padding to C+14 computes B^k phi exactly for every k<=14. In
particular p6(B^2)phi and its fourteenth moment are genuine infinite-source
quantities, not finite-source extrapolations. The dyadic trial coefficients
in the sealed finite receipt are only trial data; every energy moment is
evaluated afresh against this declared infinite operator.

Let [Rlo,Rhi] enclose <phi,H6 phi>, let r bound the corresponding polynomial
residual, and let delta bound the vector remainder above. Then the actual
Rayleigh energy lies in [Rlo,Rhi+delta], and its variance is at most
(r+delta)^2. Projection onto the orthogonal complement of phi removes any
scalar Rayleigh correction. Temple's inequality gives

    Rlo-(r+delta)^2/(17-Rhi-delta) <= E_infinity <= Rhi+delta.

The trial's norm outside G is at most e=(r+delta)/(17-Rhi-delta). Its
normalized projection into G has pure-state trace distance at most2e from
phi. Since phi has no photons above C, that projected state has photon-tail
probability at most2e^2. The photon cutoff projector commutes with rotations;
its compression to the irreducible G is scalar. The same tail bound therefore
holds for every state in G, including arbitrary mixtures. Identical-copy
factorization extends it to the full sixteen-dimensional low band. The fixed
matter-space radius-five bound from the exact finite Hermite polynomials
also applies to all states with this infinite photon field.

## Arithmetic certificate

The implementation encloses exact algebraic matrix entries at45 decimal
digits, then expands both endpoints outward by one binary64 unit. Real and
imaginary coordinates are embedded into one real symmetric sparse operator.
Exact zero entries are omitted, and assembly explicitly rejects duplicate
coordinates. Each matrix row has fewer than256 terms. With binary64
round-to-nearest arithmetic, the standard gamma_(n+4) error is below2^-43;
the matrix-vector enclosure uses2^-40 and a further relative guard16*2^-40.
It propagates input radii and coefficient radii separately. Scalar products
use math.fsum, a2^-48 product/sum guard and an additional relative guard.
Every upper endpoint is expanded with nextafter. An absolute1e-300 allowance
covers accumulated subnormal errors; the run asserts finite dimensions and
contains no overflow. This is an explicitly stated IEEE-754 arithmetic
certificate, not a formal theorem about an arbitrary BLAS implementation.

The exact normalized integer trial vector is enclosed before any matrix
operation. The result records padding, sparse row size, moment norms,
polynomial residual, Taylor bound, energy and photon-tail bounds. Independent
cap1,2,3 witnesses must give intersecting infinite-source energy intervals.
The GeV comparison uses the identical finite-cap2 electron unit without
retuning that anchor. An interval containing the anchor is a regulator
comparison, not a second prediction of the electron mass.

The transverse field strength, mass parameter, spatial trap, absolute-kernel
law and empirical electron unit are declared inputs. This result removes
the photon cap only. It supplies no spatial trap removal, longitudinal local
Gauss law, observed-particle identification or infinite-field magnetic moment.
