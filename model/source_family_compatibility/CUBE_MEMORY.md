# Certified frequency-dependent response of the full default cube

The readable-default cube at `c=1/8,U=1` has a resolved memory contribution
to its original fermionic spectral readout. At imaginary frequency `eta=1/4`,
removing that memory changes a complete original filter probability by more
than `23/1000`. The result consumes the same full 65,536-state ground,
source clock and exact composite operators as `CUBE_RESPONSE.md`. It is a
finite-source calculation, without a continuum mass or gap inference.

## Full source and certified ground

The complete half-filled, zero-spin block has 4,900 basis states. The full
ground and gap theorem from `centered_source_ground/gap.py` proves that its
ground is the unique full-Fock ground, with a gap greater than `1/128`.
The producer consumes that certificate, then proposes a finer integer
vector in this same block, with scale `2^31`. The inherited exact Rayleigh,
variance and spectral-separation verifier bounds its half-trace distance
from the true ground by

    delta=1/400000.

The finite-precision eigensolver is only a proposal. Every admitted bound
follows from integer source actions and rational norms. The exact ground
energy is enclosed within `2^-38` of the fixed rational

    Eref=-4124058733378/10^12.

This is a source-derived energy reference. The excitation generator stays
`K-E0`; replacing E0 by Eref in an intermediate solve is charged explicitly.
No trial-state expectation is silently identified with the exact ground.

For one hidden-spin flavour, all eight a creators and all eight b creators
map the ground block into the entire `N=9,Sz=1/2` block of dimension 3,920.
The opposite-spin sector has the same response by the exact full-ground
spin symmetry. The complete particle-hole involution of the literal full
source supplies the removal block, with reversed retarded time and phases
`S=diag(epsilon_x,-epsilon_x)` on the a/b columns. The physical removal
excitation energies are positive. No original charge record is discarded.

With W=[V0,V1], the true sixteen-column response obeys

    G(z)=W*(z-Lret)^(-1)W,
    Sigma(z)=zI-A-G(z)^(-1),
    A=[[h,I/2],[I/2,C]].                               (1)

Only the lower-right block of Sigma is nonzero. The actual C is the bond
correlation from `COMPOSITE_CAR.md`, evaluated with an independent even
source observable. Its trial proposal differs from the true C by at most
`(9/8)delta` in operator norm. This is the three-neighbour cube bound;
no infinite-lattice momentum is assigned to the open cube.

## Rational residuals for every numerical solve

Let p be the integer ground proposal and `nu=||p||^2/2^62`. The producer
checks `3/4<=nu<=5/4`. The unnormalized charged right-hand sides are exact
integer actions on p, divided by `2^31`. For each eta the complex sparse
solve at `z=i eta` is rounded to denominator `2^40`. The residual in

    [Eref+i eta-K] X=B

is evaluated over exact integers with a common rational denominator. Its
Frobenius norm is at most `rho=2^-27`. Hermiticity of K bounds the solution
error by `rho/eta`. Normalization and both retarded charge blocks contribute
at most `4rho/eta` to the error of G.

The removal proposal obtained by conjugating the addition proposal can
use a particle-hole transformed trial vector. Both trial vectors are
within `sqrt(2)delta` of the same exact invariant ground after a phase
choice. Each a or b response block therefore differs from its true block
by at most `2delta`, and the combined W differs by at most
`epsilon=3delta`. This conservative estimate does not assume the rounded
trial is an exact symmetry eigenvector. Thus the retained G proposal has
operator error

    gamma=(2epsilon+epsilon^2)/eta
             +(1+epsilon)^2 2^-38/eta^2+4rho/eta.      (2)

The first term is ground error, the second is the scalar energy reference,
and the third is the certified linear-solve error.

The inverse of the rational G proposal is itself proposed numerically and
rounded to denominator `2^40`. Let Q denote that dyadic inverse proposal,
q an exact row/column-sum norm bound on Q, and r an exact norm bound on
`I-Gproposal Q`. If `e=r+gamma q<1`, then

    ||Gtrue^(-1)-Q|| <=q e/(1-e).                       (3)

Every condition is checked with rational arithmetic. Adding the C error
gives the complete operator enclosure for the memory in (1). The
memory-free comparison `(zI-A)^(-1)` has a separately certified inverse
residual and the resolvent error `||delta C||/eta^2`. A Ritz value or a
rounded inverse is never accepted as an exact spectral quantity.

## Complete original records and the resolved memory effect

For a fixed local mode let `A_+=a*`, `A_-=a`, and let
`Delta=K-E0` on the two physical charge sectors. The exact finite filter is

    F_eta=eta(eta-i Delta)^(-1).

The four Kraus records

    F_eta A_+, sqrt(I-F_eta*F_eta) A_+,
    F_eta A_-, sqrt(I-F_eta*F_eta) A_-

sum to the identity effect on the original full source. On the invariant
ground their total successful probability is the original Lorentzian score

    p_eta=int eta^2/(eta^2+E^2) d(mu_a^++mu_a^-)(E)
         =-eta Im G_aa(i eta).                        (4)

Each local charge contributes half this value; the two failure probabilities
are each `(1-p_eta)/2`. No normalization by success is used. This ideal
finite instrument is accessible through the stated complete coherent
finite-instrument service or finite spectral readout approximations. Its
large native tape is not executed here; a service half-trace error must be
added to the original probability enclosure. The retarded sixteen-column
calculation is not a universal physical preparation channel for arbitrary
coefficient inputs and references.

The exact fixed enclosures admitted by fresh source solves are:

| eta | Im Sigma_bb(i eta) | original p_eta | memory-free comparison |
| --- | --- | --- | --- |
| 1/4 | [-0.0395,-0.0392] | [0.2175,0.2177] | [0.1943,0.1945] |
| 1/8 | [-0.0239,-0.0235] | [0.0667,0.0669] | [0.0569,0.0571] |

The jointly certified difference at eta=1/4 exceeds `23/1000` and is less
than `234/10000`. The imaginary memory entry changes by more than
`15/1000` between these two frequencies. These statements are admitted
from the full source with ground, clock, solve and inverse errors included.
They do not follow merely from a positive R*R or a compressed eigenvalue.

The finite cube has no translation momentum. The separate all-volume
`NODAL_MEMORY.md` bound concerns the genuine magnetic momentum of the
specified invariant ground average. Neither result proves a nonzero
original weight in a shrinking ballistic energy window. A dressed branch
can coexist with this nonzero memory; no sharp subbranch is ruled out by
the failure of the memory-free comparison.
