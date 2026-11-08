# Source configuration algebra for multipacket comparison

## 1. One source and the available symmetries

Use the ordered spin source, selected ground Ω, positive generator L,
spatial translations and regular canonical fields from
`spin_canonical_field`, `spin_canonical_energy`, and
`protected_defect_catalog`. The configuration fields Q3(F),Q2(G) strongly
commute for all real admitted tests. The source P2 Weyl group shifts Q3
by its exact coefficient m integral FG and fixes Q2. The latter identity
holds before scaling because Q2 and P2 use the same spin component.
The actual m is positive, with the retained conservative bound m>1/10.

The source circle action rotates Q3,Q2, fixes Ω, and commutes with L and
spatial translations on the configuration-time-history subspace. It is
not assumed to act on the whole canonical net or to rotate P2 into an
additional momentum operator. All uses of this action below are on
configuration operators or their vacuum distribution.

The bounded complex creator is the existing source operator

    B(F)=(x-i y) h(x²+y²),
    x=sqrt(r) Q3(F),    y=sqrt(r) Q2(F),
    h(u)=1/[sqrt(1+u)(sqrt(1+u)+1)].                        (1)

The cancellation-free form of h equals the angular charge projection
of t/(1+t²). It is real analytic for u>=0, with h(0)=1/2.
B has source charge+1, B* charge-1, and norm at most1/2. Its complete
three-record instrument retains the positive square-root complement.
No particle wavefunction or Gaussian inner product is specified here.

## 2. Absolute continuity of the actual configuration distribution

Choose finitely many linearly independent real tests F_1,...,F_d from
the common test space. Their L2 Gram matrix is positive. Dual tests G_j
in their real span satisfy integral F_i G_j=delta_ij. The original P2
Weyl groups for these tests supply strongly continuous unitaries V(a),
a in R^d, whose conjugation shifts the whole x tuple by -a and fixes the
y tuple. The nonzero source coefficient m is absorbed in the group
parameter. This is an application of the constructed canonical
pair, not an added conjugate momentum field.

Let E be the joint PVM of x,y on the common representation and let
mu(S)=<Ω,E(S)Ω>. With the stated sign convention,

    E(S+(a,0))=V(a) E(S) V(a)*.                           (2)

For a Lebesgue-null set N in R^d, Tonelli's theorem gives

    integral da mu((N+a) times R^d)=0.

The integrand is nonnegative and continuous in a by (2) and strong
continuity of V. Its value at zero is therefore zero. Thus the x
marginal of the actual vacuum distribution is Lebesgue absolutely
continuous. The configuration circle symmetry at angle pi/2 sends
(x,y) to (y,-x) and fixes this vacuum distribution. Its y marginal is
therefore absolutely continuous as well. No action on P2 is used.

For a joint Lebesgue-null set S in R^(2d), Fubini says its x section S_y
is null outside a Lebesgue-null exceptional set of y values. The y
marginal assigns that exceptional set zero weight. Consequently

    integral da mu(S+(a,0))
      =integral mu(dx dy) volume_d(S_y)=0.

The same continuity argument in (2) proves mu(S)=0. Hence

    mu is absolutely continuous with respect to Lebesgue measure
    on the complete finite configuration tuple R^(2d).             (3)

This is a vacuum-distribution theorem. It does not assert that every
state on the full canonical net has this joint measure class, that the
circle action extends to P2, or that a second independent canonical pair
has been reconstructed. It also gives no pointwise density, Gaussian
law or quantitative lower density bound.

Distinct translates of a nonzero smooth annular test are linearly
independent. In Fourier space a relation would make a finite exponential
sum vanish on an open set where the test transform is nonzero. Choose a
direction for which all projected translation centres are distinct.
Derivatives along that direction of orders0 through d-1 give an
invertible Vandermonde system. Every coefficient must vanish. Thus (3)
applies to every finite family of distinct translated packet profiles
used in the separated-creation construction.

## 3. No polynomial creation relation is lost on the vacuum

For finitely many tests as above, write B_j=B(F_j). Joint spectral
calculus identifies a polynomial P(B_1,...,B_d,B_1*,...,B_d*) with a
real-analytic function of the2d configuration coordinates. For each j,
the coordinate map in (1) has nonsingular real Jacobian at the origin:
it is (x_j,y_j)->(x_j/2,-y_j/2) to first order. Its image contains an
open neighbourhood. A nonzero polynomial in distinct commuting formal
variables z_j and conjugate(z_j) therefore cannot become identically
zero under this substitution.

The zero set of a nonzero real-analytic function on connected Euclidean
space is Lebesgue null. One can verify this locally by expanding in one
coordinate, applying the one-variable isolated-zero property off the
zero sets of its nonzero analytic coefficients, and using induction and
Fubini. A countable analytic cover supplies the global statement. For a
complex-valued function use a nonzero real or imaginary component.
Combining this fact with (3) gives

    ||P(B,B*) Ω||²>0 for every nonzero polynomial P.        (4)

All these operators are bounded. In particular every finite collection
of distinct creation monomials has a strictly positive source vacuum
Gram matrix. This is a concrete source-derived comparison algebra,
without a guessed Fock metric. Its Gram matrix may be ill-conditioned;
(4) does not supply an effective uniform eigenvalue bound as profiles
separate or their scale changes.

For a chosen finite generator set the homogeneous degree-n, positive
charge monomials give an injective map from its algebraic symmetric
power into the source Hilbert space. The coefficient space is generated
by these actual bounded operators. Since B(F) is nonlinear in F, it is
not identified with the symmetric power of the linear test-function
space. Its inner product is the actual n-point source Gram matrix and
need not equal a tensor power of the one-creation inner product.
The same construction permits mixed signs, graded by total source charge.

## 4. The complete internal charge spectrum and microscopic grading

The configuration circle action is strongly continuous and obeys
R(2pi)=I. Its Fourier projections

    P_q=(1/2pi) integral_0^(2pi) exp(-iq theta)R(theta)dtheta,
    q in Z,

are orthogonal and sum strongly to I, for example by the positive Fejer
averages. The Stone generator therefore has spectrum contained in Z.
`source_charge_scattering/CONTINUUM_PAIRS.md` realizes a nonzero vector
of each integer charge using complete separated creation histories.
Thus its spectrum on the constructed configuration-history
representation is exactly Z. This classifies internal charge values,
not multiplicities, elementary species or electromagnetic charges.

The microscopic spin operators used here are even fermionic bilinears
on the protected singly occupied code. They preserve the original
occupation excess while changing this internal spin charge. Products
of disjoint local spin preparations commute; the configuration-field
creators also commute exactly at equal time. These facts determine the
permutation algebra of the admitted preparation operators. They do not
prove a physical boson spin, a spin-statistics theorem or a Standard
Model gauge representation.

Odd-CAR occupation sectors, nonlocal strings or a different emergent
grading require their own source map, state, energy and detector
attachment. The protected-spin occupation-energy suppression in
`particle_mass_readout_gate` is retained. It excludes a specific
same-clock occupation-charge shortcut; it is not a universal exclusion
of emergent fermions in the source.

## 5. What this supplies for scattering

Equation(4) supplies nondegenerate finite multipacket Gram data and a
source-defined commutative creation algebra. The strict-local positive
charged preparations and their controlled time-smearing are constructed
in [LOCAL_CREATORS.md](LOCAL_CREATORS.md). Asymptotic scattering
requires a source-derived one-excitation dispersion and a comparison
identification whose norm and incoming/outgoing defects are controlled
at large time. Neither strict Gram positivity nor charge-sector counting
is a factorization theorem or an asymptotic particle assignment.
