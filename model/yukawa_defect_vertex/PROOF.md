# Complete Higgs intervention and vertex records

The source is the existing [quantum Higgs construction](../coherent_defect_wave/chiral_higgs/PROOF.md):
six lepton CAR modes, the five scalar states, the literal native charged-lepton
pair tensors, scalar number energy one, and primitive lepton interaction one.
The source clock is `exp(-itH)`. Quark occupations are the declared empty
spectator preparation. No family matrix or measured particle input enters.

The executable reconstructs every matrix entry and checks equality with the
owning source producer. It retains all 320 active occupations, not a vacuum
cyclic block or one of the joined sectors satisfying `N_L=N_H`. It verifies
all twelve total-current commutators and the exact rank-two gauge-singlet
ground projector. The source marginal is `rho=P_ground/2`, with energy minus
two. Its preparation uses the existing complete energy measurement on the
maximally entangled source/dual-reference state. All fourteen preparation
energies are retained. Ground records carry their original probability
`1/160`; the other thirteen energies are not discarded.

## Complete two-point functions

For every one of the fourteen distinct source energies E, calculate the
entire matrices

`C+_ij(E)=Tr[rho O_i P_E O_j^dagger]`,

`C-_ij(E)=Tr[rho O_j^dagger P_E O_i]`.

The polynomial spectral projectors exhaust all 320 dimensions, including
the 198-dimensional energy-one space. These matrices retain the coherent
reference coordinates, and zero weights are explicit. With `Delta_E=E+2`,
the fermion and scalar retarded resolvents are respectively

`G_ij(z)=sum_E [C+_ij(E)/(z-Delta_E)+C-_ij(E)/(z+Delta_E)]`,

`D_ij(z)=sum_E [C+_ij(E)/(z-Delta_E)-C-_ij(E)/(z+Delta_E)]`.

The fermion sum is exactly `sum_E(C++C-)=I6`, as required by the full CAR
algebra. For the two scalar components,

`sum_E C+=sum_E C-=(7/10) I2`,

`sum_E(C+-C-)=0`,

`Z_H=sum_E Delta_E(C++C-)=(14/5) I2`.

The executable independently obtains the last quantity from
`Tr[rho [phi_i,[H,phi_j^dagger]]]`. Thus `G(z)=I6/z+O(z^-2)` and
`D(z)=Z_H/z^2+O(z^-3)`. Both the vanishing first scalar term and the positive
full-rank second term follow from the complete records. This defines the
finite spectral convention `phi_canonical=Z_H^(-1/2) phi`. It does not assert
a Klein-Gordon limit, a pole residue at a selected energy, or the physical
Higgs kinetic coefficient.

The fermion instrument is the coherent isometry

`V_F=(sum_i c_i^dagger tensor |i dual,+> + c_i tensor |i,->)/sqrt(6)`.

Its squared norm is identity by CAR. The addition reference carries the
dual one-body representation, the removal reference the original one; both
have odd parity. Consequently each branch preserves total gauge charge
and combined parity. The internal reference indices remain quantum. Reading
only the branch and source energy gives all 28 outcome probabilities
`Tr C+/6` and `Tr C-/6`.

For the scalar probe use both component orientations with Kraus operators
`phi_i/2` and `phi_i^dagger/2`, retaining their dual references. Their positive
complement is `K0=diag(0,sqrt(3)/2,sqrt(3)/2,sqrt(3)/2,sqrt(3)/2)` on the scalar
memory. This gauge-invariant no-operation outcome has ground probability
`3/10`. All three branches and all fourteen energies give 42 records. The
unconditional fermion and scalar experiments contain 41 and 55 records,
respectively, after including every failed ground preparation.

## The three-point reconstruction

The source has the all-left Weyl convention, so the relevant fermion legs
are an anomalous pair. Define the actual equation-of-motion insertion

`K_ij={ [c_i,H],c_j}=i d/dt {c_i(t),c_j}|_(t=0)`.

The full 320-state operator identities give

`K_04=K_13=-phi_0`, `K_05=K_23=-phi_1`,

with the opposite signs on reversed pairs and zero on all other pairs.
These identities fix the scalar conjugation and the CAR signs without
projecting any intermediate charged state away.

Measure the mixed three-point traces `T_ij,k=Tr[rho K_ij phi_k^dagger]` and the
scalar Gram `C_ab=Tr[rho phi_a phi_b^dagger]`. The matrix C is positive and
invertible, so the coefficients are `A_ij=T_ij C^(-1)`. Reconstructing
`sum_k A_ij,k phi_k` returns K as an operator on every source state. The
nonzero coefficients have modulus one. Fermion normalization is already
fixed by the CAR spectral sum; substituting `phi=sqrt(14/5) phi_canonical`
therefore gives contact coefficient modulus `sqrt(14/5)`.

This is the contact vertex in the exact field equation, extracted from
complete finite three-point traces. It is not the full amputated one-particle
irreducible vertex at arbitrary external frequencies. The primitive unit
coefficient was specified by the source architecture; the calculation
derives the finite normalization and readout consequences of that law,
without claiming that idle defect repair uniquely forces the law.

## An actual charged Higgs intervention

Adjoin a three-state apparatus reference: a blank and a dual doublet. Let
`R_k^dagger=|k dual><0|` and define

`X=sum_k(phi_k^dagger tensor R_k^dagger+phi_k tensor R_k)`.

For each native Hermitian Higgs generator t the apparatus generator is
`diag(0,-t^T)`. Direct full-matrix commutators show that X preserves all twelve
total Gauss generators. This is a scalar creation/annihilation pulse;
neither the scalar cost nor the coefficient multiplying the interaction is
being varied. Its reference states permit the charged transitions excluded
by the smaller joined sectors.

Let P_k be the literal fermion pair annihilator in
`H_Y=sum_k(P_k phi_k^dagger+P_k^dagger phi_k)`. The gauge-invariant pair pointer
is `A=i(B-B^dagger)/8`, where `B=sum_k P_k tensor R_k`. The bound
`||P_k||<=2` implies `||A||<=1`; hence `(I+A)/2` and `(I-A)/2` are a complete
positive pointer measurement. Its Lüders instrument retains the source and
quantum reference. Both outcomes are reported with their original weights.

The pulse is `U(theta)=exp(-i theta X)`. The complete scalar/reference
generator has squared eigenvalues zero, one and two. Its exact polynomial
spectral formula keeps the whole null space. The executed finite-pulse
pointer expectation agrees at zero pulse with

`d/dtheta <A> = i Tr[(rho tensor |0><0|)[X,A]]`.

This nonzero scalar-insertion response is a gauge-invariant three-point
record. Its exact value is `-3/10`. The full finite-pulse expectation is

`<A>_theta=-(3/20) sin(theta) cos(sqrt(2)theta)
            -(3sqrt(2)/40) sin(sqrt(2)theta)`.

At `theta=pi/(2sqrt(2))`, the two pointer probabilities are
`1/2 +/- 3sqrt(2)/80`. The receipt retains both probabilities with their
original `1/160` preparation factor and all thirteen failed preparations.
The idle preparation reference remains retained throughout; no charge or
event weight is supplied by dropping a reference branch.

Spatial transport, a three-family construction, on-shell amputation and a
physical particle interpretation require their own joined calculation.
