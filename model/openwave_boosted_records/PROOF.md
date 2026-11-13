# Moving records on the longitudinal matter source

## A fixed source family and its reducing sector

At each spatial lattice address place one complete two-cell, 192-mode
longitudinal matter block from `openwave_common_electromagnetism`. This
replication is a declared construction. The internal two-cell Green energy
is not reinterpreted as the external lattice's three-dimensional Coulomb
potential. Every block retains its native strong action, semileptonic current,
quantum flux/reference system and self-energy convention. Their complete
physical space is the tensor product of the Gauss isometries W.

For each of the actual 48 positive Dirac excitation modes per block, with
mass m_j>0, add the positive one-body matrix

    K_j = sqrt(m_j^2+L_a)-m_j,
    L_a(p)=4 sum_i sin(a p_i/2)^2/a^2.

Its full Fock extension is dGamma(K_j). Dress every matrix element by the
product W. A hop between blocks changes their quantum reference charges,
as well as their rotor flux labels; dressing only the leaf flux is
insufficient. Positive spectral calculus proves K_j>=0 on every retained
momentum, so the complete full-Fock lower bound and unique physical vacuum
of the block source survive. Zero momentum attains every single-lepton rest
mass and the full charge-one threshold. No charged or neutral branch is
removed. The square root is nonlocal in the spatial lattice; finite-range
hopping and relativistic microcausality are not claimed.

The colourful vacuum with exactly one positive Dirac excitation is an
exact reducing sector: the semileptonic term annihilates the colourful
vacuum, and the self-energy-normalized electric term vanishes on one
charged excitation. Neutral branches also have zero electric energy.
Thus the exact sector Hamiltonian is

    E_j(p)=sqrt(m_j^2+L_a(p)).

This is a restriction of the specified interacting source family, rather
than a claim that the full Fock evolution has been numerically enumerated.
The source is distinct from the 90-mode first-order Higgs-Wilson source;
no finite-range, chirality or compiler theorem transfers between them.

## Complete native records and apparatus work

The literal Bogoliubov coefficients c=Aa+Ba* define the 96 by 48 isometry
T_+=(A; conjugate B). The charge identities and T_+*T_+=I are executed.
After the unitary spatial Fourier transform, its rows define a complete
native particle/hole probe POVM on the one-excitation sector. These probe
outcomes are not asserted to be bare many-body occupation probabilities.
An absorbing realization has Kraus operators |vacuum><phi_(x,alpha)|,
with a charge-taking apparatus/reference and one distinct classical record
for every x and alpha. Their effects sum to the identity in the sector.
The other Fock sectors and Gauss-constraint failures have explicit
complement outcomes. Their probability is zero for the stated input.

A declared preparation apparatus creates the normalized stated Fourier
packet with probability 1/2; the complementary record leaves the source
vacuum. The source preparation work is therefore half the packet energy.
The destructive probe returns the source to vacuum and changes its mean
energy by the negative of that amount. This is an explicit apparatus energy
account, without an autonomous zero-work detector or passive clock claim.
The selected source evolution parameter is the same in the analytic and
portable implementations; conversion to seconds requires the supplied unit
dictionary and does not derive a physical observer clock.

The portable kernel executes all 48 spectral coordinates at every retained
momentum and all 96 native probe rows. It computes spectral phases, the
native transform and the complete position transform in Taichi float64.
The independent NumPy calculation acts on the complete same sector and
checks the resulting vectors, without probability clipping or
renormalization. The 16, 32 and 64-site probes retain every momentum in
one periodic direction, with the other two momenta zero. The circumference
is 8pi/3 electron inverse-mass units. Electron, positron, muon and tau
preparations retain every zero-probability coordinate and the preparation
failure. Reversing the evolution sign supplies a failed propagation control.
The wider three-dimensional statement below is analytic; a three-dimensional
interacting Fock simulation is outside these executed probes.

## Regulator speed, unique minimum and refinement

For every retained three-momentum,

    |grad E_a|^2 = sum_i sin(a p_i)^2/(a^2 E_a^2) <= L_a/E_a^2 < 1.

Thus the group speed is bounded by one. The only minimum of L_a on the
Brillouin torus is p=0, giving one minimum per massive branch. This is a
positive-energy lattice-Laplacian construction; it does not remove fermion
doubling from an asserted first-order chiral source.

The elementary sine remainder gives, for |p|<=P,

    0 <= |p|^2-L_a(p) <= a^2 sum_i p_i^4/12 <= a^2 P^4/12,
    0 <= sqrt(m^2+|p|^2)-E_a(p) <= a^2 P^4/(24m).

On the same finite Fourier band, Duhamel and instrument contractivity bound
any complete same-input record law by

    TV <= min(1,s |t| a^2 P^4/(24m)),

where s is the retained preparation success probability. The comparison
uses the identical native probe and source reference. Its continuum
comparator is the periodic, band-limited positive-energy Hamiltonian with
the same Fourier amplitudes. Finite position probes are its orthonormal
cardinal-mode apparatus, not continuum pointwise position measurements.
The displayed rational upper bounds use pi<22/7 and exact analyzer masses.
Floating kernel discrepancies are separately reported; the analytic
bounds do not treat a vendor's trigonometric output as exact arithmetic.

## Boosted mass-shell and velocity checks

Let E=sqrt(m^2+|p|^2) and let (E_b,p_b) be its exact continuum Lorentz boost.
For both momenta inside the admitted Brillouin domain, write
Delta(p)=a^2|p|^4/(24m). The same-source dispersion obeys

    |E_a(p_b)-gamma(E_a(p)-beta.p)| <= Delta(p_b)+gamma Delta(p).

For the executed electron p=3/4, the continuum pair (E,p)=(5/4,3/4)
boosts to rest at beta=3/5 and gamma=5/4. The energy-shell defect and
velocity converge quadratically with spacing. In one dimension,

    |v_a-p/E| <= a^2[P^3/(6m)+P^5/(24m^3)].

Both inequalities retain the finite preferred frame. They do not define
a unitary Lorentz representation on the finite momentum grid, prove
length contraction of a localized interacting object, or supply interacting
Lorentz covariance. The rest and moving plane-wave preparations are
reported separately from the finite packets; they are not substituted for
a derived boosted-packet apparatus. No transverse photon source or universal
matter/light speed identification is asserted.
