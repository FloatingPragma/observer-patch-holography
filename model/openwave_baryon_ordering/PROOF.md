# Native baryon sector and signed intrinsic records

## Source and exact reduction

The carrier is the two-cell native15-coordinate, two-spin CAR algebra.
The displayed baryon has three particles in the actual colourful weak-doublet
module. Other number and native-coordinate sectors remain reducing complements.
The Hamiltonian is the sum of number energies, positive graph differences,
local colour Casimirs with coefficient1 and local colourful-spin and weak
Casimirs with coefficients1/4 each. A weak square means the actual doublet
current, not a particle-name assignment. Every term commutes with total colour,
total spin, total weak isospin and electric charge. These statements follow
from the Casimir commutators and scalar hopping; unequal number energies retain
charge and spin and are scalar in each fixed flavour-count sector.

The integer colour-epsilon embedding R has120 mutually orthogonal columns in
the2024-dimensional three-particle occupation space. Its Gram matrix is diagonal
with entries6 times the product of orbital multiplicity factorials. Native
CAR arithmetic gives integer Hnum=48H. The execution checks exactly

```
R^T R=diag(g),
L Hnum R=R diag(L/g) (R^T Hnum R),
```

where L is the least common multiple of the Gram entries. Thus no discarded
colour-singlet vector affects the reduced source. The normalized exact matrix
has entries (R^T Hnum R)_ij/(48 sqrt(g_i g_j)). Positivity controls every remaining
Fock sector; this certificate is a specified colour-singlet observation and
does not assert the global full-Fock ground is a baryon.

## All-eigenvalue certificate and exact preparation

Outward interval arithmetic encloses each algebraic matrix entry. The operator
error between H and its dyadic representative is at most dimension times the
largest entry error. A numerical eigensolver proposes a dyadic frame O and
dyadic energies E; neither proposal is trusted as a certificate. Exact integer
products bound eta=norm(O^T O-I) and rho=norm(HO-OE). Both Frobenius bounds imply
their operator bounds. For eta<1/2 the polar unitary U obeys norm(O-U)<=eta.
A verified K bounds H and E. Therefore

```
norm(HU-UE)<=rho+2 K eta<=2 rho+4 K eta.
```

Weyl's inequality encloses all120 ordered eigenvalues in the stored intervals.
The lowest four form an isolated band. Three spin halves and three weak-doublet
halves allow only half-integer irreducibles, each of dimension at least two.
Colour-singlet spin/isospin representations therefore have dimension at least
four. A four-dimensional isolated band is exactly one (1/2,1/2) representation;
Schur's lemma makes its Hamiltonian scalar. At fixed charge0 or1 and spin_z1/2
its ground line is unique.

A rational trial is rounded only to propose integer epsilon coefficients.
Charge and spin_z selection and the exact spin projector15I-4S^2 retain the
specified sector. The stored original occupation vector v has exact norm N,
Rayleigh quotient E and certified residual rho_v. If gamma is the distance
from E to the certified next band, the unhalved trace distance to the unique
ground line is at most2 rho_v/gamma. The exact rank-one projector vv^T/N acts
on a normalized epsilon-column seed. Its success probability is the stored
rational squared overlap, and failure is its complement. All2024 unconditional
occupation probabilities are stored as rationals and sum to success. The
failed quantum state and idle opposite-charge apparatus reference are retained.
At arbitrary source time t, the complete record trace distance from the t=0
trial records is at most min(2,2 abs(t) rho_v), by Duhamel's bound. No exact
stationarity is attributed to a rounded eigenvector.

## Intrinsic charge, mass ordering and obstruction

For each original outcome use number centre X=(x1+x2+x3)/3 and actual quark
charges2/3,-1/3. The observable at radius r sums charges at abs(x-X)=r. Its
operator norm is at most2, so the exact ground expectation differs from the
rational trial expectation by at most twice the trace bound. The neutral
state's certified signs are positive at r=1/3 and negative at r=2/3; its
coincident charge at r=0 is zero. This is a translation-invariant two-cell
graph record. It is not an isotropic physical nucleon form factor or a radius
in femtometres.

Let a,b be the light up/down free vectorlike gaps divided by the unchanged
electron gap. Their exact source inequalities give0<a<b<1. In the fixed
three-particle charge sectors the number-energy change is scalar:

```
E_p=E_equal+2(a-1)+(b-1),
E_n=E_equal+(a-1)+2(b-1),
0<E_n-E_p=b-a<1.
```

This proves the ordering and simultaneously excludes the ground-state channel
n to p+electron+antineutrino, even with a massless neutrino. No choice of the
positive vectorlike mixing parameters in this common-family construction
removes that obstruction. Multiplication by the inherited electron unit gives
conditional graph energies in GeV, whose nucleon identification fails.
Rational analyzer intervals, integer square-root enclosures and exact interval
arithmetic propagate the free mass ratios through every displayed energy
interval; decimal endpoints round outward. The empirical electron anchor is
held at its declared central value. Its experimental uncertainty and physical
model error are separate from these mathematical source intervals.
The added equal spin/weak weights and the free-gap number-energy prescription
are declared architecture. Their consequences do not supply an interacting
projection theorem for the free Weyl/Higgs source or a physical beta-decay law.
