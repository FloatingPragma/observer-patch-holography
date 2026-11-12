# Bound colour singlets from native currents

## Source and retained carrier

At each vertex of a finite graph there are all fifteen native internal
coordinates and two spin coordinates, with the full fermionic Fock space.
For one family on two vertices this is sixty fermionic modes. The vacuum is
empty. An idle apparatus can carry the opposite electric charge or conjugate
colour representation; its energy is declared zero. These references retain
all outcomes and do not modify the matter Hamiltonian. The preparation and
readout times are external services. Between them evolution is autonomous.

Write T_x^a for the second quantizations of the eight actual native colour
currents at vertex x. The native Cartan generators in the source file are
not orthogonal: the orthonormal pair is T6 and (T6+2T7)/sqrt(3).
Thus the fundamental Casimir is 4/3. The projector onto native colourful
coordinates is (3/4) times that internal Casimir. S_x is the spin-half current
on this entire colourful subspace. It acts on both quarks and conjugate
quarks, and vanishes on the three colourless native coordinates.

The declared action is

```
H = sum_f m_f N_f
    + t sum_edges,internal,spin,f (c_x-c_y)^*(c_x-c_y)
    + U sum_x,a (T_x^a)^2 + J sum_x,i (S_x^i)^2,
    m_f>0, t,U,J>=0.                                      (1)
```

All terms preserve global native colour, weak and hypercharge transformations.
The action has no local gauge links or local Gauss constraint. The positive
graph kinetic term and number-energy term are declared choices, distinct
from the homogeneous Weyl/Higgs action. Coupling constants are architectural
inputs. For the finite execution m=1, t=1/8, U=1, J=1/4, with one edge.
The full Fock inequality H>=m_min N follows term by term. It proves the empty
matter vacuum and controls every uncomputed occupation sector. It does not
prove a continuum particle interpretation or identify an internal spin with
spacetime angular momentum.

## Exact sector reduction

Occupation numbers of each native quark/conjugate species are conserved.
The three-doublet-quark sector has 24 modes and dimension choose(24,3)=2024.
Its complete colour-singlet subspace is Sym^3(C^8), of dimension120. An
orthogonal basis consists of epsilon-colour contractions of three orbital
creators. The orbital indices include both vertices, both spin states and
both weak components. The normalization squared is 6 times the product of
orbital multiplicity factorials. The exterior Cauchy decomposition gives
120 singlets, 168 adjoints and 56 decuplets; 120+168*8+56*10=2024.

The quark/conjugate-quark sector has dimension24*24=576. Its complete colour
singlet has dimension8*8=64. The conjugate colour contraction is solved from
the actual native current matrices, including the exterior-basis signs.
The complement consists of64 adjoints. These are exact reducing sectors
of (1), not a deletion of their complements. The implementation verifies
B*B=I, T_total B=0 and H B=B h directly in each full fixed-occupation sector.
Their two reduced matrices have sizes120 and64. This finite representation
shortcut replaces diagonalization on a space of dimension2^60.

## Binding inequalities

On a graph with maximum degree z the lowest one-quark excitation is bounded
below by m+4U/3+3J/4. On the two-vertex graph the uniform spatial mode attains
this value. Any two same-sign quarks have total local colour Casimir at
least4/3: two quarks on one vertex transform as bar3 or6; separate vertices
contribute8/3. Their kinetic and spin costs are nonnegative. Every two-plus-one
fragment channel therefore has threshold at least

```
3m + 8U/3 + 3J/4.                                      (2)
```

A local colour-singlet three-quark state with internal spin1/2 and weak
isospin1/2 has trial energy at most3m+3J/4+3zt. The positive spin interaction
preserves this whole subspace. The variational principle places its lowest
eigenvalue below every two-plus-one threshold by at least8U/3-3zt, when this
quantity is positive. Three separate quarks have a higher lower threshold.
For z=1 the certified gap is55/24. Charge-zero and charge-one members are
exactly degenerate by global weak flavour symmetry. Spin3/2, including
charges-1 and2, is retained as an additional multiplet.

A local quark/conjugate-quark colour singlet with internal spin0 has trial
energy at most m_a+m_b+2zt. Its two-fragment threshold is
m_a+m_b+8U/3+3J/2. Its binding gap is therefore at least8U/3+3J/2-2zt,
equal to67/24 for the executed parameters. The colourless meson sector
contains charges-1,0,1, including both neutral basis combinations and spin1
partners. No neutral singlet is silently removed. Turning off U and J gives
free graph particles with zero binding below the fragment threshold.

These bounds are relative to the same source and the same mass and reference
conventions. A fragment threshold is the sum of excitation energies of
noninteracting copies of that finite source, restricted to its exactly
conserved native particle occupation sectors. A two-plus-one channel uses
all276 two-quark and all24 one-quark states. A meson channel uses all24
quark and all24 conjugate-quark states. This is a finite spectral binding
comparison, not an asymptotic separation or scattering theorem on the
two-vertex graph. These bounds establish finite bound eigenstates. They neither assert
that isolated colour charges are asymptotic QCD particles nor establish
physical proton, pion, kaon or neutron masses and lifetimes.

## Complete preparation and records

For each stated charge and internal spin, diagonalize the full reduced
Hamiltonian and form the complete lowest-energy-band projector P. Choose a
basis seed deterministically and use the instrument {P,I-P}. Both outcomes
are retained. The successful state is P|seed>/sqrt(p); the failure probability
is1-p. The original occupation instrument has every2024 or576 outcomes,
including zeros, plus the explicit preparation failure. The receipt records
the complete successful occupation distribution after time7/10, its charge
expectations and spatial co-location probability. Outside the reducing
sector the instrument is the identity failure channel. Idle references are
retained. As the success branch is an exact energy eigenstate, its record
probabilities are invariant under autonomous time evolution.
The receipt also gives the unconditional probabilities of every
success-and-occupation record, multiplied by the preparation success
probability, and the complementary preparation failure. Those records sum
to one. Conditional successful distributions are labelled separately.

The spectral preparation service is architectural. This result supplies no
autonomous hadron formation, baryon-number violation, strange decay, unequal
up/down mass mechanism, or neutron positive-core/negative-shell distribution.
Equal masses and the exact weak symmetry force the proton-like and
neutron-like energies to coincide. Physical identification needs observables
that distinguish these members beyond their native charges and binding.
