# Two-baryon test on the signed colour vacuum

This is a complementary audit of the frozen four-site source in
`hadrons/CHAIN4_MANIFEST.json`. It is not a physical nuclear calculation:
this source has no identified nucleon spin, isospin or three-dimensional
scattering continuum. The test asks whether its two units of baryon number
even lie below twice its one-unit excitation in the same finite box.

The audit reconstructs the CAR operators without importing the producer's
matrices. Instead of taking the nullspace of the colour Casimir, it imposes
zero Cartan weight and both simple SU(3) raising currents. In a finite
unitary representation, a zero-weight highest-weight vector belongs to the
trivial representation. This gives the complete singlet dimensions
1, 20, 50, 20, 1 at quark numbers 0, 3, 6, 9, 12. Every Hamiltonian component
preserves this kernel. Its nonorthogonal Gram matrix is retained.

The source is signed nearest-neighbour hopping with t=1/8, staggered mass
equal to the inherited raw charged-register coefficient, and the sum of
prefix colour Casimirs with the three predeclared coefficients 1/2, 1, 2.
Independent rational LDL inertia checks the complete ground brackets.
The computational mass rounding is accounted for by the same explicit
operator bound 12 times its error, checked independently.

For each member, the actual vacuum has N=6. Write its sector ground energy
as E6, and the others as E9 and E12. The one- and two-baryon-number gaps are
E9−E6 and E12−E6. The two-excitation excess is therefore

    D = E12 + E6 − 2 E9.

Its strictly positive rational interval in `NUCLEAR_THRESHOLD_REVIEW.json`
rules out a positive binding energy against this same-box reference.
Subtracting an arbitrary common energy zero cannot change D. The identical
statement holds for the antibaryon sectors N=3 and N=0. This does not prove
an infinite-volume absence of bound states: the N=12 state saturates this
small box, and independent-particle scattering thresholds require their
own volume limit. Nor does the result identify a deuteron or a halo.

There is a useful exact interpretation of the previously derived
second-order virtual hopping term. Let n_x in {0,1} denote a local colour
singlet's occupancy. For C=9t²/(4κ)>0,

    −C 1_(n_x != n_y) = −C(n_x+n_y) + 2C n_x n_y.

Its density interaction is repulsive; the negative total virtual energy
shift includes a one-body contribution. A negative second-order energy
shift alone is consequently insufficient evidence for nuclear attraction.
This algebra concerns the perturbative coefficient. The finite threshold
test uses the full Hamiltonian and does not omit the higher orders.

The next substantive alternatives require actual additional source
structure: spin and flavour currents, a derived exchange interaction, and
an outgoing relative-coordinate channel. Installing a target-selected
binding term or assigning a lifetime to a discrete finite spectrum would
not repair the missing nuclear construction.
