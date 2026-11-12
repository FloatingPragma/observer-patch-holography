# Separate volume and spacing controls for the finite Gauss atom

The Hamiltonian is exactly the Dirichlet Green construction in
[openwave_field_observables](../openwave_field_observables/PROOF.md),
with nonrelativistic mass one and declared coupling squared 4 pi. Its
potential is obtained by inverting the same three-dimensional incidence
Laplacian, rather than inserting samples of -1/r. The reflection-sector
maps retain their original normalization, including the fixed reflection
plane. The ground and two second-shell sectors are computed separately.

## A domain comparison without a second field derivation

Take two centred, nested odd cubes at the same spacing, with the point
charge at their common centre. The larger cube's Green function restricted
to the smaller cube differs from its Green function by a discrete harmonic
function with nonnegative boundary values. The finite maximum principle
therefore gives G_large >= G_small pointwise on the smaller cube. The
attractive potential obeys V_large <= V_small.

Zero extension of a smaller-cube wavefunction preserves its kinetic
quadratic form, including boundary edges. It weakly decreases the
potential quadratic form. The min-max principle therefore makes every
ordered eigenvalue nonincreasing under this fixed-spacing domain
extension. The same comparison holds in each reflection sector because
zero extension intertwines reflections. It does not imply monotonicity
of energy gaps, and it does not compare different spacings.

The executed volume sequence uses spacing 1/3 and side lengths 24, 32, 40.
The independent spacing sequence holds side length 24 and uses grids
71, 95, 127, 159. Their common endpoint is evaluated once. Thus the effect
of increasing the box is distinguished from the effect of refining its
stencil. A separate small nested-cube calculation directly checks the
Green comparison before any spectral fitting.

## Finite numerical controls and their limits

For a normalized trial vector u and reported E, the residual norm
||Hu-Eu|| bounds distance from E to the spectrum of that finite numerical
matrix. It does not alone identify which eigenvalue lies nearby. For a
real symmetric Hamiltonian with nonpositive off-diagonal entries and a
connected stencil, its ground eigenvector is strictly positive. For any
strictly positive trial vector v, the exact Collatz inequalities give

    min_i (Hv)_i/v_i <= lambda_ground <= max_i (Hv)_i/v_i.

Subtract H from a sufficiently large multiple of the identity and apply
the elementary positive-matrix comparison to prove this statement. The
receipt prints the corresponding floating diagnostics and every residual.
These calculations do not claim outward-rounded evaluation of the exact
Green inverse or continuum eigenvalue enclosures.

The reference hydrogenic gap 3/8 supplies a comparison for the declared
mass and charge convention. Differences between successive grid results
are convergence diagnostics, not proved continuum uncertainty intervals.
The regulator, boundary, supplied mass/coupling and physical hydrogen
identification retain their separate roles.
