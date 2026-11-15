# Compact modes and static one-axis flux controls

`flux_controls.run()` reconstructs finite bare-operator identities without a
worker checkout or result-file dependency. The carrier and canonical ring torus
come from the portable labelled geometry. Exact integer nullspace reduction
finds two carrier eigenvalue-6 vectors vanishing at the ring ports. Extending a
vector by zero to the 2-by-2 torus preserves its Laplacian eigen-equation.
The bare mass-zero port/transit evolution is independently exponentiated and
compared with its two-dimensional breathing formula. A joint port/transit
eigenvector is stationary up to phase; the ports-only vector changes its transit
weight. These states do not identify a stable particle in the integer-record
lift or monitored instrument.

The supplied static flux is class 3 on a reference face and its antipode,
represented by five class-3 seam phases. The gauged carrier Laplacian commutes
with its chosen +120-degree parallel-transport rotation R. Its six doublets
have ground value 0.8775137836. Exact signed-permutation arithmetic gives
R cubed = -I, while the equally commuting rephased operator -R has cube +I.
Removing the gauge phases from the rotation rejects commutation. Replacing the
antipodal flux by the two adjacent face classes induced by one class-3 seam
splits the ground level by 2.059363009. The sign depends on the transport
convention; doublet transfer requires the flux-placement hypotheses.

These finite observables give neither a field-carried angular momentum density
nor a three-dimensional electron rotation representation. The Z6 algebra has
five phase sectors with fermion sign, all at odd flux class. It supplies no
exchange dynamics, flux binding or consistent full-generation fermion dictionary.
The measurement is a static algebra control with supplied phases and axis;
rows 4, 6 and 7 retain their physical-state and observation requirements.

Implementation: [flux_controls.py](../oph_openwave/flux_controls.py).
Tests: [test_flux_controls.py](../tests/test_flux_controls.py).
All local checks and rejecting controls return explicit booleans;
`full_criterion_met` is false. Workspace audit provenance is kept separately
from this portable module.
