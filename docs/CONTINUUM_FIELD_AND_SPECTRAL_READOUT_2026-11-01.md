# Continuum fields, complete ground-family spectra and energy readout

This extension is followed by the
[native detector and centered-ground construction](NATIVE_DETECTOR_AND_CENTERED_GROUND_2026-11-01.md),
which adds a complete native phase clock and resolves the charge-clock
and reference choice left open here.

The original Weyl/contact family now has a spatial continuum construction:
all eight valleys, bounded smeared CAR fields, their complete multi-time
correlations, continuous space and time translations, positive energy and
finite propagation speed. This advances beyond the previous
[fixed-spacing field construction](NATIVE_GROUND_AND_SPATIAL_RESPONSE_2026-11-01.md).
It is a common subsequential limit under an explicit bounded-contact
refinement. It does not yet prove that a nonzero interaction survives that
limit or supply the physical Standard Model and Einstein identifications.

## The continuum construction

The [continuum field net](../model/continuum_field_net/README.md) uses
`c(a)=v/a`, `v>0`, and uniformly bounded contact coupling. At each spacing
it uses the original infinite-volume ground, full CAR algebra and full
interacting dynamics. Cell averages embed physical smearing functions;
the eight alternating lattice phases retain every Weyl valley.

A new full-operator estimate bounds the contact commutator of a smeared
field by `U||f||`. Together with the cell-average kinetic bound, it gives
uniform continuity in physical time. One spacing subsequence then carries
every field correlation. Positivity and the CAR survive, and the lattice
locality estimate becomes exact graded locality outside a conservative
physical cone. The cone speed `96v` is a proof bound, not a measured light
speed or a derived Lorentz metric.

On fixed physical momentum windows, the inverse retarded response also
has a bound uniform in spacing:

\[
\|G_a(z)^{-1}\|\le |z|+vP+U_{\max}/2+
\frac{U_{\max}^2}{4\operatorname{Im}z}.
\]

This removes the earlier global ultraviolet norm divergence on these
windows. It does not remove the positive imaginary-frequency restriction.

There is a controlled special case: if the contact strength tends to
zero, the full-CAR dynamics converge to the eight free massless Weyl
fields. This statement now uses the positive-energy ground representation,
rather than the earlier empty-Fock finite-particle comparison. It gives a
free massless control, not a proof of the required interacting massive and
massless particle theory.

## The complete finite ground family

The [spin-resolved spectrum](../model/spin_resolved_spectrum/README.md)
now covers all five vectors in the exact spin-two ground space and all
sixteen linear creator modes. Every such preparation has addition support
within the same 78 positive frequencies. The line `sqrt(3)/16` has total
hidden spin `5/2`; the other 77 lines have total hidden spin `3/2`.

The reachable `N=5` space has dimension 1256 inside the complete
4368-dimensional sector. The remaining 3112 states are dark for every
single-creator probe from this ground family and its time translates.
These are total internal-spin labels in a spin-two background, not
identified Lorentz spins of particles. The dark-sector energies are not
claimed to have been computed.

For the invariant mixed ground, each fixed-spin local addition has total
weight `3/4`. After conditioning on addition, the `sqrt(3)/16` line has
weight `2/5`, and every other line has `2/3` of its previous highest-ground
down-addition weight. All remain positive. Removal has positive energy
cost `sqrt(3)/16`; it appears at the negative of that frequency in the
retarded convention. The full mixed-reference retarded response has
79 signed poles.

The invariant mixed reference has conditional mean addition energy
`(1+sqrt(3)/16)/3`, approximately `0.369417725158` model units, and
variance `67/288-sqrt(3)/72`. Removal has the deterministic positive cost
`sqrt(3)/16`, approximately `0.108253175473`. These are moments and costs
of the specified excitation distributions, not relativistic rest masses.
The earlier highest-weight down-addition probe has mean energy `1/2`;
changing the preparation changes its spectral weights, while the allowed
energy lines remain the same.

Consequently ground polarization, coherent superposition and linear probe
choice cannot rescue the rejected finite-gap charged-lepton dictionary.
They do not introduce a new energy beyond the previously certified
envelope. This statement still leaves nonlinear probes, other reference
states and continuum particles open.

## A complete energy measurement with a resolution bound

The [spectral readout](../model/spectral_readout/README.md) constructs a
positive, complete finite-resolution POVM directly from the spectral
measure. Its Fejer filters use finitely many characteristic-function
samples and retain every clock outcome. The proof separates resolution,
aliasing, native-control errors, statistical errors and the still-unknown
finite-to-infinite ground-state discrepancy.

For the existing 78-line finite source at time step one, the calculated
probability of a readout error of at least `1/2` model energy unit is:

| Clock outcomes | Probability, rounded |
| ---: | ---: |
| 8 | 0.100284517428 |
| 16 | 0.055193158232 |
| 64 | 0.020209031434 |

At 16 outcomes, the leading readout centers are zero and `3pi/8`, with
probabilities about `0.483266606026` and `0.358720295560`. These are filter
bins, not asserted particle energies. In particular a bin centered at
zero does not establish a massless particle. The higher-time samples are
source calculations; the previously executed native experiment supplies
the first harmonic only.

Exact infinite-ground sum rules also give a local signed-energy mean in
`[1/8,1/2]`, second moment in `[19/128,67/128]` and variance in
`[17/128,35/128]` at `c=1/8,U=1`. The signs label addition and removal
response frequencies, not negative physical excitation energies above
the ground. Compressed retarded resolvents converge strongly on the
common one-particle space at fixed spacing; the varying GNS generators
are not silently identified with one operator on one Hilbert space.

## Ground constraints and the remaining interaction problem

The [ground constraints](../model/continuum_ground_constraints/README.md)
prove that the density tends to one fermion per lattice site when `U/c`
tends to zero. A variational bound controls the occupation defect and
every fixed microscopic local word. Under `c=v/a` with bounded
`U(a)->U0`, each fixed smooth physical packet has limiting first energy
moment equal to its Weyl kinetic moment plus `U0/2` times its squared
norm. The scalar contribution is retained without an energy shift.

This does not prove a Gaussian macroscopic limit. A normalized continuum
packet can concentrate in a shrinking momentum region even when its
integrated lattice occupation defect tends to zero. Conversely, a
nonzero microscopic memory Gram can escape to arbitrarily high energies
while the finite-frequency limit becomes free. Explicit positive-measure
controls exhibit both possibilities. Determining which happens in the
actual nonzero-coupling source requires further uniform estimates.

The remaining physical requirements are therefore precise: surviving
continuum interaction, Lorentz and gauge structure, the appropriate
massive and massless spectral sectors, a common detector and unit map,
and the same-state modular, stress, area, curvature and coarse-graining
premises of the existing Einstein theorem stack. None is replaced by a
chosen finite-bin energy or an assumed geometric interpretation.

## Validation of this extension

The four component packages pass 83 strict tests and fresh-process write
and check replays. Independent proof and source reviews cover each
component. The integrated state passes 499 integration tests and 22
subtests, plus 136 supporting tests. All 1,414 model artifact hashes,
generated claim/premise/completion views and postdiction-ledger parents
were checked.

The map records 71 realized scoped attachments among 76 reviewed domains.
The custody audit preserves all 84 canonical premise dispositions, 346
claim classifications and existing completion-gate statuses. It adds no
physical comparison or accepted physical gate. All eight papers were
rebuilt and passed publication warning checks; visual review covered 64
changed or adjacent pages. The extension adds no new Lean theorem:
the mathematical arguments, exact witnesses and executable checks have
their distinct roles stated in the component records.
