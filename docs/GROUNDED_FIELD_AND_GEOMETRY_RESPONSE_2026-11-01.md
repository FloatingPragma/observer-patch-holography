# An interacting ground, energy density and matter response

The subsequent [native ground and spatial-response extension](NATIVE_GROUND_AND_SPATIAL_RESPONSE_2026-11-01.md)
adds finite ground preparation and charged readout, the full matrix
spectrum and spatial inverse-response control on the infinite ground.
The earlier results below retain their original proof and validation scope.

The original three-dimensional Weyl/contact action has a proved
infinite-volume ground representation with positive excitation energy.
The construction keeps its original zero chemical potential and all
particle-number sectors. Its ground energy density has a convergent
finite-volume calculation. On the certified finite cube, the same action
also gives an exact matter-dependent propagation response, which connects
the particle calculation to the existing geometry theorem interfaces.

These results extend the [finite ground and relativistic mass calculation](INTERACTING_VACUUM_AND_RELATIVISTIC_MASS_2026-11-01.md).
They close particular ground and response premises. They do not supply
the controlled interacting chiral 3+1-dimensional continuum or a match to
measured particle species.

## The original infinite-lattice ground

The [ground construction](../model/thermodynamic_ground/README.md) uses the
complete ground projector of each finite open box, extends it by an
exterior vacuum, averages its spatial translates, and takes a weak-star
subsequence. The local dynamics estimates proved for the full CAR
algebra control the boundary terms, including odd charged probes.

The limiting state is translation invariant and gauge invariant. Its GNS
time-translation generator is positive, so its spectral energies are
excitation energies above a ground. The contact interaction remains a
nonzero operator in this representation. This last statement establishes
an interacting operator law; it does not by itself establish nontrivial
scattering or a connected four-point function in the limiting state.

The density satisfies

\[
\frac14\le\omega(n_{0\uparrow}+n_{0\downarrow})\le1.
\]

No half-filled reference was inserted. A product of exact four-particle
cube trial states supplies the lower density bound; the full particle-hole
energy identity supplies the upper bound. Spatial averaging before taking
the limit is essential to retain both bounds.

The theorem constructs a state by compactness. It supplies neither a
selected ground functional nor a native preparation algorithm. Invariance
under the continuing repair and controller process remains a separate
requirement.

## A convergent energy-density readout

Let `E_ell` be the absolute ground energy of an open cube of even side
`ell`, at fixed `c>0,U>0`. The thermodynamic energy density exists and obeys

\[
\frac{E_\ell}{\ell^3}-\frac{3c}{\ell}
\le e_\infty\le\frac{E_\ell}{\ell^3}.
\]

The surface width `3c/ell` is an explicit error bound. It follows by counting
all intertile bonds and bounding their operator norms; the upper bound uses
even fixed-number trial states whose intertile hopping expectations vanish.
The constructed ground has local energy expectation `omega(h_0)=e_infinity`.
A certified finite-box energy interval passes through the same bounds with
its input uncertainty retained. No convergence rate for the ground state
is inferred from this energy-density estimate.

At the original `c=1/8,U=1`, the existing full cube ground energy is
`E_2=-sqrt(3)/4`. Combining its trial bound with the global kinetic and
particle-hole bounds gives

\[
-\frac{\sqrt3}{8}\le e_\infty\le-\frac{\sqrt3}{32},
\]

or approximately `-0.216506350946110 <= e_infinity <= -0.0541265877365274`
in model energy per site. This is an energy-density bracket, not a particle
rest mass. It can be narrowed by certified ground calculations on larger
boxes.

## Matter response on the same action

The [retarded response construction](../model/retarded_matter_response/README.md)
uses the exact ground `G` of the original cube and the normalized stationary
charged state `sqrt(2) c_0up^dagger G`. Its addition and removal spaces are
invariant sectors of the full action. Their exact resolvent retains every
intermediate state.

For the down component above `G`, the first two response moments give

\[
M_1=h+\tfrac12 I,\qquad M_2-M_1^2=\tfrac14 I.
\]

The second identity is a nonzero dynamical memory contribution, beyond a
static potential. Adding the localized up excitation changes the down
component's first moment by

\[
\operatorname{diag}(1/2,1/6,1/6,0,1/6,0,0,0).
\]

The up-component retarded response in these references remains exactly
`(z-h)^(-1)`. These statements follow from the full source action and CAR
identities, rather than from a fitted effective Hamiltonian. In the original
field frame, the response is the corresponding onsite unitary conjugate.

The Schur resolvent expresses the down response through an exact
frequency-dependent memory kernel. The certified charged gap controls a
low-frequency analytic expansion on the cube. Separately, the new
thermodynamic ground subsequence supplies local retarded-response
convergence on compact subsets of the upper half-plane. This convergence
does not license exchanging the infinite-volume limit with a real-axis
pole limit.

The local down-addition spectrum above the finite ground is fully
resolved: **78 distinct positive energy lines, with zero unresolved probe
weight**. Exact polynomial annihilation and a matching cyclic-rank witness
establish the complete support. Integer interval arithmetic isolates every
root and encloses every residue to 24 decimal places. The
[complete result](../model/retarded_matter_response/results.json) retains all
78 intervals and their normalized weights.

Some of the strongest lines are:

| Energy above `G`, in model units | Normalized probe weight |
| ---: | ---: |
| 0.0299509383421 | 0.313011377216 |
| 0.108253175473 | exactly 1/10 |
| 0.136046750488 | 0.111967695928 |
| 1.08963895147 | 0.148656969335 |
| 1.10825317547 | exactly 1/10 |

The earlier isolated symmetry contributions at `sqrt(3)/16` and
`1+sqrt(3)/16` are proved to be their entire respective line weights.
The other 76 lines carry exactly `4/5` in total. The highest visible energy
is approximately `1.83607446544239`. The lowest visible line is a property
of this local probe; it is not asserted to be the minimum of every state
in the full charged sector. These are computed source energies, with no
new native preparation, finite-shot detector execution or relativistic
rest-mass interpretation.

## Using the existing geometry theorem stack

The relevant existing interface is
[`geometric_matter_join.source_tetrad`](../model/observer_particles/geometric_matter_join.py),
which consumes a spatial principal tensor and the action-density
normalization through `intrinsic_geometry.wave_density_metric`. The
subsequent symbol and stress constructions establish their geometric
conclusions when those inputs agree across the matter sectors.

The response above supplies an actual matter-dependent effective operator
on the unchanged microscopic graph. A new fundamental metric register is
not required for this route. The remaining construction must control its
local low-energy spatial symbol and show that the relevant sectors share
the resulting metric and normalization. A Hartree potential or a temporal
Taylor expansion alone is insufficient for that identification.

The Einstein stack then requires the same state's regional modular,
stress, entropy and coarse-graining inputs. Its theorems remain available;
the work is to realize their joint premises on this ground and its
observables. No separate gravity model or assumed Einstein response has
been substituted for that attachment.

## Particle masses and physical identification

The earlier calculated finite cube excitation energies and the massive
one-dimensional continuation keep their stated domains. The latter's
anchored rest energy `exp(-pi/8)/pi` and rest mass `64 exp(-pi/8)/pi` are
properties of a chosen continuation; other continuations through the same
bare point give different masses. The current ground theorem does not
select that continuation or convert those values into masses of observed
species.

The outstanding particle calculation requires charged spectral measures
over a selected thermodynamic ground, controlled poles and residues under
the actual spatial refinement, and common spin, charge, clock and detector
identifications. Those are the remaining inputs to a physical comparison.
The ledgers retain the distinction between a finite energy line, a
relativistic rest mass on a specified limit, and a measured particle.

## Validation and ledger status

The two packages pass 76 focused tests and fresh receipt replays. The
integration suite passes 492 tests and 22 subtests; the checks covering the
reported CI failures pass another 123 tests. Independent review preserves
all 84 canonical premise dispositions, 346 claim classifications, prior
attachment rows and frozen comparison contracts. The finite map has
64 realized scoped attachments among 69 reviewed domains, with no new
physical identification or integrated completion gate accepted.

All eight paper PDFs build with no overfull boxes, unexplained underfull
boxes or reference, citation, glyph or font errors. Their 48 affected pages
were rendered and visually inspected. Claim views, source hashes, the
surface inventory, completion views and postdiction links are synchronized.
The selected GitHub issues remain open where their full acceptance
requirements have not been met.
