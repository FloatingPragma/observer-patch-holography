# Spin, charge and an executed collective mode of the same field

The complete three-dimensional Weyl/contact action has an exact onsite
reformulation as a spin-independent **pi-flux Hubbard model**. Its commuting
spin and eta SU(2) algebras constrain the entire Fock carrier, not only the
earlier dimer candidates. An exact certificate checks the original action on
all 65,536 states of the open cube, including every number sector and every
local creation/removal Ward identity. A separate source experiment prepares,
evolves and reads its collective eta-pair mode on all three coordinate
plaquettes, with every export failure and terminal outcome retained.

This is a connection to established Hubbard physics, not a new abstract
symmetry theorem: see [Yang's eta-pair construction](https://doi.org/10.1103/PhysRevLett.63.2144)
and [Yang and Zhang's SO(4) symmetry](https://doi.org/10.1142/S0217984990000933).
The contribution is the exact map from the existing OPH field action, its
complete representation and source interfaces, and their independent checks.
It advances the finite components of CP-03 through CP-06 without accepting
their physical gates or assigning a Standard Model species.

## Exact map, generators and particle/hole consequences

On any finite open rectangular box, write the actual hopping block as
`h_(x,x+e_j) = (c/2) U_j`, with `U_j=i sigma_j`, and retain the onsite
interaction `U n_up n_down`. Set

```
R_x = (i sigma_0)^x0 (i sigma_1)^x1 (i sigma_2)^x2,   d_x=R_x c_x.
R_x U_j R_(x+e_j)^* = (-1)^(sum_(k>j) x_k) I.
```

Each R has determinant one, so the onsite contact and pair creator are
unchanged. The induced exterior-power unitary acts on every occupation
sector, including the empty and full states. The scalar hopping signs are
`((-1)^(y+z), (-1)^z, 1)`. Their product around every elementary plaquette
is -1, matching the original Wilson holonomy `U_i U_j U_i* U_j*=-I`.
The finite residue certificate checks all 64 coordinate classes modulo four;
the periodic formula establishes the link identity for all integer sites.

In this frame define

```
S_a = (1/2) sum_x d_x* sigma_a d_x,
eta+ = sum_x (-1)^(x0+x1+x2) c_x,up* c_x,down*,
eta- = (eta+)*,   eta_z=(N-V)/2,   K=H-U*N/2.
```

The two SU(2) algebras commute. `[H,S_a]=0`, `[H,eta+]=U eta+`, and
`[K,eta+]=0`. In the original chart the spin matrix at site x is
`(-1)^(sum_(j!=a) x_j) sigma_a/2`; a spatially uniform bare spin rotation
is generally not this symmetry. The CAR creation/removal operators form a
spin/eta `(1/2,1/2)` tensor, checked by exact local commutators on the full
carrier. Consequently each K eigenspace decomposes into spin/eta multiplets
and has exact `Q=N-V` versus `-Q` partners, at the same centered energy.
Probe weights agree only when both reference and probe are transformed as
well. Neither this symmetry nor a dimer mass calibration fixes their energies.

The kinematical catalogue is complete for this exhibited symmetry action.
For `a=2m_S`, `b=2m_eta`, the weight multiplicity is

```
w(a,b) = binom(V,(V+b+a)/2) binom(V,(V+b-a)/2),
m(s,e) = w(2s,2e)-w(2s+2,2e)-w(2s,2e+2)+w(2s+2,2e+2).
```

Noninteger or out-of-range binomial indices give zero. Highest-weight
subtraction proves this formula and the dimension sum is exactly `4^V`.
The receipt retains every multiplet for V=1,2,4,8. It is not an energy
classification, a list of stable particles, or the full commutant of H.
Special couplings and geometries can have additional symmetries.

The centers satisfy `exp(2pi i S_z)=(-1)^N` and
`exp(2pi i eta_z)=(-1)^(N-V)`. Their product is `(-1)^V I`.
Even-site boxes carry a linear SO(4) state representation; odd-site boxes
descend only projectively. Conjugation on observables descends in both cases.
This internal symmetry is not a derivation of physical Lorentz spin.

## An exact collective solution and complete source history

For arbitrary pair angle theta and phase phi,

```
|theta,phi> = product_x [cos(theta/2)|00>_x
                  + (-1)^|x| exp(i phi) sin(theta/2)|11>_x].
```

This is a normalized superposition of the V+1 eta-tower states
`(eta+)^m|vac>`, which have exact H energies `m U` and K energy zero.
Every complete spatial hopping bond commutes with eta+ and annihilates
vacuum, so it annihilates the whole tower. H evolution therefore sends
`phi -> phi-U t`; K evolution is stationary. These are excited collective
states, not a selected ground or a faithful thermal reference. For an
eta-singlet reference, both eta ladders annihilate its support and this
collective channel has zero residue. The zero K frequency is not a
dispersing massless particle pole.

The source executes theta=pi/2, phi=0 on a complete four-site plaquette:
eight blank field bits, twenty-eight existing preparation pulses, the unchanged
routed full-field evolution word, its explicit optional number-phase word
for K, and the twenty-eight-pulse inverse preparation as analyzer. Each plane
uses the actual Weyl matrices on its two axes. The H experiment uses **2,216
native pulses**, counting 4,432 proper-code load/export service events,
8,160 first-failure Kraus outputs, and all 256 terminal addresses. K uses
sixteen additional number-phase pulses. Initialization, coherent services,
destructive final readout and their physical calibration remain supplied.

Each final site has probability `sin^2(U t/2)` of being doubly occupied;
singly occupied terminal sites are ideally absent, but all numerical outcomes remain
in the receipt. At U=t=1 the probability is approximately **0.2298488471**,
with exact ideal lower bound `529/2304` and upper bound `1/4`. The K and U=0
controls give zero. Raw occupation measurements before the analyzer do not
see the phase. The pair preparation and analyzer are even but change number
by two: the complete onsite macro is a Y rotation conjugated by CNOT, and
its full four-state matrix preserves parity. Individual encoded primitive
pulses need not preserve field parity. The macros explicitly supply a coherent
charge-phase apparatus. This does
not select a laboratory clock or break a physical superselection rule by
native repair.

The existing split schedule keeps each spatial bond's two disjoint commuting
mode hops consecutive. Each whole bond fixes the tower, and each complete
onsite layer advances its common phase. Thus the ideal split error is exactly
zero **on this tower**, for every split count. The general-input Strang bound
is reported separately and is not set to zero. Wrong staggered preparation,
an incomplete bond, or a charge-preserving local potential need not obey this
identity. Tests compare the literal source output with a separate complete
Hamiltonian exponential, including negative time and refinement.

For each distinct pulse the actual six-dimensional response, including its
export complement, is checked against independent analytic primitive columns
with numerical norm ceiling `delta=10^-11`. For L pulses, telescoping gives
the declared bound `L delta/(1-L delta)`, allowing numerical norm growth.
The complete returned complex vector and all first-failure outputs must
satisfy this budget and probability accounting; no success conditioning is
used. These are numerical admission gates, not outward certificates for
floating arithmetic or hardware. Retained raw residuals are observations;
the fixed ceilings are the acceptance bounds. Replaying executes all sources
again and checks every receipt field.

## Topology, locality and the remaining particle map

The -I plaquette is a fixed holonomy of the supplied hopping background. It
cannot be flattened to identity by onsite basis changes, but it is not a
state-dependent localized defect charge. On a connected box an onsite pair
creator `sum_x f_x c_x,up* c_x,down*` commutes with hopping precisely when
`f_y=-f_x` on every edge (at c>0). Hence f is proportional to the global
staggered profile: a pair confined to a proper region is not protected by
this identity. Independent tests check that kernel and its local perturbation
control. Arbitrary SU(2) links preserve eta on a bipartite graph; central
holonomy supplies the additional hidden spin symmetry here. Odd periodic
identifications can break the open-box generators, and a deliberately twisted
boundary is a separate action.

The supplied Hamiltonian symmetries do not override the existing native C/E
source-only conservation obstruction. Autonomous localized protection,
same-excitation native maintenance energy, physical clock normalization,
electromagnetic/Lorentz identification, a chiral limit and the measured
spectrum remain required. The earlier dimer family exclusion and failed
finite dispersion test stand. This construction supplies a genuine common
field multiplet and collective-response interface; it supplies no new
independently tested physical particle or topological defect.

```
python -m model.field_symmetry --check model/field_symmetry/results.json
python -m pytest -q -W error model/field_symmetry
```
