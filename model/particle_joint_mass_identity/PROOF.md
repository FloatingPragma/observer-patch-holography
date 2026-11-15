# Joint quantum numbers and rest coefficients

## Two unchanged sources

The calibrated construction uses exactly the 180 physical CAR modes,
360-dimensional native Nambu symbol, 56 canonical boson modes and global
boson cap two of
[conditional_particle_spectrum](../conditional_particle_spectrum/PROOF.md).
The full source Hilbert dimension is `2^180 × 1653`. Its preparation is
the native Bogoliubov vacuum times the capped oscillator vacuum. The
dual-charge and dual-colour quantum references are retained. All Yukawas,
the running Higgs background, gauge and scalar parameters, and optional
Weinberg inputs retain their empirical ancestry.

The compiled lepton calculation separately uses the neutral free source,
filled reference, canonical CAR transformation and reservoir SWAP of
[lepton_mass_compiler](../lepton_mass_compiler/PROOF.md). Its coefficients
come from the fixed native analyzer, followed by one electron energy
anchor. These two sources and their mass conventions are never merged.

## Actual commuting fermion operators

On the original two-momentum CAR coordinates use the native electric
charge Q, native colour Casimir C, spin S_z and momentum P_z. On the
Nambu doubling their matrices are

```
Q_N=diag(Q_CAR,−Q_CAR^T),
C_N=diag(C_CAR,C_CAR^T),
S_N=diag(S_CAR,−S_CAR^T),
P_N=diag(P_CAR,−P_CAR^T).
```

The source momenta are `(0,0,±1 GeV)`, so helicity is `S_N P_N` in the
declared units. The executable constructs these operators in the actual
360-dimensional source basis, verifies their mutual commutators and
their commutators with the native Nambu Hamiltonian, and checks every
named positive-energy vector. Both charge conjugates and both momenta
remain. The negative Nambu partner completes the canonical transform;
it is not another copy of the particle catalogue.

The exact internal identities are `Q^T K+KQ=0` and `[C,K]=0`. The exact
family Takagi frame is reconstructed by the parent producer, rather
than inferred from nearly degenerate numeric eigenvectors. In each
analytic mass pair the block `[[p,m],[m,−p]]` squares to
`(p²+m²)I`. This gives the exact energy attached to its projector even
when the difference from one GeV is below binary64 resolution. The
record explicitly separates this exact attachment from the numerical
commutator diagnostics.

Group the canonical positive modes by the **measured** tuple
`(Q,C,2 helicity,P_z,E)`, without putting a species name in that tuple.
This gives 84 joint eigenspaces whose ranks sum to 180. Quark colour
components share rank-three projectors; they are not separately
measured. Names in `candidate_species` describe the supplied comparison
dictionary after grouping, and do not prove a physical identification.

## Boson polarizations and colour

The existing gauge-Higgs mass matrix and scalar Hessian supply the
frequencies. A circular change of polarization basis diagonalizes the
actual spin-about-z matrix on each transverse pair. Massive vector
blocks also retain their longitudinal polarization with helicity zero.
The radial scalar has helicity zero. Photon and gluon blocks retain
only their two transverse modes; no longitudinal massless oscillator
is introduced.

The eight equal-frequency colour gauge coordinates are the declared
canonical adjoint fiber of the parent free action. The executable
constructs the eight normalized SU3 adjoint matrices from the fundamental
Gell-Mann commutators and checks their quadratic Casimir `3I8` exactly.
This specifies an orthonormal basis within that existing degenerate
fiber, rather than measuring eight classical colour labels. It is not
a derivation of a confined physical gluon state or a conversion from
an unnormalized microscopic gauge-current basis.

The polarization unitary, electric charge, colour Casimir, helicity,
momentum and the parent Hamiltonian are checked on all 56 modes. Their
joint projectors form 28 groups with total rank 56. Gluon groups retain
rank eight; the other groups have rank one. The source's canonical
free polarization and gauge structure are declared construction data.

## Complete joint instrument on the full source

There are 112 joint one-excitation projectors: 84 fermionic and 28 bosonic.
Add the source vacuum projector and the **whole** orthogonal
many-excitation complement. Their ranks are respectively 236 in total,
one, and `2^180×1653−237`. The 114 orthogonal projectors sum to identity
on the complete source. On the many-excitation complement the apparatus
records that outcome; it does not pretend to resolve its complete
interacting or many-body colour/energy spectrum.

Use a fair classical choice between the parent's complete fermionic
and bosonic probes. The ungrouped fermion Kraus operators are

```
K_F+ = sum_j a_j† tensor |dual j>/sqrt(360),
K_F− = sum_j a_j  tensor |j>/sqrt(360).
```

The ungrouped boson operators use the same formula with denominator
sqrt(116), together with the original capped no-operation effect divided
by two. The reference states are orthogonal mode coordinates and retain
their charge and colour representation. The probe records only five
pointers: fermion creation/removal, boson creation/removal and boson
no-operation. None contains a species or family name.

The full CAR anticommutator identity supplies half of the total identity.
The parent's capped-boson identity supplies the other half on all 1,653
boson states. Composing with the 114 joint projectors therefore produces
a complete instrument on every source input and every untouched
reference. A direct sum of the parent's reference spaces has dimension
474, including its unused blank fermion coordinate.

On the stated vacuum preparation the joint probabilities are rank/360
for fermionic one-excitation groups, rank/116 for bosonic groups and
1/58 for the bosonic no-operation vacuum. All other pointer/projector
records are zero and are retained. Thus 570 original records have 113
positive and 457 zero probabilities and sum to one. Actual capped
creator columns reproduce the bosonic weights; the full native CAR
frame supplies the fermionic creators. With every fermion mode occupied,
all outputs instead lie in the retained many-excitation complement.
Dropping that complement would fail the all-input instrument identity.

Within a rank-r colour group the successful source/reference output is
a coherent maximally entangled vector, not a classical colour mixture.
Dephasing its colour coordinates changes its conditional density by
unhalved trace norm `2(1−1/r)`. Multiplying by every original group
probability and adding the orthogonal output blocks gives the negative
control recorded in the receipt. Consequently colour labels cannot be
silently measured or discarded while claiming the same full law.

## What distinguishes a species, and at what resolution

Without energy, charge, colour, helicity and momentum leave the three
generation masses degenerate. They also fail to separate some neutral
vector/scalar targets in helicity-zero sectors. Energy resolves those
distinctions in the exact declared source. Colour components remain
degenerate by design. Neutral opposite-helicity Majorana readouts do
not create an extra sterile species or an independent antiparticle copy.

Every mass squared in this calibrated packet is a nonnegative rational
number. Integer-square-root enclosures therefore give rigorous rational
intervals for `m` and `sqrt(1+m²)` without rounding tiny neutrino masses
to zero. These intervals bound arithmetic error in the fixed input
packet, not empirical or physical-model uncertainty.

The executable also merges the exact projectors into fixed energy bins.
It checks that every rational energy enclosure lies wholly in its bin.
With half-bin width `10^−6 GeV`, the optional neutrino mass eigenfields
remain indistinguishable by energy at momentum one GeV. With the much
smaller declared half-bin width `10^−26 GeV`, their exact source energies
occupy distinct bins. No physical detector accuracy, source service
compilation or finite experimental exposure at that resolution is
claimed. A roundoff-sensitive numeric diagonalizer is not used to
certify the distinction.

## Compiled lepton spin refinement

On the literal sixteen-state Dirac CAR block, retain the parent sea and
reservoir SWAP. The dimensionless Hamiltonian is
`H=b† beta b+2I`, the native electric charge is `Q=2−N_b`, and
`S_z=sum_j ((−1)^j/2)b_j†b_j`. These commute exactly. Multiplying the
parent energy projectors by the diagonal charge/spin projectors gives
fifteen nonzero joint effects with ranks summing to sixteen. Their
colour Casimir is zero and momentum is zero; helicity at rest is not
defined or reported.

Both original reservoir settings and both reservoir outputs are retained.
Each family therefore has sixty spin/charge/energy records. Forgetting
spin reproduces all thirty-six original records exactly. The probe can
fail and leave a neutral pair at energy twice the same mass; it is not
renormalized to successful one-particle events.

The same parent finite-compiler Ramsey bound and positive denominator
guard give conditional GeV intervals for these three free coefficients.
Only the electron value fixes the empirical unit. The phase and family
dictionary retain the parent's retrospective ancestry, and anchor
uncertainty, finite-shot error, source-to-laboratory time and interacting
pole matching remain outside the finite construction. The calibrated
17-target table has its own supplied Yukawa and gauge inputs and is
reported separately.
