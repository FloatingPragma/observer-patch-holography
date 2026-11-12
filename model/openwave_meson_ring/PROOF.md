# Native meson emission on an outgoing ring

## Source, reference and clock

The source is the complete two-site native meson and semileptonic current of
`model/openwave_meson_weak_current`, with its colour-singlet reduction and
native family word `(2,3,10)`. Its one-quark/one-antiquark preparation, spectator
antiquark weak component, three positive family number energies, strong
couplings and contact coefficient `G=1/64` are unchanged. The outgoing charged
and neutral lepton coordinates occupy an even ring of `L` sites. Contacts are
only at the original hadron sites 0 and 1. All original 180 fermionic modes
remain in the carrier; the extension has `180+12(L-2)` modes. The source includes
all other occupation sectors and spectators. The clock parameter is the same
source Hamiltonian time convention, with no identification as seconds.

Write `K_L=q laplacian(C_L)/4`, counting the two opposite edges separately at
`L=2`. The charged free operator for family `f` is `m_f I+K_L`, and the neutral
operator is `-K_L`. All `6L` neutral modes are filled in the reference, and the
reference energy `-3qL` is subtracted. Thus a neutral hole has positive kinetic
energy `K_L`; the charged mode and hole both have dispersion
`q(1-cos(k))/2`. The three charged energies are the same native ratios
`(1,1/m,M/m)` as the source meson construction. The chosen box
`q=(1/m-2)/2-11G/4`, approximately 102.342189236, is inherited from the two-site
contact-resonance architecture and is held fixed through the volume test.
This choice is a declared resonance selection, not a physical strange-meson
calibration.

Every current is the same local native weak generator, Pauli current and
family matrix as in the two-site source. The quark source has weak coordinate
`w` in `{0,1}`. Conservation of `Q=-w-N_e=-1` proves `N_e<=1` in this preparation.
The local current also preserves the combined charged-plus-neutral occupation
of each lepton family, so the emitted hole and electron belong to the same
family. These are source invariants; no particle-pair cutoff is installed.
All failed preparations and orthogonal sectors are retained by the instrument.

## Exact reduction and injection strength

The colour-singlet, charge-minus-one, zero-spin-component sector has dimension
`24+216L^2`. Spin zero and the reflection exchanging sites 0 and 1 reduce its
prepared invariant subspace to dimension `6+36L^2`. The embedding includes both
four-spin singlets and the full emitted position, family and hadron sectors.
The program verifies the full intertwiner `H B=B h`, not only a projected
Hamiltonian. The first six columns of `B` are the no-emission source sector.
This reduction extends to arbitrary quantum references by tensoring with the
identity.

In that decomposition the complete restricted Hamiltonian is

    h = [[A, G C*], [G C, D]].

Here `D` contains the native neutral-current interaction as well as every
emitted free and strong term. The injection obeys the exact identity
`C* C=24 I_6`. One local Pauli contraction is

    T_(t,b,a;s) = sum_mu eta_mu bar_sigma_mu[t,s] bar_sigma_mu[b,a].

The exact two-spin identity is `sum_(t,b,a) conj(T_(t,b,a;s))
T_(t,b,a;s')=8 delta_(s,s')`. The three orthogonal emitted lepton families give
24; unitarity of the native quark family word removes its family rotation.
Emissions at distinct source sites have orthogonal position records and the
one-quark number operator sums to one. Spectator spin, antiquark position and
colour do not alter the identity. Consequently every normalized parent state
has the same short-time emission law

    P_pair(t)=24 G^2 t^2+O(t^3)=(3/512)t^2+O(t^3).

This is an exact volume-independent coefficient of the full source evolution.
It is neither an exponential rate nor a numerical fit.

For an initially empty emission sector, eliminating its amplitude gives the
exact six-component Volterra equation

    da/dt=-i A a-G^2 integral_0^t C* exp[-iD(t-s)] C a(s) ds,
    b(t)=-iG integral_0^t exp[-iD(t-s)] C a(s) ds.

The receipt checks its resolvent Schur complement against the complete
source inverse for `L=2,4`. No free-reservoir replacement is made in this
identity. A finite ring has a discrete spectrum and supports recurrences.

## Volume control

The infinite-ring source acts on the finite hadron space tensored with the
zero/one electron-hole sector over the integers, in the filled-reference hole
representation. The kinetic shifts, finite mass operators and current at the
two contacts are bounded, so their sum defines a bounded self-adjoint source
on this invariant sector. This statement does not assert equivalence with an
infinite-volume vacuum representation of every unprepared Fock sector.

Split the source into `D_0+T`, where `D_0` contains the strong source, on-site
energies and contact current and `T` contains only outgoing nearest-neighbour
hops. `D_0` can create/annihilate a pair only at sites 0 and 1; it cannot move
an existing outgoing coordinate. Each coordinate has off-diagonal hopping
norm at most `q/2`, so `||T||<=q<103`. Embed the even ring around the two
contacts in the integer line. A word capable of reaching its differing wrap
edge requires at least `R=(L-2)/2` hops. Interaction-picture Dyson words of
length less than `R` agree on the whole no-emission preparation domain.
Unitary factors generated by `D_0` have norm one, including the contact
creation and annihilation. Thus for any such preparation and quantum reference,

    ||(U_L(t)-U_infinity(t))P_parent||
      <=2 sum_(k>=R) (103|t|)^k/k!
      <=2 (103|t|)^R/[R! (1-103|t|/(R+1))]

when `103|t|<R+1`. The comparison uses the natural embedding of localized
outgoing states, with the outside-region outcome retained. Complete record
probabilities of common embedded/local instruments differ by at most twice
this vector bound. The outside-line complement is retained. This statement
does not identify the finite-ring Fourier instrument with an infinite-line
momentum spectral instrument. The same proof covers
the failed-preparation parent branch because it has no outgoing pair.
For `L=32,t=1/100` the exact rational upper bound is below
`2.546754e-12`. This is a certified source-volume error; it does not include
floating matrix-exponential roundoff. The receipt labels the independently
computed 36,870-state propagation as a floating diagnostic. An analytic
`L=1024,t=1` bound is supplied without claiming an execution of that carrier.

## Complete records and work

The native orbital/spin seed has preparation success probability approximately
0.248344602311. The failed branch is a retained quantum state, with its own
orthogonal flag. On success, the exact source is evolved before either a
position instrument or a separate momentum instrument is applied. The
receipt gives every occupation probability in the complete charge/spin
sector at time 720. Each reduced colour-singlet probability splits into three
equal native matching-colour records. Other colour and source sectors remain
in the instrument with probability zero for this preparation. All recorded
probabilities are multiplied by the original success probability; failure
probability is not redistributed.

The momentum instrument is the unitary two-coordinate discrete Fourier
transform in the canonical hole basis; the hole has its own conjugate
momentum coordinate, and its dispersion is the same even cosine. It is an alternative to the position
instrument, not a joint sharp position/momentum record. Momentum energies use
the same source dispersion and the exact contact diagonal averaged over its
position-basis support. Both instruments resolve the native hadron colour,
so their measured strong energy uses the full-colour Hamiltonian diagonal,
including the colour-dephasing cost. The receipt computes
`W=success*(Tr(H rho_after)-Tr(H rho_before))` for each instrument from the
full strong, kinetic and contact source. The initial nonselective spectral
preparation has zero mean work: it commutes with the parent strong Hamiltonian
and the contact is `-3 I` on the complete no-emission sector. Source evolution
conserves energy. The failed branch is kept without an additional fine readout.

## Finite-volume observations and physical scope

At source time720, pair probabilities for `L=2,4,6,8` are approximately
`0.99490795,0.24404062,0.06322543,0.02351293`. The selected two-plane-wave
bound-meson matrix element has exact magnitude `2|C_01|/(sqrt(2)L)`.
The parent contact mean is -3 and the target contact mean is `3-1/L`, so the
fixed box detunes the target by `G(1/2-1/L)`. A declared recentering control
`q_L=(1/m-2)/2-G(6-1/L)/2` equalizes these two means; it neither diagonalizes
the full interacting source nor supplies a physical linewidth.

The ring calculation derives a stable inclusive short-time coefficient and
an exact memory kernel, and rejects treating the two-site conversion peak as
a volume-stable irreversible lifetime. Physical strange-meson identity,
measured decay lifetime and an infinite-time exponential law are not claimed.
