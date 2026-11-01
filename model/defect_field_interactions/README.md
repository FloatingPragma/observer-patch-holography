# Complete particle-changing Higgs drive

This is a specified external drive on the same spatial CAR model used by
the static quasiparticle and propagation constructions: twelve ports,
ninety internal modes per port, and **all** `2^1080` occupation states.
Its finite implementation is a tensor product of complete two-mode
source gates. No global statevector of that size is materialized.

```sh
python3 -m model.defect_field_interactions.packet --check
python3 -m pytest -q -W error model/defect_field_interactions
```

The actual retained Higgs invariant tensors are evaluated at
`H=(0,v exp(i phi))`. Their `conjugate_Higgs` flags matter: the creation
pairing coefficient has phase `−phi` in the up channel and `+phi` in the
down and charged-lepton channels. Uniformly multiplying all of Delta by
one phase would be a different drive. The gauge background is held fixed;
this is an external time-dependent field, not a passive gauge change and
not a calculation of a quantized Higgs or gauge interaction.

The complete local Nambu matrix has the form
`B_phi=cos(phi)B_0+sin(phi)B_1`, with `{B_0,B_1}=0`. A canonical change only
inside degenerate energy/charge spaces decomposes this actual matrix into
42 pairs of opposite-charge equal-energy quasiparticles. Every original
energy/charge projector is unchanged. The corresponding complete CAR
Hamiltonian of each pair is

\[
H_{ab}(\phi)=E\cos\phi(n_a+n_b-1)
             +E\sin\phi(a^\dagger b^\dagger+ba).
\]

`bare_component_chart((27,78))` additionally executes the four-state
canonical map for the actual electron component. It intertwines both
the static mass and the phase-driven family; an extra commuting
quasiparticle phase fixes the driven off-diagonal convention. This map
allows a source-executed propagation factor written in the original bare
CAR coordinates to be composed with the present quasiparticle gates.
Identifying those two coordinate systems without this map would be wrong.

All four occupations `00,10,01,11` are retained. The odd sector stays
unchanged, while the even sector mixes vacuum and a neutral pair. The
six massless modes at every port are identity spectators. Thus there are
504 driven blocks and 72 spectator modes. These are a decomposition of
the full domain, not a projection onto a chosen occupation sector.

`pair_word` compiles the exact two-level rotation using the existing
source scalar and Pauli controls. Every distinct frozen calibrated
energy, every native multiplier `1+j/4`, and all four input columns are
executed and checked. A separate two-bit occupation pointer is executed
on all sixteen system/pointer input columns. The saved arrays retain
every output amplitude. Proper-code service failure channels retain the
existing exact-zero law; floating execution checks are not hardware or
interval certificates. Initial state preparation, storage and physical
readout costs are not silently included in the pulse counts.

The modes of one pair need not be adjacent in the fixed global CAR order.
`routed_pair_word` moves them together using source-compiled **fermionic**
swaps and undoes the route afterwards. These swaps put a minus sign on
double occupancy, preserving the parity signs of intervening modes. An
independent four-mode Jordan–Wigner test covers every spectator occupation
and rejects an ordinary qubit-swap replacement. Routing costs and their
contribution to the conservative global numerical bound are retained.

At `phi=0` this is exactly the previous onsite static Hamiltonian. The
spatial reference remains explicit: twelve copies contain `12 E_vac`,
and the same fixed subtraction `11 E_vac I` leaves the earlier single
global scalar. That scalar is also source-executed. Restricting the static
limit to the hard-core two-particle subspace gives the prior `U⊗U` law.
At nonzero phase the pair term leaves that subspace. It conserves electric
charge and fermion parity while changing quasiparticle number by two.

For a vacuum pair, the complete pointer gives pair probability
`p=sin²(E tau_eff) sin²(phi)`. Returning the external background to its
initial static value supplies expected excitation energy `2 E p`. This
is switching work in the specified source action, not a rest-mass
derivation from the number of compiler pulses. The full spatial example
starts with two particles on distinct spectator ports and retains every
even-number outcome through `N=1006`; small positive tails are stored as
decimal values rather than rounded to zero by floating underflow.

The [`full_domain`](../defect_field_bridge/full_domain.py) repair extension
now makes the next repair step defined on these pair-created collisions.
The package executes the actual complete source repair isometry between
a drive and its inverse, keeping all 32 coin/raw records and all 16
local occupation states. A collision flag distinguishes vacuum from a
same-port pair. Consequently the inverse-drive return probability is
`1−2p(1−p)`, compared with one without the repair. Two occupied ports
outside the seam can be tensored as spectators, turning this into an
`N=2`/`N=4` example. The collision hold is therefore explicitly a
measurement disturbance, not a claim that repair leaves the interaction
coherent.

There is also an exact full-spatial echo calculation at the declared
interval `tau=0.002`. Before repair each initially empty seam endpoint has
even occupation. The noncollision record therefore means both endpoints
are still vacuum, with probability `A` obtained from all 84 local pair
factors. The other record has probability `1−A`. The full inverse-drive
return is `A²+(1−A)²`; all other sites and their particle sectors remain
present. This follows from two orthogonal measurement projectors and does
not assume that their mixed output remains Gaussian.

Propagation, this local drive, and repair may be composed as an explicit
ordered model law on the shared spatial CAR space. The common clock
advances once for the entire physical interval. Multiplying every term
by the same `1+j/4` does not mean advancing the clock once per compiled
term. A chosen split-step product is not silently equated with continuous
simultaneous evolution. Interacting physical poles, scattering data and
quantized gauge/Higgs dynamics remain further scientific work.
