# Autonomous native charge-pair records with retained flux loops

The constrained Green minima do not by themselves define an autonomous
particle Hamiltonian. The unit-Yukawa interacting source sends a simple
one-ec state into 26 onsite output terms. Its selected particle/string
subspace is therefore not invariant. No pinned energy is substituted into
that source as a complete equation of motion.

This construction declares a finite autonomous specialization on the full
N=4 cube carrier: 30 fermionic modes and a five-state scalar at each of
64 cells, all 119 states on each of 192 native links, and the full
30-mode reference. The source is

    H=2 N_matter+N_scalar+4 N_reference
      + sum_{x=0}^3 [a_x c_x^* U_x c_{x+1}+adjoint]
      +(kappa/2) rho_Y^T G rho_Y.

The four sites are (0,0,z), with periodic z. U is the literal ec matrix
coordinate of the compressed Peter-Weyl construction, not a fitted scalar
or a replacement hopping amplitude. Both native spin entries of sigma_z
are retained. The coefficients are the actual nearest z terms of the
periodic native kinetic matrix: three i sqrt(2)/8 entries and the boundary
-i sqrt(2)/8. Other path/species hopping coefficients and Yukawa vertices
are declared zero. The number mass coefficient is two. G is the
three-dimensional cube Green operator with its zero mode removed, acting
on the actual total matter/scalar/reference hypercharge at every cell.
The unit nonblank-link penalty is absent; there is no transverse-loop
energy. These choices are an explicit Hamiltonian specialization, without
a transfer of an interacting ground state, particle spectrum or clock.

The full definition is block diagonal in reference occupations, with the
appropriate native reference current and number energy in every block.
The selected input has all reference modes filled except the spin-up ec
mode. Its number is 29 and its native hypercharge is -1. It is a weak and
colour singlet. One active ec particle has hypercharge +1 and either spin.
Other matter modes are vacuum, all scalars are blank, and off-ring links
are blank. Each of these is a specified preparation, not a quotient of
the carrier. Other preparations and references are retained in the complete
measurement complement.

## Exact invariant physical sector

The ec multiplication coordinate and its adjoint preserve exactly the
three link states with native labels 0, 46 and 105. Their outgoing
hypercharges are 0,-1,+1, and every weak/colour endpoint current vanishes.
In this ordered basis the actual coordinate is

    U = [[0,0,1],[1,0,0],[0,0,0]].

Its 116 outside rows vanish on these columns, also for U*. This is checked
against the full 119-state coordinate. A native hop preserves one ec
particle and its spin, shifts the active link according to this matrix,
and preserves every native Gauss equation. For a particle at y and ring
fluxes f_z the physical condition is

    delta_{z,y}-delta_{z,0}+f_z-f_{z-1}=0.

Enumerating all fluxes in {-1,0,+1} gives nine physical position/flux states
for each spin, hence an eighteen-dimensional invariant sector. This
retains all allowed circulation patterns in that preparation sector.
The full native column action has no output outside it. Every one of its
864 active-site Gauss equations is checked with actual matter and endpoint
matrices; other cells contain only singlet factors. The unused source
sectors remain in the full carrier and are invariant under this declared
Hamiltonian. This argument does not discard their possible different
spectra or assign them the selected preparation probability.

## Autonomous current and force

The clock is exp(-itH), with no position pinning and no intermediate
measurement. Let J_z be the derivative of H with respect to an inserted
phase on the actual link hopping. The full matrices satisfy

    d N_z/dt = i[H,N_z] = J_{z-1}-J_z,
    (d J_z/dt)_electric = i[V,J_z],
    V=(kappa/2)rho_Y^T G rho_Y.

The first identity is the native charge continuity equation. The second
is a term in the same autonomous Ehrenfest equation, not a force installed
from a classical formula. On the ring the Green energies are
0, 21/128, 31/160, 21/128 at unit kappa, independent of circulation.
They equal the associated pinned-minimum differences, while the unpinned
spectrum and motion are computed from the complete eighteen-row H.

The explicit coherent input on connected positions one and two has kinetic
coherence with the negative-energy sign. Across that bond the electric
energy increases by 19/640. The exact electric contribution to the outward
current acceleration is

    - (sqrt(2)/8)(19/640) = -19 sqrt(2)/5120.

It is directed toward the lower-energy separation. This is an exact
finite quantum current statement. It does not imply a universal classical
trajectory or monotone arrival probabilities as the coupling changes.
Coherent interference and the finite flux cutoff affect those probabilities.

Complete position/spin/flux probabilities are reported at eight times for
zero, unit and fourfold electric coupling. They sum to one; the full
complement has zero probability on the declared ideal input. Both spin
branches remain in the matrices. Mean energy includes the reference
contribution 116 and particle number contribution two. Independent
matrix-exponential and Taylor propagation agree, with an explicit Taylor
tail. Native phase differentiation and finite-time current differentiation
check the exact commutator identities. Removing circulation states leaves
a noninvariant four-state path subspace and changes the arrival law.

The construction supplies autonomous quantum motion and an electric
current-acceleration response on a declared finite gauge source. It does
not derive physical alpha, GeV units, a hydrogen spectrum, a far-field
continuum force, transverse radiation, or autonomous motion of the
unchanged interacting Peter-Weyl construction.
