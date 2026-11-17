# Weak-current transfer without empirical hadron coefficients

The source is the hadron worker's declared current-square amendment on the
complete colour-singlet three-quark carrier. This package does not reuse
the calibrated strong source or any of its state vectors. The electron is
one and the neutral is massless in this new source's declared units. The
common 432-state lepton mass-slot construction is not identified with it.

## Ground states and actual current intervals

For each charge assignment, the 14-dimensional Sz=1/2 Hamiltonian is rebuilt
from its radical entries and the new diagonal current square. Vectorlike
number energies shift a fixed charge sector by a scalar, so they do not
change its eigenvectors. A separate exact rational congruence calculates
enclosures for the lowest and second eigenvalues. The positive gap between
those enclosures proves that the first state is unique in this projection.

Let q be the rational approximation to a normalized trial vector, lambda
the midpoint of the ground interval, and d a lower bound on the separation
of lambda from every other eigenvalue. Interval arithmetic on the actual
radical matrix gives r >= ||(H-lambda)q||. The normalized q has spectral
weight outside the ground state at most r/(0.999d). With its phase chosen
to have positive overlap, the distance to the true ground vector is at most

    epsilon = 2r/(0.999d) + |qᵀq-1|.

`amended_current.py` checks 0.999<||q||<1.001 before using this bound.
It obtains distances below 1.3e-8. For any quark current J of norm at most
three, the difference between the rational-vector bilinear and the exact
ground-state bilinear is at most

    3(epsilon_proton + 1.001 epsilon_neutron).

This is added to the exact radical interval for q_protonᵀ J q_neutron.
The resulting certified total-current intervals are

    fV in [0.9997181305, 0.9997182368],
    fA in [1.6659572990, 1.6659574053],
    fA/fV in [1.6664268370, 1.6664271202],

with the exact rational endpoints in `AMENDED_CURRENT.json`. Printed
decimals here are outward-rounded descriptions, not the certificate itself.

Every Hamiltonian term is a spin scalar. On three spin-one-half quarks,
the possible total spins are one-half and three-halves. The unique Sz=1/2
ground has a certified raising norm below 0.1, excluding the sqrt(3)
raising norm of spin three-halves. It therefore has spin one-half.
Cell reflection commutes with the same Hamiltonian. Its certified
deviation from even parity is below 0.1; uniqueness excludes the odd
eigenvalue. Both ground states are cell even.

Spin covariance and equal cell parity imply the local current identities

    <p,s|J_x^0|n,up> = (fV/2) delta_(s,up),
    <p,s|J_x^j|n,up> = (fA/2) (sigma_j)_(s,up).

The actual proton spin-down partner is constructed with the total spin
lowering operator. All eight local current blocks are checked independently
against these identities, with maximum numerical residual 1.47e-13.

## Complete three-dimensional weak source

Declare one electron/neutral Dirac family, masses one and zero, nearest
neighbour kinetic p=2 sin(theta), and Wilson coefficient one. At the
neutral zero mode choose the negative-beta reference subspace explicitly;
the active covariance is I/2 and all zero-energy complements are retained.
The scalar neutral mass is zero, not an imported oscillation value.

Use the same two contact positions and the full native semileptonic
cross-current, with G=1/64. The source-independent 64-state local norm
theorem is rerun from its exact matrix. Its characteristic polynomial is

    x^28 (x-1)^17 (x+1)^3 (x-2)^3 (x-3) (x+3)^3
    (x+6) (x²+2x-9)^4.

Its norm is six. For one lepton family and three quarks, the complete
weak operator has norm at most 18. The two contact positions are orthogonal
blocks of each quark's position, so they do not add an extra factor two.
Projection onto the colour-singlet carrier cannot increase the norm.
This bound retains every lepton Fock sector.

The proton-channel spin polynomial in `PROOF.md` is multiplied by fV²,
with g=fA/fV. Its two-cell form factor is unchanged. The interval bounds
prove fV>0.999 and 1.66<g<1.67. Its spin weight is below 16, so the same
selected first-order norm bound 4Gt applies. At t=1/100000,

    delta = (18Gt)²/2,
    probability error <= (8Gt+delta)delta.

The staggered-momentum pairing argument of `PROOF.md` supplies positive
directional sums. The exact electron-negative margin is

    (Gt)² (0.999)² (1.66)(0.66) 2/(7L)
    - 16(Gt)²(14t)²/12 - (8Gt+delta)delta.

The corresponding antineutrino margin replaces 0.66 with 2.66. Both are
strictly positive at L=4 and L=8. Dividing by 120 keeps the actual uniform
hadron-input preparation and its 119/120 failure. Exact rational margins
and all twelve finite cube/control records are in `AMENDED_BETA.json`.

The control labelled zero chirality is an equal classical mixture of left
and right contacts, with its fair coin retained. It is not a coherent sum
of the contacts. The unchanged w=0 and equal-mixture controls have zero
velocity contrast; the right-current control reverses the left-current
contrast.

## An open energy shell, with its limit stated

The source gap has the exact interval

    1.11317360314155 < Q < 1.11317360317066.

Thus Q exceeds the electron-plus-massless-neutral threshold. In the
continuous momentum family, take electron momentum zero and neutral
momentum (0,q,0). Its energy is

    E_nu(q) = sqrt(4 sin²(q) + (1-cos(q))²).

On (0,pi/2) this is continuous and strictly increasing from zero to sqrt(5),
so exactly one q has E_nu(q)=Q-1. Its derivative is nonzero. The contact
form factor is one and its spin weight is (fV²+3fA²)/2>0. The implicit
function and coarea formulas then give positive free transition spectral
density on a neighborhood of this regular five-dimensional energy shell.

This statement belongs to the continuous momentum extension. The executed
small periodic cubes are not asserted to contain an exactly on-shell pair.
Their short-time conversions include off-shell weight. Neither the positive
free density nor the strict short-time sign is a theorem about the full
interacting survival law, irreversible release, or a neutron lifetime.

## Remaining physical and implementation obligations

No empirical hadron mass, axial value or neutrino mixing parameter enters
this amended beta candidate. Its predicted Q and gA still disagree with
their physical targets. The hadron worker has executed the added native
diagonal, but the new ground filter and full weak source have not been
compiled and executed natively. The complete Fock evolution is bounded
analytically at short time, not simulated in its enormous full dimension.
The common lepton/colour clock, physical masses, long-time limit and
physical lifetime remain open. These distinctions prevent promotion of
this constructive subproblem repair into a completed physical beta row.
