# Three-dimensional native beta angular records

The source uses the complete 120-dimensional reducing colour-singlet baryon
sector of `openwave_axial_state`, its literal local weak currents and its
four empirical calibrations. Quark spectators retain their reducing vacuum
sectors. The lepton carrier is a periodic cubic array with 24 canonical
fermionic coordinates per site, all charged and calibrated Dirac-neutral
families, and kinetic vector p_j=2 sin(k_j). Every valley and helicity is
retained. The two hadron cells couple to lepton sites (0,0,0) and (1,0,0).
This spatial embedding is a declared architecture. It does not supply a
centre-of-mass kinetic action for the composite hadrons.

The reference fills the negative bands of this same three-dimensional
massive Dirac operator. Its active covariance at momentum p is
(I+sigma.p/E)/2. Opposite momenta give local covariance I/2 and zero weak
current means. The producer reconstructs the offsite covariance and the
complete native current-column norm by Wick's identity

    <J K> = Tr[J(I-rho)K rho].

The literal quark currents, all three weak generators, three families and
four Pauli-current components enter this contraction. It includes every
first-pair channel; the full action includes all subsequent multipair
sectors. The one-quark, one-family contact norm is six by the exact local
polynomial in `openwave_common_mass_source`. Summing three quarks and three
lepton families proves ||V||<=54, independent of the cubic volume.

## Angular polynomial from the native current

The actual native proton-ground blocks are I for the temporal current and
g sigma_j for the spatial currents, with g in the certified interval
(1.2752,1.2754). Each contact supplies half this block. No effective axial
transition matrix defines the source. For spin-up neutron input, proton
spin s and lepton spin matrix, the summed-contact amplitude is

    M_s = delta_(s,up) I - g sum_j sigma_j(s,up) sigma_j.

Let a=p_e/E_e and b=p_antineutral/E_antineutral. The left spin densities
are (I-sigma.a)/2 and (I-sigma.b)/2. Direct symbolic Pauli contraction gives

    sum_s Tr[D_e M_s D_antineutral M_s*]
    = [1+3g² +(1-g²)a.b +2g(1-g)a_z +2g(1+g)b_z]/2.

The executable checks the polynomial symbolically and independently against
the actual 120-state native current. Its normalized coefficients are

    a0=(1-g²)/(1+3g²), A=2g(1-g)/(1+3g²),
    B=2g(1+g)/(1+3g²), D=0.

In the conventional neutron sign convention lambda=-g, these are the
zero-recoil allowed-current coefficients. The finite two-contact source
multiplies this polynomial by cos²[(k_e,x-k_h,x)/2], where k_h denotes the
canonical hole coordinate and the physical antineutral momentum is -p_h.
It also retains the exact finite-time sinc² factor of the strong energy
and both emitted energies. Thus the finite cube instrument has its actual
form factor; equating it with a pointlike continuum angular distribution
would require a transfer theorem.

The named PDG 2025 neutron table provides zero-recoil a0 and A comparisons
and the B comparison. The axial calibration has shared experimental ancestry
with that table. These comparisons are conditional checks, without an
independent-fit or statistical-significance claim. The code reports the
differences, including the a0 and B discrepancies, rather than declaring
all three measured coefficients matched.

## A strict parity-odd record for the full Fock source

Select the actual proton band and the light charged family, sum all three
neutral masses, and resolve both momentum z signs. The prepared neutron
has native spin up. The density formula sums every outgoing spin. Under
reflection of the electron z momentum, the contact form factor and the
energy mismatch are unchanged. Summing the hole z momentum cancels its odd
terms exactly, including for massive neutrals. The difference between the
positive and negative electron z records consequently has negative sign
from 2g(1-g)a_z. This is a three-dimensional momentum statistic.

For L=4 or L=8 at least L/4 positive z modes have p_z>=sqrt(2), and
E_e<=sqrt(13)<4. Thus a_z>1/3 for at least N/4 electron modes, N=L³.
The sum of the contact form factor over all hole modes is N/2. Unitarity
sums the neutral family weights to one. At zero time-integral mismatch,
the negative probability contrast is at least

    (G t)² g(g-1)/12.

Every selected mismatch has magnitude below nine. The inequality
1-sinc²(x)<=x²/3 and the selected spin-weight bound below 16 give the
conservative finite-time correction 16(Gt)²(9t)²/12. This argument uses
massive bands directly; no massless zero-mode replacement occurs.

Duhamel's identity on the entire admitted Fock space gives

    ||psi_exact-psi_free-psi_first|| <= delta=(54Gt)²/2.

The selected event is orthogonal to the free input, and its first-order
norm is below 4Gt. A bounded parity observable supported on that event
therefore has probability error at most (8Gt+delta)delta. At G=1/64 and
t=1/10000 the rational lower margin in the receipt is positive for both
cubes. The exact full-source parity contrast is strictly negative. No
pair-sector truncation is used to establish this statement.

## Complete records, units and work

The rank-one neutron preparation on the uniform 120-state baryon input has
success 1/120 and retained failure 119/120. Nine joint momentum-sign
outcomes plus the entire complementary quantum source form a complete
terminal instrument. Record intervals are jointly dependent bounds on this
one instrument, not independently normalizable postselections. The receipt
also gives complete neutral-family totals and electron radial-energy bins.
Position and momentum are separate instruments. Group velocity is
v_j=4 sin(k_j)cos(k_j)/E, which is not the momentum sign on every valley.

The complete autonomous source conserves mean energy. Preparation has zero
nonselective mean work because its strong filter is spectral and every
lepton current mean vanishes. The terminal instrument commutes with H0;
its retained interaction-readout work obeys |W|<=2G||V||<=27/16.
The electron anchor converts free energies to GeV. The displayed positive
strong/free-channel endpoint is not a full interacting pole endpoint.
Neither a physical weak clock nor a neutron lifetime is established here.
