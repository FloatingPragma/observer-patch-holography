# Native polarized beta records on an outgoing ring

## One declared source

The strong action and actual quark weak currents are those of
`openwave_axial_state`. The four empirical inputs are proton, neutron,
charged-pion masses and the axial magnitude. The vector current equals one;
the actual ground axial magnitude lies strictly between1.2752 and1.2754.
The source uses its full colour carrier, including the edge colour cost.
The complete120-state three-quark colour-singlet space is reducing. No
proton/neutron-only transition matrix defines the interaction.

Extend the192-coordinate common native source of
`openwave_common_mass_source` by24 Dirac-lepton coordinates per additional
outgoing cell. Quarks occupy the original two cells; other native species
and quark families are retained. Their vacuum sectors are reducing because
the semileptonic current preserves their number and family. The admitted
three-quark preparation therefore has24 native quark coordinates and the
full lepton Fock space on24L coordinates. This source extension changes
the lepton kinetic operator to `p_k=2 sin(2pi k/L)`. Every lattice valley,
charged and neutral mass, antiparticle and helicity is retained.
The momentum scale2 is a declared source coefficient, not a fit to a
physical beta spectrum.

The charged masses are the native analyzer ratios in electron units.
The neutral tensor is exactly the empirical calibration construction in
`openwave_neutrino_calibration`: its six stated oscillation inputs, chosen
lightest mass and normal ordering remain inputs. The tensor is installed
in the native active/sterile pairing coordinates and contracted against
the unchanged native weak current. No current form factor is fitted in
this extension. No physical weak coupling or lifetime calibration enters
`G=1/64`.

The reference is the full massive Dirac Gaussian vacuum. For a momentum p,
spin sign s and mass m, its active positive/negative weights are
`u^2=(1-sp/E)/2` and `v^2=(1+sp/E)/2`, with `E=sqrt(p^2+m^2)`.
Opposite momenta make the local active covariance I/2 and the local weak
current mean zero. Off-site covariance is not diagonal and is recomputed
explicitly. The filled-neutral scalar source used in
`openwave_dressed_baryon_lines` is a different reference. Its contact energy
shifts are not transferred to this Gaussian reference.

The full interaction remains the local native semileptonic current at
cells0 and1. Its exact64-state one-quark/one-family Fock norm is6. Summing
three quarks and three lepton families gives the bound54 on the entire
admitted lepton Fock space. This includes every multipair state; no pair
cutoff is imposed on the source evolution.

## Complete first-current column

For every hadron eigenstate, both contact cells, all four Pauli-current
components, both spins and all momenta, the producer contracts the actual
native current with the massive positive/negative bands. Its24 family
channels comprise nine electron-antineutral, nine neutral-positron, three
electron-positron and three neutral-antineutral channels. The neutral
current is diagonal in the neutral mass frame because the family matrix
is unitary. Both charged-current orientations and all neutral currents are
included. The local vacuum mean vanishes, so these pairs exhaust V applied
to the prepared free vacuum; higher pairs enter the exact subsequent
unitary and its error bound.

The time integral is evaluated with the complete hadron energy and pair
energy. It is `t exp(-i DeltaE t/2) sinc(DeltaE t/2)`. The receipt retains
all channel totals, every hadron spectral weight and all momentum records
for the actual proton/light-electron/neutral terminal instrument. The
remaining full source is a retained quantum outcome. These are first-order
numerical distributions with explicit exact-source error; they are not
asserted to be exact full-unitary distributions.

An independent Gaussian calculation verifies completeness. On the two
contact cells, build the complete native active covariance rho from all
momentum and family bands. For native lepton matrices J and K,
`<JL JK> = Tr[J(I-rho)K rho]` because their means and active-active anomalous
covariance vanish. Contract this matrix with every actual quark-current
column. Its squared norm agrees with the sum of all24 pair channels,
including the off-site covariance. The zero-momentum common-source
covariance is not reused on the ring.

## Exact parity contrast

Prepare the exact strong neutron spin-up state. The terminal observable
is +1 for an actual strong proton with a light electron of positive
Dirac momentum, -1 for negative Dirac momentum, and zero for zero Dirac momentum
and the entire source complement. All neutral masses and both helicities
are summed. Spin is the fixed native spin label; this is a one-dimensional
momentum-parity statistic, not the physical three-dimensional neutron coefficient.
The group velocity is dE/dk and differs from the sign of p on some lattice
valleys. No spatial outgoing-velocity interpretation is assigned to these bins.

Reflection and rotation imply that each local ground proton current is
half the global native current. With axial magnitude g, its spin tensor is
`delta_sp,up delta_se,sn - g sum_a sigma_a(sp,up) sigma_a(se,sn)`.
After summing proton spin, the three nonzero squared spin coefficients are
`(1-g)^2`, `4g^2` and `(1+g)^2`; the fourth is zero. These identities are
checked against the actual120-state matrices. The contact form factor is
`cos^2((theta_h-theta_e)/2)`.

In the zero-neutral-mass limit with the massive-reference zero modes
retained at half weight, set `a=p_e/sqrt(1+p_e^2)` and `b=sign(p_h)`.
Simultaneously reflecting both momenta gives the spin-weight difference
`2g[(1-g)a-(1+g)b]`. Summing hole momenta, the term in b is nonpositive
because its surviving form-factor term is proportional to
`sin(theta_e) sum_h sin(theta_h) sign(sin(theta_h))`.
For positive electron Dirac momentum, the remaining contribution is strictly
negative whenever g>1. At least L/4 electron modes have `a>4/5`, and the
sum of the contact form factor over hole momenta is L/2. Therefore the
negative first-order contrast has magnitude at least

    (G t)^2 g(g-1)/5 times the minimum sinc^2 weight.

The selected energy mismatch has magnitude below7, so this weight exceeds
`1-(7t/2)^2/3`. Exact neutral squared masses imply `mu_max<10^-7` in electron
units. For L divisible by8, the smallest nonzero momentum exceeds8/L.
The receipt bounds the massive-neutral correction by
`24(Gt)^2[mu_max^2/(4p_min^2)+t mu_max]`. The unitary family sum is complete;
no empirical mixing angle is needed to prove this sign.

For the exact full interaction-picture evolution, two applications of
Duhamel's identity and unitarity give

    ||psi_exact-psi_free-psi_first|| <= (54 G t)^2/2 = delta.

The terminal event is orthogonal to the initial state. Its first-order
norm is at most4Gt, since the total native spin weight is below16 and the
family rows are normalized. Its odd-probability error is consequently at
most `(8Gt+delta)delta`, rather than an error linear in the full state
norm. At t=1/10000 the negative lower bound survives this full-Fock error
for L=8,16,32,64. The exact positive-minus-negative momentum probability is below
`-1.2684e-13` conditional on preparation. The original-input bound is
smaller by120. This analytic bound uses certified axial and exact mass
intervals; floating spectral sums are diagnostic values inside that
proved domain.

The scalar kinetic and filled-neutral control commutes with ordinary
spatial reflection. Its even spin-polarized preparation gives zero
positive/negative momentum contrast for all times. The chiral kinetic/reference
extension is explicit; the scalar control is not silently reinterpreted.

## Records, calibration and work

The input is uniform on the120-state native baryon carrier and has the
specified full lepton Gaussian reference. A rank-one neutron filter
succeeds with probability1/120. Its119/120 failed quantum branch remains
in the record. Nonselective preparation has zero mean work: the strong
filter is spectral and every lepton current mean vanishes. The full
Hamiltonian is autonomous and conserves its mean energy.

The terminal proton/momentum instrument commutes with H0. Its full-current
readout work is retained, with absolute bound `2G||V||<=27/16`. Positive momentum,
negative momentum, zero-momentum and the complementary quantum outcome form one
complete instrument. The receipt gives joint-dependent probability
intervals; they are not separately postselected normalized distributions.

The four finite rings provide an executed volume comparison. This proves
neither an infinite-volume angular limit nor an irreversible asymptotic
release. The source time and contact coefficient have no physical weak
clock calibration here. The source neutron/proton masses and axial input
are conditional calibrated strong properties; full relativistic nucleon
poles, the physical three-dimensional beta spectrum and the neutron
lifetime are not established by this finite statistic.
