# Complete finite baryon charged-current conversion

## Common source and reducing sectors

The hadronic action is exactly `model/openwave_calibrated_hadrons`, with
the same electron unit, mass calibration and positive native current
terms. Its complete three-quark occupation sector has dimension 2024;
the exact colour-singlet embedding has dimension 120. There are two
cells, two spins and all fifteen native internal coordinates, hence sixty
fermionic modes for this declared one-generation member. The eight
active lepton modes are the native neutral/charged weak doublet. The
other native occupations are retained as invariant spectators.

Fill the four neutral modes and leave charged leptons empty. Total lepton
occupation four is conserved, giving all `binomial(8,4)=70` states. The
charged one-body operator is `I+qL/2`; the neutral operator is `-qL/2`,
with the filled-reference constant `2q`, where
`L=[[1,-1],[-1,1]]`. Neutral holes have positive kinetic energy. Their
rest energy is zero in this member. The reference filling, including its
zero modes, is an explicit preparation choice.

The action adds the literal native cross-current

`G sum_(x,mu) eta_mu (Q_plus^mu L_minus^mu + Q_minus^mu L_plus^mu + 2 Q_3^mu L_3^mu)`,

with `G=1/64`, `eta=(1,-1,-1,-1)` and
`bar_sigma=(I,-sigma)`. Using `sigma` for both spatial bilinears gives
the identical product. All spin components, neutral-current terms and
lepton occupation states occur. There is no selected neutron/proton
transition coefficient in the action. The native colour-current
intertwining established for the form factors applies to every term.

The product colour-singlet/lepton space has dimension 8400. Electric
charge conservation gives `N_e=N_up-1`; its charge-zero sector has
dimension `40+40*16+20*36=1400`. In particular, two charged electrons and
two neutral holes are allowed when all three quarks are up type. They
must be retained even though their asymptotic threshold is closed.
Fixing twice the total spin projection to one gives 398 states. Every
discarded entry of the corresponding off-diagonal Hamiltonian block is
zero. These exact conserved-quantity reductions leave the complete
current-reachable complement of the prepared state inside the calculation.

## Preparation, box and thresholds

An explicit normalized colour-singlet orbital seed is filtered onto the
even-parity neutron ground state with spin up. The filter is rank one;
the opposite parity, other spin/charge states and all other Fock states
belong to its complement. The success probability is approximately one
third. The full failed vector and its unchanged flag are retained and
evolved. The seed energy and the mean work of the complete filter are
computed with the full interacting Hamiltonian.

The box coefficient is declared as
`q=(delta-1)/2`, where `delta=(m_n-m_p)/m_e` is the fitted strong-source
splitting. This makes the free two-odd-mode electron/hole energy equal
to the strong nucleon gap. Contact terms and all their detunings remain
in the Hamiltonian. Choosing this box uses the calibrated gap; it is
neither a natural decay-volume prediction nor an empirical lifetime fit.

The one-electron endpoint is `0.00078233355931 GeV`. It follows from the
three hadron calibration inputs and the unchanged electron anchor. The
quoted endpoint belongs to the strong source plus free outgoing leptons.
The complete weak source has a filled-neutral-reference contraction
`V_reference=2 I_3`. Its neutron and proton expectation shifts are
respectively `-G m_e` and `+G m_e`, about minus and plus7.984keV.
These shifts are retained in the evolution. The fitted strong masses
are therefore not asserted to equal the physical poles of the weak
extension, and the free-channel endpoint is not an interacting endpoint
certificate. Expectation energies alone also do not define a decay pole.
The
lightest charge-two baryon plus two electron rest energies lies above
the neutron. Virtual two-electron amplitudes are nevertheless present.
The distinct vectorlike/balanced source whose nucleons were about
1.7 MeV retains its failed mass and beta-threshold tests. No transfer
identifies that failed source with this calibrated one.

## Complete records and conservation

Separate complete experiments measure the lepton positions or their
two-cell momenta, together with native hadronic occupation records.
The momentum transformation is the fourth exterior power of the actual
eight-mode Fourier matrix, including fermionic signs. It preserves
charge and spin projection. These are alternative instruments; no joint
position-and-momentum value is assigned.

Every stored product-label probability expands into native colour
occupation probabilities through the exact colour-embedding weights.
The receipt contains every product index, lepton mask and nonzero native
colour mask/weight. Its original probability is preparation success
times the stored conditional probability times that colour weight.
Both complete experiments also contain preparation failure and the
orthogonal Fock complement. The latter has zero prepared-state weight
by the proved invariant sectors. Neither emission nor success is
silently renormalized in the original record law.

The full Hamiltonian conserves energy during evolution. Fine colour
readout breaks a coherent singlet and has nonzero work. Position and
momentum instruments have different work: the calculation uses the
actual full occupation diagonal, or its lepton-Fourier transform,
rather than the lower coherent colour-singlet energy. The failure
state evolves with the same source and preserves its original norm.
The general finite Lüders-instrument theorem supplies complete
positivity, normalization and external-reference preservation for
these complete orthogonal projectors. The concrete instances are
executable algebra checks, not additional Lean theorems.

At evolution parameter 512, floating evaluation gives approximately
`0.96855747` conditional exactly-one-electron probability,
`0.96246817` conditional proton-band probability and `0.00020450`
conditional two-electron probability. Multiplication by preparation
success gives about `0.32285249` original exactly-one-electron
probability. Subsequent return exhibits finite coherent recurrence.
The source parameter orders this declared experiment; it is not an
external universal clock or an identified physical lifetime.

Independent matrix-exponential and spectral propagation agree at
parameter 16 within the stored residual. Zero contact gives zero
emission. The record laws preserve normalization and energy. Numerical
probabilities are floating evaluations, without outward enclosures of
every arithmetic error.

## Physical boundary

This source supplies an allowed finite nucleon conversion using one
calibrated mass action and its actual native weak current. It supplies
no outgoing continuum, irreversible neutron lifetime or observed
angular distribution. The same current's certified axial ratio near
`5/3` fails the held-out magnitude near `1.2753`. A mass fit and a
finite conversion signal therefore do not establish a complete
physical neutron-decay model.
