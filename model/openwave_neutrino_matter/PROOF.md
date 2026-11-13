# Matter oscillations from the native forward current

## Declared forward source and its exact current

The three neutral Dirac masses and family tensor are those of
`openwave_neutrino_calibration`. The medium consists of fresh electron
ports, each with the exact positive-energy rest-electron active weight
1/sqrt(2) and unpolarized spin state I/2. Density, composition, length and
number of ports are declared apparatus parameters. The forward-only
positive-frequency interaction is a declared Hamiltonian: its coefficients
are the forward matrix elements of the literal native weak current. No
matter potential or probability table defines those coefficients.

This member is not asserted equivalent to the unfiltered momentum-changing
contact source. In particular it supplies no cross-section or suppression
bound for all channels outside the forward port architecture. The full
forward source is a number-conserving quartic polynomial
sum M_(is,jr) a_i^dagger a_j b_s^dagger b_r. Its complete one-neutral,
one-electron sector has dimension twelve; all other occupation sectors
are retained. The neutral coordinate includes three masses and both
positive-energy helicities. Other native species retain free spectator
vacua. A fixed electron input is a medium preparation, not a discarded
outcome of scattering.

Let A_(mu,a)=sigma_mu tensor T_a tensor I_family, with T_a formed from
the actual fifteen-state native raising current. In the neutral flavour
spin basis, the connected forward matrix element is

    sum_(mu,a) eta_mu [
      (nu_i* A nu_j)(e_s* A e_r)
      -(nu_i* A e_r)(e_s* A nu_j)].

The two ordered current products give a factor two; the rest-electron
active weights give one half. The second term is the canonical fermion
exchange term. Exact symbolic contraction yields characteristic polynomial
x^9(x-1)(x+1)^2, hence norm one. Tracing the electron spin with I/2 gives

    Mbar0=(P_e/2-I/4) tensor I_spin.

Thus the electron-flavour term arises from the native exchange contraction.
The common neutral-current term is retained, including its relation to
sterile outputs. The pure weak-isospin current is the source here; no
unproved full electroweak neutral-current completion is inferred.

At momentum p, the active weights are u_(i,s), with the complementary
sterile weights v_(i,s), u²+v²=1. Let Q=(U tensor I_spin)diag(u). The actual
massive current is M=(Q* tensor I) M0 (Q tensor I), and
Mbar=Q* Mbar0 Q. Consequently ||M||<=1, ||Mbar||<=1/4 and
||M-Mbar tensor I||<=5/4. The family matrix, both helicities and every
sterile factor are retained. Charge conjugating the neutral forward
bilinear reverses its normal-ordered one-body sign and transposes its
matrix. In the dual helicity basis the antineutral matrix is minus the partial
transpose on its neutral indices; the electron indices are unchanged.
Its unpolarized mean is -Mbar*, with U replaced by U*. The exact uncompressed
antineutral polynomial is (2x-1)^7(2x+1)^5/4096, so its norm is one half.
The common norm-one fluctuation bound also covers this sector. Conjugating
the entire joint matrix would incorrectly transpose the electron medium.

## Complete finite medium and the coherent comparison

Set T to the declared flight parameter and alpha=G_Weyl n_e T. The complete
unitary for one of N equal intervals is

    exp[-i(H_free T tensor I +alpha M)/N].

For each initial and final electron spin, its matrix block divided by
sqrt(2) is a Kraus operator. Four such operators exactly complete one
medium step. Their N-fold products retain every electron-spin record;
tracing them is an explicit completely positive channel, not a selected
survival branch. The executable evaluates this complete channel by powers
of its 36-dimensional superoperator. That is a finite-algebra shortcut to
the full instrument with 4^N record words.

The comparison unitary has Hbar T=H_free T+alpha Mbar. In its interaction
picture the fluctuation W has mean zero in the electron input and norm
at most B=5|alpha|/4. The reduced-channel Dyson series has zero first-order
term. Its higher terms have trace-norm sum bounded by
exp(2B/N)-1-2B/N <= 2B² exp(2B/N)/N². Telescoping N channels, using trace
norm contraction for arbitrary external quantum references, proves

    TV(complete neutral records, coherent-mean records)
      <= 2B² exp(2B/N)/N.

The bound is conservative by a possible factor two and covers any neutral
terminal instrument. It does not compare the complete environment record
with a nonexistent coherent-medium environment. The finite medium has its
own explicitly retained environmental records. The numerical channel
comparison is diagnostic; the analytic bound does not assume floating
matrix exponentiation is exact.

## Physical units and relativistic comparison

The weak unit dictionary in `openwave_lepton_widths` supplies
G_Weyl=2sqrt(2)GF, with GF's muon-decay calibration ancestry. Number density
in cm^-3 converts to energy cubed using hbar c, and the declared baseline
converts to T=L/(hbar c). The chosen material convention is
n_e=rho Y_e N_A with one gram per mole of nucleon number. It is an apparatus
scenario rather than a density determination. The native potential is

    V_CC=G_Weyl n_e/2=sqrt(2)GF n_e,
    V_common=-G_Weyl n_e/4.

The active relativistic comparison has
H=diag(m_i²-m_1²)/(2p)+V_CC U* P_e U. The common term is a phase only on
the ideal active block; it is retained in the complete six-state source.
For antineutrals the sign reverses and U conjugates. Stable exact dispersion
differences use (m_i²-m_1²)/(E_i+E_1).

Let chi=(m_max/(2p))². On the prepared negative-helicity block,
|1-u_i|<=chi, and the suppressed readout/preparation weights are at most
chi. The mean operator has no helicity offdiagonal block. The dispersion
phase error is bounded by 2chi times the largest relativistic phase;
the active potential change is bounded by |alpha|chi. Triangle inequalities
for preparation, evolution and the complete detector give the conservative
record bound 20chi(1+max|phase|+|alpha|). All physical scenarios here have
chi below 2e-21. The finite-medium bound includes its helicity-changing
fluctuations independently.

No measured matter-event likelihood is fitted. The atmospheric two-flavour
resonance scale is a diagnostic of the derived coefficient; the executed
records use all three flavours and the calibrated complex phase.

## Preparations, readout and work

A native active-coordinate preparation port of flavour mu and specified
helicity creates amplitude U_(mu,i)*u_i. Its failure probability is computed
as sum_i |U_(mu,i)|²v_i², without subtracting nearly equal floats. Its retained
reference has energy p and carries the Dirac number on failure. The
preparation work is the success-weighted sum of E_i-p, evaluated by the
stable quotient m_i²/(E_i+p).

At the terminal apparatus, each active flavour port has rows U_(alpha,i)u_i
for both helicities. The sterile diagonal operator diag(v) completes the
instrument. Active absorption transfers Dirac number to an explicit
zero-energy reference port. The output work includes the removed neutral
energy; the sterile quantum branch retains its energy. This apparatus does
not assign simultaneous sharp position and momentum.

Every medium unitary conserves its source energy between switches. Electron
spin reference ports have the same rest energy on all branches. External
switching work and the final detector work are retained; the absolute
interaction switching bound per interval is 2G_Weyl n_e. The construction
supplies a finite observer medium with declared clocks and lengths. It does
not derive an autonomous beam release, a macroscopic material state from
first principles, or a measured detector event rate.
