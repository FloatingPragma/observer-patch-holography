# The same family coefficient in a free relativistic field

The configured coherent family program supplies the Hermitian positive
coefficient `Y=C_ex^2/3`. The full-CAR construction in
[CAR_KOIDE.md](CAR_KOIDE.md) realizes that coefficient with complete quantum
archives. Its free kinetic attachment uses the existing native charged-Higgs
line and the same spin/family ordering. The resulting continuum statement
is a free Dirac theorem, with all lattice valleys retained. Dynamical gauge
interactions, a source-selected Higgs vacuum and the physical clock are
separate requirements.

## Exact mass shell of the native quadratic action

Let `D=P_L C_neutral P_R` be the actual directed charged-lepton Higgs tensor
on the fifteen-component register. The source verifies `rank(D)=1`,
`D^2=0`, and unit singular normalization. Its symmetric all-left kernel is

\[
 K=D\otimes\overline Y+(D\otimes\overline Y)^T.
\]

The native ninety-mode fermion convention has

\[
 h(p)=-c\,\sigma\cdot p\otimes I_{45},\quad
 \Delta=\epsilon\otimes\overline K,\quad
 H(p)=\begin{pmatrix}h(p)&\Delta\\\Delta^\dagger&-h(-p)^T\end{pmatrix},
 \qquad \epsilon=i\sigma_2.
\]

For every spatial axis, `sigma_j epsilon+epsilon sigma_j^T=0`.
Together with the three-dimensional Clifford identity, this gives exactly

\[
 H(p)^2=\operatorname{diag}
 \big(c^2|p|^2I+\Delta\Delta^\dagger,
      c^2|p|^2I+\Delta^\dagger\Delta\big).
\]

The six nonzero singular coefficients of `K` are the three singular values
of `Y`, each twice. The other thirty-nine internal coefficients vanish in
this charged-only action. A neutral Higgs norm `w_H` therefore gives

\[
 E_j(p)^2=c^2|p|^2+m_j^2,\qquad m_j=w_H\sigma_j(Y).
\]

The charged block has twelve positive-energy Nambu rows. Nambu doubling
does not add species; after canonical quantization these are three free
Dirac fields per valley with two spin states and conjugate charged excitations. The
same native charge action gives opposite charges on creation and conjugate
creation components. This is a field identification within the configured
model, without identification of its three fields as observed leptons.

The constant family diagonalization and the native hole-spin conjugation
put each charged block into Dirac form. One convention is
`alpha_j=tau_3 tensor sigma_j`, `beta=tau_1 tensor I`,
`gamma^0=beta`, `gamma^j=beta alpha_j`. The sixteen exact anticommutators
give `2 diag(1,-1,-1,-1)`. The generators
`[gamma^mu,gamma^nu]/4` therefore supply the spinor Lorentz action, and
the mass term is a Lorentz scalar. This proves the free transformation
structure as well as the dispersion identity.

`relativistic_mass.py` checks the exact spin identities and reuses
`common_source.native_record_yukawa_field.dispersion_certificate` on the
complete 180-dimensional native symbol. Numerical controls check two
nonempirical phase choices and all eight valley orientations. They are
execution checks of the exact formula, not outward roundoff certificates.

## Controlled lattice limit and reference state

For lattice spacing `a`, use the same mass kernel and
`p_j -> sin(k_j)/a`. Every corner `k_b=pi b`, `b in {0,1}^3`, has a branch
`k=k_b+a p`. Its limiting kinetic momentum is `((-1)^b_j p_j)_j`.
All eight branches have the same masses; their chiral current assignments
retain the mirror obstruction in the parent construction.

For `|p_j|<=P`, the elementary remainder bound
`|sin(a p_j)/a-p_j|<=a^2|p_j|^3/6` proves

\[
 \|H_a(k_b+a p)-H_b(p)\|\le\delta_a,
 \qquad \delta_a={c a^2P^3\over2}.
\]

Eigenvalue errors are at most `delta_a`, and Duhamel's identity gives
unitary error at most `|t| delta_a` in the same declared clock. No additional
factor of capacity or lattice spacing multiplies the rest masses.

Assume the three charged masses have a proved positive floor `m_*`.
Choose the filled negative-energy Gaussian reference of this quadratic
action. After a common, spacing-independent family diagonalization, each
Dirac block has `H^2=E^2 I` and negative projector `(I-H/E)/2`.
Since `E_a,E>=m_*` and `|E_a-E|<=delta_a`,

\[
 \|P_{-,a}-P_-\|\le\delta_a/m_*,\qquad
 \|e^{-itH_a}P_{-,a}-e^{-itH}P_-\|
 \le\delta_a(|t|+1/m_*).
\]

These estimates give convergence of the charged two-point functions
smeared with fixed compact-momentum test functions; Wick's rule gives
convergence of each fixed Gaussian correlation. The limiting free Dirac
equation is hyperbolic with characteristic speed `c`. Its causal
anticommutator is obtained by applying the Dirac operator to the
Klein-Gordon causal propagator and has light-cone support. This supplies
the free relativistic field structure in this chosen reference. It is not
a proof that the native repair process selects that vacuum.

A momentum packet concentrated near one corner defines a prepared valley
response. It does not remove the other corners from the theory. The complete
local free detector has the parent's original whole-cell weights, and a
sharp single-valley projection is not covariant under local gauge rotations.

## Operational mass and physical identification

In this free model, the finite functional-calculus instrument in
[`spectral_readout`](../spectral_readout/README.md) supplies complete energy
effects for a provided controlled propagator. The charge-paired construction
in [`LOCAL_SPECTRAL_READOUT.md`](../source_family_compatibility/LOCAL_SPECTRAL_READOUT.md)
also applies to normalized Bogoliubov modes by CAR and Fourier orthogonality;
its unrelated source-specific collar constants are not imported. Under
these controlled-evolution and reference-preparation services, the finite
clock reads `E_j(p)` with its resolution and alias bounds. The numerical
symbol checks here execute no full-Fock energy tape. The family-resolved rest energy is `m_j`;
equivalently `E_j(p)^2-c^2|p|^2=m_j^2`. The phase score in
[PHASE_AND_PREMISES.md](PHASE_AND_PREMISES.md) reads its family coefficients
without a mass label. Canonical quantization and the source charge map then
identify them as free charged field mass parameters. The conserved free
family labels distinguish these branches; the mass is not the threshold of
the entire electric-charge sector, which includes other source components.

The weak-feedback limit and the spatial limit use specified pulse clocks.
Their error budgets can be composed when the clock, ancillary services and
kinetic evolution are implemented together. Fixed analyzer gate overhead
is not automatically negligible in the native scheduler's attempt clock.
This statement supplies no uniform physical-capacity service schedule.

The theorem does not establish an interacting chiral continuum, measured
weak decay widths, gauge-dressed charged states or experimental mass
matching. Turning on the finite gauge links changes the reference and
spectral problem. The mass-response theorem in
[MASS_RESPONSE_BRIDGE.md](MASS_RESPONSE_BRIDGE.md) states the estimates
required to carry these action coefficients into an interacting observable;
its correction bounds must be supplied by that same theory.
