# Regional entropy, modular motion and the native response

## One fixed interacting source

Use the radius-one, 27-site bounded radial source of
[`radial_interacting_core`](../radial_interacting_core/PROOF.md), at spacing and
velocity one, mass 64, and centered contact coupling \(U=1/65536\). Its full
54-mode Fock Hamiltonian is \(H\). Write \(\rho_G=G/2\) for its actual ground
doublet mixture and \(\rho_C=(C_-+C_+)/2\) for the equally weighted actual
lowest states in the neighboring particle-number sectors. The two charged
states have the same energy \(E_G+g\), with \(g>0\). These are spectral
projectors of the full Hamiltonian, including its entire virtual complement.

At the origin define
\[
 D=(n_\uparrow+n_\downarrow-1)^2
   =|00\rangle\langle00|+|11\rangle\langle11|,
 \quad W_0={U\over2}D-{U\over4}I.
\]
The parent proves \(n_G=\operatorname{Tr}\rho_GD<1/1000\) and
\(n_C=\operatorname{Tr}\rho_CD>999/1000\). Put \(d=2^{54}\),
\(\epsilon=1/64\), and, for \(1/4\le q\le3/4\),
\[
 a_C(q)={q-\epsilon/2-(1-\epsilon)n_G\over n_C-n_G},
 \quad a_G(q)=1-\epsilon-a_C(q),\qquad
 \sigma_q=a_G(q)\rho_G+a_C(q)\rho_C+\epsilon I/d.                 \tag{1}
\]
Both coefficients exceed \(6/25\): their numerators are bounded below by
\((1-\epsilon)(999/1000)-3/4+\epsilon/2>6/25\), and their common denominator
is at most one. Thus (1) is a faithful stationary state of this same \(H\),
with \(\sigma_q\ge2^{-60}I\), and \(\operatorname{Tr}\sigma_qD=q\).
No thermodynamic or zero-temperature limit is used. Selecting the spectral
states and calibrating their mixture from \(n_G,n_C\) are supplied preparation
inputs; this existence construction does not execute that preparation.

Number conservation removes the local matrix elements between different
occupation sectors. The source's particle-hole symmetry exchanges empty and
double occupation. Its commuting antiunitary
\(\Theta=R D_{\rm stagger}(i\sigma_2)K\) restricts at the origin to spin
time reversal, so the one-particle local block is a scalar matrix. The ground
mixture, each unique neighboring charged state under \(\Theta\), their
particle-hole paired mixture, and \(I/d\) have these invariances. Therefore
the **actual origin reduction** of (1) is
\[
 \tau_q={q\over2}\Pi_{\rm even}+{1-q\over2}\Pi_{\rm odd}.          \tag{2}
\]
The same argument applies to the faithful Gibbs state of this full \(H\),
with its own source-determined parity expectation. Formula (1) supplies the
explicit interval of regional states used below.

The actual ground mixture also has a faithful central-site reduction, without
the identity admixture. In the origin/complement block matrix the coupling
satisfies \(V^*V=(3/2)I_2\): each of the six distinct neighbors contributes
\(I_2/4\). The complement \(A\) has \(\|A\|\le64+2=66\), so
\(\|A^{-1}V\xi\|\ge\|V\xi\|/66\). Its scalar normalized central zero
weight therefore obeys
\[
 \alpha\le\alpha_+={1\over1+(3/2)/66^2}={2904\over2905}.
\]
The parent's full-state trace-distance bound gives an effect error at most
\(t=27/16252874\). Its free-ground parity is \((1-\alpha^2)/2\), whence
\[
 {1\over4096}<{1-\alpha_+^2\over2}-t\le n_G<{1\over1000}.
\]
Thus the actual origin density has minimum eigenvalue at least
\(((1-\alpha_+^2)/2-t)/2>0\). The full ground state has rank two and is
not faithful on the full Fock algebra. Its exact \(n_G\) is not numerically
computed or prepared by this certificate.

The full-source energy and local contact satisfy
\[
 {d\over dq}\operatorname{Tr}(H\sigma_q)={g\over n_C-n_G}>0,
 \qquad {d\over dq}\operatorname{Tr}(W_0\sigma_q)={U\over2}.       \tag{3}
\]
The identity admixture cancels from the derivative. These are finite source
energies, without an identification with a conserved Lorentzian stress tensor.

## Exact regional modular motion and its native word

On the complete four-dimensional site algebra, including odd fermionic
operators, let \(\kappa_q=\log((1-q)/q)\). Then
\[
 K_q=-\log\tau_q=\log{2\over1-q}\,I+\kappa_qD,
 \qquad \tau_q^{is}A\tau_q^{-is}
   =e^{-is\kappa_qD}A e^{is\kappa_qD}.                          \tag{4}
\]
Since \(D=(I+Z_\uparrow Z_\downarrow)/2\), the source word is an identity
phase of angle \(-s\kappa_q/2\), followed by CNOT, a Z pulse of the same
angle, and inverse CNOT. The native CNOT is its actual three-pulse H,J,H
word, including the source scalar phase. The resulting eight pulses act on
all four occupations and extend as the identity on any remote field and
reference. The certificate executes the full four-column reference isometry,
retaining every native complement failure output; success is neither
postselected nor renormalized.

For the local contact's Heisenberg flow
\(\alpha_t(A)=e^{itW_0}Ae^{-itW_0}\), (4) is exactly
\(\alpha_{-2\kappa_q s/U}(A)\). At \(q=1/2\) it is the identity. The time
factor changes sign above \(1/2\). This is the local contact clock, not the
full interacting Hamiltonian flow, whose hopping terms leave the site
algebra. The full Gibbs state's modular factor is instead \(-\beta\) for
the full Hamiltonian; (1) is generally not that Gibbs state.

All nonrational controls have explicit rational enclosures. With
\(z=1-2q\), \(|z|\le1/2\),
\[
 \kappa_q=2\sum_{j=0}^{N}{z^{2j+1}\over2j+1}+R_N,
 \qquad |R_N|\le {2|z|^{2N+3}\over(2N+3)(1-z^2)}.               \tag{5}
\]
The two variable-angle errors sum to at most \(|s|\) times the logarithm
error. Fixed pi angles and all twelve rounded source ports have separate
complete-isometry bounds. For exact rational integrated port commands of
precision \(b\), the refinable isometry error is at most
\(8B2^{-b}+|s||R_N|\), where \(B\) is the source sum of generator norm
majorants. Complete-channel unhalved trace distance is at most twice this
bound, for arbitrary untouched references. These rational commands have a
different interface from the executed binary64 controls; physical control
precision and absence or compensation of Hamiltonian drift are apparatus
inputs.

For the larger calibration domain \(1/4096\le q\le4095/4096\), containing
the actual \(n_G\), range reduction writes the absolute log odds as
\(k\log2+\log y\), with \(1\le y\le2\) and \(0\le k\le11\).
Both logarithms use (5) with \(z=(y-1)/(y+1)\le1/3\), giving a uniform
rational algorithm. The executed examples cover both domain endpoints,
\(q=1/4,1/2,3/4\), and modular times in \([-1/2,1/2]\); the exact ground
value requires source calibration. The refinable commands explicitly
approximate all twelve port coefficients of every pulse, including the
fixed pi factors, with a separate logarithm remainder. Refining both series
order and port precision sends the bound to zero.

## A shared entropy, energy and response relation

Put \(h(q)=-q\log q-(1-q)\log(1-q)\). Equation (2) gives
\[
 S(\tau_q)=\log2+h(q),\quad S'(q)=\kappa_q,
 \quad S''(q)=-{1\over q(1-q)}.
\]
For any two members of (1), the exact first-law remainder is
\[
 S(\tau_q)-S(\tau_r)=\kappa_r(q-r)
   -D(\tau_q\Vert\tau_r),\qquad
 2(q-r)^2\le D(\tau_q\Vert\tau_r)\le{8\over3}(q-r)^2.           \tag{6}
\]
The relative entropy is precisely the Bernoulli divergence. Its Hessian on
the interval lies in \([4,16/3]\), proving the inequalities by the integral
Taylor remainder. At \(r=1/4\), \(\kappa_r=\log3>0\), and (3) yields
\(\Delta S=\beta_{\rm eff}\Delta E-D\), with
\(\beta_{\rm eff}=\log3\,(n_C-n_G)/g>0\). This is an exact relation along
the specified affine source family, without a global Gibbs claim.

The fixed parity instrument of
[`matter_geometry_response`](../matter_geometry_response/PROOF.md) has
\[
 K_0={7(I-D)+D\over\sqrt{50}},\quad
 K_1={i((I-D)+7D)\over\sqrt{50}},\quad
 p_1={1\over50}+{24q\over25},\quad
 \theta={3\over25}(1-2q).                                      \tag{7}
\]
The parity-to-probe encoding and its execution are the separate
[`radial_geometry_join`](../radial_geometry_join/PROOF.md) construction;
the entropy identities use only the exact two Kraus operators and their
complete output instrument. This package has no dependency on that join's
receipt. Its native modular word is independently executed.

The same \(q\) in (1), (3), and (7) controls the all-seam native mean law
\(M=I-[L+\theta(D_{(0,1)}-D_{(0,5)})]/64\), its retained pair clocks, and
its uniformly isolated rank-three Gram matrix. Consequently
\[
 \Delta S=-{25\kappa_r\over6}\Delta\theta-D(\tau_q\Vert\tau_r). \tag{8}
\]
Fresh independent source preparations and independent native work are inputs
to iteration of this mean law. The one-attempt instrument identities hold
for arbitrary global field correlations; neither fresh-state iteration nor
source preparation follows from them.

For any four-state density matrix \(\omega\), including coherences,
\[
 S(\omega)-\kappa_r\operatorname{Tr}(D\omega)
 =S(\tau_r)-\kappa_r r-D(\omega\Vert\tau_r).                    \tag{9}
\]
Thus \(\tau_r\) uniquely maximizes the regional entropy at fixed parity
expectation \(r\), equivalently at the fixed native response statistic
\(\theta_r\). The unconstrained functional
\(S(\omega)+(25\kappa_r/6)\theta(\omega)\) has the same unique maximizer.
The coefficient is derived from (2), and \(\theta\) is the source kinetic
bias, not an area or a volume. In particular \(q\) and \(1-q\) have equal
regional entropy and opposite bias. No entropy definition of area is made.
Sampling the full pointer records estimates \(q\) with inverse slope
\(25/24\); on the interval \(|\kappa'|\le16/3\), which supplies a
controlled modular-calibration error in addition to (5).

## Complete records and measurement work

For any field state \(\rho\), the unread instrument (7) is
\[
 \Phi(\rho)={7\over25}\rho+{18\over25}\Delta_D(\rho)
 ={16\over25}\rho+{9\over25}(I-2D)\rho(I-2D).                   \tag{10}
\]
It preserves parity expectation and the one-site state (2) for the
symmetry-constrained family. It can
change the global field's entropy and energy because that state need not
commute with \(D\). Let \(R\) purify the full input and \(B\) be the
committed pointer bit. Each branch has one Kraus operator and a pure
conditional state on field plus \(R\). Therefore
\[
 S(B,\mathrm{field})_{\rm out}-S(\rho)
   =H(B)-I(R:B).                                               \tag{11}
\]
This retains the classical record; the unmeasured coherent pointer has a
different entropy. For a regional input \(\tau_q\), (7) is a binary
symmetric observation of its parity with crossover probability \(1/50\):
\[
 I(D:B)=h(p_1)-h(1/50),\qquad
 S(B,\mathrm{site})_{\rm out}-S(\tau_q)=h(1/50).                \tag{12}
\]
The last identity cannot be substituted for (11) on the full correlated
source. Tracing the remote field changes the reference information.

The native scheduler retains its eight-bit word, placement coin, selected
seam or idle label, committed current and terminal native work state. Given
seam, output and current, the integer input is recovered by adding the
current at the first endpoint and subtracting it at the second. Thus for
fixed pointer and coin word this map is injective, including zero-current
and idle records. For an independent native work input with any classical
law \(X\), the complete classical/quantum output entropy is
\[
 S(\mathrm{out})-S(\rho)-H(X)=H(B)-I(R:B)+9\log2.              \tag{13}
\]
The local version is \(h(1/50)+9\log2\). Current entropy is included in
the reversible record map; adding it separately would count it twice.
The exact ideal proper-code source has zero complement failures. Equation
(13) concerns that ideal instrument, rather than deleting failures from a
rounded implementation. The native modular certificate retains its own
rounded-control failures, and the separate parity-join certificate bounds
the rounded implementation of (7).

Every local mass and centered contact term commutes with \(D\), while the
hoppings incident at the origin change local parity. Let their sum be
\(H_\star\). Equation (10) proves, without a projected-core replacement,
\[
 W_{\rm read}=\operatorname{Tr}H(\Phi\rho-\rho)
             =-{18\over25}\operatorname{Tr}(H_\star\rho).       \tag{14}
\]
This is work performed by the controlled readout. Its calculation neither
assigns heat to a reservoir nor includes the apparatus preparation/reset
cost. For the full faithful Gibbs state \(\gamma_\beta\) of this same \(H\),
\[
 \beta W_{\rm read}
  -[S(\Phi\gamma_\beta)-S(\gamma_\beta)]
  =D(\Phi\gamma_\beta\Vert\gamma_\beta)\ge0.                    \tag{15}
\]
The random-unitary form of (10) gives the entropy difference nonnegative,
hence work nonnegative at \(\beta>0\). For (1), which need not be Gibbs,
the relative-entropy identity with \(-\log\sigma_q\) holds but cannot be
replaced by \(\beta H\). Pure ground inputs satisfy passivity separately.

## Domain and Einstein boundary

The construction gives a faithful regional state, its native modular word,
an exact constrained entropy relation, and complete-record work/entropy
identities on the same finite radial matter source and its controlled native
response. The supplied spectral-state preparation, radial mass texture,
probe service and addressed controls are explicit inputs. It provides no
geometric area functional, fixed-volume small-ball entropy stationarity,
Lorentzian conserved stress tensor, regional geometric modular inclusion,
continuum interacting limit, or Einstein dynamics. The signed kinetic bias
in (8) and (9) cannot fill the area slot of the conditional Einstein theorem.
No canonical physical premise or acceptance gate is discharged by that
substitution.
