# Mass coefficients and charged response

The native coherent analyzer supplies the family coefficient
\(Y=C_{\rm ex}^{\dagger}C_{\rm ex}/3\), with its original scalar normalization.
The [same-current Higgs map](../common_source/native_matter_attachment.py)
admits this coefficient in the charged-lepton quadratic action. Canonical
kinetic normalization gives the free rest-energy parameters
\(m_i^{(0)}=v_H\sigma_i(Y)\). Squaring these singular values would change the
declared action. A free quadratic response with a specified vacuum is a
spectral realization of these parameters. An interacting charged response
requires a matching theorem for the selected state, clock and observable.

## Finite stability theorem

Let \(M_0\) be a nonsingular complex family mass matrix. Let the positive
Hermitian kinetic matrices and a mass correction obey

\[
 \|Z_L-I\|\leq\epsilon_L<1,\qquad
 \|Z_R-I\|\leq\epsilon_R<1,\qquad
 \|E\|\leq\eta\sigma_{\min}(M_0),\quad 0\leq\eta<1.                 \tag{1}
\]

For constant forms multiplying the complete Lorentz-covariant Weyl kinetic
terms, canonical normalization gives

\[
 M=Z_L^{-1/2}(M_0+E)Z_R^{-1/2}.                                  \tag{2}
\]

Possible transposes from an all-left field convention do not change the
singular-value estimates. Whitening only the time derivative does not prove
that the spatial kinetic terms have the same limiting speed.

For every ordered singular value,

\[
 a\,\sigma_i(M_0)\leq\sigma_i(M)\leq b\,\sigma_i(M_0),\qquad
 a=\frac{1-\eta}{\sqrt{(1+\epsilon_L)(1+\epsilon_R)}},\quad
 b=\frac{1+\eta}{\sqrt{(1-\epsilon_L)(1-\epsilon_R)}}.                \tag{3}
\]

The singular-value min-max principle gives
\(|\sigma_i(M_0+E)-\sigma_i(M_0)|\leq\|E\|\). Since
\(\sigma_i(M_0)\geq\sigma_{\min}(M_0)\), this gives the factors
\(1\pm\eta\). For invertible \(A,B\), the same principle and its application
to \(X=A^{-1}(AXB)B^{-1}\) give

\[
 \sigma_{\min}(A)\sigma_{\min}(B)\sigma_i(X)
 \leq\sigma_i(AXB)\leq\|A\|\|B\|\sigma_i(X).
\]

The extremal eigenvalues of \(Z_L,Z_R\) prove (3), without commutation,
normality or a family-diagonal self-energy assumption. If only
\(\|E\|\leq d\) is known, the lower numerator is
\(\max(0,\sigma_i(M_0)-d)\) and the upper numerator is
\(\sigma_i(M_0)+d\). This absolute statement covers zero coefficients; it
provides no small relative error for a small or vanishing coefficient.

For three positive coefficients define
\(Q(m)=\sum_i m_i/(\sum_i\sqrt{m_i})^2\). Equation (3) implies

\[
 \frac{a}{b}Q(m^{(0)})\leq Q(m)\leq\frac{b}{a}Q(m^{(0)}),\qquad
 \frac{a}{b}\frac{m_i^{(0)}}{m_j^{(0)}}
 \leq\frac{m_i}{m_j}\leq\frac{b}{a}\frac{m_i^{(0)}}{m_j^{(0)}}.      \tag{4}
\]

Intersect the first interval with \([1/3,1]\). A near-balanced input interval
\([q_-,q_+]\) becomes \([(a/b)q_-,(b/a)q_+]\). Common scale factors cancel.
The code encloses the square roots by integer-square tests on dyadic
rationals. It does not use floating-point logarithms or eigenvalues to admit
the error multipliers.

A distinct matching hypothesis may state that three specified observable
parameters \(m_i^{\rm obs}\) obey
\((1-\rho)\sigma_i(M)\leq m_i^{\rm obs}\leq(1+\rho)\sigma_i(M)\).
Then replace \(a,b\) by \((1-\rho)a,(1+\rho)b\) in (4). The helper records
\(\rho\) as a supplied bound. No physical value of \(\rho\), \(E\), or
\(Z_{L,R}\) is computed here. In particular, (3) does not turn a generic
energy-dependent complex self-energy into a positive constant matrix.

## An exact native-coefficient countercontrol

At phase zero the [coherent source](../standard_model_source_map/KOIDE_SOURCE.md)
has

\[
 C=a_0I+b_0(R+R^2),\quad a_0=\sqrt{10/72},\quad b_0=\sqrt{5/72},
 \quad R^3=I,\qquad P_h=(I+R+R^2)/3.
\]

The actual family matrix is
\(Y=\ell(I-P_h)+\ell z^2P_h\), where
\(\ell=(a_0-b_0)^2/3>0\), \(z=4+3\sqrt2\), and
\(\operatorname{tr}Y=5/18\). Its three singular values give \(Q=2/3\).
For any rational \(0<e<1\), take the permitted family-coefficient change

\[
 E_e=((1+e)^2-1)YP_h,\qquad Y_e=Y+E_e.                            \tag{5}
\]

It uses the same charged Higgs invariant line, internal representation and
neutral-background electric charge. The full gauge statement is simultaneous
Higgs covariance. It is an allowed correction, not a calculation of an
interaction in the source. Its root ratios are \(1,1,(1+e)z\), so

\[
 Q(Y_e)-\frac23
 =\frac{e\big((2+e)z^2-8z\big)}{3(2+(1+e)z)^2}>0.                \tag{6}
\]

Here \(z^2=8z+2\) proves strict positivity directly. The certificate takes
\(e=1/10000\), reconstructs the actual coefficient and projector, and keeps
the original scale. Its relative norm is
\(\|E_e\|/\ell=(2e+e^2)z^2<1\).

Kinetic whitening alone also permits a violation. Choose
\(Z_L=I+((1+e)^4-1)P_h\), \(Z_R=I\), and \(E=0\).
The heavy mass coefficient is divided by \((1+e)^2\). For sufficiently
small positive \(e\), its root ratio is \(z/(1+e)>1\), and its Koide ratio
is strictly below \(2/3\). This is an exact positive kinetic-form control.
Neither control identifies an actual source loop correction.

For arbitrary positive masses the differential identity is

\[
 d\log Q=\sum_i\left(\frac{m_i}{\sum_jm_j}
       -\frac{\sqrt{m_i}}{\sum_j\sqrt{m_j}}\right)d\log m_i.        \tag{7}
\]

Its coefficients sum to zero. Family-dependent corrections can change the
ratio even when the common normalization is known exactly.

## Matching to charged particles

The free action and an interacting charged theory have different spectral
questions. In a free theory with separately conserved family numbers, each
family sector can have its own threshold. Within the charged-lepton
subsystem, electric charge alone includes the lightest family. The full
charged-only M15 action also has massless charged spectators outside that
subsystem, so its complete electric-charge threshold may be lower. Neither
threshold supplies three independent charged-lepton masses.

For a stable electric charge with massless electromagnetic radiation, a sharp
mass-shell eigenvector is not a required mass definition. Under the
asymptotic-flux and bounded-fluctuation hypotheses in
[Buchholz's Gauss-law theorem](https://doi.org/10.1016/0370-2693(86)91110-X),
a charged state cannot be a mass-operator eigenvector. A lower charged
energy-momentum boundary with the specified infrared dressing is compatible
with this infraparticle behavior. The theorem does not establish that the
configured finite link source satisfies its continuum hypotheses.

An unstable charged excitation requires a resonance or inclusive line-shape
prescription. A lower threshold in the full electric-charge sector is not
its resonance mass. [Kniehl and Sirlin](https://arxiv.org/abs/0801.0669)
derive gauge-independent pole mass and width definitions for unstable
fermions in perturbative gauge theory, with flavor mixing omitted in their
analysis. Applying a complex-pole prescription requires its analytic and
infrared conditions; no such continuation follows from finite Hermitian
whitening. The electron threshold and muon/tau resonance questions therefore
cannot be replaced by a demand for three isolated real spectral atoms.

[Sumino, Eq. (2)](https://arxiv.org/abs/0812.2090) records the one-loop
modified-minimal-subtraction conversion

\[
 m_i^{\rm pole}=\bar m_i(\mu)
 \left[1+\frac\alpha\pi
 \left(1+\frac34\log\frac{\mu^2}{\bar m_i(\mu)^2}\right)\right].    \tag{8}
\]

The family-dependent logarithm changes Koide's relation; a common
multiplicative term cancels from it. The proposed cancellation uses an
additional family gauge model. It is not an identity of the native repair
instrument or its twelve supplied currents. Equation (8) is a perturbative
matching formula, not a nonperturbative sharp-electron-pole assertion or an
OPH correction computation.

The [same-regulator kinetic attachment](../standard_model_source_map/KINETIC_ATTACHMENT.md)
retains eight sine valleys with both chiralities. Its mirror-current
obstruction and the noncommutation of a bare valley projection with local
gauge charge are unaffected by (3). A quantum raising link cannot be
replaced by a unit background without a state and approximation bound.

A measured-mass comparison consequently needs the following actual data:
the common selected gauge/fermion state and clock; a charged observable with
its infrared prescription and original record weights; a continuum
kinetic/response matching bound, including mirror treatment; and
self-energy, kinetic and observation errors for the stated threshold or
resonance definition. Existing complete coefficient inverses and local
spectral instruments estimate their own finite response. They supply none
of these interacting matching estimates by renormalizing away the other
records. The finite bound (3) converts such estimates into mass-ratio and
Koide error bars when those estimates are proved.

## Verification

`test_mass_response_bridge.py` checks exact saturation, noncommuting left and
right kinetic forms, nonnormal mass perturbations, zero-coefficient absolute
bounds, rational square-root enclosures, and the complete exact native
coefficient in (5). It independently checks the source characteristic
polynomials and the whitening-only countercontrol. JSON certificates are
fresh values, with no empirical masses, physical assignments or computed
interacting self-energies.
