# A localized charged packet close to its calibrated rest energy

## One source and one compact preparation

Use the complete periodic kinetic/analyzer source and its exact native
charged frame from `openwave_lepton_mass_shell`. The physical carrier has
90 fermion modes per site, and the analyzer extension has 324000. Every
character, valley and massless spectator is retained. The chosen reference
fills the negative bands and leaves zero modes empty. There is no claim that
this reference is unique or that the full charge-one sector has a positive
mass floor. The selected electron-compatible character has the exact native
coefficient \(m=e_1>0\), charge minus one, and four-row Hamiltonian

\[
h=-\tau_3\otimes\boldsymbol\sigma\cdot D+m\beta,
\quad \beta=\tau_1\otimes1,\quad
D_j=(T_j-T_j^*)/(2i),\quad h^2=m^2+K^2,
\quad K=-\tau_3\otimes\boldsymbol\sigma\cdot D.
\]

The periodic side is \(N\). Choose an integer \(L\ge4\), with \(N\ge L+5\),
and set

\[
g_j=\sqrt{\frac8{3(L+1)}}\sin^2\frac{\pi j}{L+1}
\quad(1\le j\le L),\qquad g_j=0\quad\hbox{otherwise}.
\]

The product \(f(x)=g_{x_1}g_{x_2}g_{x_3}|\beta+,+z\rangle\) has norm one
and support in the cube \(C=\{1,\ldots,L\}^3\). This is a compact
preparation on the same lattice. No spatial trap or Hamiltonian change is
introduced. A source-defined fermion-addition instrument acting on the
filled reference has success vector \(P_+f\), where
\(P_+=(1+h/E)/2\), \(E=\sqrt{m^2+K^2}\). Its original success is
\(w=\|P_+f\|^2\); the entire quantum failure output is retained. A charged
reservoir and opposite idle reference preserve total charge and parity.
All negative and zero modes are present in the apparatus domain.

The conditional wavefunction \(\psi=P_+f/\sqrt w\) is not exactly compact:
the positive-band filter has a tail, bounded above. Its position effects
are positive-band-filtered native source effects, with the full source
complement retained. No unfiltered local electron detector is identified.

## Exact finite moments and rest-energy bounds

Write \(q=\pi/(L+1)\). Direct finite shifts, including both boundary collars,
give

\[
\begin{split}
d_2=\|Dg\|^2&=\frac{\sin^2(2q)}3+
                 \frac{4\sin^4q}{3(L+1)},\\
d_4=\|D^2g\|^2&=\frac8{3(L+1)}\left[
 \frac L8\sin^4(2q)-\frac12\sin^2(2q)\cos(2q)\sin^2q
 +\frac14\sin^4q\right].
\end{split}
\]

For example, the second difference on \(1\le j\le L\) is
\(-\tfrac12\sin^2(2q)\cos(2qj)\) times the normalizing constant, with
an additional \(\sin^2q/4\) at \(j=1,L\). The four outside collar values
are retained. These give the displayed fourth moment. Product structure
and Pauli anticommutation then give the exact full-dimensional moments

\[
a_2=\langle f,K^2 f\rangle=3d_2,\qquad
 a_4=\langle f,K^4 f\rangle=3d_4+6d_2^2.
\]

They are independent of \(N\) once the collar does not overlap. For
\(L\ge4\), \(\cos(2q)>0\). Dropping its negative contribution and using
\(\sin q\le q\), \(\sin(2q)\le2q\) gives explicit rational upper bounds
when \(\pi\) is enclosed by the alternating Machin arctangent series.
The executable emits these exact fractions. The small-box calculation
independently checks both moments by a complete discrete Fourier transform.

The \(\beta+\) preparation has vanishing expectation of \(K F(K^2)\).
Consequently

\[
w=\frac{1+m\langle E^{-1}\rangle_f}{2},\qquad
1-w\le\eta:=\frac{a_2}{4m^2}.
\]

Assume \(\eta<1\). Since \(f\) vanishes outside \(C\),
\((1-Q_C)P_+f=-(1-Q_C)P_-f\), so

\[
\|(1-Q_C)\psi\|^2\le\frac\eta{1-\eta}.
\]

The source excitation energy above the same filled reference satisfies

\[
0\le\langle E\rangle_\psi-m
 =\frac{\langle K^2/E\rangle_f}{2w}
 \le d:=\frac{a_2}{2m(1-\eta)},\qquad
\langle(E-m)^2\rangle_\psi
 \le b^2:=\frac{a_4}{4m^2(1-\eta)}.
\]

Every inequality remains valid when the native mass is replaced by its
certified lower endpoint in a denominator. These are bounds on a specified
localized packet's energy relative to its own source ground. No binding or
reference energy from a different Hamiltonian is subtracted. Markov's
inequality also bounds the complete energy detector's probability above
\(m(1+10^{-6})\) by \(d/(10^{-6}m)\).

The executed certificate takes \(L=2^{30}\), \(N=4(L+1)\). Its mean energy
exceeds the branch rest energy by at most \(3.014\times10^{-9}\) relatively,
about \(0.003014\) parts per million. The sole energy anchor remains
\(m=0.00051099895069\) GeV. This calibration produces a localized state
with a mean energy interval immediately above that anchor; it is not a
prediction of the electron's absolute mass. The carrier and compiler cost
are printed in full. Its enormous full-volume evolution is not numerically
enumerated; exact finite shift identities and rational inequalities certify
this member. Complete smaller boxes execute the native mass, source
propagator and every position/energy record.

## Persistence, volume and disturbances

The same free source evolves \(\psi_t=e^{-iht}\psi\). Its occupations and
energy law are conserved. Spectral calculus gives

\[
\|\psi_t-e^{-imt}\psi\|\le |t|b,\qquad
\|(1-Q_C)\psi_t\|^2\le
 \left(\sqrt{\frac\eta{1-\eta}}+|t|b\right)^2.
\]

This is persistence of spatial concentration for a finite time, with an
explicit spreading bound. No soliton or shape-minimization premise is used.
For the stated packet and \(|t|\le10^9\), the conditional outside-cube
probability is at most \(1.103\times10^{-7}\). Original success/failure
weights are retained. Doubling \(N\) at fixed \(L\) leaves both local
moment bounds unchanged; the two finite source sizes and their resources
are recorded separately. A smaller-width control has its own bounds.

If the actual prepared vector differs by at most \(\epsilon_0\), and a
Hermitian perturbation of the complete source has total operator norm at
most \(\delta\), Duhamel's identity changes the right-hand amplitude by
at most \(\epsilon_0+|t|\delta\). The executable supplies a finite example
with \(\epsilon_0=10^{-5}\), \(\delta=10^{-14}\). These are total source
norm bounds, not bounds per site; an extensive local perturbation does not
inherit them. All source/reference components remain in the comparison.
No interacting particle stability follows without that interaction bound.

The kinetic compiler's finite full-Fock theorem transfers every complete
instrument uniformly over retained references and ancillary outputs.
For the physical target energy, the 90 modes per site give a kinetic
Fock norm at most 135 per site, and the pairing norm is at most six. Thus
\(0\le H_{\rm phys}-E_G\le282N^3\). Choose the unhalved full-output
trace error at most \(d/(10^6\,282N^3)\), as well as at most \(10^{-12}\).
Then the additional mean physical excitation energy is at most \(10^{-6}d\).
This is conditional only on the initial preparation success; every final
ancillary output is retained, and the energy observable is tensored with
identity on those outputs. The printed compiled interval also keeps its
rest-edge lower endpoint: since \(m<1\), \(E(p)<2\), \(w\le1\) and
\(a_2\ge\sin^2(2q)\), the target excess is at least
\(\sin^2(2q)/4\). Its rational lower bound exceeds the additional
compiler energy error in every stated member. It is not the full enlarged-source energy, and
the ancillary vacuum is not asserted to be its ground. A fixed trace error
alone cannot bound a tiny mean-energy error on a growing carrier. The
printed cycle count uses this stricter budget, the full \(N^3\) volume and
time \(10^9\), and exceeds \(10^{52}\). This is a finite resource bound,
without an efficient implementation claim.

## Autonomous phase records on this preparation

Use an idle opposite-charge reference, with zero Hamiltonian as in the
lepton source. A pointer initially in \(|+\rangle\) controls the success
Kraus operator

\[
K_s=|0\rangle\langle0|\otimes\sqrt w\,1+
 |1\rangle\langle1|\otimes c^\dagger(f)d^\dagger_{\rm ref}.
\]

It is a contraction on the full source/reference carrier. Its failure
operator \(\sqrt{1-K_s^*K_s}\) completes the instrument. On the chosen
reference vacuum and source ground its success probability is \(w\), and
the normalized output is the equal coherent superposition of the source
ground/reference vacuum and this same packet/reference occupied. Each arm
has zero total charge and even total fermion parity. The coefficient \(w\)
is the specified finite spectral expectation above; its preparation and
positive-band services are declared architecture. There is no claim that a
primitive repair sequence supplies them.

Let \(C(t)=\langle\psi,e^{-iht}\psi\rangle\). Project onto the two coherent
return-arm combinations and their full orthogonal complement. Choosing
one of the two quadratures with probability one half gives

\[
p_{a,\pm}=\frac w8\bigl(1+|C|^2\pm2c_a\bigr),\quad
p_{a,\perp}=\frac w4(1-|C|^2),\quad
p_{\rm fail}=1-w,
\quad c_X=\Re C,\quad c_Y=-\Im C.
\]

These seven original-input records sum to one, including every nonreturn
and preparation failure. After the single preparation the declared quadratic target evolves
without a continuing periodic drive. This autonomous statement concerns
\(H_{\rm phys}\). Its finite analyzer phase-cycle implementation uses
continuing control operations and supplies a bounded approximation to
sampled target records, rather than an uncontrolled microscopic clock. Positivity of \(E-m\) gives

\[
|C(t)-e^{-imt}|\le |t|d,\qquad |C(t)|\ge1-|t|d.
\]

For \(|t|d<1\), differentiating the same spectral expectation gives a
uniform instantaneous phase-frequency bound

\[
\left| -\frac{d}{dt}\arg C(t)-m\right|
 \le\frac d{1-|t|d}.
\]

The stated packet retains visibility greater than \(0.9997728\) through
about 12000 rest cycles, with relative frequency error below
\(3.015\times10^{-9}\). The target phase and its energy use the same packet,
Hamiltonian, ground and reference. Trace closeness of its sampled compiler
outputs does not transfer an instantaneous frequency derivative between
control pulses. A reference energy \(\mu\) would instead
give \(C_\mu=e^{-i\mu t}C\); it is not silently subtracted from a measured
beat. The zero reference cost is an explicit inherited construction choice.

There is no passive onset theorem: dephasing the initial arms removes the
quadrature contrast. The source's native differential scalar drive
\(v(t)Q_{\rm source}\) shifts this charge-minus-one packet by minus
\(v(t)\), giving \(C_v=e^{i\int v}C\). Thus this phase is not rigid against
that differential drive. A common scalar gauge transformation of source
and compensating reference cancels, which is a different intervention.
Applying the actual full-generator inverse after evolution gives a unit
return; the inverse uses both spatial inversion and charged-leg parity,
not leg parity alone. These controls separate autonomous prepared phase
coherence from passive self-starting or nonlinear drive resistance.
