# A controlled observable limit of the local carrier walk

The source is the declared local twelve-slot walk W=SR on cubic-routed
icosahedral cells. R reflects the uniform local event vector and S is the
holding or neighboring-port involution. The exact Bloch matrix W(k) has
dimension 144: twelve ports and twelve event labels per cell. Source and
event padding, native export complements and reset outcomes are retained by
the compiled instrument. They are not assigned a conditional success law.

The source's cell Laplacian H(k) obeys

\[
A^*S(k)A=I-H(k)/12,
\qquad \lambda_0(k)=|k|^2/8+e(k),
\qquad |e(k)|\le3\|k\|_\infty^4/8
\]

for \(\|k\|_\infty\le1/2\). The inherited rational Schur bound establishes
this inequality on all twelve internal components. Its zero-momentum
eigenvalues are 0, 4, 6 and 8 with multiplicities 1, 3, 2 and 6. The local
denominator is twelve; the global-event denominator is a different source.

## Source-derived continuum generator

Let P be the orthogonal projector onto ker(W(0)−I). It has rank 26. It is
computed exactly as a rational polynomial in C=(W(0)+W(0)*)/2:

\[
P=\prod_{a\in\{-1,2/3,1/2,1/3\}}\frac{C-aI}{1-a}.
\]

Every asserted projector identity is checked over the rationals. Write
\(W'_a(0)=iD_aR\), where \(D_a\) is the signed external-edge derivative.
The three effective Hermitian matrices are \(H_a=P D_aR P\).

Let e be the uniform vector on all 144 internal/event labels, with components
1/12. Put \(g_a=H_ae\). Exact arithmetic gives

\[
H_a=g_ae^*+eg_a^*,\quad e^*g_a=0,\quad
g_a^*g_b=\delta_{ab}/48.
\]

Thus \(b_a=\sqrt{48}g_a\) and e form an orthonormal four-dimensional space.
The other 22 vectors in P are stationary under every effective generator.
At momentum p, the generator on P is

\[
K(p)=i\sum_a p_a H_a.
\]

Its two nonzero eigenvalues are \(\pm i|p|/\sqrt{48}\); the other 24 are
zero. Both acoustic signs and the stationary modes are part of the limit.
The matrices above are derived from the actual carrier permutation and
reflection, rather than installed as an execution input.

Use the shared clock interface
`clock/clock_contract.py::effective_time(n,h,tau0)`, which declares
\(t=n\tau_0h\), with cell spacing h and one positive scale \(\tau_0\).
In the basis \((e,b_1,b_2,b_3)\), the effective fields satisfy

\[
\partial_t\phi=\frac{1}{\sqrt{48}\tau_0}\nabla\cdot\chi,
\qquad
\partial_t\chi=\frac{1}{\sqrt{48}\tau_0}\nabla\phi.
\]

Consequently \(\partial_t^2\phi=\Delta\phi/(48\tau_0^2)\). The transverse
part of \(\chi\) is stationary. These equations do not supply Maxwell radiation.
The vector \(\chi\) is an amplitude coordinate, rather than a native edge-current record or electromagnetic field. Their conserved density is \(|\phi|^2+|\chi|^2\), with the additional stationary
components included when prepared. Its probability flux is \( -2c\,\mathrm{Re}(\bar\phi\chi)\), where \(c=1/(\sqrt{48}\tau_0)\). A shared physical-second attachment is
an additional hypothesis of the clock construction.

## Native displacement readout

For axis a, the matrix \(D_aS(0)\) is diagonal: it is +1 on the two positive
outgoing external endpoints with moving event label 11, −1 on the two negative
endpoints, and zero elsewhere. After the source transport, the terminal
signed displacement has observable \(J_a=-D_aS(0)\). Transport from a positive
outgoing port to its paired negative port increments that cell coordinate. The actual
source/event record determines this label through the native inverse recoder.

Exact rational arithmetic verifies

\[
PJ_aP=-H_a.
\]

Consequently the first moment of this native signed displacement on the
continuum amplitude state is \(-2\operatorname{Re}(\bar\phi\chi_a)/\sqrt{48}\).
Multiplying a displacement by cell spacing h and dividing by the common
step duration \(\tau_0h\) gives the probability-flux coefficient above.
This identifies a bilinear amplitude observable with a native readout.
It does not identify the vector amplitude itself with an individual edge
record, a material charge current or an electromagnetic field. The complete
terminal-instrument bound controls its source-to-limit measurement error:
if that total variation bound is \(\epsilon\), the signed mean differs by
at most \(2\epsilon\), and its velocity mean by at most
\(2\epsilon/\tau_0\).

## Uniform complete-state bound

Let Fourier support satisfy \(|p|_\infty\le K\), with \(hK\le1/2\).
Initial states lie in P and may be entangled with an arbitrary unchanged
reference. This includes a bandlimited scalar envelope multiplying e.
It does not include an arbitrary one-cell delta packet.

Put \(Q=I-P\), \(W_1(p)=\sum_a p_aW'_a(0)\), and

\[
X(p)=-[(W(0)-I)|_Q]^{-1}QW_1(p)P,\qquad J_h=P+hX.
\]

The smallest distance from 1 to another eigenvalue of W(0) is
\(\sqrt{2/3}\). Since the external ends partition the ports,

\[
\|W_1(p)\|\le K,\quad
\|X\|\le\sqrt{3/2}K\le5K/4,\quad
\|K(p)\|=|p|/\sqrt{48}\le K/4.
\]

The second-order remainder of each external phase obeys
\(\|W(hp)-W(0)-hW_1(p)\|\le h^2K^2/2\). The defining equation for X is

\[
W(0)X+W_1P=X+K(p).
\]

Expanding one step and using the integral remainders of the unitary
exponential yields

\[
\begin{aligned}
\|W(hp)J_h-J_he^{hK(p)}\|
&\le h^2K^2\left(\frac54+\frac12+\frac1{32}
                  +\frac5{16}+\frac58hK\right)\\
&\le\frac{77}{32}h^2K^2.
\end{aligned}
\]

Telescoping n steps uses unitarity and retains every output component:

\[
\boxed{\|W(hp)^nP-Pe^{nhK(p)}\|
\le\frac52hK+\frac{77}{32}nh^2K^2.}
\]

The same bound holds uniformly on a Fourier direct sum or direct integral,
and after tensoring an arbitrary reference. For normalized states, their
trace distance is bounded by the displayed vector bound, capped at one.
Every complete terminal instrument therefore has total variation distance
at most this bound. No optical amplitude or stationary mode is discarded
from the actual evolved state. At fixed t, the bound is
\(h[(5/2)K+(77/32)(t/\tau_0)K^2]\), which tends to zero.

On a periodic unit cube, take finitely many momenta \(p=2\pi m\). The
Fourier transform gives the actual finite carrier state. The fixed internal
embedding by e and the three b vectors maps the continuum fields to native
port/event amplitudes. Coarse cell-position records compare the carrier
probabilities with the sampled continuum density, including all prepared
stationary components. The bandwidth condition prevents sampling aliases.
Fine port/event records use the corresponding fixed component projectors.

To compare with continuum probability integrated over centered cells, rather
than sampled density, average a grid translation by q over
\([-h/2,h/2]^3\). Fourier sampling is isometric for every such translation.
Its state-vector change is at most \(K\sum_a|q_a|\), so the total variation
change in any cell/component probability is bounded by the same quantity.
The average of \(\sum_a|q_a|\) is \(3h/4\). Thus continuum cell-integrated
position probabilities obey the preceding terminal bound plus \(3hK/4\).
The comparison keeps a common finite outcome space of detector cells.

A source instrument between evolution intervals can leave P and broaden
Fourier support. An arbitrary sharp position measurement therefore does
not satisfy this theorem's preparation premise for the next interval.
The exact full-state same-W instrument executes that continuation, including
the flat modes. A multi-interval continuum claim needs a bound for each
intervening operation's failure to preserve the admitted domain.

## Scalar spectral and optical bounds

For a scalar envelope multiplying e, uniform internal preparation requires
no spectral filtering. With \(r=hK\), its optical weight is at most

\[
\eta=\frac{r^2}{(4-r-r^2/2)^2}.
\]

If \(\theta=\arccos(1-\lambda_0/12)\), then

\[
\left|\frac{\theta^2}{h^2}-\frac{|p|^2}{48}\right|
\le h^2K^4\left[\frac1{16}
 +\frac1{576(1-r^2/48)}\right].
\]

For p=0 the phase error is zero. For other p, dividing by the sum of the
two nonnegative frequencies and using \(\sqrt{48}<7\) gives the uniform
n-step phase bound

\[
\epsilon_\theta=7nh^3K^3\left[\frac1{16}
 +\frac1{576(1-r^2/48)}\right].
\]

The unmeasured return amplitude is the complete spectral average of
\(\cos(n\theta_j)\). Its difference from the scalar continuum cosine is
at most \(2\eta+\epsilon_\theta\). The X-probe total variation bound is
\(\eta+\epsilon_\theta/2\); the Y probe has probabilities one half exactly.
Both acoustic branches enter the cosine. For the normalized modewise low-band embedding, choose each eigenvector's
phase so that its overlap with the uniform internal vector is nonnegative.
That overlap is at least \(\sqrt{1-\eta}\). Keeping the same normalized
Fourier amplitudes, including any reference, gives initial trace distance at
most \(\sqrt\eta\). The trace-distance comparison after evolution is
therefore at most \(\sqrt\eta+\epsilon_\theta\). The vector-norm bound
is instead \(\sqrt{2-2\sqrt{1-\eta}}+\epsilon_\theta\), and the executable
uses its rational upper bound \(2\sqrt\eta+\epsilon_\theta\). This
separate diagnostic does not project the full-state theorem onto a prepared
eigenbranch.

## Execution and scope

The executable checks the exact projector and gradient-vector identities,
the corrector equation, complete 144-component operator errors and rejecting
domain/sign controls. The Bloch action is independently compared with the
actual graph at development sizes 2, 3 and 4. Formal refinement bounds on
larger sizes are evaluations of rational inequalities, not acceptance runs
or enumerations of those global carriers.

The native compiler supplies the size-two member and its complete apparatus.
A uniform native-pulse implementation of every refined member is a separate
construction obligation. Its physical errors cannot be inferred from this
ideal source theorem. Source-word, preparation, reset and detector bounds
compose additively with the complete-observable bound on their admitted
domains. Supplied control rates and service latencies retain their ancestry.

For the frozen size-two member, `SOURCE_ERROR.json` certifies an operator
error at most \(6.5592856182\,10^{-11}\) per source word and source-preparation
error at most \(1.1113827752\,10^{-14}\). The executable composes n times the
word bound with the preparation and explicitly supplied instrument and reset
bounds. It rejects an attempt to apply that certificate to another size.
At fixed effective time, n grows as 1/h. A physical implementation limit
therefore requires per-step error \(\delta_W(h)=o(h)\), together with vanishing
preparation, detector and total timeout errors. A fixed finite-precision word
bound alone does not establish that convergence.

This result concerns the free one-load dynamics. A common many-body extension,
autonomous vector field, interacting continuum limit, particle species,
masses and force measurements require their corresponding same-source
construction and records.
