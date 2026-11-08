# Dynamical memory of the selected spin source

The source is the cubic spin-one-half antiferromagnet

\[
 H_J=J\sum_{\langle xy\rangle}(\mathbf S_x\cdot\mathbf S_y-1/4),
 \qquad J>0.
\]

The infinite reference is its infinitesimal staggered-field-selected ordered
ground. The [bounded-history construction](../spin_canonical_field/README.md)
uses its positive ground generator and the single scaling
\(X=px\), \(T=Jpt\). It supplies a regular canonical field algebra, positive
energy and a finite propagation cone. The following construction extracts its
exact response memory through bounded operations. It also evaluates that
memory and a noncentral time history on the finite selected source, with
certified errors. No Lorentz transformation, sharp dispersion or scattering
amplitude is inserted into the generator.

## 1. Exact projected dynamics, including the complement

Let \(K\geq0\) be a selfadjoint ground generator and let
\(a\in\operatorname{Dom}K\) be a unit source-created vector. Write

\[
 \Pi=|a\rangle\langle a|,\quad Q=1-\Pi,\quad
 \epsilon=\langle a,Ka\rangle,\quad r=QKa,\quad D=QKQ|_{Q\mathcal H}.
 \tag{1}
\]

These are fixed by the source and its prepared vector. The bounded
off-diagonal operator \(|a\rangle\langle r|+|r\rangle\langle a|\)
can be subtracted from \(K\). The remaining selfadjoint operator is block
diagonal on \(\mathbb Ca\oplus Q\mathcal H\), since \(a\in\operatorname{Dom}K\).
Thus \(D\) is selfadjoint on \(Q\operatorname{Dom}K\) and positive. This also
specifies the compression when \(K\) is unbounded.

For \(c(t)=\langle a,e^{-itK}a\rangle\), block elimination gives

\[
 \dot c(t)=-i\epsilon c(t)-\int_0^t\kappa(t-s)c(s)\,ds,
 \qquad \kappa(t)=\langle r,e^{-itD}r\rangle,\quad c(0)=1. \tag{2}
\]

Indeed the complementary vector has initial value zero and is
\(-i\int_0^t e^{-i(t-s)D}r c(s)\,ds\). Substituting it in the first block
proves (2); negative times follow by reversal. The complementary dynamics
is retained exactly. Its positive spectral measure has total mass

\[
 \kappa(0)=\|r\|^2=\|Ka\|^2-\epsilon^2. \tag{3}
\]

For real \(z>0\), define the full response and its memory by

\[
 G(z)=\langle a,(K+z)^{-1}a\rangle,\qquad
 \Sigma(z)=\langle r,(D+z)^{-1}r\rangle.
\]

The positive block Schur complement gives

\[
 G(z)=\frac1{z+\epsilon-\Sigma(z)},\qquad
 \Sigma(z)=z+\epsilon-\frac1{G(z)},\qquad
 0\leq\Sigma(z)\leq\frac{\|r\|^2}{z}.                 \tag{4}
\]

For every integer \(j\geq0\),
\((-1)^j\Sigma^{(j)}(z)=j!\langle r,(D+z)^{-j-1}r\rangle\geq0\).
Hence the memory has a positive measure and is completely monotone on the
positive real axis. It vanishes at one positive \(z\) exactly when \(r=0\),
which means the prepared response is one energy eigenvector. Equal mean
energies do not imply equal response: the measures \(\delta_2\) and
\((\delta_1+\delta_3)/2\) have the same mean, while their memories at
\(z=1/4\) are zero and \(4/9\), respectively.

Equation (2) is a linear equation for a response amplitude even when the
underlying source is interacting. Nonzero memory is not a proof of
interacting scattering; a composite probe of a free theory can have it.
For an unbounded Stone field, finite first energy moment does not imply
\(a\in\operatorname{Dom}K\). Section 4 supplies a bounded, operational
regularization rather than dropping this domain requirement.

## 2. Complete readout of the source memory

Consider a finite source with \(0\leq K\leq E_{\max}\) and a normalized
admitted probe vector \(a\). A complete excitation-energy measurement scored
by \(z/(E+z)\) has mean \(s_z=zG(z)\). All detector outcomes remain in
the distribution. A local encoded \(Z\) is a unitary probe, so it needs
no successful-branch renormalization. A probabilistic probe would require
its independently bounded preparation weight throughout the calculation.

Suppose the mean energy is known to error \(\delta_\epsilon\), and the
complete bounded score has error \(\delta_s\). Put
\(g_*=1/(E_{\max}+z)\), \(\delta_G=\delta_s/z<g_*\). Since
\(G(z)\geq g_*\), equation (4) gives

\[
 |\widehat\Sigma-\Sigma|
 \leq\delta_\epsilon+
 \frac{\delta_G}{g_*(g_*-\delta_G)}.                       \tag{5}
\]

For interval inputs \(\epsilon\in[e_-,e_+]\),
\(G\in[g_-,g_+]\), \(g_->0\), the sharper direct interval is

\[
 \max(0,z+e_- -g_-^{-1})\leq\Sigma(z)
 \leq\max(0,z+e_+ -g_+^{-1}).                              \tag{6}
\]

The energy score is Lipschitz with constant \(1/z\). A complete clock
whose label error exceeds \(d_E\) with probability at most \(\ell\)
therefore contributes at most \(d_E/z+\ell\), plus its stated complete
instrument, input and sampling errors. The no-alias bound and summed clock
tail are supplied by [the complete shell readout](../spin_shell_readout/PROOF.md).
A clock bin width alone is insufficient. The source's controlled evolution
and energy-readout service are the declared coherent observer interface;
an autonomous bare repair implementation is not inferred from a matrix.

The finite certificate uses the actual 12-edge open cube, \(J=1\),
staggered field \(B=1/4\), and its unique ground \(\Omega_B\).
The existing integer witness verifies every one of its 256 eigenvalues
and a norm bound on its ground vector before any response is used.
The global bit flip commutes with its Marshall-frame Hamiltonian and
anticommutes with the local \(Z_0\). Positivity and uniqueness of the
ground consequently give \(\langle Z_0\rangle=0\).
The vector \(a=Z_0\Omega_B\) has norm one and contains no ground component.

The actual generator is \(K_B=H_B-E_0\), not its zero-field replacement.
Using the exact integer Hamiltonian and all 256 spectral coefficients,
`finite.py` certifies

\[
 \epsilon\in[1.825861195052,1.825861195070],\qquad
 \|r\|^2\in[0.436013120972,0.436013121316].
\]

At \(z=1/4\) the result is

\[
 G(z)\in[0.532773987262,0.532773987604],\qquad
 \Sigma(z)\in[0.198892679494,0.198892680711].
\]

The finite response has a strictly nonzero complement. Its dynamics is
therefore not the single frequency obtained by replacing its spectrum by
\(\epsilon\). No finite gap or open-cube character is interpreted as a
continuum mass or small momentum.

For clarity, the mean is computed independently as
\(\langle Z_0g,H_B Z_0g\rangle-\langle g,H_Bg\rangle\),
using the normalized integer ground proposal \(g\). The variance is
\(\langle Z_0g,H_B^2 Z_0g\rangle-
 \langle Z_0g,H_B Z_0g\rangle^2\), so the energy origin cancels.
With \(\|H_B\|\leq14\) and ground-vector error \(e_g\), their errors are
at most \(56e_g\) and \(6\cdot14^2 e_g\). The response uses the complete
polar spectral certificate: if its shifted-operator error is \(d\) and
Gram defect is \(\gamma\), its error is at most
\(d/z^2+[\gamma(2+\gamma)+2e_g]/z\). Every displayed interval is rounded
outward after rational arithmetic.

## 3. A noncentral bounded time history

For two source unitaries \(U,V\), define their time-dependent group
commutator and its variance on the reference by

\[
 C(t)=\alpha_t(U)V\alpha_t(U)^*V^*,\qquad
 d(t)=1-|\langle\Omega,C(t)\Omega\rangle|^2.              \tag{7}
\]

This is exactly \(\|(C(t)-\langle C(t)\rangle)\Omega\|^2\).
If \(d(t)>0\), the commutator does not act as a scalar phase even on
the reference. Both real and imaginary expectations are complete two-outcome
Hadamard experiments. Their probabilities are
\((1\pm\operatorname{Re}\langle C\rangle)/2\) and
\((1\pm\operatorname{Im}\langle C\rangle)/2\); no outcomes are discarded.
Controlled source evolution, inverse evolution and local rotations are
the admitted coherent interface of
[the native spin observer](../native_spin_observer/README.md).

On the same selected cube take \(U=Z_0\), \(V=Z_1\), and source time
\(t=1/16\). Both are selfadjoint unitaries. The exact target is

\[
 C(t)=e^{itH_B}Z_0e^{-itH_B}Z_1
       e^{itH_B}Z_0e^{-itH_B}Z_1.                          \tag{8}
\]

For each of its four evolutions use the order-18 Taylor polynomial.
The unitary real-time integral remainder is
\(R=(14|t|)^{19}/19!\). The product error is at most
\(\delta=4R(1+R)^3\); no exponential factor is needed in the single
remainder. The integer Hamiltonian, time and ground vector make the
polynomial expectation an exact complex rational number. Its error
against the true selected-ground expectation is at most
\(\delta+2e_g\). Rational square-root intervals then certify

\[
 d(1/16)\in[0.006279146488,0.006279146491].                \tag{9}
\]

The receipt includes all four interferometer probabilities. The same
operator error gives a complete unhalved trace-norm error
\(2\delta+\delta^2\) for the controlled polynomial diagnostic against
the target isometry, including an arbitrary untouched reference. A physical
unitary compilation has its own synthesis error added; the polynomial
certificate is not itself claimed to be an executed unitary circuit.

For a canonical Weyl algebra with linear symplectic time evolution,
\(\alpha_t(W(f))=e^{i\theta_t(f)}W(S_tf)\), the Weyl relation makes
every such group commutator a scalar phase. Thus a positive limit of (7)
for the actual scaled Weyl probes would exclude linear Weyl dynamics
regardless of whether the state is Gaussian. The finite Pauli probes in
(8) are not these scaled continuum Weyl probes. Equation (9) establishes
noncentral source dynamics; it does not assert that this noncentral part
survives the specified field scaling, or that noncentrality alone proves
scattering.

## 4. The memory belongs to the same bounded-history continuum

Let \(A_p=\sin(sQ_p(F))\), or \(\sin(sP_p(G))\), with fixed real
\(s\ne0\), and use \(L_p=K_{\rm source}/(Jp)\).
These operators of norm at most one are differences of the source Weyl
unitaries. On the common subsequence defining the actual field histories,
their spectral measures \(\mu_p\) converge weakly to that of
\(A\Omega\), where \(A=\sin(sQ(F))\) or \(\sin(sP(G))\).
This follows from convergence and continuity of their time characteristic
functions, including their total mass at time zero.

The inherited energy theorem proves \(Q(F)\Omega\ne0\) and
\(P(G)\Omega\ne0\) for the admitted nonzero annular tests. Since
\(\sin(sA_0)\Omega/s\to A_0\Omega\) for either Stone field \(A_0\),
some positive rational \(s\) gives \(A\Omega\ne0\). This supplies a
nonzero bounded probe without an assumed pole or a selected energy tube.
For any such \(s\) and any fixed \(\lambda>0\), injectivity of
\((K+\lambda)^{-1}\) gives

\[
 n_\lambda=\int(E+\lambda)^{-2}\,d\mu(E)>0.
\]

Then the actual normalized filtered response is

\[
 b_\lambda=
 \frac{(K+\lambda)^{-1}A\Omega}{\sqrt{n_\lambda}}
 \in\operatorname{Dom}K.                                \tag{10}
\]

It is the normalized successful branch of the source resolvent filter,
whose original probability is \(\lambda^2n_\lambda>0\). All complementary
records retain total probability \(1-\lambda^2n_\lambda\), as constructed in
[the relaxed response](../spin_relaxed_response/README.md). At finite
source scale the regularizer is \(Jp\lambda\); the clock is unchanged.
The source norm and moments converge to

\[
 \begin{split}
 n_\lambda&=\int\frac{d\mu(E)}{(E+\lambda)^2},\\
 \epsilon_\lambda&=
 \frac1{n_\lambda}\int\frac{E\,d\mu(E)}{(E+\lambda)^2},\\
 h_{2,\lambda}&=
 \frac1{n_\lambda}\int\frac{E^2\,d\mu(E)}{(E+\lambda)^2},\\
 G_\lambda(z)&=
 \frac1{n_\lambda}\int\frac{d\mu(E)}{(E+\lambda)^2(E+z)}.
 \end{split}                                                   \tag{11}
\]

Every integrand is bounded and continuous on \([0,\infty)\).
Consequently these limits require neither uniform integrability of the
unfiltered first moment nor an unbounded-field projection identification.
For any fixed verified lower bound on \(n_\lambda\), interval arithmetic
propagates the unnormalized readout errors through (11). Applying (1)--(4)
to (10) gives the actual continuum response memory, as an ordered limit
of the same source, state prescription and clock. It does not substitute
a Hilbert space with an imposed free dispersion.

For each fixed filtered response, the values \(G_\lambda(z)\), \(z>0\),
uniquely determine its entire probability measure. To see this, write
\((E+z)^{-1}=\int_0^\infty e^{-t(E+z)}dt\). Equality of all positive-axis
resolvents gives equality of the Laplace transforms of the bounded
continuous functions \(\int e^{-tE}d\nu(E)\), and therefore equality of
those functions. Uniqueness of the finite-measure Laplace transform gives
equality of the original spectral measures. This is an identification
from the complete response, without selecting eigenvalues by analogy
with a desired particle list. It supplies no finite-shot reconstruction
of an arbitrarily narrow spectral feature.

Likewise, every expectation in (7) for fixed scaled Weyl controls and
times is a bounded-history expectation. Its limit therefore exists on
that same subsequence and its variance is the limit of the source
variances. A strictly positive verified limit would establish noncentral
continuum time evolution on those fields. No nonzero lower bound on this
limit is derived by the finite cube experiment.

These statements are constructive for every fixed finite source and
filter accuracy through the admitted coherent controls. They supply no
effective thermodynamic box-size or field-selection modulus. They also
do not evaluate the fixed-\(\lambda\) continuum measures in (11), permit
interchange with \(\lambda\downarrow0\), or establish a sharp pole.
The ordered-source current/force residual in
[the shell diagnostic](../spin_shell_diagnostic/README.md) is the relevant
source correlation for that spectral problem. The finite causal-poset
theorems constrain operational order after their source placement and
clock hypotheses; they do not remove a nonzero memory or prove boost
covariance of these histories. Lorentz transformations of the full field
algebra and nonzero continuum scattering require additional source proofs.

## 5. A uniform source estimate and its exact limitation

There is a quantitative bounded probe whose filter probability is uniformly
positive in the actual scaling limit. The source-parity theorem in
[the defect catalogue](../protected_defect_catalog/PROOF.md) supplies,
for a real test \(F\) of unit \(L^2\) norm with Fourier support
\(1/2\leq|q|\leq1\), the odd bounded probe

\[
 s=2/m_0,\qquad A=\sin(sQ(F)),\qquad 0<m_0\leq m.
\]

That theorem uses the same Weyl relation, positive-energy representation
and source symmetry to prove

\[
 \|A\Omega\|^2\geq n_0:=\frac{25m_0^2}{2304}.           \tag{12}
\]

The complete instrument retains \(\sin(sQ_p(F))\) and
\(\cos(sQ_p(F))\); the first branch is odd and the second is even under
the exact source parity. Both are linear combinations of opposite Weyl
controls, and their effects sum to the identity. The inherited
bounded-history limit gives
\(\liminf_p\|A_p\Omega\|^2\geq n_0\).
The odd vectors are orthogonal to the even reference. Orthogonality to all
of \(\ker K\) requires a separate ground-space statement; this bound
alone is not a positive-energy-window certificate.

The source energy-kick estimate supplies

\[
\limsup_p\|L_p^{1/2}A_p\Omega\|^2
 \leq h_0:=3s^2=12/m_0^2.                                \tag{13}
\]

For any positive measure of mass \(n>0\) and first moment \(h\), Jensen's
inequality for \((E+\lambda)^{-2}\) gives

\[
 n_\lambda\geq\frac{n^3}{(h+\lambda n)^2}.               \tag{14}
\]

The right side increases with \(n\) and decreases with \(h\). Equations
(12)--(14) therefore give the actual-source lower bound

\[
 \liminf_p n_{p,\lambda}
 \geq\frac{n_0^3}{(h_0+\lambda n_0)^2}>0.                \tag{15}
\]

For all sufficiently small \(p\), an explicit uniform lower constant is
obtained by replacing \(n_0\) by \(n_0/2\) and \(h_0\) by \(2h_0\)
in (14). Call it \(b_\lambda>0\). Since \(\|A_p\|\leq1\),
\(E^2/(E+\lambda)^2\leq1\) and (4) imply

\[
 \lambda^2 n_{p,\lambda}\geq\lambda^2 b_\lambda,
 \qquad 0\leq\Sigma_{p,\lambda}(z)
 \leq\frac1{z b_\lambda}.                                \tag{16}
\]

This is a uniform memory bound and a nonvanishing original preparation
probability on the actual source. It supplies no effective threshold in
\(p\). Its upper bound is of order one, not a vanishing force estimate.
Replacing it by a vanishing estimate would require additional information
about the selected ground correlations.

An explicit local-field countercontrol shows why positivity, locality,
canonical commutators and even equality of the leading first moments do
not force a single Lorentz cone. This control is not a state or a limit
of the OPH spin source. Take two independent free scalar fields on
\(\mathbb R^3\), with canonical pairs \((q_i,\pi_i)\), and Hamiltonian

\[
 H=\frac12\sum_{i=1}^2\int
     [a\pi_i^2+b_i|\nabla q_i|^2],d^3x,
 \quad m=\kappa=1/4,\quad a=12,\quad(b_1,b_2)=(1/2,3/2).
\]

Define \(Q=\sqrt{m/2}(q_1+q_2)\),
\(P=\sqrt{m/2}(\pi_1+\pi_2)\). Then \([Q(F),P(G)]=im\int FG\).
The two positive frequencies are \(c_i|q|\),
\(c_1^2=6\), \(c_2^2=18\). In the product vacuum, their summed
quadrature measures have coefficients

\[
 \begin{array}{c|cc}
 &Q&P\\ \hline
 \text{mass}&\dfrac{ma}{4|q|}(c_1^{-1}+c_2^{-1})
            &\dfrac{m|q|}{4a}(c_1+c_2)\\[4pt]
 \text{first energy moment}&6\kappa&\kappa|q|^2/2\\[4pt]
 \text{inverse energy moment}&1/(6|q|^2)&1/96.
 \end{array}                                               \tag{17}
\]

These obey all the quoted canonical-field moment upper bounds and the
allowed \(m^2/2\leq\kappa\leq1/2\). The first moments equal the source
form, so there is no escaping-energy defect in this control. It has a
positive generator, continuous spatial rotations and translations, and
local commutation outside speed \(\sqrt{18}\). This also lies within
the inherited conservative \(\ell^1\) cone of speed 72.
Nevertheless its \(P\)-response has two different energy/momentum ratios,
with exact ratio-moment defect \(7-4\sqrt3>0\). No single Lorentz speed
preserves both nonzero cones under boosts: preservation of
\(E^2-c_i^2|q|^2=0\) requires that same boost speed to equal \(c_i\).
Its group commutators are scalar and its scattering is zero despite the
nonzero projected memory. Thus neither memory alone nor the quoted
kinematic and moment conclusions establish the desired interacting
relativistic dynamics. Equation (16), together with the source-specific
filtered force/current correlations, is the proved uniform estimate.

## 6. Removing the selection torque on the same finite cube

The nonzero finite memory is not solely an effect of the applied staggered
field. In the Marshall frame the exact source matrices obey
\(H_B=H_0-(B/2)\sum_{x=0}^7X_x\). The code forms \(H_0\) by removing
this term from the source's integer matrix; it does not change any bond.
The complete 256-vector zero-field witness is an untrusted rounded proposal.
Its exact Gram defect, full Frobenius spectral residual and polar comparison
certify every eigenvalue and an isolated lowest eigenvalue. This proves
uniqueness of the finite zero-field ground directly; positivity of all
its basis coordinates is not required.

For the same local \(Z_0\) probe and reference-subtracted zero-field clock,
the certified values are

\[
 \begin{split}
 E_0&\in[-7.8200893743749,-7.8200893743747],\\
 \operatorname{gap}(H_0)&\in[0.820089374367,0.820089374383],\\
 \|QK_0Z_0\Omega_0\|^2
  &\in[0.608562479916,0.608562480007],\\
 \Sigma_0(1/4)&\in[0.314093786867,0.314093787470].
 \end{split}                                               \tag{18}
\]

Here \(Q\) in the third line projects off the normalized probe vector,
as in (1); it is not a spin-charge operator. Every level contributes to
the response. The exact spin-current identity at zero field contains no
selection-field torque, while the response memory in (18) is positive.

There is also a quantitative same-box field-removal limit. Write
\(\delta=4B\), so \(\|H_B-H_0\|\leq\delta\), and let the certified
zero-field gap lower bound be \(g\). If \(2\delta<g\), Weyl's inequality
keeps the lowest level simple. Choosing ground phases to have nonnegative
overlap gives

\[
 \|\Omega_B-\Omega_0\|\leq e_B:=\frac{2\delta}{g-\delta}.
\]

Indeed the component orthogonal to \(\Omega_0\) is bounded by
\(\delta/(g-\delta)\), by applying the zero-field resolvent to the
ground eigenvalue equation; the unit-vector distance is at most twice
that component. Since \(|E_B-E_0|\leq\delta\),
\(\|K_B-K_0\|\leq2\delta\). Resolvent comparison and the norm-one
probe therefore give the bounds below. Each zero-field bond is minus a
singlet projector, up to the fixed Marshall conjugation, so
\(-12I\leq H_0\leq0\). For the mean-energy bound compare the two
expectations of \(H_0\), with norm at most 12, and separately subtract
the two ground energies. This contributes \(24e_B\) and \(2\delta\).

\[
 |G_B(z)-G_0(z)|\leq\frac{2\delta}{z^2}+\frac{2e_B}{z},
 \qquad |\epsilon_B-\epsilon_0|\leq2\delta+24e_B.       \tag{19}
\]

Both resolvents are at least \(g_z=(24+2\delta+z)^{-1}\), so

\[
 |\Sigma_B(z)-\Sigma_0(z)|\leq2\delta+24e_B+
 \frac{2\delta/z^2+2e_B/z}{g_z^2}.                       \tag{20}
\]

At \(B=2^{-40}\), \(z=1/4\), the rational bound is less than half
the positive lower endpoint in (18). It decreases to zero with \(B\).
This is an actual fixed-volume source limit, not a moment countermodel.
Taking \(B\downarrow0\) on one cube is different from first taking the
thermodynamic limit and then selecting an ordered ground. Consequently
(18)--(20) do not fix the ordered infinite source's low-momentum memory
or replace the missing ground-specific correlation estimate.

## 7. Memory inside each conserved source charge

The [source charge construction](../protected_defect_catalog/PROOF.md)
supplies the same selected cube's exact conserved internal charge
\(Q_s=\sum_x\epsilon_xX_x/2\), ground charge zero, and
\(Z_0=A_++A_-\) with \([Q_s,A_\pm]=\pm A_\pm\).
The native full-Fock attachment is replayed, including all 6,976 branch
columns. No convenient charge is chosen after observing the energies.

Let \(w_\pm=\|A_\pm\Omega_B\|^2\) and
\(a_\pm=A_\pm\Omega_B/\sqrt{w_\pm}\).
Both weights have certified lower bound greater than \(4/25\).
The complete preparation instrument has exactly the two Kraus operators
\(A_q=Z_0(I-q\epsilon_0X_0)/2\), \(q=\pm1\). It measures only the
local \(X_0\), then applies \(Z_0\), and retains both outcomes. These
are admitted local pointer and rotation operations on the same source;
their exact matrices obey \(\sum_q A_q^*A_q=I\), including arbitrary
references. The program does not claim to execute that complete native tape.
On the charge-zero ground the two branches prepare precisely the stated
coherent vectors. A separate 256-record local-X measurement gives the
charge weights but dephases within a charge sector. Grouping its labels
does not undo that backaction, so it is not the preparation used here.
Since the charge commutes with \(K_B\), its two coherent response
measures are orthogonal under all time evolution. In particular

\[
 G(z)=w_+G_+(z)+w_-G_-(z),\qquad
 \epsilon=w_+\epsilon_++w_-\epsilon_-,
\]

\[
 \operatorname{Var}_{Z_0\Omega_B}(K_B)
 =w_+\operatorname{Var}_{a_+}(K_B)
  +w_-\operatorname{Var}_{a_-}(K_B)
  +w_+w_-(\epsilon_+-\epsilon_-)^2.                       \tag{21}
\]

Applying (1)--(6) to the same generator in each channel gives

| Source charge | Original weight | Conditional memory strength | \(\Sigma_q(1/4)\) |
|---|---|---|---|
| \(-1\) | \([0.835475070930,0.835475070931]\) | \([0.338256843285,0.338256844037]\) | \([0.169680040792,0.169680042091]\) |
| \(+1\) | \([0.164524929069,0.164524929070]\) | \([0.921555893760,0.921555895454]\) | \([0.314469337260,0.314469338280]\) |

Thus the full response memory is not solely a mixture of two different
charge means. Each exact charge channel has its own strictly positive
complement memory on the actual selected source and clock.

The interval calculation normalizes the integer proposal image explicitly.
If the original ground-vector error is \(e_g\), its error as a normalized
branch is at most \(2e_g/\|A_qg\|\), because \(\|A_q\|\leq1\).
The proposal norm is an exact rational squared norm, bounded away from
zero. Mean, variance and resolvent bounds therefore use this branch error
instead of the unconditioned ground error. All 256 spectral coefficients
enter each certificate, even though the exact vector has a sharp charge.
The original complete scored probability is \(w_q zG_q(z)\); its
complement retains the other charge and all remaining records. Conditional
normalization does not erase the herald weight.

The charge is an internal source-spin quantum number. A positive finite
variance is neither a continuum decay width nor a physical particle mass.
The table gives two measured source channels, without identifying them
as Standard Model species or assigning physical electromagnetic charge.
