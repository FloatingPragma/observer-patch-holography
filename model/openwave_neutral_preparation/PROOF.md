# Local-product preparation of a nearly stationary neutral band

## Source and local input services

Use exactly the periodic eight-cell Hamiltonian
\(H=H_{\rm loc}+2N_{\rm nonblank}+2^{-34}T\) of
`openwave_neutral_stationary`. Its physical raw-lepton-seven, baryon-zero
band has atomic projector \(P=\sum_yP_y\), rank \(r=1024\), interacting
projector \(E\), and \(\|E-P\|\leq\eta<1.900\times10^{-6}\).
The complete local source, link carrier and coefficients are unchanged.

The neutral local line and both ground spin lines are explicit algebraic
projectors of actual sixteen-dimensional physical blocks. For square-free
annihilator \(q\) and selected root \(e\),

\[
 P_e=\frac{1}{q'(e)}\sum_{i=0}^{\deg q-1}H_{\rm loc}^{i}
                      \sum_{j=i+1}^{\deg q}q_j e^{j-i-1}.
                                                               \tag{1}
\]

`local.py` evaluates (1) with the native Gram metric, selects an explicit
normalized basis seed \(b\) with \(p=\langle b,P_eb\rangle>0\), and emits
the exact coefficients of \(v=P_eb/\sqrt p\). Since
\(\langle b,v\rangle=\sqrt p\), the Householder service

\[
 U=I-\frac{2|b-v\rangle\langle b-v|}{2(1-\sqrt p)}
                                                               \tag{2}
\]

maps \(b\) to \(v\) unitarily and extends by identity outside their span.
It preserves the local gauge-singlet sector, raw charges and fermion parity.
All full local Fock/scalar complements are present. The local blocks and
selected roots are known; (2) does not require an unknown interacting
many-cell eigenvector.

Choose the neutral cell uniformly among eight and each of the other seven
ground spin lines uniformly. Eight local services (2), twenty-four supplied
blank link registers and ten retained fair classical bits prepare the
reduced source density

\[
 \rho_0=\frac18\sum_y\frac{P_y}{128}=P/1024.              \tag{3}
\]

The eight native seed states and blank links are declared reservoirs, not
outputs of an arbitrary-state reset theorem. An occupied gauge-neutral
fermionic reference balances the odd source parity and is retained during
the source clock. Each local service has the same parity on its input and
output. No external charged reference is required. The recipe archive is
retained, so the joint source/archive state is not asserted stationary or
decoupled. Equation (3) describes its reduced source state.

## A complete band filter without a full-dimension penalty

The isolated band lies in explicit rational \([a,b]\) with its entire
same-sector complement above \(c>b\); \(U\) bounds \(\|H\|\). Define

\[
 x(z)=\frac{2z-U-c}{U-c},\quad z_a=-x(a),\quad z_b=-x(b)>1,
 \quad f(z)=\frac{(-1)^mT_m(x(z))}{T_m(z_a)},
 \quad M(z)=1-(1-f(z))^K.                                \tag{4}
\]

The joint physical/raw-charge projector is implicit; outside it set \(M=0\).
All eight momenta are retained. The exact Chebyshev growth and positive
flattening estimates from the band-filter theorem apply to (4). The rational
search in `filtering.py` gives \(m=5091\), \(K=2\), degree 10182 and

\[
 \|M-E\|\leq\varepsilon=10^{-16},\qquad \|M\|\leq1.
                                                               \tag{5}
\]

On the band \(M\in[1-\varepsilon,1]\); on the complement its absolute
value is at most \(\varepsilon\). Thus on \({\rm Ran}(P)\),

\[
 A=PM^2P,\qquad q_0 I\leq A\leq I,
 \qquad q_0=(1-\eta^2)(1-\varepsilon)^2.                  \tag{6}
\]

The plain complete filter \((M,\sqrt{I-M^2})\) has success at least \(q_0\)
on (3), rather than a rank-over-full-carrier probability. Since
\(EPE=E-E(I-P)E\) and \(\|E(I-P)E\|\leq\eta^2\), its reduced conditional
source law has all-time stationarity defect at most
\((4\eta^2+8\varepsilon)/(1-\eta^2)\).

## One polar correction gives a uniform band density

Let \(\alpha=1-q_0\). The positive inverse-square-root series gives

\[
 A^{-1/2}=\sum_{n=0}^{\infty}\frac{\binom{2n}{n}}{4^n}(I-A)^n,
 \qquad B_1(A)=\frac{3I-A}{2}\leq A^{-1/2}.
\]

All coefficients are positive and at most one. Therefore

\[
 0\leq I-A^{1/2}B_1(A)\leq\tau I,
 \qquad\tau=\frac{\alpha^2}{1-\alpha}.                 \tag{7}
\]

The complete success operator, extended by zero outside \(P\), is

\[
 L=\sqrt{q_0}\,M P B_1(A)P,
 \qquad L^*L=q_0B_1(A)AB_1(A)\leq q_0P\leq I.           \tag{8}
\]

Its failure operator is \(\sqrt{I-L^*L}\). This defines a trace-preserving
instrument, including arbitrary untouched references and the retained recipe
archive. The lower bound applies to every input vector in \(P\), so

\[
 q_0(1-\tau)^2\leq\Pr(L)\leq q_0,
 \qquad\Pr(L)>0.99999999999639.
                                                               \tag{9}
\]

The decimal lower bound in the receipt is a display; its rational endpoint
is the certificate. No successful branch is dropped from the complete law.

Let \(F\) project onto \({\rm Ran}(MP)\). Equation (6) makes this range
rank 1024. Since \(\|(I-E)MP\|\leq\varepsilon\), equal-rank subspace
geometry gives \(\|F-E\|\leq\varepsilon/\sqrt{q_0}\). The polar map from
\(P\) onto \(F\) is an isometry; (7) bounds every eigenvalue of the
successful density in that frame. For the reduced conditional density
\(\rho_L=L(P/r)L^*/\operatorname{Tr}[L(P/r)L^*]\),

\[
 \left\|\rho_L-\frac Er\right\|_1
 \leq\frac{4\tau}{(1-\tau)^2}+\frac{2\varepsilon}{q_0}
 <2.000001\times10^{-16}.                              \tag{10}
\]

This is unhalved trace norm. Since \(E/r\) commutes exactly with \(H\),

\[
 \|e^{-itH}\rho_Le^{itH}-\rho_L\|_1
 <4.000002\times10^{-16}\quad\text{for every real }t.    \tag{11}
\]

No dephasing time or unresolved tiny many-body gap enters (11). The finite
polynomial polar correction replaces an exact spectral dephasing service.
The statement is about the specified reduced density; an individual recipe
state, its purification or its archived correlations may evolve.

For one successful vector evaluation, (8) uses three evaluations of the
degree-10182 polynomial and a bounded number of known atomic-projector
services. This gives 30546 source-H applications as a formal operator cost.
The linear combination and failure square root are declared finite
Stinespring services. Neither this count nor the eight local Householder
maps supplies a native repair gate compiler, hardware runtime or free work.

## Every record and quantum output

The exact plane/link symmetries of `openwave_neutral_dispersion` commute with
\(P,H,M,A,L\). Each pure recipe is their simultaneous eigenvector. These
symmetries flip the translation characters, so a complete momentum
measurement after success has probability \(1/8\) at each momentum for
every recipe. The original-input zero-momentum success is at least the
lower bound in (9) divided by eight. Selection of that outcome retains the
other seven outcomes.

Resolve the full-H band versus its complement, then the eight atomic
location projectors versus their complete profile complement, then the six
gauge-invariant local occupation outcomes. Let \(q_j\) be success for
recipe \(j\), \(h_j\) its conditional band weight and \(w_j\) the probability
of each atomic location conditional on a momentum and band. Translation
covariance makes the eight location probabilities equal, while the exact
plane/link symmetry makes the same weights apply to all momenta. Then

\[
 h_j\geq1-(\varepsilon/q_0)^2,\qquad
 (1-\eta^2)/8\leq w_j\leq1/8.
\]

The six local probabilities \(p_\gamma\) are the exact native algebraic
values from the neutral projector. The complete original-input law is

\[
 \frac{1-q_j}{1024},\quad
 \frac{q_j(1-h_j)}{1024\,8},\quad
 \frac{q_jh_j(1-8w_j)}{1024\,8},\quad
 \frac{q_jh_jw_jp_\gamma}{1024\,8}.                     \tag{12}
\]

There are 1024 recipe records, eight momentum records, eight cell records
and six occupation records, together with the displayed complements:
410624 complete records in the factorized representation. The symbolic
sum is exactly one. Position and momentum are sequential interventions;
they are not simultaneous sharp properties. The occupation readout can
change energy, and the quantum degeneracies and all references are retained.

The exact source, actual ground and conditional electron-energy conversion
are unchanged. This preparation theorem supplies neither motion on the
exactly flat source nor a physical neutrino/dark-matter interpretation.
