# Complete observation histories and the necessary phase strength

Use the selected open Wilson source, Gauss-polynomial density, idle
conjugate reference and Hamiltonian clock of the
[primary contract](../particle_common_contract/PROOF.md). Write
\(a=2^{-n}\), \(R=1/a\), \(L=2R^2\), and
\(K=H-E_0\ge0\). The selected mean satisfies
\(\operatorname{Tr}\rho K\le\epsilon=a/16\). Ground subtraction occurs
only in the estimates and does not change the source clock or records.

This construction proves a restriction on an actual family of complete
observations. It does not establish collapse or survival of every source
observable. In particular the original local phase probes are outside
the vanishing-strength domain proved below.

## 1. The normalized local density on the actual carrier

For the retained classical origin \(X\), use the physical-radius-one tent
\(f_R\) from the [local phase construction](../particle_primary_local_fluctuation/PROOF.md).
Define
\[
 Q_X={1\over30R^3}\sum_y f_R(y-X)(N_y-30I).
 \tag{1}
\]
Since \(\sum_y f_R=R^3\) and \(\|N_y-30I\|\le30\),
\(\|Q_X\|\le1\). Number phases preserve the link cutoff and every
Gauss constraint. Only apparatus coefficients wrap on the address list;
the source has no wrapping bond.

For an actual hopping bond \(T_b\), including its \(1/a\), put
\(d_{Xb}=\Delta_b f_R/(30R^3)\). The literal source identity is
\[
 [Q_X,[H,Q_X]]=-\sum_b d_{Xb}^2T_b.
 \tag{2}
\]
Electric, magnetic and onsite terms commute with \(Q_X\). The full-Fock
bond bound \(\|T_b\|\le30/a\) includes all fifteen internal rows.
The exact one-axis squared-gradient sum is
\(D_R=2(2R^2+1)^2/(9R^3)\). For \(R\ge2\), absent wrapping bonds
only reduce the coefficient sum, so
\[
 \|[Q_X,[H,Q_X]]\|\le C_n={D_R\over10aR^6}
 ={a^4(2+a^2)^2\over45}\le{9\over80}a^4.
 \tag{3}
\]
At member zero the support is one corner with three actual bonds, giving
\(C_0=1/10\), also below the uniform bound. These are full-carrier
operator estimates, valid on every origin block and arbitrary reference.

## 2. Complete coherent circuits

A phase call of strength \(r\) is
\[
 V_X(r)=\exp(i\pi r Q_X\otimes A),\qquad A=A^*,\quad\|A\|\le1,
 \tag{4}
\]
with \(A\) acting on apparatus registers. Rational \(r\) and rational
control eigenvalues give finite roots-of-unity phase gates on the actual
occupation basis. This includes deterministic signs, controlled signs
and controlled no-operation branches. Direct-sum control of the origin
has the same bound. Apparatus-only gates can change the control basis
between calls.

For the joint nonnegative energy \(K\otimes I\), the double commutator
of the generator in (4) is
\(\pi^2r^2[Q_X,[H,Q_X]]\otimes A^2\). The semibounded energy-form
lemma in [common phase histories](../particle_primary_current_tail/PROOF.md#2-semibounded-energy-controls-the-unitary-form-domain)
and \(\pi^2<10\) therefore give, on every joint vector,
\[
 \|(K^{1/2}\otimes I)V_X(r)\psi\|
 \le\|(K^{1/2}\otimes I)\psi\|+{3\over4}|r|a^2\|\psi\|.
 \tag{5}
\]
Arbitrary apparatus-only channels and complete pointer measurements
preserve the source energy mean. Source evolution, controlled source
evolution and ideal source energy measurement also preserve it. Thus a
complete coherent protocol with phase-call strength bounds \(r_j\),
including coherent or classical feedback implemented by these controls,
obeys
\[
 S=\sum_j|r_j|,\qquad
 \operatorname{Tr}\rho_{\rm out}K
 \le e_n(S):=(\sqrt\epsilon+3a^2S/4)^2.
 \tag{6}
\]
The proof is the Hilbert--Schmidt triangle inequality applied to (5),
with all apparatus and reference registers retained. A call's strength
is the largest absolute controlled strength in that call. No branch
probability is divided out. A complete cosine/sine readout is included:
prepare a balanced control, apply (4) with \(A=Z\), then read it in the
appropriate basis. Its two source Kraus operators are
\(\cos(\pi rQ_X)\) and \(i\sin(\pi rQ_X)\). Reading the sign basis
instead retains the two signed phase operations. Both instruments have
all outcomes. Temporal filtering by coherent source evolutions and such
phase calls is within (6), including its failure outputs.

Measure source energy initially and finally on this complete protocol.
The initial readout preserves the stationary source marginal, while its
reference coherences can change. If the original transfer exceeds
\(\eta>0\), final ground-relative energy exceeds \(\eta\), because
initial ground-relative energy is nonnegative. Hence
\[
 \Pr_{\rm original}(E_f-E_i\ge\eta)
 \le\Pr_{\rm original}(E_f-E_0\ge\eta)
 \le\min\{1,e_n(S)/\eta\}.
 \tag{7}
\]
This uses the same finite nonnegative-energy principle as
`Lean/QFT/GroundRecordTail.lean`, with the final mean supplied by (6).
For every strength family satisfying \(a^2S_n\to0\), (7) tends to
zero at every fixed positive threshold. A declared fixed SI/GeV clock
conversion preserves that statement. A capacity-dependent clock change
would define a different represented interval.

## 3. Full-reference waiting-time indistinguishability

Let a protocol prefix have joint source energy mean at most \(e\), and
compare insertion of source waiting time \(t\) with no wait. Purifying
the entire prefix and removing the irrelevant ground-energy phase gives
\[
 \|(e^{-itK}\otimes I-I)\Psi\|^2\le2|t|e.
 \tag{8}
\]
The inequality follows pointwise from
\(|e^{-itx}-1|^2\le2|t|x\) for \(x\ge0\).
Half trace distance is at most \(\sqrt{2|t|e}\). Every subsequent
complete channel contracts that distance, including retained pointer
and idle-reference outputs. For finitely many waiting-time changes,
a hybrid comparison bounds the final half trace distance by
\[
 \min\left\{1,\sum_j\sqrt{2|\Delta t_j|}\,
       (\sqrt\epsilon+3a^2S/4)\right\}.
 \tag{9}
\]
Thus bounded fixed uncontrolled waiting-time comparisons become indistinguishable
under the same subcritical strength condition. This is a complete-output
statement, without stationarity of the joint finite purification. It
does not make all pointer records deterministic; zero-energy records
and degeneracy can survive.

Controlled source evolution remains an allowed energy-preserving operation
in (6), but its query-time changes are outside (8)--(9). On a ground input,
\(e^{-itH\otimes A}\) includes the observable apparatus phase
\(e^{-itE_0A}\). The unchanged source energy zero does not permit deleting
that phase. Even a one-dimensional source supplies this countercontrol.

Equations (6)--(7) also give a necessary observation strength. An original
probability at least \(p>0\) above \(\eta>0\) requires
\[
 S\ge{4\over3a^2}(\sqrt{p\eta}-\sqrt\epsilon)_+.
 \tag{10}
\]
If each normalized phase call has \(|r|\le1\), this lower-bounds its
number of calls. It is not a lower bound on arbitrary native hardware
gates. Receipts enclose the square roots outward by exact rationals.

## 4. Positive finite work does not remove the obstruction

For a uniform retained origin and balanced signs, let the phase strength
be \(r\). In the sine-chord domain \(|r|\le30R^4\), its original
signed source work is enclosed by
\[
 {2r^2D_R\over900R^6}\,{h_n\over a}
 \le W_n(r)\le{9\over16}r^2a^4,
 \qquad h_n=-a\langle H_{\rm hop}\rangle_{\rho_n}/V>0.
 \tag{11}
\]
The lower bound uses the actual selected hopping density, not an anchor
expectation, and follows from the original local-phase address average
and \(1-\cos(\pi x)\ge2x^2\) for \(|x|\le1\). The upper bound
uses \(1-\cos x\le x^2/2\) and (3). Work is positive for a nonzero
fixed strength at every finite member but tends to zero under bounded
strength. The complete-history obstruction therefore does not rely on
a vanishing finite matrix element.

The original local Weyl unitary has precisely
\[
 r_{\rm Weyl}=30/a^2,\qquad
 \pi r_{\rm Weyl}Q_X=\pi a\sum_y f_R(y-X)(N_y-30I).
 \tag{12}
\]
Its bound from (5) is \(45/2\), agreeing with the existing source form
lemma. Its positive finite work has a uniform lower bound. It is outside
the condition \(a^2S_n\to0\), so neither (7) nor (9) proves its
nonelastic survival or collapse. That apparatus needs an actual selected
energy-weighted tail or a fixed-time contrast. Passivity and the source
polynomial constrain spectral slopes and work, not those missing current
correlations. An operator-norm replacement is excluded by the separate
actual physical-anchor ultraviolet bound; that anchor is not the selected
state. No physical particle, rest mass or interacting continuum follows
from the apparatus restriction proved here.
