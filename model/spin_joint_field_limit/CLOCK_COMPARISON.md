# Selection rotation and the clock of the joint family

Use the physical spin frame of the actual finite bipartite source,

\[
 H_0=J\sum_{xy}(S_x\cdot S_y-1/4),\quad
 Z=\sum_x\epsilon_xS_x^1,\quad H_B=H_0-BZ,
\]

with J,B,p>0. The time variable T denotes source time T/(Jp). A common
scalar reference can be subtracted from both Hamiltonians. Scalars cancel
in the Heisenberg relations below; a coherently controlled clock retains
its relative scalar phase. The direct joint family uses the prepared
H_B ground and this same H_B clock.

## Exact selection rotation

For real finite tests define

\[
 Q_p(F)=p^2\sum_x\epsilon_xF(px)S_x^3,\qquad
 P_p(G)=-p\sum_xG(px)S_x^2,\qquad D=-{B\over Jp}Z.
\]

The exact spin commutators give

\[
 i[D,Q_p(F)]={B\over J}P_p(F),\qquad
 i[D,P_p(G)]=-{B\over Jp^2}Q_p(G).                    \tag{1}
\]

These identities hold on the full finite source space. Let
\(\theta=BT/(Jp)\). The selection-only automorphism
\(\beta_T(A)=e^{iTD}Ae^{-iTD}\) therefore satisfies

\[
 \beta_T(Q_p(F))=\cos\theta\,Q_p(F)+p\sin\theta\,P_p(F),
\quad
 \beta_T(P_p(G))=\cos\theta\,P_p(G)-{\sin\theta\over p}Q_p(G).
                                                               \tag{2}
\]

Equation (2) is a finite source rotation, with no continuum equation of
motion assumed. In the family p=n^-2, B=Jn^-5, its angle is T/n^3.
The Q-to-P coefficient is O(n^-5), the P-to-Q coefficient is O(n^-1),
and 1-cos(theta) is O(n^-6) at fixed T.

For any normalized finite source vector psi put

\[
 q_F=\|Q_p(F)\psi\|,\quad r_F=\|P_p(F)\psi\|,\quad
 c_p(F,G)={p^3\over2}\sum_x|F(px)G(px)|.
\]

The operator commutator norm is bounded by c_p(F,G). Duhamel's formula
for the two exponentials in (2), followed by integration of the
equal-time commutator bound, gives

\[
\begin{split}
 \|(\beta_T(e^{iQ_p(F)})-e^{iQ_p(F)})\psi\|
 &\le(1-\cos\theta)q_F
       +p|\sin\theta|\bigl(r_F+c_p(F,F)/2\bigr),\\
 \|(\beta_T(e^{iP_p(G)})-e^{iP_p(G)})\psi\|
 &\le(1-\cos\theta)r_G
       +{|\sin\theta|\over p}\bigl(q_G+c_p(G,G)/2\bigr).
                                                               \tag{3}
\end{split}
\]

For example, in the first line the integrand contains P_p(F)e^{iuQ_p(F)}psi.
Its norm is at most r_F+u c_p(F,F); integration over 0<=u<=1 explains
the factor 1/2. This uses finite matrices and makes no assumption about
a common domain of continuum field products. The same argument passes
to finite absolutely summable tests by norm approximation.

For an equal-time word applied to the source ground, a Q_p(F) norm is
at most its vacuum norm plus the sum of c_p(F,G_j) over preceding P
gates. A P_p(G) norm has the analogous sum over Q gates. This follows
by commuting the field through each bounded unitary and using its
operator commutator bound. The joint family's uniform vacuum moments
therefore make (3) tend to zero on every fixed equal-time word.

This statement also controls complete equal-time instruments assembled
from these gates, finite coherent settings and ancillary operations.
On each coherent setting use the largest applicable bound in (3).
The field-norm increment at a preceding controlled gate is bounded by
the largest c_p over its settings. Ancillary operations commute with
the source fields. A telescoping sum bounds the difference of the full
Stinespring output vectors, including an untouched initial reference.
Half trace distance and the total variation of all records are at most
that sum. Failure and complementary records are included; no successful
branch is renormalized. This result concerns the specified ground input
and its equal-time preparation circuit, rather than an all-input channel
norm bound for unbounded source fields.

## The additional estimate required to remove the field from the clock

Let alpha_B and alpha_0 be the two scaled Heisenberg evolutions on the
same finite source and the same initial state. For t>=0,

\[
 \alpha_B^t(A)-\alpha_0^t(A)
 =\int_0^t\alpha_0^{t-s}
       \bigl(i[D,\alpha_B^s(A)]\bigr)\,ds.             \tag{4}
\]

Consequently, on any vector psi,

\[
 \|(\alpha_B^t(A)-\alpha_0^t(A))\psi\|
 \le\int_0^t
 \|[D,\alpha_B^s(A)]e^{-i(t-s)L_0}\psi\|\,ds,         \tag{5}
\]

where L_0=H_0/(Jp) with any common scalar reference. Negative times
obey the corresponding oriented-integral formula. Applying (5) to
each factor and its actual preceding mixed-history vector gives a
sufficient criterion for equality of all fixed finite bounded-history
correlations: the finitely many integrals must vanish, including their
adjoint and complete-instrument factors. This is a condition on the
actual source commutators at nonzero time. Equation (3) proves only the
equal-time case and supplies no such uniform mixed-history estimate.
The Heisenberg criterion does not compare two differently recentered
coherently controlled clock unitaries.

The inherited all-input energy transport bound does not fill this gap.
For a smooth staggered rotation supported in a scaled region of fixed
size, the selection angle is O(B/(Jp)) and there are O(p^-3) edges.
Its single-axis energy kick is consequently bounded by

\[
 \delta^2=O\bigl(B^2/(Jp^5)\bigr),\qquad
 {\delta\over\sqrt{Jp}}=O\bigl(B/(Jp^3)\bigr).        \tag{6}
\]

The existing overlap estimate for histories of energy O(Jp) and time
O(1/(Jp)) therefore vanishes under B=o(Jp^3). The concentration
construction uses Jp^3/B->0. In the joint family B=Jp^(5/2), the
bound (6) fails to tend to zero. This failure of a sufficient estimate
does not establish that the two limits differ. It excludes using that
estimate to assert their equality. The direct H_B construction is
independent of an identification with the ordered zero-field prescription.

## Complete native histories with the retained field

A finite H_B clock needs no field removal. Let a finite box have N
sites and e bonds, set A=Je+BN/2, and let the absolute physical durations
of k evolution slots be t_j=|T_j|/(Jp). These are upper bounds on every
classically adaptive branch. With M product slices in each slot, apply all actual
exchange gates and then all selection rotations. The product-formula
operator error summed over the slots is at most

\[
       \eta_{\rm split}={A^2\over2M}\sum_j t_j^2.       \tag{7}
\]

Each slice costs 62e+11N native pulses, including the exchange scalar
phase. Thus P=kM(62e+11N) is a pulse upper bound. The existing dyadic
compiler supplies complete isometry error at most 29P2^-b. Setting
M>=A^2 sum_j t_j^2/eta and choosing b with 29P2^-b<=eta/2 makes their
sum at most eta. Zero-duration slots require no gates. No source
eigenvector or diagonalization is used by this exact rational planner.
Coherently controlled durations require their own controlled-gate pulse
count; the count above does not include that compilation.

Interleaved complete instruments preserve the accumulated half trace
error. Their own compilation errors are charged separately. If the
prepared input has error delta, the complete output, including all
failures and initial independent reference systems, has error at most
delta+eta+eta_instruments relative to the same finite H_B history.
Positive-field marginal and source locality errors can be added using
their existing complete local-history bounds. This circuit keeps the
declared field in both the prepared state and every clock slot.

`clock_word` emits the literal source gates lazily on an encoded spin
register. `clock_plan` and `cube_clock_plan` count them without allocating
large matrices or tapes. The short two-spin execution audits the complete
native implementation with a maximally entangled reference. It is not an
execution of a growing-volume continuum schedule. These clock statements
provide no sharp particle dispersion, Lorentz covariance or scattering
identification.
