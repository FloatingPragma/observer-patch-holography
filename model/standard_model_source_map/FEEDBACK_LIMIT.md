# Correlated record feedback and its static action limit

The complete source record law admits a controlled weak-pulse limit with
the same four flavor coefficients as
[`native_flavor_orbit`](../native_flavor_orbit/PROOF.md). The work register,
record-conditioned pulse and coherent clock are specified controller
services. This result neither selects those services from passive repair
nor selects a physical vacuum or particle pole. The finite-step channel
contains noise; the theorem bounds its limit without dropping its records.

## Source law and stationary coefficient

Use the twelve-port memory, reset port2, geometric rank-one projectors
\(P_m\), and half-turns from
[`native_family_record_memory.py`](../common_source/native_family_record_memory.py).
The independent complete event flags have probabilities

\[
 w_{00}=(1-p)(1-q),\quad w_{01}=(1-p)q,\quad
 w_{10}=p(1-q),\quad w_{11}=pq,\qquad 0\le p<1,\quad0\le q\le1. \tag{1}
\]

For each old memory \(m\), the actual destination is \(F(m,C,E)\): a
non-C event resets to2; a C event applies its half-turn and then the E
half-turn if E occurs; a forbidden core destination resets to2. All twelve
inputs and all four events occur in the instrument, including zero weights.
For duty \(D\in\{1,C,E,CE\}\), put \(V_{m,r}=\chi_D(r)X\otimes P_{F(m,r)}\).
Here \(X\) is the Pauli matrix of the declared two-state logical register,
not a derived physical particle-hole identification.

With stationary law \(\pi\), the unconditional mean is

\[
 \bar V=X\otimes Z_D,\qquad
 Z_D=\sum_{m,r}\pi_m w_r\chi_D(r)P_{F(m,r)}.             \tag{2}
\]

The source code checks (2) against the independently retained
`native_flavor_orbit.algebra.coefficient`. Its traces are \(1,p,q,pq\);
no division by duty probability is made. In particular,
\(Z_C=Y-(1-p)P_2\), rather than \(pY\). Writing \(a=p(1-q),b=pq\),

\[
 Y=\frac{(1-a)P_2+a(1-a)P_3+abP_{10}}{1-a^2+ab},
 \qquad P_9=P_{10}.                                    \tag{3}
\]

The baseline unchanged C/E law has \(p=1/12,q=1/6\) and128 raw outcomes.
The configured current-gated source uses the actual thirteen-outcome
current instruments of
[`native_gauge_family_records.py`](../common_source/native_gauge_family_records.py):
\(p=\sin^2(\theta\sqrt{c_L})/12\),
\(q=\sin^2(\theta\sqrt{c_R})/6\), with the source Higgs-domain Casimirs.
Its21632 native/current fine labels per old memory reduce to512 labels
after grouping the current outcomes into two binary jumps. Scalar coarse
effects on those invariant input domains prove (1), even for entangled
auxiliary legs. The domain, leg preparation and supplied angle are part of
the source specification. The theorem does not apply these probabilities
to an arbitrary target-field state outside those domains.

## Reduced evolution, with a common Hamiltonian

The argument holds for any bounded Hermitian record Hamiltonians
\(\|V_{m,r}\|\le L\) and any common selfadjoint \(H_0\) for which their
unitary evolutions exist. In finite dimension this requires no extra domain
condition. Require that each reset-record Hamiltonian is independent of
the old memory. All four duties obey this condition, including E-duty on
non-C events. At step size \(h=t/N\), apply

\[
 K_{m,r}=\sqrt{w_r}\,|F(m,r)\rangle\langle m|
       \otimes e^{-ih(H_0+V_{m,r})}.                    \tag{4}
\]

Old memory is an explicit classical record. Equation (4) dephases a
coherent memory input; no unmeasured memory coherence is assumed. For each
old input, \(\sum_rK_{m,r}^*K_{m,r}=|m\rangle\langle m|\otimes I\).
Every history and every disabled pulse is retained.

Equip classical-quantum blocks with \(\|x\|=\sum_m\|x_m\|_1\), allowing
an arbitrary reference space. The zero-pulse transition has
\(T=(1-p)R+pS\), where R resets and S is stochastic. If \(\sum_m\eta_m=0\),
\(\|T\eta\|\le p\|\eta\|\). Set
\(\Pi x=\pi\otimes\sum_mx_m\) and \(\eta=(I-\Pi)x\).

Conjugate one step by the common \(e^{ihH_0}\). Duhamel gives
\(\|\operatorname{Ad}_{U_{m,r}}-\operatorname{Ad}_{U_0}\|_\diamond\le2hL\).
The reset terms cancel on zero-sum blocks, since their pulse is independent
of old memory. Therefore

\[
 \|\eta_{j+1}\|\le p\|\eta_j\|+4hL,\qquad
 \|\text{marginal coupling from }\eta_j\|_1
       \le2hpL\|\eta_j\|.                             \tag{5}
\]

In the interaction picture, each commutator generator has diamond norm
at most \(2L\). Its Dyson remainder after the first integral is bounded by
\(2h^2L^2\): the double integral has area \(h^2/2\), and the exact
intervening unitary channels have norm1. The stationary mean of the first
integral is exactly the first integral for \(\bar V\), even when
\([H_0,V_{m,r}]\ne0\). Comparing the two remainders costs at most
\(4h^2L^2\). No factor \(\|H_0\|\) enters this estimate.

Sum (5), use \(\|\eta_0\|\le2\), and telescope against the exact target
unitary \(\bar U(t)=e^{-it(H_0+\bar V)}\). The reduced channel error,
in **unhalved** trace distance uniformly over input and reference, is

\[
 \delta_N\le
 \frac{4pLt+4(1+p)L^2t^2}{N(1-p)}.                    \tag{6}
\]

A stationary **product** memory preparation removes the linear transient
term in this reduced bound. A stationary memory marginal correlated with
the work/reference does not. The proof uses contraction of completely
positive channels throughout, so the reduced estimate is reference-safe.

For the logical service L=1. Equation (6) extends the existing always-duty
bound to all four source duties, with the exact means (2). It supplies no
finite-N equality between a noisy memory channel and its mean unitary.

## The complete archive comparison

For fixed initial memory m, let \(r\) denote the entire event path,
\(p_r\) its original probability, and \(U_r\) its work unitary. Define the
comparison instrument to keep the **same path** and final memory, applying
\(\bar U\) on the work register on every path. This comparison archive is
not a unitary system and is not replaced by a stationary product archive.

For purified work/reference input \(\psi\), the reduced error (6) implies

\[
 \sum_rp_r\{1-|\langle\bar U\psi,U_r\psi\rangle|^2\}
 \le\delta_N/2.
\]

This follows by measuring the target pure state in the reduced output.
The unhalved distance of the complete flagged outputs is consequently

\[
 \sum_r2p_r\sqrt{1-|\langle\bar U\psi,U_r\psi\rangle|^2}
       \le\sqrt{2\delta_N}.                           \tag{7}
\]

Purifying each initial CQ block separately and taking convex combinations
proves the same bound for arbitrary initial memory/work/reference
correlations. Those correlations are retained in the ideal archive. The
archive estimate uses the uniform fixed-memory (6), even if the reduced
comparison started in a stationary product.

The source fine records can be restored without increasing (7). Conditional
on each coarse flag, its scalar effect makes the normalized auxiliary
instrument trace preserving on the admitted domains. Apply this same
conditional CPTP refinement to the auxiliary legs and the archived coarse
path in both comparisons, keeping every fine outcome. Work pulses and
auxiliary-leg operations commute. This also covers auxiliary/reference
correlations; it does not identify their state with a blank or stationary
register. Unspecified adaptive access to the archive changes the controller
and is outside (4). No path is conditioned on successful feedback.

## Lift to the configured odd-CAR action

The existing Higgs-line map
[`native_record_yukawa_field.assemble_kernel`](../common_source/native_record_yukawa_field.py)
is real-linear in its family matrices:

\[
 K_s(P)=D_s\otimes\overline P+(D_s\otimes\overline P)^T,
 \qquad \Delta_s(P)=\epsilon_{\rm spin}\otimes\overline{K_s(P)}.
                                                               \tag{8}
\]

Here the neutral unit-Higgs directed tensors \(D_s\) are partial
isometries with disjoint initial/final spaces, of ranks3,3,1 for the
up, down and charged-lepton action channels. This is the source tensor
normalization verified by `higgs_tensor_certificate`, not an empirical
sector normalization. On the ninety physical odd-CAR modes define the
bounded pairing operator

\[
 V_s(P)=\tfrac12\sum_{i,j}\Delta_s(P)_{ij}a_i^*a_j^*+\mathrm{h.c.}
                                                               \tag{9}
\]

An antisymmetric matrix has unitary skew blocks with singular values
duplicated in pairs. CAR rotation makes (9) a sum of commuting disjoint
two-mode pair terms. Each has eigenvalues \(\pm d,0,0\), proving
\(\|V_s(P)\|=\|\Delta_s(P)\|_*/2=\|K_s(P)\|_*\).
For rank-one P the two disjoint blocks of (8) have2 rank(D_s) unit
singular values. Thus the three sector bounds are6,6,2, and a simultaneous
controller has \(L\le14\) per cell. The mean of (9) is exactly
\(V_s(Z_D)\). Complex conjugation/transpose in (8) is essential to this
all-left convention; Nambu rows are not extra CAR flavors.

This lift explicitly **configures** record-conditioned physical CAR
Hamiltonians (9). The logical six-state pulse alone does not implement
them. The auxiliary current legs that generate the scalar probabilities
are separate from an arbitrary target field state. With independent sector
memories at V finite cells, a simultaneous all-reset event gives

\[
 p_{\rm joint}=1-\prod_{x,s}(1-p_{x,s}),\qquad
 L\le14V|v_H|.                                        \tag{10}
\]

Use these constants in (6)-(7), with the common kinetic Hamiltonian included
in \(H_0\). This is a quantitative finite-volume action limit. The global
reset estimate (10) is not uniform in source growth. A neutral fixed Higgs
background preserves electric charge; invariance under all twelve gauge
currents requires transforming the Higgs tensor as well. No dynamical gauge
vacuum is substituted for that specified background.

If each implemented step has complete-output diamond error at most
\(\varepsilon_j\), composition adds \(\sum_j\varepsilon_j\) to (7).
An input-state error adds its unhalved trace distance. The theorem therefore
requires a coherent pulse service with \(h=t/N\), sufficient attempt and
archive capacity, and vanishing total compilation error. It does not assert
execution of an arbitrarily long native tape, a clock selected by passive
memory, or a volume-uniform continuum construction.

Here t counts the specified target coherent-pulse time. The fixed-angle
auxiliary current gate and native repair attempts have their own execution
costs each round. Equation (6) does not equate t with their total elapsed
time or infer arbitrarily fast memory from those fixed gates. A common
implementation clock must account for gate overhead and target evolution
during it, or supply the corresponding parallel service explicitly.

## Certificate and physical boundary

`feedback_limit.py` checks the complete48-row coarse memory and all four
exact source means, computes rational outward error budgets, and supplies
finite path Kraus controls with noncommuting H0. Tests independently check
all-input completeness, reference-entangled complete histories, singular
event factors, the physical pairing norm on an explicit Fock space, and
the source Higgs tensor ranks. They do not replace the source state by a
trace or isolated ground.

The limit preserves the existing real matrices Z_D, their singular spectra
and their zero charged-current CP invariant. A common Higgs amplitude
multiplies their action singular values; it does not square signed
circulant roots. No Koide relation or phase is generated by averaging.
Physical Higgs selection, source clock calibration, interacting poles and
the continuum/mirror-removal requirements are distinct from this finite
complete-record action construction.
