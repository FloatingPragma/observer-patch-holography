# A finite autonomous controller for native contrast events

The construction is a time-independent Hamiltonian that executes a specified
finite program of the [native contrast law](../contrast_memory_spherical_v1/law.py), local
event counters and the [eleven-bit seam reader](../contrast_memory_spherical_v1/reader.py).
It produces simulation-generated observer records. The clock, its couplings,
the program, blank apparatus and the initial adjacent pair are supplied
architecture and preparation. The result concerns this finite controlled
process; it supplies no physical proper-time calibration, autonomous emitter,
continuum limit or absorbing halt.

The circuit-to-Hamiltonian construction follows the Feynman clock framework
described by Costales, Gunning and Dorlas,
[Efficiency of Feynman's quantum computer](https://arxiv.org/abs/2309.09331).
The engineered path uses the spin-chain transfer construction of Christandl,
Datta, Ekert and Landahl,
[Perfect state transfer in quantum spin networks](https://arxiv.org/abs/quant-ph/0309131).
The normalization, full data-space conjugacy, prefix instruments and finite
observation interval are derived below. Runtime estimates for a uniform clock
or an adiabatic construction do not apply to this engineered path.

## 1. Native source and finite program

For one seam, let \(J\) embed the nine load states into the blank sector of
the complete 99-state load and memory space. The exact identities in
`contrast_memory_spherical_v1/law.py` give

\[
 h_\kappa=U G_\kappa U^*,\qquad
 h_\kappa J=J a_\kappa,\qquad
 a_\kappa=\tfrac12(I-T)+\kappa(n_a-n_b)^2.
\]

Here \(U\) is the full Hermitian echo unitary, including the complementary
memory states. Functional calculus gives, for every real event action
\(\theta\),

\[
 e^{-i\theta h_\kappa}
   =Ue^{-i\theta G_\kappa}U^*,\qquad
 e^{-i\theta h_\kappa}J=Je^{-i\theta a_\kappa}.
\]

These are operator identities on all nine source inputs, including doubled
loads, and remain identities after tensoring with an arbitrary reference.
The controller uses this complete native event as a supplied gate in one
transition. It does not assign separate clock transitions or intermediate
observations to the three factors of the echo.

For each port of the cell, a reference load register and its memory supply
the contrast-only term. Its full generator is \(\kappa U G_1U^*\), with
\(G_1=0\oplus D^2\otimes I\). On blank memory and reference load zero its
exponential reduces to \(e^{-i\theta\kappa n_v^2}\), and keeps that reference
load zero. Including all twelve such terms gives the degree-six convention
of the [primary specification](../contrast_memory_spherical_v1/SPEC.md).

In particular, this program does not identify a native event with the
shifted gate in `contrast_memory_spherical_v1/sectors.py:SeamEvents`. The latter has
generator \(\tfrac12(I-T)-2\kappa n_an_b\); the difference is the local
operator \(\kappa(n_a^2+n_b^2)\). The reference phases and native seam
phases are retained explicitly.

Order the thirty cell seams as `CELL_SEAMS`, followed by the twelve
reference terms in port order. Write the resulting reduced generators as
\(A_1,\ldots,A_{42}\). One period is the chronological sequence

\[
 e^{-i\tau A_1/2},\ldots,e^{-i\tau A_{41}/2},
 e^{-i\tau A_{42}},
 e^{-i\tau A_{41}/2},\ldots,e^{-i\tau A_1/2},
 \qquad \tau=\tfrac18.
\]

Let \(S\) denote this product in its acting order. The declared program
applies eight periods, four stages of a fresh reader on seam zero, eight
periods, and four stages of a second fresh reader on seam one. There are
\(16\cdot83=1328\) source operations and eight reader stages, giving
\(c=1336\) nonpadding transitions. It appends \(1336\) identity
transitions, so the path has \(L=2672\) edges and \(L+1\) vertices.

The reduced source space contains every total-two configuration on twelve
load qutrits: 66 distinct pairs and twelve doubled loads, dimension 78.
This whole space is invariant under the source gates. The reported adjacent
basis preparation is supplied; it is not the finite-fuel emitter of the
primary model. The operator and instrument theorems also hold for any
source-reference density operator supported on this 78-dimensional space.

The finite product \(S\) is the specified evolution. Equality with
\(\exp(-i\tau\sum_jA_j)\) is not asserted. A continuous-source theorem
requires a separate approximation bound before transfer to this product.
Two readings in one run form a sequential instrument with the backaction
of the first reading; the primary model's fresh-run reading bounds do not
remove that backaction.

## 2. Counters, carrier and support

Each seam event increments the counters of its two port observers; each
reference event increments its one owning port's counter. Reader stages
and padding do not increment source-event counters. Each increment is a
cyclic shift on a declared finite alphabet, and is a unitary on every
counter input. Its inverse decrements the counter on a backward clock
transition.

Starting from zero, each of ports zero through ten receives 192 increments;
port eleven receives 176. The latter reference is the central, single
full-action gate of each period. Counter dimensions are consequently
193 for eleven ports and 177 for port eleven. No prefix of the declared
program wraps these counters. The value records the source events in the
current coherent prefix, rather than a monotonically growing archive of
every traversal of the rail.

The full data space includes twelve source qutrits, thirty seam memories,
twelve reference qutrits with twelve reference memories, the twelve
counters and two eleven-bit readers. Its dimension is

\[
 d=3^{24}\,11^{42}\,193^{11}\,177\,2^{22}.
\]

The mathematical clock is a single particle with \(L+1\) possible sites,
so the declared clock and data carrier has dimension \((L+1)d\). A unary
implementation with one bit per clock site allocates \(L+1\) bits, and
the path theorem applies to its invariant one-excitation sector. It is
not a spectrum theorem for the entire \(2^{L+1}\)-dimensional unary space
or for states with several clock excitations.

A seam gate touches its two load qutrits, private seam memory and the two
associated counters. A reference gate touches the owning port, its
reference apparatus and its counter. Reader gates act on the chosen seam
predicate and the allocated reader block; the predicate itself depends
on both endpoint loads. Each clock transition couples its two adjacent
clock sites to this fixed gate support. Full gates on that support are
declared interaction primitives. An elementary two-register synthesis and
a routing of the clock path exclusively through the original cell seams
are additional implementation tasks, not consequences of support size.
The rail connectivity is supplied apparatus connectivity.

## 3. Complete clock conjugacy and spectrum

Let \(V_1,\ldots,V_L\) be the full data unitaries just specified, including
counter shifts, the full reader permutations and the padding identities.
Define

\[
 W_0=I,\qquad W_k=V_k\cdots V_1,\qquad
 D=\sum_{k=0}^{L}|k\rangle\langle k|\otimes W_k,
\]

and

\[
 C_L=\sum_{k=0}^{L-1}J_k
   (|k+1\rangle\langle k|+|k\rangle\langle k+1|),\qquad
 J_k=\frac{\sqrt{(k+1)(L-k)}}{L}.
\]

The autonomous Hamiltonian is

\[
 H=\sum_{k=0}^{L-1}J_k
  (|k+1\rangle\langle k|\otimes V_{k+1}
  +|k\rangle\langle k+1|\otimes V_{k+1}^*)
  =D(C_L\otimes I_d)D^*.
\]

The equality follows block by block from
\(W_{k+1}W_k^*=V_{k+1}\). Since each \(W_k\) is unitary on the full data
space, this is a conjugacy on that entire space, including dirty private
memories, nonzero counters and all reader words. A reducing source space
is needed for the source interpretation, not for this conjugacy.

Identify the clock basis with the spin-\(L/2\) representation. Then
\(C_L=2J_x/L\), so

\[
 \operatorname{spec}(H)=\{-1+2r/L:0\le r\le L\},\qquad
 \operatorname{mult}(-1+2r/L)=d,\qquad \|H\|=1.
\]

The largest individual coupling is
\(\sqrt{\lfloor(L+1)^2/4\rfloor}/L\), slightly larger than one half.
The exact norm bound comes from the spin spectrum, rather than bounding
every edge strength by one half. No independent always-on native source
Hamiltonian is added to \(H\); adding one generally changes the conjugacy.

## 4. All prefixes and arbitrary-reference instruments

For any data vector, including an entangled untouched reference, the
propagator on a clock initially at zero is

\[
 e^{-itH}(|0\rangle\otimes|\psi\rangle)
  =\sum_{k=0}^{L}\alpha_k(t)|k\rangle\otimes W_k|\psi\rangle,
\quad
 \alpha_k(t)=(-i)^k\sqrt{\binom Lk}
       \sin(t/L)^k\cos(t/L)^{L-k}.
\]

The real powers in this expression retain their signs at all times.
The prefix distribution is the complete normalized binomial law

\[
 p_k(t)=\binom Lk p(t)^k(1-p(t))^{L-k},\qquad
 p(t)=\sin^2(t/L).
\]

Tracing out the clock yields the channel
\(\mathcal E_t(\rho)=\sum_k p_k(t)W_k\rho W_k^*\). This follows from
orthogonality of the clock states, without a measurement during execution.
All prefixes, including unfinished readers and padding, contribute.

Let \(J_{\rm prep}\) append blank memories, zero references and counters,
and fresh reader words to an arbitrary permitted source input. For public
reader label \(z\), let \(Q_z\) be its actual computational projector on
the apparatus. The unnormalized public instrument, retaining the source
and an untouched reference, is

\[
 \Phi_{t,z}(\rho)=\sum_{k=0}^{L}p_k(t)\,
 \operatorname{Tr}_{\rm apparatus}
 [Q_z W_kJ_{\rm prep}\rho J_{\rm prep}^*W_k^*Q_z].
\]

Equivalently, choose a complete private apparatus basis \(\eta\) and
use Kraus operators
\(\sqrt{p_k(t)}\langle z,\eta|W_kJ_{\rm prep}\).
The sum runs over every private word; it cannot replace a retained failure
label by an identity channel on the source. Mismatch or unavailable public
labels can merge private states that distinguish the source predicate.
Tracing those states retains their dephasing effect. All formulas tensor
with the identity on an arbitrary reference and preserve total trace when
summed over the complete public outcome set.

For nominal completed readers let \(P_b^{(e)}\) denote the projector onto
predicate value \(b\in\{0,1\}\) on seam \(e\). Their public records are
decision one and public triplet \(7b\). The completed sequential source
instrument has Kraus operators

\[
 K_{b_1,b_2}=P_{b_2}^{(e_1)}S^8P_{b_1}^{(e_0)}S^8,
 \qquad \sum_{b_1,b_2}K_{b_1,b_2}^*K_{b_1,b_2}=I.
\]

This proves the native source-to-archive mapping for arbitrary references,
including coherences within either predicate sector. No-arrival,
unfinished, mismatch, unavailable and agreed-false records remain distinct
whenever their actual decision and validity bits distinguish them. A zero
public triplet alone is insufficient to identify readiness.

## 5. Completion throughout a finite observation interval

For \(k\ge c\), identity padding gives \(W_k=W_c\). Put
\(q(t)=\sum_{k<c}p_k(t)\). After discarding the clock, and with any fixed
terminal instrument, the channel is

\[
 \Phi_t=(1-q(t))\Phi_{\rm complete}+q(t)\Phi_{\rm unfinished,t}.
\]

Both terms are trace-preserving channels when their complete public
outcome labels and residual quantum outputs are retained. Thus

\[
 \tfrac12\|\Phi_t-\Phi_{\rm complete}\|_\diamond\le q(t).
\]

In particular, every terminal public probability law differs in total
variation by at most \(q(t)\), uniformly for all permitted source states
and untouched references. This is an unconditional bound; it does not
discard unfinished mass. The comparison discards the clock on both
sides and does not identify their full clock distributions.

On fresh reader words, the second reader's decision bits stay zero until
its final comparison stage. That stage writes one of decisions one, two
or three for the nominal, mismatch and unavailable reader variants.
Consequently the observable predicate “second decision is nonzero” equals
\(k\ge c\) for these variants. It marks completion, not successful
agreement. An intervention that omits the final comparison does not
inherit this observable readiness statement.

For every

\[
 t\in[\pi L/3,\;2\pi L/3],
\]

\(p(t)\ge3/4\). Binomial stochastic monotonicity therefore gives the
exact finite bound

\[
 q(t)\le\frac{1}{4^L}\sum_{k=0}^{c-1}\binom Lk3^k,
 \qquad c=L/2.
\]

For completeness, if \(X\) is binomial with success probability at least
\(3/4\), applying the exponential Markov bound with parameter
\(\log3\) yields

\[
 \Pr(X<L/2)\le
 3^{L/2}(1-p+p/3)^L
 \le(\sqrt3/2)^L\le e^{-L/8}.
\]

The last inequality follows from
\(\log(4/3)=\log(1+1/3)\ge(1/3)/(1+1/3)=1/4\).
With \(L=2672\), the looser displayed ceiling is \(e^{-334}\).
This proves a whole interval, without inference from a sampled time grid.
Binary floating execution does not certify this tiny error numerically;
it is an exact theorem for the declared Hamiltonian. Numerical residuals
and comparisons retain their separate finite precision tolerances.

There is also a conditional implementation bound. Suppose every realized
source gate is unitary and differs from its prescribed full-domain gate
by at most \(\delta\) in operator norm, while the reader gates and clock
couplings are exact. Telescoping each prefix product gives
\(\|\widetilde W_k-W_k\|\le1328\delta\). Orthogonality of the clock
states gives the same bound on the corresponding prepared-clock
isometries. Contractivity under the terminal instrument therefore bounds
the half-diamond error by \(1328\delta\), uniformly in the clock time and
the untouched reference. The complete-channel comparison is bounded by
\(q(t)+1328\delta\), capped at one. For \(\delta=10^{-6}\), the added
ceiling is \(0.001328\). This assumes a bound on the actual unitary gates;
a floating-point residual on selected columns does not establish that
physical implementation hypothesis. Errors in the clock couplings need
their own bound.

## 6. Transfer, return, energy and speed

At \(t_*=\pi L/2\), the only nonzero clock amplitude is
\(\alpha_L=(-i)^L\). The completed program is implemented with unit
probability. At \(t_{\rm return}=\pi L\),

\[
 e^{-it_{\rm return}H}=(-1)^L I.
\]

Here \(L\) is even, so the return is exactly the identity. In uninterrupted
closed evolution the data, counters and both internal reader archives
return to their prepared state. There is no permanent completion record
or absorbing halt. A terminal measurement changes the state; the exact
unmeasured return formula does not describe that altered trajectory or
an external archive storing the measurement result. The observation
interval concerns any one terminal snapshot, not repeated monitoring.

For the prepared clock at zero, regardless of the data state,

\[
 \langle H\rangle=0,\qquad
 (\Delta H)^2=J_0^2=1/L,\qquad
 \langle H\rangle-E_{\min}=1.
\]

Energy relative to the controller ground is positive. A complete
computational measurement of the clock removes all neighboring-clock
coherences, so its postmeasurement mean controller energy is also zero.
No equal-energy conclusion for a pointer-only measurement follows from
this argument. The spectrum and conserved energy are those of \(H\),
not of the primary always-on native Hamiltonian. Source event actions
\(\theta\) are parameters in the gate matrices; they are not elapsed
times between clock transitions. Gate synthesis, state preparation and
external terminal readout costs are not a thermodynamic work ledger.

Let \(X=\sum_k k|k\rangle\langle k|\otimes I\). Since \(D\) commutes
with \(X\), the spin representation gives
\(\|i[H,X]\|=1\). Therefore
\(|d\langle X\rangle/dt|\le1\) for every state. For the prepared clock,
\(\langle X\rangle=L\sin^2(t/L)\), in agreement with this bound.
This is a speed bound in the declared rail coordinate and energy units,
not a limiting speed in effective spacetime. Scaling the entire
Hamiltonian by a positive factor scales energies and inverse times
together.

## 7. Executable verification domain

The native reduction identities and full reader permutations establish
the full operator domain. Numerical execution may use the invariant
78-state source, deterministic prefix counter labels and the reachable
reader branches, provided its representation preserves their coherent
joint isometry. Tracing or measuring those branches before later source
operations must reproduce the specified instrument, not suppress its
backaction. Doubled loads, dirty apparatus checks and every failed control
remain explicit domains or outcomes; they are not certified by testing
only the prepared adjacent basis state.

A complete receipt can store the full clock law and the prefix-dependent
source and reader laws in factored form. Such storage must reconstruct
all declared outcomes without selecting completed prefixes or dropping
zero and failure labels. Agreement and the contrast control are separate
checks: the reader with a common flipped source predicate can agree on
an incorrect value. The finite clock theorem supplies execution and an
archive observation interval; empirical particle or field conclusions
require their own measured records and source-domain arguments.
