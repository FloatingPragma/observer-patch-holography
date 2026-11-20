# The binding theorem of the A-contrast-memory-spherical-v1 model

This document proves that a pair of loads on a seam stays a pair under the contrast-memory law
(`A-contrast-memory-spherical-v1`, a candidate native repair law) on every federation of complete cells, with a
bound independent of the number of cells; gives the exact constants of one cell; and bounds the emitter and the
reader. The law, carrier and notation are those of [SPEC.md](SPEC.md). The constants are measured by
[theorem.py](theorem.py) and stored in [RESULTS.json](RESULTS.json).

## 1. Setting

$G = (V, E)$ is a federation: complete twelve-port cells, every internal seam kept, external seams forming a
partial matching on ports, so every port has degree $d_v \in \{5, 6\}$. On loads $n_v \in \{0, 1, 2\}$,

$$H = K + \kappa C,\qquad K = \tfrac12 \sum_{e} (I - T_e),\qquad C = 6 \sum_v n_v^2 - 2 \sum_{\{a,b\} \in E} n_a n_b .$$

Every summand of $K$ is a projector, so $K \ge 0$. $H$ conserves the multiset of loads. The sector of two loads
contains every pair of distinct ports and every doubled load; its dimension is $V(V+1)/2$. Its values of $C$ are

| Configuration | $C$ |
| --- | ---: |
| two loads on the two ends of a seam (a pair) | 10 |
| two loads on ports that share no seam | 12 |
| one doubled load | 24 |

since a binary configuration has $C = 6 \cdot 2 - 2\,[\text{adjacent}]$ and a doubled load has $C = 6 \cdot 4$. The
swaps conserve the multiset of loads, so the doubled loads form a reducing subspace; they stay part of the
complement below and of every outcome.

Let $P$ project onto the span of the pairs on seams and $Q = I - P$ in the sector of two loads. Write

$$H = \begin{pmatrix} A & B^* \\ B & D \end{pmatrix},\qquad A = PHP,\ B = QHP,\ D = QHQ .$$

## 2. The three inequalities

**Theorem 1.** On every federation of complete cells with a matching of external seams, and every $\kappa \ge 0$,

$$PHP \le (10\kappa + 6)\,P,\qquad QHQ \ge 12\kappa\,Q,\qquad \|QHP\| \le \sqrt{10}.$$

*Proof.* $C$ is diagonal, so $A = 10\kappa P + PKP$, $D = QKQ + \kappa\,QCQ$ and $B = QKP$.

*First inequality.* The diagonal of $PKP$ at the pair on $(a, b)$ is $(d_a + d_b - 2)/2 \le 5$. An off-diagonal
element of $PKP$ joins the pairs on $(a, b)$ and $(a, t)$ when one swap moves the load from $b$ to $t$; then $t$ is
adjacent to $a$ and to $b$, so $a, b, t$ span a face, and the element is $-1/2$. Two seams are so joined exactly when
they are two sides of one face. All faces lie inside cells, so an external seam has no such neighbour, and inside a
cell the seams with this relation form the line graph of the dodecahedron (faces of the cell are vertices of the
dual, seams are its edges). The adjacency matrix of a line graph is $M^{\top} M - 2I \ge -2I$ with $M$ the incidence
matrix, so the off-diagonal part of $PKP$, which is $-\tfrac12$ times that adjacency, is at most $1$. Hence
$PKP \le 6P$ and $A \le (10\kappa + 6)P$.

*Second inequality.* $QKQ \ge 0$ as a compression of $K \ge 0$, and $QCQ \ge 12\,Q$ by the table. Hence $D \ge 12\kappa Q$.
The condition $\kappa \ge 0$ is needed here: the uniform doubled-load state lies in $Q$ and has energy
$24\kappa$, which is below $12\kappa$ for negative $\kappa$.

*Third inequality.* $B$ has entries $-1/2$. Its column at the pair on $(a, b)$ has one entry for every swap that moves
a load of the seam to a neighbour that leaves the two loads apart: at most $(d_a - 1 - m) + (d_b - 1 - m)$ entries
with $m$ the number of common neighbours of $a$ and $b$. For an internal seam $m = 2$ and the count is at most 6; for an
external seam $m = 0$ and the count is at most 10. The largest absolute column sum is therefore at most 5. The line
of $B$ at two separate loads on $x, y$ has one entry for every common neighbour $t$ of $x$ and $y$ and every choice of
which load came from $t$: at most $2 \cdot 2 = 4$ entries, because a nonadjacent pair of ports has at most two common
neighbours (inside a cell by the geometry of the icosahedron; across cells because the external seams form a
matching, so at most one common neighbour lies in either cell and none in a third). Doubled loads have no entry.
The largest absolute line sum is at most 2. By the Schur test $\|B\| \le \sqrt{5 \cdot 2} = \sqrt{10}$. $\square$

## 3. Persistence at every value of the evolution parameter

**Theorem 2.** Let $A \le a P$, $D \ge c\,Q$, $\|B\| \le \beta$ with $\gamma = c - a > 0$, $m = \operatorname{rank} P$, and let $\Pi$
be the spectral projector of $H$ onto its lowest $m$ levels. Then

$$\|Q\,\Pi\| \le \beta/\gamma,\qquad \|(I - \Pi)\,P\| \le \beta/\gamma,$$

and for every normalized state $\psi$ with $P\psi = \psi$ and every real $t$,

$$\|Q\, e^{-itH} \psi\| \le \frac{2\beta}{\gamma},\qquad
\Pr(\text{no pair on a seam};\, t) \le \min\Big(1, \frac{4\beta^2}{\gamma^2}\Big).$$

*Proof.* By min-max with the trial spaces $\operatorname{ran} P$ and $\operatorname{ran} Q$, $\lambda_m(H) \le a < c \le \lambda_{m+1}(H)$, so the
gap of the ordered blocks is not closed by $B$. Put $H_{\rm low} = \Pi H \Pi$ (levels at most $a$) and $X = Q\Pi$. From
$QH\Pi = X H_{\rm low}$ and $QH = QHP + QHQ$,

$$D X - X H_{\rm low} = -B\,P\,\Pi .$$

The spectra of $D$ on $\operatorname{ran} Q$ and of $H_{\rm low}$ on $\operatorname{ran}\Pi$ are separated by $\gamma$, with $D$ above. The
solution of this Sylvester equation is

$$X = -\int_0^\infty e^{-s(D - a)}\, B P \Pi\, e^{s(H_{\rm low} - a)}\, ds ,$$

and $\|e^{-s(D - a)}\| \le e^{-s(c - a)}$, $\|e^{s(H_{\rm low} - a)}\| \le 1$, so $\|X\| \le \beta/\gamma$. In the same way
$Y = (I - \Pi)P$ solves $H_{\rm high} Y - Y A = (I - \Pi)\,Q\,B$ with $H_{\rm high} \ge c$ and $A \le a$, so $\|Y\| \le \beta/\gamma$.

For $\psi = P\psi$ write $\psi = \Pi\psi + (I - \Pi)P\psi$. The first part stays in $\operatorname{ran}\Pi$ under $e^{-itH}$, so
$\|Q e^{-itH}\Pi\psi\| \le \|Q\Pi\| \le \beta/\gamma$; the second part has norm at most $\beta/\gamma$. The sum is at most
$2\beta/\gamma$, and the probability of the outcome "no pair", which is $\|Q e^{-itH}\psi\|^2$, is at most $4\beta^2/\gamma^2$.
Mixed states supported in $\operatorname{ran} P$ follow by convexity. $\square$

**Corollary.** On every federation of complete cells with a matching of external seams, Theorems 1 and 2 give
$a = 10\kappa + 6$, $c = 12\kappa$, $\beta = \sqrt{10}$, $\gamma = 2\kappa - 6$, and for $\kappa > 3$

$$\Pr(\text{no pair on a seam};\, t) \le \frac{40}{(2\kappa - 6)^2}$$

at every value of the evolution parameter, hence at every reading of every observer. The bound does not depend on
the number of cells or on the routing of the external seams. At $\kappa = 128$ it is $40/250^2 = 0.00064$.

The dressed band $\operatorname{ran}\Pi$ is invariant under the law. The bare outcome "no pair" also contains its virtual
dressing, so it need not vanish or decrease along the evolution: the pair may separate briefly and is read on
other seams.

**Measured constants.** The proof uses the bounds $a - 10\kappa \le 6$ and $\beta \le \sqrt{10}$. On the federations of the
model the largest level of the adjacent block of the exchange, $\max \operatorname{spec}(PKP)$, and the largest singular value of
the coupling, $\|QKP\|$, are

| Federation | $\max \operatorname{spec}(PKP)$ (bound 6) | $\|QKP\|$ (bound 3.1623) |
| --- | ---: | ---: |
| 20 cells | 5.7008 | 2.0821 |
| 80 cells | 5.7000 | 2.0811 |
| 320 cells | 5.7000 | 2.0822 |

With these constants the bound is $4\beta^2/(2\kappa - 5.70)^2 = 17.3/(2\kappa - 5.70)^2$: 0.163 at $\kappa = 8$ and 0.00028 at
$\kappa = 128$ on 20 cells.

## 4. One cell

On one cell without degree term every port has five seams, the pair costs 8, separate loads 10 and a doubled load 20
($C = \sum_e (n_a - n_b)^2$). There are 30 pairs on seams, 30 pairs at distance two, 6 antipodal pairs and 12 doubled
loads. The exact constants are

$$a = 8\kappa + 5,\qquad c = \min\big(10\kappa + 4 - \sqrt6,\ 20\kappa\big),\qquad \beta = 2.$$

- $\max\operatorname{spec}(PKP) = 5$: the diagonal is 4 and the line graph of the dodecahedron has smallest adjacency level $-2$.
- $\|QKP\| = 2$: the largest absolute line and column sums of $QKP$ are both 2, and the uniform vector on the 30 pairs
  attains the value.
- The bottom of the block of separate loads is $4 - \sqrt6$. The two orbits of nonadjacent pairs (distance two and
  antipodal) give the quotient $\begin{pmatrix} 3 & -\sqrt5 \\ -\sqrt5 & 5 \end{pmatrix}$, whose lower level $4 - \sqrt6$ has
  a positive lift; since the block has nonpositive off-diagonal elements, the Perron argument makes it the bottom.
- Doubled loads cost $20\kappa$ and the exchange on them is nonnegative.

The degree term of the law adds $\kappa \sum_v n_v^2$ on one cell: $2\kappa$ on every binary pair and $4\kappa$ on a doubled
load. Then $a = 10\kappa + 5$, $c = \min(12\kappa + 4 - \sqrt6,\ 24\kappa)$, and the gap $\gamma = 2\kappa - 1 - \sqrt6$ (for
$12\kappa + 4 - \sqrt6 \le 24\kappa$) is the same. The bound is $16/\gamma^2$ capped at
one: 0.1016 at $\kappa = 8$ and 0.00025 at $\kappa = 128$.

## 5. The emitter

Add one fuel state $|F_e\rangle$ at the seam $e = (a, b)$ with energy $E_f = \langle e|H|e\rangle = 10\kappa + (d_a + d_b - 2)/2$, and
during the pulse the term $g(|F_e\rangle\langle e| + |e\rangle\langle F_e|)$ for the duration $\tau = \pi/(2g)$, with $H$ acting throughout.

**Theorem 3.** The emitted state $\psi(\tau)$ satisfies

$$\big\|\psi(\tau) + i\, e^{-iE_f\tau}\, |e\rangle\big\| \le \frac{\sqrt{d_a + d_b - 2}}{2g}.$$

*Proof.* On $\operatorname{span}\{|F_e\rangle, |e\rangle\}$ the comparison dynamics with diagonal $E_f$ and coupling $g$ is an exact
resonant rotation, $\phi(s) = e^{-iE_f s}(\cos(gs)\,|F_e\rangle - i \sin(gs)\,|e\rangle)$, which reaches $-i e^{-iE_f\tau}|e\rangle$ at
$s = \tau$. The full generator differs from the comparison generator on $\phi(s)$ only through
$r = (H - E_f)\,|e\rangle$, the off-diagonal column of $H$ at $e$: $d_a + d_b - 2$ allowed swaps with amplitude $-1/2$ each, so
$\|r\|^2 = (d_a + d_b - 2)/4$. Duhamel's formula gives
$\|\psi(\tau) - \phi(\tau)\| \le \|r\| \int_0^\tau |\sin(gs)|\,ds = \|r\|/g$. $\square$

After the pulse the converter is off and the law preserves the error in norm. Residual fuel reads "no pair". With
Theorem 2,

$$\Pr(\text{no pair, doubled load or fuel};\, t) \le \min\Big(1, \Big[\frac{2\beta}{\gamma} + \frac{\sqrt{d_a + d_b - 2}}{2g}\Big]^2\Big).$$

On the federations $d_a + d_b - 2 \le 10$. At $\kappa = 128$, $g = 64$, $\beta = \sqrt{10}$, $\gamma = 250$ the bound is the exact
rational

$$10\Big(\frac{2}{250} + \frac{1}{128}\Big)^2 = \frac{64009}{25600000} = 0.0025003515625$$

on every seam of every federation, and the weight outside the dressed band, $(\beta/\gamma + \epsilon)^2$, is at most
$35721/25600000$.

## 6. The reader with the law acting

The reader of a seam (SPEC.md, section 7) runs as four stages $G_j$, each a reversible involution on the eleven
bits, with Hamiltonian $L_j = \frac{\pi}{2\tau_j}(I - G_j)$, so $e^{-i\tau_j L_j} = G_j$ on every word. Its only control
on the loads is the diagonal projector $\mathbf 1\{n_a = n_b = 1\}$, so every stage commutes with the diagonal part
$V = \kappa C$ (and the fuel energy) of the law, and stages on different seams commute. Write $H = V + K$. The comparison
evolution omits $K$ and is exactly $e^{-iV\tau_{\rm read}} \prod_j G_j$, the ideal reader at the start of the window.
Duhamel's formula gives $\|U_{\rm actual} - U_{\rm ref}\| \le \tau_{\rm read}\|K\|$. In the sector of two loads every line of $K$ has
nonnegative diagonal at most $\Delta = 6$ and off-diagonal absolute sum equal to it, so $\|K\| \le 12$, and the total
variation between the public laws is at most $12\,\tau_{\rm read}$, for any number of seams read in parallel. At
$\tau_{\rm read} = 1/16384$ this is $3/4096$, and the probability that the public records show no agreed pair after the
emitter at $\kappa = 128$ is at most $64009/25600000 + 3/4096 = 0.0032327734375$ at every reading. After one complete
read round the ideal reader acts as $\rho \mapsto \sum_e P_e \rho P_e + Q\rho Q$, each branch on a seam is a pair covered by
Theorem 2, and later "no pair" stays below $0.0025003515625 + 0.00064 + 0.000732421875 = 0.0038727734375$.

The finite readout of one cell propagates all 161,792 states (78 two-load states, the fuel and 2048 reader words)
with the law acting in every stage: the public total-variation error is $2.09 \times 10^{-7}$ against the bound
$\tau_{\rm read}\|K\| = 4.2 \times 10^{-4}$ on the cell.

## 7. Scope

The theorem concerns the sum of the seam terms with a continuous evolution parameter. As local events the same
binding holds when each port observer serves each of its seams once per round, and it needs small event angles
when seams act independently at random ([events.py](events.py)). The bound controls which seams the two loads are
read on; it states no speed of propagation, and no continuum limit.

## 8. Pair presence and the state after a reading

Two observers can ask whether a pair was reported anywhere on the cell. In one instrument they keep
only the aggregate answer. In the other they also keep the public words from individual seams. Both
give the same first presence probabilities. Their disturbance of a later reading differs.

The construction uses the cell, native blank-memory reduction, degree references, fuel converter and
eleven-bit reader of SPEC.md, sections 1 to 7. Its supplied apparatus is a reversible tree connecting
the thirty reader blocks and a choice of which registers to copy. Port and seam labels address hardware.
The measured events below are functions solely of copied observer records. The source has 78 total-two
states and one residual fuel state. No fuel or doubled-load outcome is discarded.

**Archive circuit.** On fresh words apply the four stages of the primary reader to every seam. Encode
each public word (decision and triplet) into four flags by XOR into a fresh leaf. Agreed true sets the
first flag, mismatch the second, unavailable the third, and an unread or malformed word the fourth.
Agreed false sets none. Pair the leaves in the carrier's seam order, promoting an unpaired node
unchanged. At each parent XOR the bitwise OR of its two children into its target. Copy the root into a
fresh public register. The address switch also copies each seam's five public bits into its own archive
register. Run the tree, leaf encoders and all four reader stages backwards. Each operation is an
involution on its full register domain. On the ideal fresh preparation the scratch registers return to
zero for every source basis state. Linearity gives the same clean return on entangled inputs. Only
orthogonal record labels are copied. The local words must remain coherent until the copy and inverse
operations finish. Releasing a seam address to any additional archive changes the instrument.

One read allocates 330 reader bits, 236 flag bits at 30 leaves and 29 internal nodes, four root archive
bits and 150 seam archive bits. Both settings allocate all 720 bits. The campaign uses two fresh sets.
The tree has depth five. Including the reverse gates, the read has 21 layers, with one layer for the
archive copies. Its tree wiring, reversible control gates, blank workspace and pulse schedule are
apparatus choices. Actions and interaction strengths are in the construction's exchange units.

**Complete instrument.** For a source basis label $s$ let $f(s)$ be the actual final archive string.
The circuit acts as $|s,0\rangle\mapsto |s,f(s),0_{\rm scratch}\rangle$, up to the known source phase
described below. The Kraus projector for record $z$ is $\sum_{s:f(s)=z}|s\rangle\langle s|$.
This includes the complete fault alphabet. A mismatch on one seam and a positive report elsewhere set
both flags. Agreement with a common flipped predicate can give a positive aggregate on every input.
States whose archive strings coincide keep their relative coherence. These identities tensor with
every untouched reference.

For nominal words write $P_e$ for each of the thirty rank-one seam projectors, $P=\sum_eP_e$ and
$Q=I-P$ of rank 49, including fuel. Stripping off the common diagonal source phase, the two conditional
maps for a positive answer are
$$X\longmapsto PXP,\qquad X\longmapsto\sum_e P_eXP_e.$$
The second formula traces the retained seam address. The negative map is $X\mapsto QXQ$ for both.
This realizes the projective instrument characterized by the
[finite instrument theorem](../finite_instrument_calibration/PROOF.md) for the aggregate circuit. Its projectors sum to the identity, and every state supported on a single answer block is invariant
after removing the known phase. The resolved circuit passes all basis return tests, while the
uniform superposition of the thirty seam states has return probability $1/30$ in its positive branch.
These are mathematical probe states for the channel identity. The experiment uses the local emitter.

For the comparison read, the half-diamond distance between the two instruments with only the
first aggregate visible and the residual source retained is exactly $29/30$. On the positive block, averaging the thirty
unitaries that multiply seam $j$ by $\exp(2\pi i mj/30)$, $0\le m<30$, deletes exactly the
off-diagonal seam entries. The identity term has weight $1/30$, so the triangle inequality on the
remaining channels gives the ceiling $29/30$, also with a reference. The negative block agrees
and the positive input projection is trace nonincreasing. The uniform positive superposition
attains the ceiling: its pure output minus the resolved mixture has eigenvalues $29/30$ once and
$-1/30$ twenty-nine times. This optimized distance concerns a mathematical input. The emitted
history below gives a witness on public records alone.

**Finite service time.** Give each read total duration $w=1/16384$. During a layer of duration
$w/21$, implement each involution $G$ by the pulse $\pi(I-G)/(2w/21)$. Gates in each layer commute,
including readers that share a port, because their source controls are diagonal and their targets
are disjoint. The largest norm of a single prescribed pulse is $21\pi/w$.
Write the source as $V+K$, with $V=\kappa C$ and the uncoupled fuel energy. All apparatus gates
commute with $V$. The comparison read therefore executes the complete circuit and accumulates
$\exp(-iwV)$. With exchange also acting, Duhamel's inequality gives an operator error at most
$w\|K\|\le12w$ per read by section 6. Both unitaries act on the same enlarged workspace.
Two reads give at most $24w$ in total variation on the complete joint record, for each instrument,
by telescoping and contraction under measurement. Comparing two experiments costs at most $48w$.
This includes the probability of any physically reachable fault or workspace leakage absent from
the ideal service law, without conditioning on a successful comparison. The executable table gives
the diagonal-source comparison law. The pulse bound encloses the always-on law.

The first aggregate has the same marginal distribution in both always-on experiments, even if
exchange leaves the workspace entangled. Up to the copy layer the joint state is identical.
During that layer the optional seam archive pulse acts only on the seam words and their archive.
It commutes with the root copy on disjoint registers and with the source Hamiltonian. Its unitary
therefore cancels under the partial trace that defines the root marginal. Every subsequent gate
leaves that archived root untouched. The separation of the joint laws thus tests a difference in
the later reading given the same first record statistics.

**Observer histories.** The [protocol](presence_protocol.json) fixes the primary emitter on seam
$(0,1)$ at $g=64$, four units of free action before each reading, and a final nominal archive with
all seam words. These actions are supplied schedule parameters. At the first read the two choices
retain the same public aggregate. For comparison, forget only the first seam address after computing
the full resolved two-read distribution. At $\kappa=8$ the comparison calculation gives first
nominal presence probability $0.992329$ for either choice and total variation $0.440513$ between
the joint public laws. The common always-on first probability is within $12w$ of the exact
comparison probability.
It is $0.196148$ for exchange alone and $0.249676$ for reversed contrast. The first-seam mismatch
and unavailable controls each give $0.260336$ at $\kappa=8$. A common flipped value on all seams
gives $0.446650$. The [complete receipt](presence_results.json) contains the archive alphabets,
all joint probabilities on their reachable support and the residual fuel. A raw seam address or
source coordinate has no independent identification as physical space. A comparison with
laboratory data requires a matching pair of observation procedures.

**Certified public event.** Ask for presence in the first aggregate and a positive final record on
the emitting seam. The comparison law gives probabilities $0.0182640$ for the aggregate circuit
and $0.1029264$ for the resolved circuit with its first address hidden. For a separate error
certificate let $a$ be the exact emitted state after the first free interval and $b$ the emitting seam
basis state propagated backwards through the second free interval. For a first positive record
the intermediate read phase is common to all thirty seam states. The final read phase cancels in
the record probabilities. Put $z_e=\overline b_e a_e$. The two event probabilities are
$|\sum_e z_e|^2$ and $\sum_e|z_e|^2$. The twelve doubled states are unreachable from this emitter because
every native term preserves the multiset after the converter acts, and the converter only joins
fuel and the binary emitting pair. Fuel is included in both propagations.

The [integer certificate](presence_certificate.py) builds twice the binary-plus-fuel Hamiltonian
from the allowed swaps. It propagates complex integer pairs in units of $2^{-48}$. Each substep
has exponential argument of norm at most $1/2$, from the maximum absolute row sum of the symmetric
matrix. The degree-16 Taylor remainder has operator norm at most
$2(1/2)^{17}/17!$. In each recurrence step flooring the real and imaginary components costs at most
$\sqrt{2d}\,2^{-48}$ in vector norm for $d=67$. Since the preceding term error contracts by at least
a factor of two in the recurrence, summing the 16 term errors costs at most
$32\sqrt{2d}\,2^{-48}$ per substep. The program rounds this bound upwards. With norm below two for
every approximate input, twice the Taylor remainder suffices. The accumulated error bound below one
verifies that norm hypothesis by induction from a normalized exact input. Unitarity carries the
incoming error unchanged. The emitter uses the rational action $(355/113)/128$ and adds the norm
bound for its difference from $\pi/128$, enclosing $\pi$ by alternating arctangent series in
Machin's identity.

For input vector errors $\epsilon_a,\epsilon_b$, the scalar amplitude error and the Euclidean
error of the vector $(z_e)$ are each at most $\eta=\epsilon_a+\epsilon_b+\epsilon_a\epsilon_b$
by Cauchy-Schwarz. The probability error is at most $\eta(2+\eta)$ for each circuit.
The rational bounds put this error below $3\times10^{-7}$. After subtracting both numerical errors
and $48w$ from the event difference, the total variation between the two **always-on public
histories exceeds 0.0817**. This lower bound uses the emitted preparation and a single unconditional
joint record. The full-law floating-point estimate and the certified event have different scopes.

## 9. Binding energy and the motion of a cluster

A common relativistic dispersion would give \(E(p)=M c^2+p^2/(2M)+O(p^4)\) near rest.
The ratios of rest energies and inverse curvatures would agree. A common rescaling of clock or momentum
cannot alter this test. On the ring, \(k\) labels the phase of a prepared translation character. Identifying
it with a physical momentum requires an independent length and propagation measurement. The ring uses the
local law of section 1 with four fixed contrast references per port. The ring calculation in
[extras.py](extras.py) omits those references. They add \(4\kappa N\) on binary configurations of \(N\) loads.
The ring is a separate carrier choice. The following single-load identity also holds on the official federations.

A single nonzero load of value \(q=1,2\) evolves by \(L/2+6\kappa q^2 I\), as follows directly by moving
the load through one native swap. Under the identification of the two position bases, every free interval
differs by a scalar phase. Projecting the same position sets, coupling the same position-controlled meters
or adapting the next such service to a recorded outcome gives identical joint record laws by induction.
This includes the two successive cell readings on 20, 80 and 320 cells in the [counts](material_inertia/records/).
The two rest offsets have ratio four and the transport is identical. Adding a common chemical potential
\(\mu N\) leaves \(E_2=2E_1+12\kappa\) for the doubled load. This particular stable load therefore also fails
the common-mass interpretation whenever \(\kappa>0\) and \(E_1>0\).

### All cluster sizes on the ring

On binary configurations a homogeneous, reflection-symmetric, load-conserving Hermitian two-port term has
one real hopping matrix element and diagonal entries at \(00,01=10,11\). On the chain the sum, relative to
vacuum, has diagonal \(cN-U b\) and off-diagonal \(-t\) for one legal hop, where \(b\) counts adjacent occupied
bonds. Here \(t>0\), up to the usual staggered phase convention on a line. The native law gives
\(2t=1,\ c=6\kappa+1,\ U=1+2\kappa\). A common per-load shift changes only \(c\). A bound cluster at rest
in the attractive easy-axis family has \(\Delta=U/(2t)=\cosh\eta>1\).
Its infinite-chain energy is

\[
 E_N(k)=(c-U)N+
 2t\,\frac{\sinh\eta}{\sinh(N\eta)}\bigl(\cosh(N\eta)-\cos k\bigr).
\]

This is the ferromagnetic XXZ droplet formula in the source normalization. For an explicit substitution,
order the occupied ports \(x_1<\cdots<x_N\), write \(g_j=x_{j+1}-x_j-1\), and use the amplitude
\(e^{ikx_1}\prod_{j=1}^{N-1}a_j^{g_j}\). Put \(a_0=e^{ik}, a_N=1\) and
\(a_j=[e^{ik}\sinh((N-j)\eta)+\sinh(j\eta)]/\sinh(N\eta)\). At a touching pair the two forbidden hops
and the contact energy cancel by \(a_{j-1}+a_{j+1}=2\Delta a_j\). Evaluating at consecutive occupied
ports gives the displayed energy. The internal coefficients have modulus below one. At zero momentum the
relative-coordinate wavefunction is strictly positive and is the ground droplet. The thermodynamic
ground-energy result and Bethe construction are given by
[Nachtergaele, Spitzer and Starr](https://doi.org/10.1007/s00023-006-0304-6).

Thus the kinetic mass ratio \(R_N=E_1''(0)/E_N''(0)=\sinh(N\eta)/\sinh\eta\) is greater than \(N\) for
every \(N>1\). The binding energy \(N E_1(0)-E_N(0)\) is positive because
\(\tanh(N\eta/2)<N\tanh(\eta/2)\). For any positive constituent rest energy the rest ratio is therefore
**below** \(N\). This rules out a common relativistic mass for these bound clusters for the entire stated
two-port family. A site energy linear in load cannot reverse the inequalities. The gapless threshold at
\(\Delta=1\) has no normalizable two-load bound state at rest. These statements concern the line and its
finite-ring tests. They establish no obstruction for all OPH realizations.

The ring instrument compares the nonzero shifts \(E_1(k)-E_1(0)\) and \(E_N(k)-E_N(0)\) at the same
enrolled translation character. The ratio equals the curvature ratio for the infinite-chain cosine bands
above. On a finite ring it is a finite difference. In particular the observed ratios for the three-port
and controlled-memory laws do not give a zero-momentum derivative without a separate limit.

The predicted kinetic ratios satisfy \(R_1=1,\ R_2=2\Delta,\ R_{n+1}=R_2R_n-R_{n-1}\), with no separate
species coefficients. In particular \(R_3/R_2=R_2-1/R_2\).

### Predicting a position distribution from a binding measurement

The binding reader compares a compact pair at \(k=0\) with two independently evolving uniform singles.
The interferometer switches which registers evolve, then recombines the paths. Both active arms carry the
same total load, so every common per-load energy shift cancels. The isolated one-load dimer has level
splitting \(2t=1\), which sets the energy unit. In that unit \(B_2=(\Delta-1)^2/\Delta\).
The positive attractive root determines \(\Delta\) from \(B_2\). For the finite ring, take the positive
relative amplitude \(z^{r-1}+z^{L-r-1}\), where \(z=1/\Delta\). Its only residual is at contact. The
componentwise ratios of the finite Hamiltonian applied to the trial vector bound its lowest energy,
giving \(0\le B_{2,L}-B_2\le(1-z^2)z^{L-3}/(1+z^{L-2})\). The [reader](material_inertia/prediction.py)
first uses the upper bound one, then tightens it using only the inferred lower bound on \(\Delta\).

For a translated compact cluster, the cosine band predicts the Fourier inverse of
\(\exp[-i(1-\cos k)t/R_N]\) for its translated return amplitudes, up to a common phase. The terminal
position instrument has \(L\) occupied-run outcomes and a complementary outcome for every configuration
which is not one run. The source replay propagates every eigenmode in every translation sector. A second
law adds a positive three-port term \(s(|110\rangle-|011\rangle)(\langle110|-\langle011|)\).
For a pair at \(k=0\) its action vanishes on the entire sector, so the two laws have identical binding
inboxes at every wait and on every ring in the experiment. In the two-body recurrence it replaces \(U\)
by \(U-2s(1-\cos k)\) in the binding boundary condition. In native units the pair curvature becomes
\(1/(2\Delta)+2s(1-\Delta^{-2})\). Setting this equal to \(1/[E_2(0)/E_1(0)]\) determines \(s\) from the
source. Agreement at pair curvature is a design choice. Three- and four-load motion independently tests it.

A finite-volume error bound uses only the binding inbox and the native family. Periodize the ordered
droplet amplitude by summing over the \(N\) cyclic choices of its first occupied port, with \(e^{ikL}=1\).
For one cut the infinite-chain recurrence fails only when the two ends touch across that cut. All \(L-N\)
holes then lie in interior gaps. Write \(\rho=2e^{-\eta}\), bounding every \(|a_j|\).
At such a configuration the residual of that term is at most \((\Delta+\rho)\rho^{L-N}\) in dimer units.
At each of the \(L\) compact configurations one term has unit modulus; all other terms together have
modulus at most \((N-1)\rho^{L-N}\). Dividing the residual norm by the resulting lower bound on the
state norm gives

\[
 \epsilon_N\le\sqrt{\binom{L}{N}/L}\,
 \frac{N(\Delta+\rho)\rho^{L-N}}{1-(N-1)\rho^{L-N}} .
\]

The bound applies with positive denominator and \(\rho<1\). In each character the compact orbit has
dimension one. In its complement there are at least two occupied runs, so the native contrast costs at
least \(4\kappa N+4\kappa\), and the exchange is positive. Interlacing places the second eigenvalue above
this threshold. The infinite band lies at most at the compact diagonal \(4\kappa N+\Delta\) for \(N\ge2\).
When \(\epsilon_N<\Delta-2\) the certified eigenvalue is the finite ground line. The compact-to-complement
coupling has norm at most one by the two end hops of magnitude \(1/2\). The ground-state block equation
therefore bounds the non-ground spectral weight of the compact input by \(1/[1+(\Delta-2)^2]\).
These estimates use \(N<L/2\).

Let \(\delta\phi\) bound the error in relative band phase across the characters, including two wrap
errors, the clock calibration and the inferred \(R_N\) interval. The normalized state differs from the
ideal compact Fourier flight in norm by at most
\(\sqrt{\delta\phi^2+4/[1+(\Delta-2)^2]}\): decompose the input into its ground component and its
orthogonal complement separately in each character. Total variation for the complete position read is
bounded by this norm plus the finite-counter rounding error. For the enrolled 19-port flights this
binding-only bound is below 0.096 for every size and wait. The observed discrepancy is at most 0.00214.
The five-load low-momentum frequency is unresolved at the enrolled precision, although the separate
large-momentum and direct position instruments resolve its motion. [Sealed forecast and held-out
kinematics](material_inertia/records/challenge/) and [position measurements](material_inertia/records/flight/)

### What the written memory can change

Let \(W\) be the official nine-to-ninety contrast write. A Hermitian operation confined to the ten
memory labels preserves the written code exactly when its tensor product with the identity on repaired
loads commutes with \(WW^\dagger\). The induced load action is the compression by \(W\). Solving the
commutator equations over the rationals using \(\sqrt2 W\) gives complex dimension 18 for the label
commutant and dimension 9 for its image. The latter is \(\mathbb C\oplus M_2\oplus M_2\). The common
scalar acts on the three equal-load states. The same two-by-two block acts on the orientations of
\(01\) and \(12\), and the other block acts on \(02\). The rational calculation verifies this form on
every basis solution and verifies rank nine, proving both inclusions. The vacuum and doubled-background
inputs constrain the equal-load label even in a test intended for the two-load sector. Thus this entire
class preserves the multiset. Its homogeneous reflection-symmetric binary restriction with positive
hopping and \(U>2t\) is covered by the ring obstruction above. The proof assumes exact preservation on each individual write,
memory action and unwrite. [Exact algebra](material_inertia/memory_algebra.py)

One output predicate suffices to break multiset conservation. Let \(m_d\) be the normalized sum of the
two coin labels at contrast \(d\). Apply the positive square of \(m_0-(m_2+m_{-2})/\sqrt2\) conditioned
on repaired output \(11\), then unwrite. On the nine blank inputs the logical term is the square of
\(|11\rangle-(|20\rangle+|02\rangle)/\sqrt2\). The identity \(GW=WV\) is checked exactly on all
columns, so the induced unitary returns the private memory to blank at every duration. Removing the
output predicate gives squared generator leakage one on the vacuum. This construction uses the native
registers with a changed memory generator. Its control blocks only the conversion between distinct and
doubled loads and retains the diagonal costs and hopping between doubled endpoints.

On the 20-cell federation the doubled preparation produces up to 797 agreed seam-pair outcomes per
1,048,576 simulated preparations. Starting from the emitting seam gives up to 316 doubled-location
outcomes. Both matched controls give zero throughout. The original public seam predicate is used with
an ideal switched read. The supplementary doubled-location meter resolves part of its no-pair outcome.
All 28,920 two-load configurations participate. The controlled memory gives a pair rest ratio of
1.707 and kinetic ratio of 24.94 on the 16-port ring. The same registered comparison as for the
native law puts both two- and three-load ratios outside tolerance on both the 12- and 16-port rings.
[Memory and detector results](material_inertia/records/memory/)
