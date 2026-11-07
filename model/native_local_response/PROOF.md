# Local response from complete native repair records

## Source and observation domain

The local source is the admitted weighted binary H66 program with the
literal C/E endpoint instrument. The connected source is the admitted fair
`port_pair` federation at the evaluated levels0,1,2. Its vertex weights
are all one. The weighted connected family in `weighted_source_family`
has the same local operator proof but is a separately declared extension;
this result does not admit it or attach its stationary law to fair dynamics.

The proof uses finite products of the original endpoint kernels. A source
configuration is a binary occupation mask, with its original global total
preserved. Embedding the fixed-number sector into the binary tensor space
is a mathematical representation of those same configurations. A restriction
to a vertex set keeps every local number sector; it does not select a local
total-two sector or add a physical source register.

Selectors are fixed schedules or independent draws from a fixed seam law.
They do not read source loads, FIFO values or archived record values.
The source starts with its ordinary empty FIFO and retains the actual
archive and FIFO updates. Memory either exports its oldest record or has
capacity for every attempted command. The locality adapter rejects a
possible intrinsic FIFO halt and every record-adaptive selector. The
original source laws and admission rule are unchanged.

Specify a final retained vertex set A and the command instances whose
complete native raw records are observed. The observation map retains
those records, their global command positions, final A and an arbitrary
untouched quantum reference R. It traces the other raw record values and
the terminal exterior source. Input-independent selector metadata and
elapsed global attempts may also be retained with the same labeling.
The root set used in spatial bounds contains A and **both endpoints of
every observed state-dependent record event**. Merely intersecting A is
insufficient. A global archive observation is a different output map.
Nothing forbids an observer from reading the source or scheduler
specifications, or from choosing that different output map.

## Exact backward dependency cone

For a deterministic seam schedule, initialize the backward support to A.
Read the commands in reverse order. Include an observed event and add
both endpoints. Include an unread event when its endpoints meet the
present support and then add both endpoints. Otherwise discard that
event only after summing its complete raw instrument. This constructs a
conservative support cone, not a minimal-support theorem.

For every native raw slot h on e=(a,b), its Kraus map has the form
\(K_{e,h}=L_{e,h}\otimes I_{e^c}\). Its domain, current or contrast,
placement and scalar coin weight depend only on a,b. The complete25-slot
fixed-seam instrument satisfies \(\sum_hL_{e,h}^*L_{e,h}=I_e\), including
zero law weights, hold and reserved raw addresses. The executable binding
compares all750 fixed raw maps on all66 H66 basis columns with the
unmodified interpreter. Four endpoint patterns on an actual intercell
seam supply the separate fair-federation attachment. The maps are scalar
partial injections with no hidden phases, so domain/image/weight equality
determines every complex matrix unit.

If e is disjoint from the support S of a pulled-back local record/source
effect, the unread Heisenberg channel gives

\[
 \sum_hK_{e,h}^*(O_S\otimes I)K_{e,h}=O_S\otimes I. \tag{1}
\]

An observed branch or an intersecting event can enlarge support only by e.
Induction over the reverse schedule proves that the entire observed
classical-quantum output is a channel of the reduced initial cone/reference
state. Equivalently, if C is the computed initial cone, there is a channel
\(\Lambda\) made from the retained original endpoint instruments such that

\[
 \mathcal O\mathcal I(\rho_{VR})
 =\Lambda(\operatorname{Tr}_{V\setminus C}\rho_{VR}). \tag{2}
\]

Equation(2) includes arbitrary initial source/exterior/reference
correlations and all complex coherences. Its exact finite checks compare
an independently summed native expanded execution, the full sparse mask
instrument and the cone-restricted instrument, retaining all625 observed
two-record slots from15625 original three-command histories. Real and
imaginary operator components are checked separately. The finite theorem
`ptraceFst_local_kraus` and its entry identity `liftLeft_conj_apply` in
`Locality/NoSignalling.lean` provide the unread-channel partial-trace step.
The owning native disjoint-record proof supplies the compatible complete
endpoint commutation identity; it does not erase chronological records.

Remote records exist throughout actual execution. Their values can affect
a controller using the single global FIFO, even across disjoint seams.
The actual `native_causal_dataflow` witness gives receiver probabilities
1/4 and0 for two source inputs with identical receiver data. Its controller
and storage ancestors are real dependencies. That program is rejected by
this spatial adapter. If a future protocol retains or reads such messages,
their source events must join its dependency cone; no graph-distance
bound here applies after silently deleting them.

## Finite-volume local multi-time comparison

Let B contain all observation roots. Compare the actual process with a
control that uses the same global selector choices and clock, applies
the original endpoint instrument on seams inside B, and acts as identity
on B for crossing or exterior choices. Remote raw values are absent only
from the declared observation map. They are not relabeled as local
successes. The comparison is an exact local channel on the reduced B
state, with internal rates inherited from the original full selector;
it is not a renormalized smaller native universe.

Couple the two processes using the same input-independent selector tape.
On every tape whose backward cone stays in B, equation(2) gives the same
complete local output channel. Let p_escape be the probability of the
complementary tapes. Each complementary output has trace p_escape for
every normalized input, so the trace-norm triangle inequality gives

\[
 D(\mathcal O\mathcal I(\rho_{VR}),
       \mathcal I_B(\rho_{BR}))\le p_{\rm escape}. \tag{3}
\]

This is also a half-diamond bound. It includes every selected local record,
all local back-action and arbitrary R. It is an upper bound on output
distinguishability, not an equality. If another initial BR state differs
by joint half trace distance epsilon, the bound becomes
\(\min(1,\epsilon+p_{\rm escape})\). A diagonal unreferenced canonical
limit does not imply that joint premise for an arbitrarily correlated R.

The backward cone under independent seam probabilities q_e is an exact
finite subset process: e adds its other endpoint when it intersects the
present support. This probability concerns the conservative process started
with all observation roots at the final time. It can exceed the event-timed
support obtained by adding record roots only when their events are reached,
and neither support is claimed minimal. Escape from B is absorbing. `escape` computes this
probability with rational arithmetic and an explicit subset-state cap.
No endpoint Hilbert matrix is constructed.

For a ball of radius r around a root set S, escape requires at least r+1
time-ordered edge choices. If the graph has maximum degree Delta and
\(q_{\max}=\max_eq_e\), then at T global repair attempts

\[
 p_{\rm escape}\le
 \min\{1,|S|\binom{T}{r+1}(\Delta q_{\max})^{r+1}\}. \tag{4}
\]

Each first-infected vertex has an earlier parent, so an escaping ancestry
contains an ordered path. There are at most \(|S|\Delta^{r+1}\) prefixes
of the required length and \(\binom{T}{r+1}\) possible time choices.
The union bound proves(4). The sharper `walk_bound` sums all first-exit
walk lengths: with weighted adjacency W on B and
\(b_v=\sum_{u\notin B}q_{vu}\), its coefficients are
\(c_j=1_S^TW_B^{j-1}b\), and
\(p_{\rm escape}\le\min(1,\sum_{j=1}^T\binom Tj c_j)\).
The recurrence uses graph vectors. Walks can overcount actual ancestry;
the exact subset process supplies an independent finite control.

These rates count actual repair choices. The explicit global bounds and
the following federation calculation use the source's repair probability
one. A separate hold probability can be handled by weighting the edge
rates by that probability; it is not silently omitted from the clock.

## A four-state halo bound on the actual federation

Take S to be every port of one native cell and B its one-edge halo. The
retained record events for this cell bound have both endpoints in S;
retaining crossing-record values requires their exterior endpoints as
additional roots and the corresponding larger-root calculation. The
literal graph has30 internal seams and3 intercell seams per cell. Its
halo contains three external ports in three distinct neighboring cells.
Each has its unique link to S and five links outside B; no two halo ports
are linked. These ownership and adjacency statements are checked for
every cell at each evaluated level0,1,2. The total number of seams is
\(M=630\,4^\ell\) and the maximum vertex degree is6.
The source attachment in this package is restricted to those three
evaluated levels. Equations(5)–(8) apply to any finite M whose graph has
the stated halo properties; applying them at other formal refinement
levels requires that geometry premise separately. The420 finite cell
checks are not a proof about unevaluated refinement levels.

Before first escape, only the subset of these three halo ports can be
added to the backward cone. If n have been added, its next attempt adds
another with probability \((3-n)/M\), escapes with probability \(5n/M\),
and otherwise holds. Thus the exact surviving four-state column kernel is

\[
 Q_{nn}=1-\frac{3+4n}{M},\qquad
 Q_{n+1,n}=\frac{3-n}{M}\quad(0\le n<3). \tag{5}
\]

Its total escape probability after T attempts is
\(1-1^TQ^Te_0\). An escape has one of15 ordered two-edge witnesses:
one of the three cell-to-halo seams precedes one of its five outward
seams. Therefore

\[
 p_{\rm escape}(T)=1-1^TQ^Te_0
 \le\min\{1,15\binom T2/M^2\}. \tag{6}
\]

The15 events overlap; their sum is a union bound. This improves the
any-cell-boundary-hit bound for the local output protected by the halo.
For level0 and T=30, the exact cone escape is less than0.015 and the
ordered-pair bound is29/1764. The proof bounds the complete local
multi-time instrument, not only mean occupation.

## The cell's internal repair clock

Stop immediately after the kth choice of an internal cell seam. No waiting
interval is added after that final record. Its probability per global
attempt is a=30/M. The elapsed global attempt T is retained as an
input-independent clock record. No conditioning on a rare T or on source
outcome values is used.

A forward selector word has k independent geometric failure blocks,
each followed by one internal event. Reversing the word puts the last
internal event first and the initial failure block last. The probability
of a failure word is \(a\prod q_e\), invariant under reversal. Hence
the reversed failure blocks have the same independent distribution.
Every internal event has endpoints in the entire-cell root S and acts
as identity on the cone. The final pre-first-internal failure block is
included, so there are exactly k failure blocks in the reverse scan.

The probability of a surviving noninternal self-loop at state n is
\(1-(33+4n)/M\); the birth probability is \((3-n)/M\). Summing one
geometric block gives the substochastic matrix

\[
 R=30(D-B_0)^{-1},\qquad
 D_{nn}=33+4n,\quad (B_0)_{n+1,n}=3-n. \tag{7}
\]

Consequently the exact cone escape probability through the kth internal
repair is \(1-1^TR^ke_0\). The factor M cancels. For k=1 it is
677/50061; this is less than1/60. For a simple general bound, let L count
noninternal choices before the kth internal choice. Its negative-binomial
factorial moment gives
\(\mathbb E\binom L2=k(k+1)(1-a)^2/(2a^2)\). Conditional noninternal
choices have probabilities \(q_e/(1-a)\). Summing the15 ordered
boundary/outward pairs yields

\[
 1-1^TR^ke_0\le\min\{1,k(k+1)/120\}. \tag{8}
\]

This is an unconditional bound for the complete local output including
T. It does not say that T has a capacity-independent distribution, or
that complete observer outputs at different capacities are identical.
It gives no vanishing interaction error as capacity grows at fixed k.
For a generic window, a finite cutoff T gives the alternative bound
\(p_{\rm escape}(T)+\Pr[\operatorname{Binomial}(T,a)<k]\), capped at1,
using matching explicit timeout outputs.

## Bounded record response

For any retained-record/source/reference score of norm at most one, let
\(C_k,C'_k\) be its expectations in the actual and halo comparison
processes at the internal repair count k. Equation(3) gives
\(|C_k-C'_k|\le2(1-1^TR^ke_0)\). The score may depend on the complete
selected local record history, including T. For a declared analysis
parameter \(0\le q<1\), absolute convergence gives the Abel-response bound

\[
 \left|(1-q)\sum_{k\ge0}q^k(C_k-C'_k)\right|
 \le2\left[1-(1-q)1^T(I-qR)^{-1}e_0\right]
 \le\min\left\{2,\frac{q}{30(1-q)^2}\right\}. \tag{9}
\]

The four-by-four rational triangular solve checks this formula directly.
Here q is an analysis weight. The infinite series describes the mathematical
kernel family, not an admitted unbounded execution. It identifies neither
an energy, a physical time, a particle pole nor a scattering observable.

## Execution and proof boundaries

Native admission permits at most \(2^{32}\) attempts. The expanded host
evaluator has64 steps and65536 branches as separate ceilings. A literal
four-command fixed schedule returns `RESOURCE_STOP` after three steps;
this is retained and tested. Equations(1)–(9) state finite-kernel identities
and analytic bounds. They do not call an unevaluated or stopped run a
completed observer experiment. The exact subset/four-state evaluators
bound T and k by256, and the general walk evaluator bounds T by64.
Their mathematical inequalities hold for every finite T or k under the
stated source conditions. Simple polynomial bounds need no astronomical
rational power at large declared capacity.

The source attachment and all-domain proofs here are mathematical prose
with exact executable controls. No additional Lean theorem is claimed.
The result supplies controlled local response on the source graph. It
does not identify that graph with physical three-dimensional space,
choose a Lorentzian or Hamiltonian clock, derive an interacting continuum
theory, or assign measured particle masses.
