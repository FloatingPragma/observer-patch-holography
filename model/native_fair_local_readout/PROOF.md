# Local readout from complete fair repair records

A finite connected set of binary ports has a native record-only occupation
readout with an explicit erasure outcome. The instrument acts on every local
number sector, including coherent superpositions compatible with the global
conserved total. It retains the terminal disturbance, complete raw records
and an arbitrary untouched quantum reference. Its error is a probability of
coverage or constant-input erasure, rather than a renormalized success law.

## 1. Source and output

Use the declared fair federation with $V=240\,4^\ell$ ports,
$M=630\,4^\ell$ seams and a declared global conserved occupation $q$.
Every attempt selects a seam independently and uniformly, chooses either
placement coin fairly, and applies the original E repair with repair
probability one. This is a literal `native_program` repair command; it does
not insert a Hamiltonian or coherent gate. The per-seam attempt clock is
$t=T/M$. Population exchange time for the same uniform attempts is $T/(2M)$.

Let $B$ be a connected set of $b\geq2$ ports and choose a deterministic
spanning tree $\mathcal T$ of its induced seam graph. Observe **every** repair
record incident to $B$, including boundary seams, in chronological order.
The selector, ordered endpoints, coin, original raw address, contrast and
current are retained. The terminal quantum state on $B$ and an untouched
reference $R$ are part of the output. Any remote source ports that undergo
repairs are source systems, rather than untouched references. Retaining
additional native records or the full terminal source gives the same
statements below. The source's FIFO has one slot and exports the oldest
record; no archive entry is discarded by the source. A specified observation
map may trace remote values and exact tick labels.

The raw alphabet has $1+24M$ addresses per uniform attempt. A fixed-seam
attempt has 25. Zero hold, C-law and unused E addresses are retained as zero
maps. The decoder reports `ZERO_RAW_WORD` when an incident zero address is
queried or a complete tree word has no feasible binary input. The latter
includes contradictory sequences of individually valid E slots, with or
without the optional global-number refinement. Their complete raw records
and zero quantum blocks are retained. An unnormalized zero-map query is
not a positive-probability run.

`source_plan` binds all incident seams against the literal graph through
level two. `rooted_plan` binds a finite certified routed window and checks
that every neighbour of each retained whole cell occurs in the outer
window. The returned source addresses and incident edges therefore include
the boundary. The all-level original-cell internal tree is also admitted.
No unchecked limit of floating port routing is inferred from these finite
attachments.

## 2. Reconstruction with moving boundary occupations

Write $n_i^0\in\{0,1\}$ for an initial local occupation and $s_i$ for its
accumulated displacement. At an observed seam $(a,b)$, the record gives
incoming contrast $d=n_a-n_b$ and current $j$. Before updating the shifts,
the first occurrence of an internal tree edge gives

$$n_a^0-n_b^0=d-s_a+s_b. \tag{1}$$

Update $s_a\leftarrow s_a-j$ and $s_b\leftarrow s_b+j$ for endpoints in
$B$. Boundary records participate in these updates. Their outside initial
occupations are not supplied to the decoder. Repeated internal records are
also retained; only the first hit is needed to acquire each tree equation.

When every tree edge has occurred, integrate (1) from a root to obtain
$n_i^0=c+u_i$, with root value $c\in\{0,1\}$ and $u_{\rm root}=0$.
If some $u_i\ne0$, exactly one choice of $c$ can be binary at all ports.
The decoder returns that compact initial mask and its current-updated
terminal mask. If all $u_i=0$, it returns `CONSTANT_ERASURE`; this is the
initial empty/full local subspace. Boundary records can resolve some such
histories, but this decoder deliberately retains the uniform erasure rule.
If a tree edge has not occurred, it returns `COVERAGE_ERASURE`.

These conclusions hold on every nonzero native history. They follow by
induction from the actual binary endpoint update. `decode` validates all
retained typed raw records and binary consistency for a reported nonconstant
mask. Completeness of the incident record list requires the source-bound
history adapter; arbitrary omitted input records are not authenticated by
the decoder alone. The tests exhibit an omitted-boundary-current example
that changes an inferred initial mask.

The reconstruction extends the complete connected-tree calculation in
`native_sector_readout`: a known global total can resolve its whole-graph
empty/full ambiguity, whereas a proper local window has no supplied local
total. The decoder does not read hidden source occupations. A supplied
input mask is used only in independent finite verification controls.

## 3. Full quantum instrument and exact erasure distance

Each fine native history operator $K_h$ is a nonnegative scalar times a
partial injection of occupation basis states. Tensoring it with $I_R$
retains arbitrary reference coherences. For a decoded history with initial
local mask $x$, equation (1) implies

$$K_h=K_hP_x^B. \tag{2}$$

For a coverage-complete constant history,
$K_h=K_h(P_{\emptyset}^B+P_{\rm full}^B)$. It need not separately measure
these two projectors. Empty/full coherence can survive, with correlated
exterior occupations ensuring the same global $q$.

Define $\mathcal D$ as the actual complete instrument with its decoder flag.
Define the mathematical comparator $\mathcal O$ by an initial projective
measurement $\{P_x^B\}$, followed by the **same** native instrument, with
true initial label $x$ appended. Its terminal backaction is that of this
specified initial measurement and the original repair history. It is not
an extra implemented source service. Both channels share the full history
label; $\mathcal D$ uses a distinct erasure flag on failures.

On every successful block (2) makes the two outputs identical. On a failure
block, the actual erasure flag and oracle true-label flags have orthogonal
classical support. Both positive blocks have equal trace: $K_h^*K_h$ is
diagonal, so the initial projective measurement changes no fine-history
probability. Their half trace distance equals that common trace, regardless
of off-diagonal source/reference terms. Thus for every input density
$\rho_{\rm source,R}$,

$$\tfrac12\|\mathcal D(\rho)-\mathcal O(\rho)\|_1
 =\gamma=p_{\rm cov}+(1-p_{\rm cov})\eta,\qquad
\eta=\operatorname{Tr}[(P_{\emptyset}^B+P_{\rm full}^B)\rho_B]. \tag{3}$$

The selector word is independent of the input. Its total instrument is
trace preserving for each fixed selector word, and the all-zero recovered
differences characterize precisely the initial constant subspace on
coverage. This proves the last equality, including correlations with $R$.
Tracing other output systems gives an upper bound by $\gamma$; equality
need not persist after such further coarse graining. Equation (3) is a
state-dependent complete-channel comparison, not a state-independent
small diamond error: the worst-case input can have $\eta=1$.

For a fixed local tree executed once with no intervening boundary repairs,
$p_{\rm cov}=0$. Its erasure error is exactly $\eta$ in the full labelled
comparison. A sector-conditioned preparation supplies its actual initial
sector weights; it does not make $\eta$ small for every input.

## 4. Coverage and canonical local number mixture

Put $k=b-1$. Independent uniform selectors give

$$p_{\rm cov}=1-\sum_{j=0}^k(-1)^j{ k\choose j}
 (1-j/M)^T
 \leq\min\{1,k(1-1/M)^T\}
 \leq\min\{1,k2^{-\lfloor T/M\rfloor}\}. \tag{4}$$

The exact formula is inclusion-exclusion for missing chosen tree seams.
The first bound is a union bound. For $M\geq2$,
$(1-1/M)^{-M}\geq1+M/(M-1)>2$ by Bernoulli's inequality; $M=1$ is
immediate. A necessary single-edge hit also gives
$p_{\rm cov}\geq\max\{0,1-T/M\}$. The code computes exact finite
inclusion-exclusion under its displayed rational arithmetic budget and
checks both inequalities. Larger budgets use the certified rational bounds.
No conditioning on coverage, a rare duration, or a selected output is used.

If the global source marginal is the fair canonical fixed-$q$ state, its
local number mixture gives

$$\eta=\frac{{V-b\choose q}+{V-b\choose q-b}}{{V\choose q}}
 =\frac{(V-q)_b+(q)_b}{(V)_b}
 \leq(1-q/V)^b+(q/V)^b. \tag{5}$$

Out-of-range binomials are zero. The falling-factorial formula avoids
allocating a global carrier. Only the source marginal enters (5); an
arbitrary extension correlated with $R$ has the same erasure weight and is
covered by (3). This does not assert closeness of that joint state to a
product reference. The declared density $q/V$ and canonical preparation
are separate premises.

As a finite certificate, a routed 348-port radius-four detector inside a
radius-five incident window, at level six with $q=V/3$, has $k=347$.
Taking $T=64M=165150720$ admitted uniform E attempts gives
$p_{\rm cov}\leq347/2^{64}$ and
$\eta\leq(2/3)^{348}+(1/3)^{348}<10^{-60}$.
The total error in (3) is below $1/10^{16}$. This is an admitted mathematical
instrument and certified finite routing attachment; the expanded host does
not execute that many branches or level-six geometry.

### Declared global-total feasibility

The global conserved total is public source specification. It can improve
readout without supplying an unobserved local total. Empty $B$ is feasible
exactly when $q\leq V-b$, and full $B$ is feasible exactly when $q\geq b$.
After complete tree coverage, the decoder checks these inequalities. If
only one constant input is feasible, its label is returned, including the
terminal mask obtained from all incident currents. If neither is feasible,
a constant record word has zero support in that declared global sector.
Valid original raw words with no feasible input in the declared sector
are retained as `ZERO_RAW_WORD`, with zero quantum block and their complete
records. They are not rejected as malformed descriptions.
For a nonconstant inferred mask it also checks
$\max(0,q-(V-b))\leq |x|\leq\min(q,b)$.

The unresolved constant weight is therefore

$$\eta_{\rm amb}=
 \mathbf1_{\{b\leq q\leq V-b\}}\,
 \operatorname{Tr}[(P_{\emptyset}^B+P_{\rm full}^B)\rho_B]. \tag{6}$$

For $0<q<V$, equation (3) holds with $\eta_{\rm amb}$ in place of $\eta$.
For $q=0$ or $q=V$, the complete global occupation is known from the
specification; the local label is resolved even without coverage and the
oracle distance is zero. Raw records and terminal state are unchanged by
this classical postprocessing. These statements use no inference of a
generic local number from the global one. In the canonical example above,
both constant inputs are feasible and the displayed bound is unchanged.
For a twelve-port original cell and declared global $q<12$, the fixed
eleven-repair tree therefore reads every feasible local number sector
exactly, including the empty cell. At $q=2$ these are local numbers zero,
one and two. This admitted readout length is independent of source level;
the host geometry limitation is separate.

`decode` accepts optional exact `population,total` arguments.
`raw_query(...,use_global_total=True)` obtains both values directly from
the authenticated source programme, rather than trusting a separate local
number label. The default erasure decoder and its bounds are also retained.
The complete-source controls include globally sparse and nearly full
sectors with identical quantum outputs before and after this label
refinement.

## 5. Source admission, query execution and one clock

The fixed-tree programme uses $b-1$ literal fixed E commands. Original-cell
internal seams are admitted at every declared level. Above level two the
source validator rejects fixed intercell addresses; a multicell tree must
therefore not be silently submitted as a high-level fixed word. The
uniform alternative uses one existing repeat command and the actual FIFO
and resource calculation. Its $T$ attempts obey the same $2^{32}$ ceiling
as every admitted native programme, together with the $2^{40}$ retained-record
bit ceiling. `uniform_ceiling` computes the smaller actual one-word limit
and checks admission at that value and rejection immediately above it.
The formal finite-kernel statements
hold at arbitrary finite $T$ with the stated graph, but that quantifier is
not a claim of admission beyond this ceiling.

The host's expanded evaluator is limited to geometry through level two,
64 steps and 65536 branches. It preserves `RESOURCE_STOP` and
`RESOURCE_STOP_UNEVALUATED_CARRIER`. `raw_query` evaluates one bounded
unnormalized complete word through the owning continuation kernel. It
checks the source programme, complete incident-edge attachment, typed raw
slots and full FIFO/archive. It is neither a normalized sampler nor a
replacement for a stopped expanded execution. The receipt replays actual
one-step uniform and eleven-step fixed-word controls.

For fixed positive per-seam duration, $T=Mt$ grows with capacity. The
single-programme ceiling therefore obstructs an all-capacity admitted
coverage guarantee. At level 200, every admitted uniform word satisfies
$p_{\rm cov}\geq1-2^{32}/(630\,4^{200})$. No implicit sequence of host
restarts or enlarged counter is assumed.

The decoder needs chronological incident records and the declared attempt
budget, but no exact individual timestamps. It can therefore be applied
after an explicitly proved binned local-history observation map. That
map's root support must include all endpoints of observed boundary seams,
not just $B$. All counted preparation attempts and any binned-clock
reservation error are charged by their owning clock theorem. Full exact
native timestamps are not identified with continuous Poisson timestamps.

## 6. Measurement-algebra boundary

For bare repair/hold commands, each complete history operator is monomial;
its effect $K_h^*K_h$ is diagonal in occupation. Record-adaptive selection
preserves this property because every fixed record history chooses one
such product. Hence equal occupation diagonals give identical complete
record laws. This is an operational restriction on that record menu. It
does not identify the terminal quantum states or exclude additional
admitted coherent services.

The separate proper-code source-word and occupation-register services have
specified input registers, loading, failure and pulse semantics. Applying
them to a growing carrier requires an authenticated local encoding and a
complete operation intertwiner, including all number and failure sectors
and its clock cost. A basis chart alone is not an implementation of that
service. This package neither assumes arbitrary local gates nor proves
that every such service attachment is impossible. Its constructive
readout uses only the literal repair grammar.

The finite locality identities in
`Lean/ObserverPatchHolography/Locality/NoSignalling.lean` and the owning
`native_sector_readout` partial-injection instrument proof apply to the
source attachment above. The tensor-reference and complete-output
statements retain their hypotheses. No physical dimension, relativistic
clock, quantum field vacuum or observed particle spectrum follows from
these occupation records.
