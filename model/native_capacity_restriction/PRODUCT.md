# Exact dynamic restriction of independent native products

The native interpreter admits finite independent products of H66 or H78
cells. Appending unused factors supplies an exact dynamic continuation of
its repair/hold programs, including the biased H66 program. This result
uses the existing independent-product primitive. It introduces no
coupling, endpoint law or field Hamiltonian.

The domain is the native program grammar of
`model/native_program/interpreter.py`: repair and hold commands, all
admitted fixed or record-dependent selectors, and the declared FIFO
export/halt memory. Coherent `source_word` services are a different grammar
and are outside this theorem. In particular, the native acquisition,
preparation and readout subprograms of `record_clock_member` can be lifted;
its compiled coherent proper-code services and their composition are not
certified as a growing observer here. Such an extension needs its own
factor-local encoding and service proof.

## Actual program construction

Let `P` be an admitted local program or an admitted independent-product
program acting on the first `r` cells, each of dimension `d=66` or `78`.
For `n>=r`, `lift_program(P,n)` replaces only the carrier by the `n`-cell
independent product. Each repair in a local program receives the literal
address `cell=0`; existing product addresses remain unchanged. Parameters,
selectors, record-law tags, repetitions and memory policy are identical.
The actual admission function recomputes and validates the resource
declaration. The theorem applies when both programs are admitted.

Write an index as `i=a+d^r b`, where `a` labels the retained cells and `b`
the appended cells. The interpreter's literal index update at addressed
cell `j<r` is

```text
i -> i + (base[(i//d^j) mod d] - (i//d^j) mod d) d^j.
```

It changes only `a`. Its input-domain condition, coin bias and raw record
also depend only on the addressed cell and declared command. Thus each
complete raw branch is exactly

\[
 K_h^{(n)}=K_h^{(r)}\otimes I_B,                                \tag{1}
\]

with the tensor ordering understood through the displayed index map.
Holds give the identity; every impossible or reserved address gives zero.
The existing one-factor spectator test checks a special case. Equation
(1), together with the actual record controller, proves the full history
statement below.

## Complete history and arbitrary correlated inputs

Let `X` be any operator on retained cells, appended cells and an untouched
reference `R`. Its factor restriction is

\[
 [\operatorname{Tr}_B X]_{a,u;c,v}
   =\sum_b X_{a+d^r b,u;c+d^r b,v}.                            \tag{2}
\]

The same complete chronological record word gives the same memory,
exported archive and last-current value in both programs. Every admitted
selector therefore makes the same weighted choice on that word. Applying
(1) and (2) proves the branch identity by induction over commands:

\[
 (\operatorname{Tr}_B\otimes\mathrm{id}_{T,R})
       \mathcal I_{P_n}(X)
   =\mathcal I_P(\operatorname{Tr}_B X).                       \tag{3}
\]

Here `T` contains the complete output tape, including zero, hold, abort
and stop slots. No factorization of the input is assumed. Correlations or
entanglement with appended factors and `R` are covered. Off-diagonal
appended-factor matrix units have zero partial trace on both sides;
the corresponding zero comparison branch is retained. The proof is
complex linear, and the rational evaluator checks real and imaginary
parts separately.

The finite stack supplies the reusable tensor algebra in
[`NoSignalling.lean`](../../Lean/ObserverPatchHolography/Locality/NoSignalling.lean).
Its `liftLeft_conj_apply` expands a lifted Kraus conjugation into matrix
entries, as used in (1)–(3). Its `ptraceFst_local_kraus` traces the acted
first factor and proves preservation of the complementary spectator
marginal after an unread complete local Kraus channel. Equation (3)
instead traces appended factors and retains the active instrument; the
no-signalling theorem supplies a complementary corollary. The actual
native index map supplies the factor-local attachment, and the retained
history induction supplies the observer composition here. Those two
source-specific steps are prose-checked and executable-tested; they are
not claimed as a separate formal Lean instantiation.

Both executions have the same native attempt count and output-block
count. Intrinsic FIFO timeouts occur at the same step. Equal evaluator
attempt or expanded-branch limits also stop them together, since the
interpreter counts structural raw branches even when their operators have
zero weight. This includes record-adaptive selectors. A new controller
that explicitly reads the appended factor count would be a different
program and would not satisfy the hypotheses.

The known source program identifiers differ. The comparison first verifies
that the source identifier is exactly the admitted lift of `P`, then
maps that identifier to `P`'s identifier while preserving complete
chronology, FIFO/archive placement, stop statuses and runtime counters.
The declared allocation accounts and direct specification queries remain
distinguishable. Equation (3) identifies a restricted process, not the
entire larger universe state. After the stated comparison the channels
are exactly equal, so their diamond distance is zero.

Successive restrictions compose. For admitted factor counts `r<=s<=n`,
`Tr_(n-r)=Tr_(s-r) o Tr_(n-s)`, and applying (3) at either stage yields
the same retained complete process. This follows directly by grouping
the finite sums in (2). It supplies a dynamic counterpart to the distinct
static canonical restriction of the fair federation.

## Evaluator and physical boundaries

The full interpreter has an early carrier guard: products above 64 cells
return `RESOURCE_STOP_UNEVALUATED_CARRIER` before executing a command.
For a local program and its 65-cell lift, this differs from an otherwise
completed local execution. Equation (3) therefore describes executed
histories only when the programs enter the common evaluated carrier
domain, with equal limits. Its raw-operator and finite-history proof
extends to larger admitted product presentations, without asserting that
the bounded host evaluator runs them. The fixed request resource envelope
must also pass admission; this theorem does not remove it.
Wall-time exhaustion, process failure or memory exhaustion depending on the
represented input size are not synchronized by the proof. The executed
comparison concerns returned interpreter histories under the stated common
limits, not equal hardware resource use or guaranteed completion on a host.

The extension preserves the active weighted source law and its complete
native records. It does not make the unused cells a connected physical
space, produce interactions, select a global reference or identify a
cosmological capacity. Nor does an arbitrary correlated extension imply
that a global reference's modular flow restricts to the active reference's
modular flow. The operational repair/record comparison (3) is the proved
statement.

The connected fair federation is a different source branch. On seam
`(0,5)`, one actual weighted repair records coin one with probability
`25/33`; a fair federation repair records it with probability `1/2`,
independently of the incoming state. Thus a direct single-attempt
comparison preserving the coin record has complete-output total
variation at least `17/66`. The executable control retains all 25 raw
slots. This excludes that direct comparison, not an arbitrary encoded
multi-attempt simulation or a separately justified observation map.

## Verification

The receipt compares actual weighted H66/H78 executions on 2, 3 and 64
factors with correlated active/appended/reference inputs, retaining both
real and imaginary operator parts. It includes adaptive selectors,
complete FIFO exports, intrinsic timeouts and evaluator attempt/branch
stops, and records the 65-factor countercontrol. Tests also address more
than one retained active factor, exact nested restrictions, nontrivial
exploratory and uniform selectors, a traced-zero off-diagonal input,
literal identifier distinction and rejection of coherent service words.
The general proof is (1)–(3); finite cases are independent executable
checks of its actual interpreter attachment.
