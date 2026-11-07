# Full-number native configuration readout

The source specification supplies its conserved total occupation. Reading
that integer, the repair law, scheduler and memory implementation is allowed.
The construction uses the existing binary endpoint E instrument and complete
raw records. It introduces no coherent register or gate service.

## Connected-tree readout

Let a connected graph have V binary ports, strictly positive endpoint coin
probabilities, and conserved total m in {0,...,V}. Choose a spanning tree,
ordered so that each edge adds one vertex to the previously visited set.
Perform one fixed-seam E repair on each tree edge, with repair probability
one. The raw instrument retains its hold slot, all C slots of zero law weight,
and all E slots, including impossible labels: 25 addresses per attempt.

For an ordered seam (a,b), the actual record gives the incoming contrast
d=n_a-n_b, the placement coin c and current j. The update is
n_a'=n_a-j, n_b'=n_b+j. Track the accumulated displacement s, initially zero.
Before updating s, the original difference is exactly

```math
\Delta_{ab}=d-s_a+s_b,\qquad s_a\leftarrow s_a-j,\quad s_b\leftarrow s_b+j.
```

Propagate these differences along the tree, assigning offset o_0=0. Every
initial configuration compatible with the word has n_i=b+o_i. Its known
total fixes b=(m-sum_i o_i)/V. The decoder checks binary initial and terminal
values and replays every endpoint transition, rejecting inconsistent words.
The sparse computation uses O(V) integers rather than binomial(V,m) columns.
It applies to the empty and full sectors as well as every intermediate m.

Each raw E map is a positive weighted partial injection in the binary basis.
Composition preserves that property. Record decoding proves that every
nonzero complete word on the known sector has the form

```math
K_h=\sqrt{w_h}|y_h\rangle\langle x_h|,\qquad
w_h=\prod_{t=1}^{V-1}p_{e_t}(c_t).
```

Every input column has exactly one live raw outcome for each coin at each
step, and the two coin probabilities sum to one. Consequently the complete
word instrument is trace preserving. Coarse-graining only by decoded initial
configuration gives precisely the effects |x><x|. There are 25^(V-1) formal
words and binomial(V,m) 2^(V-1) nonzero words. Zero outcomes are retained.
For an arbitrary source/reference operator, the exact fine output is

```math
(K_h\otimes I_R)\rho(K_h^*\otimes I_R)
=w_h|y_h\rangle\langle y_h|\otimes\rho^R_{x_hx_h}.
```

This identity proves the reference extension without assuming a classical
input or discarding back-action. Each nonzero outcome prepares a pure basis
state whose terminal label is a record function. It is an outcome-labelled
preparation, not preparation of a predetermined target or a Lüders update.
The retained archive may remain correlated with the reference.
The executable sparse operator API represents exact rational real components;
complex operators use its two real/imaginary calls and linear extension,
matching the original interpreter's domain.

## A local number mixture and its precise ambiguity

If the total is unknown, the same offset calculation permits the two root
values b=0 or b=1. Both can give binary configurations only if every offset
is zero. In that case the only inputs are the all-empty and all-full states;
their currents vanish and they remain unchanged. Every other complete word
supports at most one binary input, with its total read from the decoded mask.
The coarse effects on the full binary direct sum are therefore all singleton
projectors except that the two constant configurations share one projector.

The constant branch is proportional to
P_const=|empty><empty|+|full><full|. It retains their coherence and arbitrary
reference correlations. Replacing it by two classical labels would change
the instrument. More generally, every repair or hold has identical record
statistics on these two states: zero contrast and zero current, with the
same coin law. Record-adaptive selectors cannot distinguish them either,
by induction on the identical record prefixes. No repair-only record menu
measures their relative phase or distinguishes their occupation without a
total specification or additional access.

This direct-sum statement is a local binary comparison, not admission of an
extra global coherent carrier. The actual federation fixes its global total.
A local restriction retains all local number sectors. A known global total
does not generally determine a smaller window's total. The shared empty/full
local outcome is retained unless its known number constraints exclude one.

## Final readback of an observed repair history

Append the E tree after a body of native repairs. Tree decoding gives the
body's terminal configuration without supplying its initial configuration.
For each body record in reverse order, undo the current:

```math
n_a^{\rm before}=n_a^{\rm after}+j,\qquad
n_b^{\rm before}=n_b^{\rm after}-j.
```

The helper `reverse_currents` checks the raw C/E decoding, placement coin,
binary endpoints and complete native update. It reconstructs every mask in
chronological order from the final mask and the complete body records. Thus
occupation products or displacement functions of these masks become record
functions on this one experiment. No initial census is an input. Reading a
source specification does not supply its hidden initial configuration.

This helper authenticates endpoint semantics. The caller must additionally
authenticate the body schedule, original selector probabilities, program
digest and FIFO/archive. The final tree measures and changes the body-terminal
state. It does not commute with arbitrary quantum observables or provide an
undisturbed measurement of an unrelated process. The full body-plus-detector
word must satisfy the original source resource bounds.

## Actual source attachment and execution domains

The local weighted H66 member supplies the V=12,m=2 control. For the admitted
fair connected federation at evaluated levels 0,1,2, V=240,960,3840. Their
full spanning trees use 239,959,3839 literal commands, within the original
4096-node admission cap and attempt limit. The total is the admitted carrier's
integer m. Fixed intercell addresses beyond the evaluated levels are not
admitted by that interpreter; the general graph theorem does not change this.
Appending a body leaves at most 3857,3137,257 additional literal nodes at
these three levels, respectively, if the body is not compressed by repeats.

The builder uses one FIFO slot and `export_oldest`: every complete word has
V-2 exported records and one retained record. Storage counts are read from
the admitted source program. Intrinsic FIFO halts are distinct and rejected
by this full-tree adapter. The generic expanded evaluator retains its
65536-branch and 64-step limits. Each attached full tree reaches its branch
ceiling after three attempts, including zero blocks; it returns RESOURCE_STOP.
That output is not a completed detector.

`continue_raw` is an exact lazy query of one unnormalized full word through
the original interpreter continuation kernel. It checks equality to the
record-decoded partial map and retains the actual FIFO/archive. It neither
enumerates all words nor implements a normalized outcome sampler.
`authenticate_trace` binds the exact program digest, slots, records, input
and terminal masks, weight and storage partition to that independent query.
Testing sparse columns on each supplied sector does not enumerate the global
carrier; the matrix-unit identity supplies the all-input proof.

## A bounded repair-only reset obstruction

For any binary input x, every selected repair seam has a placement coin that
leaves x unchanged. If its endpoints are equal both coins do so. Let p_min>0
be the minimum coin probability on the finite source graph, and let r be its
repair probability. The conditional no-motion probability at every attempted
repair is at least alpha=1-r+r p_min, independently of selector choice, C/E
law or previously observed records. Holds and termination also preserve x.
Induction gives a retained no-motion event of probability at least alpha^T
for any controller with at most T repair attempts. This includes the actual
last-current and exploratory selectors; no readout of a hidden state is used.

On a nontrivial fixed-number sector, choose x different from any proposed
fixed target z. The final projector onto x has probability at least alpha^T
for the actual complete process and zero for deterministic preparation of z.
The output half trace distance from that target preparation is therefore at
least alpha^T, even after tracing the archive. Exact finite repair-only reset
of all inputs is impossible. For the actual fair federation p_min=1/2;
for the weighted H66 member p_min=8/33.

The same argument applies after a specified tree outcome. If its terminal
configuration y differs from its decoded initial x, a bounded repair-only
continuation leaves y unchanged with positive probability. Exact restoration
to x on that outcome is therefore impossible, even if the observer knows
both labels. This statement does not classify approximate restoration.

This lower bound decreases with T. It does not rule out approximate reset,
heralded preparation, an unbounded stopping protocol or declared coherent
controls. The restricted native grammar has no command taking an arbitrary
decoded configuration and deterministically relocating its particles. Such
a controller or a new preparation theorem is not silently included here.
The readout supplies no spatial field, physical clock or particle mass.
