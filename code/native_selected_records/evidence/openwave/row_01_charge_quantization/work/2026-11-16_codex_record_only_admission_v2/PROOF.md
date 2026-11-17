# Unknown-count admission with a retained reservoir

Let n=12B, B≥2 fixed before hidden first-cell preparation in[-B,B]^12,
with the remaining core initially zero. Then |Q|≤n and V_initial≤12B².
The graph is connected. Port12 is a supplied initially zero reference.

On a seam with initial difference d, the first repair current and the first
opposite-coin current obey d=2j_first+j_opposite. Intermediate identical coins
have zero current. Subtract the previously observed probe disturbance, integrate
the tree contrasts and fix the constant with initial port12. The reconstructed
current state includes every disturbance; no total or seed is needed. With
cap64 a tree edge lacks an opposite coin with probability2^-63.

Native repair preserves Q and the initial min/max box, while
V_before−V_after=(d²−(|d| mod2)²)/2. For a core span greater than one, choose
a maximum M and a shortest path to a value at most M−2. Every interior is
M or M−1. An M already supplies the next maximum; an M−1 can exchange with
its adjacent M with probability1/2 per native attempt. The final contrast
has magnitude at least two, so either coin decreases V by at least two.
There are at most6B² strict decreases, each with at most n−1 requested fair
exchanges. This proves finite adaptive selection apart from explicit timeouts.
Swap/hold preserves V and fails that strict-decrease request. It may still
inherit unit selection already performed by the shared native census; those
cases receive no control-discrimination credit.

The only adjacent-level classes with |Q|≤n are absolute {0,1}, {-1,0} or
all zero. If nonzero, start a breadth-first path at the empty target leaf and
stop at the first occupied core port. Every interior is vacant by shortestness.
Reversing the path moves one real same-sign unit through vacuum to the leaf.
The last core vertex is center0, so it ends empty; the other target leaf was
never addressed. The core retains Q−sign(Q), and all its loads remain0 or
sign(Q). At most n requested fair moves suffice. In particular Q=n is allowed:
the center's unit moves directly to the leaf, leaving n−1 reservoir units.
Zero terminates as no-unit; it is not a failed one-particle preparation.

The complete inconclusive probability is bounded by

    eta_B ≤ [2(n−1)+6B²(n−1)+n] / 2^64.

Each capped failure has its probability conditional on the entire preceding
record, so the union bound does not require independent movement counts.
This is an all-input theorem on the declared preparation domain. The finite
grid tests a disclosed subset, not all(B+1)^12 signed initial vectors.

The controller chooses one retained target unit, but neither an initial total
one nor a final global one-unit state is required. The remaining units are
present in the same source. Only target edges are addressed during the quantum
prefix; for each such edge, the literal native function depends solely on its
two target coordinates. Consequently every fine current branch equals the
one-target branch tensor identity on the fixed reservoir and marker state.
The source compiler enumerates these global columns for every successful
history. With the declared coherent memory, the extension is linear and also
holds for untouched references within that admitted target sector. It is not
claimed for arbitrary target-reservoir exchange or an interacting background.

For negative target charge, native(−q,c)=−native(q,1−c) on the target. The
actual negative arc permutation is the positive permutation conjugated by
flipping each active coin coordinate. Fair source choice is invariant under
that permutation. The complete source mean and cooling probabilities therefore
agree, while current signs and the actual negative source load remain retained.
There is no physical sign reversal or unrecorded replacement of the source.

Every extracted state starts the new prefix at the observed target leaf,
with the full recorded background. The native contrast/filter/counter emits
all complements. Prefix probability is conditional on admitted nonzero output;
the complete classical failure bound and zero outcomes must remain separate.
This proof and execution remove a total-one entry promise for this descendant.
They do not transfer the ancestor's complete224 quantum result, derive the
coherent services, or establish an electric/topological or empirical mass claim.
