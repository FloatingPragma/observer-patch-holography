# Arbitrary-load source coupling and a selected-unit response

All source values here are integers. The finite graph is connected. Every
body attempt independently chooses an edge uniformly and a fair placement
coin. Both copies use that same choice. All claims below concern this
declared joint process, not two copies with independent clocks or coins.

## 1. Exact local closure, including arbitrary backgrounds

On an addressed pair let S_X=X_a+X_b, S_Y=Y_a+Y_b and s=2c−1. The unchanged
native function has endpoint outputs

    X'_a = (S_X+s r_X)/2,  X'_b = (S_X−s r_X)/2,
    r_X = S_X mod 2 ∈ {0,1},

and the same formulas for Y. These equal the existing current formula even
for negative inputs, because |X_a−X_b| and S_X have the same parity. Put
q=X−Y, S_q=q_a+q_b=S_X−S_Y. Subtracting the actual output equations gives

    q'_a = (S_q+s(r_X−r_Y))/2,
    q'_b = (S_q−s(r_X−r_Y))/2,
    k = j_X−j_Y = q_a−q'_a.                         (1)

When S_q is even, r_X=r_Y, and both discrepancy endpoints equal S_q/2.
When S_q is odd, r_X−r_Y is either +1 or −1. Multiplying a conditionally
fair s by that past-measurable sign leaves it conditionally fair. Thus,
given the entire past and current addressed edge, the two equally weighted
outcomes (q',k) are exactly the two native outcomes on q. They need not
match under the same numerical coin label: odd background parity can swap
those labels. Equality of unordered two-coin branches is the claim.

Because each edge is independently uniform, conditioning only on q also
gives the original native transition kernel. This establishes a Markov
closure for q even though the label permutation depends on X,Y. It is
valid for every integer endpoint quadruple, with no density, sign or
vacuum-level restriction. It covers the joint emitted difference current,
not merely the population transition matrix. The new run exhausts all
6,561 quadruples in −4…4 and checks the complete branch multisets.

## 2. Synchronization and its finite expectation

The existing balancing theorem applies to q by (1). For completeness, put
V(q)=Σ_i q_i². Each native edge repair decreases V by

    ((q_a−q_b)² − ((q_a−q_b) mod 2)²)/2 ≥ 0.

If the global maximum exceeds the minimum by at least two, choose a
shortest path from a maximum M to a site at most M−2, stopping at the
first such site. The interior has values M or M−1. Permitted native
exchanges move an M along that path until a contrast of at least two is
repaired. This strictly decreases V by at least two in at most n−1
specified edge/coin attempts. If the total of q is zero, its only
adjacent-level vector is zero. It is absorbing, since equal sources remain
equal under every common choice.

Let L=n−1, p=(2m)^−L and K=V(q_0)/2. K is an integer because Σq_i=0 and
q_i²≡q_i mod2. From every nonzero q the probability of a strict decrease
in the next L attempts is at least p. There are at most K such decreases.
The waiting time in L-blocks for each decrease is dominated by a geometric
variable of mean1/p, conditional on its past. Therefore

    P(tau<∞)=1,            E[tau] ≤ L K/p.          (2)

This is a deliberately loose universal bound. Finite sampled histories
are not the proof. For the 48 selected-unit native histories, the observed
first synchronizations range from1 to285 attempts; the full4096-attempt
records remain. All48 matched swap/hold responses remain distinct.

## 3. Source response from currents, without a field evolution

Orient B with +1 at an edge's first endpoint and −1 at its second. The
two native current tapes give k=j_X−j_Y. With F_0=0, exact continuity is

    q_(t+1)=q_t−B_e k_t,     F_(t+1)=F_t+k_t e.

At tau from (2), B F_tau=q_0 on every path. The nontrivial dynamical
statement is the finite synchronization; continuity is identified as an
identity, not independent quantization evidence.

By (1), E[k|q,e]=(q_a−q_b)/2. Hence with L_G=BBᵀ,

    E[q_(t+1)|q_t]=(I−L_G/(2m))q_t = T q_t,
    E[F_tau]=Bᵀ L_G^+ q_0.                       (3)

To justify the sum, stop the discrepancy currents at tau (all later ones
are zero). The discrepancy's endpoint extrema never leave their initial
interval by (1), so each |k| is bounded by that finite initial range.
Equation (2) then gives E[Σ_(t<tau)|k_t|]<∞. We may interchange expectation
and sum. On the neutral subspace, T has spectral radius below one:
0≤T≤I because the graph Laplacian eigenvalues lie in[0,2m], and its only
eigenvalue1 is the connected graph's constant mode. Thus

    (1/(2m)) Σ_(t≥0) T^t q_0 = L_G^+ q_0.

The potential L_G^+q is a derived representation of the mean recorded currents;
it is never an evolving source variable. The executable independently
checks the exact rational stopped-reward equation on all132 primitive
sources and7,920 branches. Linearity in (3) is valid for arbitrary neutral
integer q, not just those primitive sources. Uniform shifts of either
background change no current; the same-total comparison fixes the
otherwise missing difference offset.

## 4. Same-source selection and intervention

The existing all-integer balancing theorem selects values b,b+1 for
totalQ=nb+k,0≤k<n. A nonconstant settled configuration has a mixed edge:
otherwise connectivity would make it uniform. A native repair on that
edge holds or exchanges the two values, each with probability1/2. The
nonzero emitted current therefore has magnitude one, and the relative
source produced by performing that exchange on X while leaving Y is
q=±(e_a−e_b). Repeated fair attempts succeed almost surely; capped failures
remain explicit. No edge or particle position is supplied to that
intervention policy: it chooses the first mixed edge of a native
current-only reading.

The experiment first prepares equal sources, then executes selection,
reads the selected configuration, accounts for changes caused by that
reading, performs this native intervention, and executes the response.
These are stages of one declared source protocol. All48 native preparations
settle; none of48 selection swap/hold controls settle. An actual double
source and a binary source with the same total are distinguished by the
same reader before selection. The two law-sensitive facts—selection and
response synchronization—are tested separately as well as in sequence.

At higher density these are unit excitations above b, not a proof that
every absolute occupancy is0 or1. A uniform background has no relative
excitation to choose and would produce an explicit uniform outcome.
The experiment includes b=2 and high-density lower-level holes, so a
zero vacuum is neither assumed nor inferred. Native integer states and
currents remain supplied. The sharp unit and its primitive integral
relative response class have no equal sharp half decomposition; this
last arithmetic follows from their integral carrier and is not credited
as a separate emergence mechanism or independent topology theorem.

## 5. Operational observer and controls

For one native edge, the first current j0 and the first later opposite-coin
current j1 satisfy initial contrast d=2j0+j1. Intervening repetitions of
the first coin have zero current. Each tree-edge probe corrects this
contrast by the incidence of prior probe currents. Integrating the tree
recovers the pre-probe source modulo one common constant. All statements
refer to actual native probe records; the controller never sees a hidden
load or a numeric Q. The observer then transports that reading backwards
through the preceding current tape. The two relative readings determine
q by subtracting the mean of their difference, using the protocol's
equal-total promise. This is not an absolute background measurement.

For the selection-stage reading, the probe's own accumulated current also
gives its post-probe endpoint; that is the endpoint used to choose the
intervention. Since native probes preserve the settled class, selection
does not have to be undone. Terminal response probes can alter or even
synchronize the sources, but backwards transport reconstructs the body
before those probes. A probe timeout yields no body verdict. The forced
cap1 timeout is retained as an explicit detector check.

Common swap/hold applies one identical permutation to both source vectors.
It conserves each multiset and the discrepancy multiset. Consequently it
cannot remove an initial double during selection or synchronize distinct
sources during response. It preserves the graph's symmetry, the integer
alphabet and the observer apparatus. This is a relevant altered law.
The current-corruption mutation is separate integrity evidence.

## 6. What this theorem adds to the row

This closes a specific missing join: the dynamically selected unit is now
the source of its own measured finite Gauss response, on the same expanded
integer candidate and with arbitrary common backgrounds. It removes the
old restriction of the paired-source theorem to binary common occupancies
and needs no numerical initial total or absolute zero in its reader.

It does not derive an electromagnetic force or continuum field, replace
the supplied integer alphabet, or establish a separate configuration
winding number. The old relative homology is still isomorphic to load
bookkeeping; it gains no additional credit here. The constructive
multi-load route offered in PR723 is the relevant route for adjudication.
No coherent memory is used. The separated-pair energy anomaly belongs to
the distinct old logarithmic energy candidate and is not resolved or
renamed by (1)–(3). This packet asserts no energy, spin, G6, or row award.
