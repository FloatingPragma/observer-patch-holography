# Native charge, spatial degree and the received pointer

This candidate addresses the literal static row: integer topological charge,
with a unit unable to split into two equal halves. It uses the carrier's
unchanged integer endpoint-repair process. Its source, conserved reference,
spatial dictionary and complete readout are specified in `CONTRACT.json`.
It does not identify an electron by importing a mass, spin or force result.
Those are separate rows. The internal full-row decision requires two
independent reviews of the joined evidence.

## 1. Source and exact spatial topology

Orient each native seam e=(a,b), a<b. Write B_e=delta_a-delta_b. The actual
source event transfers the integer

    j = [n_a-n_b -(2c-1)(|n_a-n_b| mod 2)]/2,
    n' = n - j B_e,        F' = F + j e.

G is the actual retained initial occupation and F initially vanishes. Thus
BF=G-n for every finite native history, including all holds. This is an
identity of the actual source and its emitted current records. The one-load
and two-load nonnegative sectors are closed because the new endpoint loads
are floor((n_a+n_b)/2)+c((n_a+n_b) mod 2) and its complement. `run.py`
compares this independent implementation against all 5,400 actual branches.

The carrier has an oriented icosahedral surface with primal face boundary A.
Its dual graph has 20 vertices and 30 edges, with one pentagonal dual cell
around each of the 12 ports. `topology.construct` reconstructs this surface
from the actual seam list, checks BA=0, and follows each oriented pentagon.

A retained current F gives a continuous circle path on every dual edge:

    z_e(t) = exp(2 pi i F_e t),       0 <= t <= 1.

Every dual vertex has the same value one. The oriented boundary of the dual
cell at v has degree (BF)_v. Degree-zero cells can be filled: lift their
closed boundary to a real phase with equal endpoint values, then contract
that lift across a disk. A nonzero-degree cell cannot be filled without a
puncture or a zero of the circle-valued field. Native charge is therefore
the signed spatial defect q=n-G=-BF, not a newly attached numeric label.

There is also an integral uniqueness statement. Two edge cochains differing
by A h are dual coboundaries and have the same spatial homotopy class.
The quotient map

    Z^30 / im_Z(A)  --B-->  {q in Z^12 : sum q = 0}

is an isomorphism. The executable verifies a spanning-tree incidence minor
has determinant +/-1. This proves that integer charges have integer current
representatives and that the kernel of B is freely parametrized by its 19
chord coordinates. It then verifies a 19 by 19 face-boundary minor in those
coordinates has determinant +/-1. Hence the face boundaries generate the
entire integer kernel, with no torsion or extra fractional sectors. Numerical
rank agreement would not prove this; the emitted determinants are exact.

For a dual-vertex integer function h, the common-vertex homotopy is

    H_s(e,t) = exp(2 pi i [s h_tail +
                         (F_e+s(h_head-h_tail)) t]).

Shared endpoints agree at every s. Every closed spatial degree is invariant.
This states the permitted geometric deformation, rather than declaring a
current count to be topological. Native changes move adjacent defects by
the actual signed current. A loop surrounding a port region has degree
sum_v(BF)_v there because all internal dual edges cancel. Moving the loop
through degree-zero cells preserves its value.

`phase_degree` separately evaluates the actual glued, deformed circle curves
on source histories. It uses a fixed 64 subdivisions and the declared
|F_e|<=3 diagnostic domain; its derivative bound is less than pi/2 per
sample. Large-turn aliasing outside that domain is rejected. The integral
topology theorem does not depend on sampling or that numerical cap.

## 2. Primitive charge and the no-half theorem

The primitive generator of pi_1(S1) is a complete turn. A primitive actual
unit current on a seam produces degrees +1 and -1 at its two adjacent dual
cells. The integral isomorphism fixes the faithful dictionary up to a common
orientation reversal. Multiplying all degrees by k with |k|>1 is a covering,
not a lattice isomorphism; it misses primitive classes. The executable
challenges a doubled winding encoder against the actual pointer, and a
sign-only encoder against measured degrees +/-2.

If a unit defect separates into two disjoint daughter regions, each daughter
boundary is a continuous closed circle map and therefore has an integer
degree. Degree is additive. Equal daughters would satisfy 2d=1 (or -1),
which has no integer solution. This excludes half-integer sharp charges;
an ensemble average of 1/2 is not a half-charge outcome. In the one-load
source, every native branch also preserves a whole load directly. A real
half-current mutation conserves the real continuity equation but gives an
antipodal endpoint on its dual edge, breaking shared-vertex gluing. It fails
the topology and record-alphabet gates. Conversely an actual native double
load splits into two whole units. The latter control prevents a vacuous
claim that the source cannot split anything at all.

The native integer alphabet and circle's primitive generator are explicit
architectural choices. No hypercharge weights, fitted e, or empirical
normalization are used to install the tested outcomes. The older Z6 voltage
with spectator-dependent updates is a different memory and is not identified
with F. Nor is this result transferred without proof to a coherent field law.

## 3. Closing the actual source-to-observer join

Finite counters retain (BF)_v modulo 6. At total L=1 or 2, every regional
degree lies in [-L,L]; reduction modulo 6 is injective on that interval.
All 6,228 preparation/configuration pairs obey the bounds, including the
overapproximation of configurations not dynamically reachable from a given G.
The proof covers every region by subset-sum bounds, not only the displayed
regions. Unrestricted charges are excluded: F=6e is an explicit alias of
zero modulo 6 with nonzero spatial degree.

Before readout, append the already frozen native AU5 thirty-seam E sweep.
Let n be the actual source and n_hat=G-B F_received the reconstruction at
its start. Their discrepancy delta remains constant during a faithfully
received sweep because the same currents update both. A matching contrast
on edge (a,b) says delta_a=delta_b. The swept graph is connected, so delta
is constant; the equal conserved total makes that constant zero. All closing
contrasts agree if and only if the observer reconstructed the actual state.
Consequently the reconstructed winding equals the actual spatial degree.

The sweep itself moves loads. The reported charge refers to latch L1 after
its thirty events, not to the earlier L0. Both latches and differing examples
are retained. Faithful receipt of the closing records is a stated premise
of this simulation experiment; matching contrasts do not authenticate a
channel that replaces the entire closing tape. Lost closed cycles with
B Delta F=0 cannot be detected by this theorem, but they are homotopically
trivial on this carrier by the integral isomorphism above. They do not alter
the regional charge, although they can alter an elapsed-attempt record.

## 4. Native pointer and complete laws

The auxiliary pointer is a second actual twelve-port carrier with one load
in an unknown initial column. Its fixed controller discovers a nonzero
native current, routes from that received record to port zero, then requests
6|w| moves round the actual cycle 0,1,8,6,2,10,0 in the sign of w. Each
requested hop has actual move/hold branches and a bounded retry count.
Every hold, discovery/routing/measurement timeout and moving record remains
in the emitted tape. A failure reports no charge. Successful paths close
and their signed moving edge count divided by six is their winding.

The added `raw_pointer_degree` reads only actual measurement currents and
spatial seams; it does not use the commanded degree or decoded target.
It agrees with the spatial source degree on every successful joined example.
The frozen observer-only facade ignores the supplied pointer final-load
snapshot, tries the unchanged AU5 schema on synthetic final columns, and
unions all compatible unknown initial columns. Source executor snapshots
are kept in the audit receipt but are absent from admission. An arbitrary
malicious final snapshot leaves the admitted observable unchanged.

For each initial pointer column x and declared retry R, let h=2^-R.
If y_k are the endpoints reached by its five incident discovery seams,

    P_reset(x)=sum_(k=0)^4 h^k (1-h)^(1+d(y_k,0)),
    P_complete(x,w)=P_reset(x) (1-h)^(6|w|).

Failure is the complementary retained probability. An unknown initial
classical distribution mixes these twelve complete laws without conditioning
on success. The producer's exact joint-law receipt contains 768 complete
instruments, each covering all 2^32 source coin words by integer recursion,
then every pointer branch. Its independent opposing audit recomputes the
complete maps. This row consumes those frozen laws, verifies their full
normalization and binds their source/kernel/pointer to the same dictionary;
it does not relabel this consumption as a new full-law run. Additional
successful and failed actual tapes are replayed through the native executor
and observer-only facade in the row's `EXAMPLES.json`.

## 5. Exact precedent, scope and verdict

The pinned upstream `MODELS.md` gives this row a static regime. Its M5
validated artifact measures seeded winding for vacuum, one positive defect,
one negative defect and the separate contours of a pair, with tolerance0.10.
The actual script disables relaxation, skips the combined pair contour and
is archived against the retired engine. Its cell does not require an
electron-mass or spin measurement. This construction uses different native
fields and exact integer topology; agreement with the precedent's gate does
not claim the models share a field or confer an official platform award.

The declared defect is relative to retained G. On a closed carrier its signed
charges sum to zero; a transmitted native unit has a target defect and a
reference deficit. The full pair contour consequently has zero degree.
This is the candidate's literal charge observable. It does not establish
that all quantities conventionally called electric charge elsewhere in the
repository coincide with this one. That requires a source transfer.

The constructive result offered for independent full-row review is the
single chain native event -> retained current -> spatial homotopy class ->
closing contrast -> finite counter -> native pointer, with integral signed
outcomes and equal-half exclusion in the stated source and intervention
domain. The observer/source distinction and construction choices are part
of this claim. The internal count remains unchanged until both independent
review records accept this exact full criterion and evidence.
