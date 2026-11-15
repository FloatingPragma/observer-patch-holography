# Native charge as a spatial homotopy class

The existing spatial adapter is `../spatial/native_dual_winding.py`. Its source
is the literal integer endpoint repair on the actual twelve-port icosahedron.
Its state is the native load vector and the committed signed current tape.
The primary preparation is one known unit at a carrier port and zero current
tape: `G=n₀` and `F₀=0`. Every subsequent event updates `n'=n−jB_e` and
`F'=F+j e`. Thus `BF=G−n` on every complete history. The source reference G
is the actual initial occupation, not an arbitrary offset or reset. All later
records retain their actual coins, seams and signed currents, including holds.

The [independent executable](audit_spatial_lift.py) verifies the producer's
132 initial/arrival choices, every coin history of their prescribed schedules,
and all 7,920 subsequent seam/coin branches. It uses a separate
floor/ceiling implementation of the endpoint rule and compares it to the actual
native function. No acceptance configuration or continuum approximation is used.

## The integral spatial connection

Write `A:C₂→C₁` for the oriented triangular-face boundary and `B:C₁→C₀`
for the primal current divergence. The actual carrier has 12 vertices,
30 edges and 20 triangular faces. The independent construction orients its
faces outward using exact arithmetic in Q(φ), rather than the producer's
orientation propagation. It checks `BA=0` and that every edge belongs to two
oppositely oriented faces. This is the oriented icosahedral sphere.

The dual graph has 20 vertices, 30 edges and a pentagonal dual face around
each native port. An integer current chain F is a dual-edge integer cochain.
An oriented dual-edge representative is the closed circle path

`z_e(s)=exp(2πi F_e s), 0≤s≤1`.

Every dual vertex has value one. The winding around the dual face of v is
exactly `(BF)_v`. A zero-winding face admits a continuous disk extension:
lift its boundary phase to a real function with equal endpoint values and
contract that real lift radially. A nonzero-winding face cannot be filled
continuously; it is a spatial vortex puncture. No independently prepared
texture field is appended to the carrier. The circle field is a geometric
realization of its actual retained current cochain.

More strongly, its homotopy class is canonical. For the dual graph,

`H¹(G_dual;Z) = Z³⁰ / image(A)`.

The map induced by B is an isomorphism onto the neutral charge lattice:

`[F] ↦ BF ∈ {q∈Z¹² : Σ_v q_v=0}`.

This statement needs an integral proof, not only matching matrix ranks.
The executable chooses a primal spanning tree. Its reduced incidence minor
has determinant one, proving surjectivity onto the neutral integer lattice.
Every integer cycle is uniquely parametrized by its 19 chord coefficients.
The face-boundary matrix restricted to those chords has Smith factors
`1,1,…,1` (nineteen entries), proving `ker B = image(A)` over the integers.
Consequently two current chains with the same native charge differ only by a
dual coboundary and define homotopic circle maps. Detailed phase interpolation,
choice of initial connecting path and contractible native loop detours do not
alter the spatial homotopy class.

The source relation is therefore a genuine equality between the actual
regional current charge and a spatial homotopy invariant, with the sign stated:
the winding is `G−n`, while excess occupation is its negative. It is not a claim
that native integrality was derived from unconstrained real-valued fields.
The integer source alphabet and the primitive circle generator are explicit.

## Native updates and non-halving

An actual current j changes just the two incident dual-face windings by
`+j` and `−j`. All other windings are unchanged. A loop surrounding a region
has winding equal to the sum of its internal face windings, because internal
dual edges cancel. Moving that boundary without passing a nonzero defect
therefore preserves its degree. Native transfer moves the endpoints of the
spatial defect; opposite endpoints can annihilate. No lifetime or fixed-core
assumption is needed for this static charge criterion.

Every daughter boundary is a continuous circle map, hence has integer degree.
Two equal daughter degrees cannot sum to one. This is not a sign-only meter:
the producer's double defect has winding two and an actual native event splits
it into two unit defects. Replacing one current by one half preserves a real
continuity equation but destroys the common-vertex gluing of the proposed
circle map. Its endpoint is the antipode, not the initial phase. Fractional
expectation values remain possible; they are not half-integer sharp degrees.

## A literal-native winding pointer

The executed load domain is finite and closed: for each known initial reference
G, one unit occupies one of twelve ports. Every native branch moves that unit
or holds it. Across reference choices the disturbance classes comprise 132
ordered dipoles plus zero, but G remains retained: zero disturbance alone does
not identify the initial port. Each regional winding is in `{-1,0,+1}`. Its residue modulo six identifies it exactly through
`0→0, 1→+1, 5→−1`. With the known preparation, a six-state observer updates
`r←r+(1_R(a)−1_R(b))j mod6` from the actual current record. The remaining
three residues are explicit domain-failure outputs. All 540,672 seed/region
decodings are independently checked, and the native branch check proves the
one-load domain remains closed at arbitrary finite time.

Thus the full signed tape is needed only to retain a chosen field
representative. The spatial class itself can be kept in 12 finite residue
counters, or one of the 133 charge configurations across G, together with the
known G. A canonical
connecting path between its two nonzero ports recovers a representative;
the integral isomorphism proves that this compression preserves the complete
spatial homotopy class. It does not replace an unknown-state preparation by
knowledge: the initial unit, G and zero counter are declared together.

No analogue phase gate is required to read this class. Realize the target
circle as a six-vertex cycle. On a dual edge with lifted current F_e, emit
`6|F_e|` signed unit transitions of that cycle. The resulting closed path has
winding F_e. Reading a dual boundary concatenates these paths in its actual
orientation. Count the signed wraps in the retained transition record.

This pointer can itself run the literal native law. The actual carrier has
the six-cycle `0→1→8→6→2→10→0`. A separate pointer carrier contains one known
unit, initially at port zero. For each requested oriented pointer hop, schedule
the next cycle seam. Its two actual coin outcomes are a move and a hold, each
with probability one half. Retry up to 16 times, retaining every hold, coin,
current and a distinct timeout outcome. No unknown pointer state is read: the
known preparation and complete native records determine it.

For K requested pointer moves, successful completion has probability
`(1−2⁻¹⁶)^K`. The first-timeout probabilities at every stage sum with it to one.
The executable checks every position, both orientations and both coins of this
actual native pointer, and retains its entire one-attempt transition table.
The maximum tested regional readout requires K moves as stated in the receipt;
its finite maximum attempts and exact failure probability are retained. A
timeout is not a charge-zero outcome. A successful readout is not silently
renormalized into an unconditional result.

The prepared auxiliary unit, retained history memory, controller and pointer
schedule are explicit observer apparatus. The pointer does not generate the
source's topological class: the integral isomorphism above already identifies
that class before the readout. It operationally displays the spatial winding
using the same endpoint repair rule.

The compressed class may instead display its regional value with zero or six
signed pointer hops. This is equivalent to the full representative readout
by the proved integral isomorphism, and bounds its pointer work independently
of the source history length. The companion native-pointer implementation
can supply reset of an unknown one-load pointer before those commands; its
initialization and timeout terms must accompany any use of that stronger
service.

## The six-class and physical boundaries

On unrestricted source domains, reading only `F mod6` loses information
needed by this observable. A chain
with six units on one edge has all-zero residues and windings `+6,−6`, whereas
the zero chain has no defect. This control lies outside the restricted
one-load/reference domain and does not refute its injective finite counter.
For unrestricted charges the signed current tape supplies the lift. The
existing spectator-sensitive active-memory voltage
`3 sorting_parity + 2j 1{selected seam and spectator}` is a different memory
construction and is not silently identified with this current lift.

One may regard a lifted tape as a piecewise oriented path through the sixth
roots of unity and take its sixfold closed path. That explains the six-state
pointer convention; endpoint phases alone do not specify the lift. No new
independently evolving electromagnetic field, force, energy or electric SI
normalization is claimed. The derived spatial degree and complete native
readout establish the mathematical carrier/topology connection in this declared
observation domain. Whether it earns the platform's charge-quantization cell
is a maintainer judgement on the source and observable, not an award made by
this audit.
