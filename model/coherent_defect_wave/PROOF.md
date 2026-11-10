# Reversible coupling of a regional defect and a carrier wave

## Source, sector and clock

This construction joins the [finite boundary-flux register](../coherent_seam_repair/PROTECTION.md)
to the [coherent seam walk](../coherent_seam_repair/PROOF.md). It uses the same
twelve-port carrier, five-level internal code, parity register and sixty
oriented seam addresses. The architecture fixes the five-face region around
port 5 and the exterior and boundary link bits. The initial negative link
is `(0,1)`. The region contains one central face flux and the exterior
contains its compensating flux.

Let `G` be the ten link words reachable by the permitted one-defect transport
rule. These are ten states of a thirty-bit link memory, not ten particle
species. The joint Hilbert space is

`H = C^G tensor (C^12 tensor C^5 tensor C^2 tensor C^60)`, `dim H=72000`.

The ten-dimensional link sector is invariant under each transport event.
The complex-state and complete-instrument framework is declared. A
thirty-bit register, link-dependent seam repair, wave-controlled link moves,
address reflection and the specified preparation and readout services are
construction choices. Their uniqueness is not required. Their consequences
are derived on this sector; arbitrary link configurations, all occupation
sectors and physical particle assignments are outside it.

The sorted interior edges are `(0,5),(1,5),(4,5),(5,9),(5,11)`. One clock tick
is the same fixed product of gates on every invocation. No external
time-varying link background or feedback selection is substituted during
the tick.

## One autonomous finite interaction

For each link word `g`, assign `G_e=Z^(g_e)` on every carrier edge, where
`Z=diag(I3,-I2)`. The existing reversible parity repair and address reflection
give the 7200-dimensional unitary `W_g=Q_g R`. Define

`B = sum_g |g><g| tensor W_g`.

This evolves the wave coherently conditional on the retained field register.
It does not measure or discard that register.

For an interior edge `e`, `T_e` flips its link bit if exactly one adjacent
region face has flux and otherwise leaves the word unchanged. The two face
bits swap on a move. Applying the same event again reverses the move, so
`T_e²=I` and `T_e*=T_e` on the ten-word sector. The no-operation case is
part of the permutation, not an omitted branch.

Let `Pi_e` project the wave onto either endpoint of `e`, acting identically
on code, parity and address. Define

`C_e = I tensor (I-Pi_e) + T_e tensor Pi_e`.

The two wave projections are orthogonal, so `C_e*=C_e` and `C_e²=I`.
This is a coherent conditional move: a superposition of endpoint and
nonendpoint wave states can become correlated with different link words.
Idle gate records identify the requested operation and tick; they do not
measure the link word or commit a moved/held branch flag. Such a measurement
would dephase this superposition and define a different channel. The link
register retains the branch information coherently. Optional detector
histories retain every actual outcome.
The fixed joint tick is

`U = C_4 C_3 C_2 C_1 C_0 B`,

with the rightmost gate applied first. Its inverse is

`U* = B* C_0 C_1 C_2 C_3 C_4`.

Thus a repeated tick is unitary on the same finite source, with a completely
specified inverse. The two factors do different work: `B` makes wave
propagation depend on link values; `C_e` makes allowed defect motion depend
on the wave. The reciprocal dependence is tested by complete observables
below. A controlled static background alone would omit the second factor.

The endpoint controls are local in the carrier incidence. The uniform
coherent reflection of all sixty addresses is an additional supplied
operation; this finite construction makes no geometrically local
continuum or finite propagation-speed claim for that service.

## Conserved sector and charge

Every permitted interior move leaves all boundary and exterior link bits
unchanged and preserves exactly one region defect. Equivalently, the mod-two
identity

`sum_(faces in region) b_f = sum_(edges in boundary) g_e`

fixes the odd regional flux. `B` leaves the link register unchanged, so the
entire tick preserves this sector, including coherent superpositions.
Allowing the boundary edge `(0,1)` to flip can remove the initial pair; the
protection theorem is relative to the stated operation class.

Every central link commutes with `Q5=diag(2,2,2,-3,-3)`, hence every `W_g`
commutes with the corresponding wave charge. The conditional moves are
identity on the internal code and also commute with it. Therefore

`[U,I_link tensor Q5]=0`.

This is an internal conserved charge. It does not identify electric charge,
physical spin, statistics or a Standard Model multiplet.

## Reciprocal complete readouts

Write `A|p,j>` for a wave at source port `p`, internal basis state `j`, blank
parity and the uniform sixty-address state. The binary detectors described
below retain both their yes and no outcomes. Face-position projectors depend
only on face fluxes, rather than a gauge-dependent raw link label.

### Link state changes the wave record

Use the same wave preparation

`psi = A(|0,3>+|5,3>)/sqrt(2)`

with either initial link word `g0=1` or `g1=3`. The first has only the
boundary edge `(0,1)` negative; the second also flips interior edge `(0,5)`.
Their region defects occupy faces `(0,1,5)` and `(0,5,11)`, respectively.

Read the wave basis projector at port 0, internal 3, parity 0 and oriented
address `(0,5)` after one complete tick. Because `psi` is in the image of A,
the first reflection fixes it. On the indicated seam the parity-zero output
amplitude is proportional to the sum of the two endpoint amplitudes for
`g0`, and their difference for `g1`. The resulting yes probabilities are

`p_yes(g0)=1/120`, `p_yes(g1)=0`.

Every `C_e` commutes with this wave-basis projector. These probabilities
therefore hold after the full interacting tick, not only after a static
wave prefix. The respective complements are `119/120` and `1`. Replacing
all controlled `W_g` by the same `W_g0` removes this particular contrast.
That control does not assert that the conditional moves lose all possible
effects on wave coherence.

### Wave preparation changes the defect record

Fix `g0` and compare the two wave preparations

`W_g0* A|5,3>` and `W_g0* A|6,3>`.

The inverse is part of the declared reversible source. The first part B of
the next tick returns these preparations to ports 5 and 6. Port 5 activates
all five conditional moves; port 6 activates none. Their final region
defects therefore occupy `(4,5,11)` and `(0,1,5)`, respectively. The complete
detector for face `(4,5,11)` has yes probabilities one and zero. The raw
final words are `4194307` and `1`, but their different face positions, rather
than those presentation labels, are the observable distinction.

Replacing all conditional moves by identity removes this contrast: the
field remains `g0` for both preparations. The inverse-wave preparation is
explicit and is not described as a natural vacuum or cosmological input.

### Coherent register correlation

The product preparation

`|g0> tensor W_g0* A(|5,3>+|6,3>)/sqrt(2)`

evolves in one tick to

`(|4194307> tensor A|5,3> + |1> tensor A|6,3>)/sqrt(2)`.

The two field vectors and the two wave vectors are orthogonal. The reduced
link state has eigenvalues `1/2,1/2` and purity `1/2`, proving entanglement
for the declared register factorization. This tensor-factor witness is not
by itself a gauge-invariant notion of physical entanglement. The reciprocal
face-position and wave-pointer laws above provide the operational
interaction witnesses.

The [conditional-frame return apparatus](GAUSS.md#gauge-invariant-recombination-readout)
also supplies an invariant interference test in the positive Gauss sector.
After coherent inverse transport, the complete return probabilities are
`(1,0)`; inserting the complete regional-face measurement first changes
them to `(1/2,1/2)`. All ten face/return histories and both final outcomes
are retained, and every branch preserves the constraint. This apparatus
is explicitly controlled by the link presentation without measuring it.

## A faithful reference for this same interacting tick

Apply the complete filter to the joint U, rather than importing a fixed-link
reference:

`K_+=(2I+U)/3`, `K_-=sqrt(2)(I-U)/3`.

Unitarity gives `K_+*K_+ + K_-*K_-=I`. Starting from the declared maximally
mixed state `I/72000`, a successful record prepares

`rho_U = [(5/4)I+(U+U*)/2]/Z`,
`Z=(5/4)72000+Re tr U`.

This state commutes with U. The numerator's eigenvalues are
`5/4+cos(omega)`, at least `1/4`, so the state is faithful on the specified
joint sector.

The trace is computable without a dense joint matrix. For each port p let
`F_p` count the link words fixed by the complete conditional-move product
at that port. Exact enumeration gives

`(F_0,...,F_11)=(6,6,10,10,6,0,10,10,10,6,10,6)`.

Their sum is 90. For each fixed `(g,p)`, ten seam addresses incident to p
have diagonal `Q_g=1/2` and fifty have diagonal one. Multiplication by the
reflection diagonal `-29/30` and the ten code/parity values gives the trace
contribution `-1595/3`. Therefore

`tr U=-47850`, `Z=42150`, `rho_U >= I/168600`.

The success probability is `281/1080`. A capped protocol with k independent
fresh preparations has terminal failure probability `(799/1080)^k` and
retains every failure record. The input mixture, controlled-U filter and
restart service are explicit architecture. Reading the retained preparation
archive changes the conditional state; stationarity concerns the reduced
joint field/wave algebra. This preparation is not a derivation of A3's
information-projection dynamics or of Born probabilities from classical
repair records.

## Accepted frame changes and the interior Gauss constraint

For fixed implemented code frames `V_a`, set `S=direct_sum_a V_a` on the
wave, with identity on parity and address. Re-present every field branch
using links `G'_ba(g)=V_b G_ba(g) V_a*`. Then

`B'=(I_link tensor S) B (I_link tensor S*)`.

The endpoint projectors commute with S, so every `C_e` is unchanged under
this conjugation. The full U and its filter/reference consequently
transform together. The apparatus uses the same inverse-frame/read/restore
wrapper as the [finite agreement theorem](../coherent_seam_repair/PREMISES.md).
Induction on branch histories proves equality of complete accepted records,
including failures, adaptive settings and arbitrary untouched references.
The transformed-link orbit is an accepted presentation domain; arbitrary
continuous local frames need not preserve the central-only link alphabet.
These are operator-level statements for the specified frame wrapper. The
sparse joint interpreter executes the central-link presentation and the
interior transformation below; it is not an arbitrary continuously
recharted-link interpreter.

One allowed central vertex transformation is the port-5 star. It flips all
five interior link bits and preserves the fixed boundary and face fluxes.
Combining that link permutation with Z on the wave at port 5 defines a
Hermitian involution Gamma. Local link covariance proves `[B,Gamma]=0`;
the syndrome-only move decision and endpoint projectors give
`[C_e,Gamma]=0`, hence `[U,Gamma]=0`. Its positive eigenspace has dimension
36000. The [interior Gauss construction](GAUSS.md) specifies its preparation,
reference and complete invariant readouts. This one finite central
constraint does not establish the complete Gauss law of a dynamical
continuous gauge field or a growing observer federation.

## Physical scope

The joint program supplies reversible finite feedback between a protected
regional flux and a coherent carrier wave, with complete operational
contrasts, one fixed clock and its own faithful reference. It supplies no
energy-minimizing field action, relativistic spatial refinement, physical
spin/statistics, scattering amplitude, mass pole or GeV mass. The five-level
code, fixed boundary, interaction grammar and gate order retain their
declared ancestry. Neither a finite walk phase nor the stationary density
is substituted for a physical particle energy.

The [executable receipt](results.json) checks the finite identities and
countercontrols on the literal carrier. The operator arguments above supply
the stated all-state and finite-history conclusions on this finite sector.

## Finite horizon and exact storage

The dimension 72000 counts the admitted joint quantum basis states. Its
logarithm is not the total memory of the executable observer or its classical
exact-arithmetic emulator. The literal link representation uses thirty bits;
wave coordinates, apparatus flags, retained records, program tables and
interpreter counters require additional storage.

Every matrix entry of U has denominator dividing 60: Q has denominator two,
R denominator thirty, and the conditional moves are permutations. For an
input with Gaussian-rational amplitudes `a_i/(D sqrt(n0))`, with integer
real and imaginary components of `a_i` bounded by M, after t idle ticks a
common amplitude denominator is `D 60^t sqrt(n0)`. Unitarity bounds the
integer numerators by `60^t sqrt(2d) M`, where `d=72000`. A dense exact
representation therefore needs
`O(d[log(1+M)+t log(60)+log(d)] + log(D n0))` bits for amplitudes, in addition
to register labels and records. Sparse storage can approach dense support.

The declared Gauss and return projectors have rational coefficients; the
filter adds the fixed algebraic number `sqrt(2)`. Every specified finite
protocol thus admits exact finite arithmetic with counted outcomes, but
unbounded time or an unbounded archive is not implemented by a fixed-size
machine. A finite execution horizon, retry cap and record budget must be
chosen for each run. This accounting supplies neither a physical capacity
law nor a cosmological capacity in nats.
