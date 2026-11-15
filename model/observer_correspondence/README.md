# One process, overlapping observer records and an effective transfer

The [quantum clock extension](QUANTUM_CLOCK.md) gives an ordinary
completely positive effective instrument for local quantum observations,
including an untouched reference. Its source-derived split-wait construction
executes coherent fields and record-controlled phase interventions on the
same twelve ports, with all timeouts and explicit approximation budgets.
It removes the transported-cone requirement for this local comparison;
the exact full-population construction below retains its distinct domain.
The separate [full quantum correspondence](../quantum_observer_correspondence/README.md)
extends that transported-cone approach to the entire total-two operator
domain with exact global decoding. The CP local construction and that
exact global representation have different comparison contracts.

This package is a finite integration test. It compares the **complete same
experiment** under direct native endpoint execution and a derived effective
representation: overlapping observer fields, both local incident clocks,
current transport, record-controlled feedback, every timeout, and the full
terminal source populations. It uses one twelve-port carrier. It does not
establish physical spacetime, a quantum field theory or the identity of the
observed Universe with OPH.

The [existing incident-field theorem](../observer_time_fields/README.md)
preserves one region's complete local histories. The
[relative-field construction](../relative_observer_geometry/README.md)
gives current reconstruction and overlap agreement, using different
executed protocols. Here their joint readouts and feedback execute under one
fixed contract. This is an exact finite change of representation, with new
joint-law checks, not a new general transfer theorem.

## The single contract

| Object | Declared realization |
| --- | --- |
| Source | The C instrument on the 66 binary total-two configurations of the native twelve-port carrier. Only diagonal preparation and interventions execute here. |
| Initial state | Uniform distribution on the sector. The bits at ports `(0,1,5)` are disclosed; all eight labels, including the impossible `111`, stay in the output menu. |
| Observers | A accesses `(0,1)`; B accesses `(0,5)`. Each receives its initial bits and every incident local current increment, including zero. They overlap at port 0. |
| Scheduled evolution | Independent uniform choice among all thirty native seams and a fair native placement coin. Remote transitions continue on the same source. |
| Local clocks | An observer ticks when an attempted seam meets its region. Zero-current attempts count. These are event counts, not seconds or modular parameters. |
| Protocol | Request two incidents on the union `(0,1,5)`, with at most two native attempts per interval. After the first incident, apply one additional native repair on `(0,1)` if reconstructed `n0=1`, otherwise `(0,5)`. |
| Feedback cost | Its current is retained and it ticks both observers. There is no reset after a read or a completed clock event. |
| Output | Initial bits, ordered scheduled/feedback marks, each local clock increment, reconstructed local fields, completion/timeout status and terminal populations. The coin and remote waiting trajectory are coarse-grained. |
| Services | Faithful initial disclosure, complete incident records, a joint comparison cut and the attempt-cap/timeout service are supplied. No physical communication delay or simultaneity is inferred. |

The direct oracle follows actual endpoint branches through each allowed remote
wait. An independent matrix recursion uses the effective maps below. Their
5,768 output words agree exactly, including 5,079 zero-probability words and
all terminal populations. There are 689 positive words. The timeout probability
is `369/625`; successful outcomes are never renormalized. Impossible menu words
can carry inconsistent formal field values, but always have zero state and
probability and are not described as physical trajectories.

Each observer's decoder receives only its local initial pair and incident
increments. At every positive-weight common prefix its fields agree with
native continuity, and the two port-0 values agree. An occupation's first entry
at port 5 is additionally decoded in B's incident clock, with initial occupancy,
no-entry completion and timeout all retained. This is current transport in the
same stationary preparation, not a separately prepared signal or a speed law.

## Exact correspondence, including finite caps

All matrices act on columns. Normalize the native remote and incident counts
by 60, writing `P=A+B` and `B=sum_m B_m`, where a mark retains both local tick
bits and the three local current increments. The union meets twelve seams,
so `A 1=(3/5)1` and `B 1=(2/5)1`. Define

```text
D = I-A,       p = 2/5,       G = D/p,
K_m = B_m D^-1,               F_m = D^-1 B_m.
```

The code stores `60D` as its integer matrix named `D`; formulas here use the
normalized `D`. Since `K_m=G F_m G^-1`, source and effective states correspond
by `u=G e`. Moreover `1^T G=1^T`. Remote repairs leave all union fields fixed,
so their projectors commute with G. Every selective internal feedback map
also commutes with G because its ports and every remote seam are disjoint.
The uniform preparation restricted to a disclosed initial field label is
fixed by G. Thus no extra initial hidden state or reset is inserted.

For cap h, the actual scheduled marked map and timeout map are

```text
K_m,h = B_m (I+A+...+A^(h-1)),       Q_h=A^h.
F_m,h = G^-1 K_m,h G,                Q_effective,h=A^h.
```

Conjugation cancels between every successive read and allowed feedback map.
The terminal decoder G then gives the complete unnormalized source state on
each word, including words ending in a timeout. This is an identity of maps;
induction proves it for every finite adaptive protocol in this declared
commuting intervention family. The retained execution exercises the stated
two-interval protocol. It does not certify arbitrary coherent preparations,
quantum references, boundary-straddling operations or inaccessible records.

Finite-cap transformed maps have 2,316 negative entries in this instance.
They preserve the **transported cone** `G^-1 R_+^66`, with decoder G; they are
not claimed to be ordinary entrywise-positive stochastic maps on arbitrary
effective vectors. Individual current-marked maps are also not asserted to
be symmetric or positive semidefinite. The checks compare decoded full states,
so a matching scalar probability cannot conceal an incorrect disturbance.

The remote probability per attempt is `3/5`, independent of state. Therefore
the probability that either of two intervals times out is
`1-[1-(3/5)^h]^2`, even with the declared intervening feedback. The same bound
couples the capped and completed output experiments without postselection.

## Which action and spatial dictionary this attaches

For the **unread completed** transfer `F=D^-1 B`, use the inner product with
metric G. Every actual native seam count is checked to equal `I+SWAP`, with
the swap a symmetric involution. Hence its normalized map is positive
semidefinite and bounded by identity. Consequently
`0<=A< I`, `B>=0`, `D-B=I-P>=0`, and `GF=B/p` is symmetric positive
semidefinite. The native triangle `(0,1,5)` contributes at least `3I/60` to B;
its exact count matrix satisfies `(C-3I)(C-6I)=0`. All other incident terms
are positive. Since `D<=I`, the G-self-adjoint transfer obeys

```text
1/20 <= F <= 1,       H=-log(F),       exp(-nH)=F^n.
```

H is defined as the matrix logarithm of the given transfer; it is the
effective Euclidean generator per completed union incident tick by definition.
The marked insertions and decoder are part of the history dictionary.
The native repairs are not the unitary continuation `exp(-isH)`.

The existing [C history action](../native_action/README.md) also attaches at
every native event boundary: `A0=pi(I-U)/2` and `exp(-iA0)=U` in its declared
20-dimensional local dilation. All 3,960 state/seam/coin branches are checked
against its actual C Kraus numerators. `A0` is likewise defined from the given dilation `U`. This is a gate-action
representation of the same source events, not a physical interpolation
between them.

The source's selected rank-three slow-response Gram is checked and restricted
to `(0,1,5)`. These three native neighbors have diagonal Gram entries 1 and
off-diagonal entries `1/sqrt(5)`. This supplies the same-source finite spatial
dictionary for the two regional restrictions and their allowed repairs.
Its incidence structure and spectral selection are premises. Coordinates
are not inferred from this two-event sample, and Gram rank is not an
established dimension of physical space.

## Controls and open physical obligations

Conditional on a union incident, the tick pairs `(A,B)` have probabilities
`P(1,0)=P(0,1)=1/4`, `P(1,1)=1/2`, and `P(0,0)=0`. Each marginal is `3/4`.
Independent clocks with these same marginals change the joint law by total
variation `1/8`; including first-interval cap-two timeout gives `2/25`.
Thus matching separately reported clocks is insufficient. Further tests
reject dropped timeouts, corrupted terminal populations with unchanged
probabilities, false nonincident labels and reuse of the stationary boundary
for an arbitrary basis preparation.

The [continuity theorem](../../Lean/Screen/RegionalContinuity.lean) supplies
the local current identity; positive exchange and finite instrument dilation
are instantiated above. No homomorphic quantum local-net retraction,
Lorentzian microcausality, modular-time identification, renormalized continuum
field theory, physical action normalization, matter correspondence or
comparison with measured spacetime is supplied. A physical equivalence claim
must specify shared accessible observations and interventions and bound their
joint output discrepancy on the same declared domain. This finite certificate
provides one exact diagonal instance of that criterion; it is not evidence
that its physical identification holds.

```sh
python -m model.observer_correspondence --check model/observer_correspondence/results.json
python -m model.observer_correspondence.quantum_experiment --check model/observer_correspondence/quantum_results.json
python -m pytest -q -W error model/observer_correspondence
```
