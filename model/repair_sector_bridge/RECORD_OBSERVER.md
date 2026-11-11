# A bounded observer of native pair separation

The native repair record supplies the coin and signed current on a known
addressed seam. For a declared known binary preparation, these records
determine every later load exactly:

    n_a' = n_a - j,    n_b' = n_b + j.

`record_observer.Observer` receives a preparation label, a capacity and the
ordered addressed records. It never receives the source's working loads.
Its seven-bit tracked-state index determines the distance between the two
occupied ports. Exhaustive comparison covers all 66 initial states and all
60 seam/coin branches. Correct initialization and the update identity then
prove correctness for every finite history by induction. Two observers
given the same preparation label and complete records consequently agree.
This is a classical record compiler, not a measurement of an arbitrary
unknown quantum state.

The preparation label and addressed uniform-seam scheduler are declared
interfaces. The fixed two-event observer has its own scheduler and does not
acquire this protocol by a relabelling. An unknown initial state is not
reconstructed by assumption: two different preparations can yield the same
first current record and different pair distances. Missing, duplicate,
out-of-order and incompatible records halt inference. No-op records remain
on the tape. The source protocol halts when the declared tape capacity is
full; the last distance does not describe unobserved subsequent motion.

## Exact relaxation of the recorded observable

Under independent uniform seam choice and the fair native coin, the actual
recorded distance has the three-state kernel in [PROOF.md](PROOF.md).
Its one-attempt Dobrushin coefficient is exactly `14/15`, giving

    TV(mu K**t, pi) <= (14/15)**t,
    pi = (5/11, 5/11, 1/11).

Direct rational matrix powers give a sharper finite answer. The largest
total-variation error over all initial distance classes first falls below
`1/100` at **39 native attempts**, where it is approximately
`0.009811214413`. The corresponding maximum at attempt 38 exceeds `1/100`.
Convexity covers every initial distribution. This certifies the admitted
readout without Monte Carlo equilibration or a fitted temperature. It
does not select a gauge-sector Gibbs law or identify physical time.

## Record and memory budget

A seam has 30 addresses. Its binary native letter is one of
`(0,0), (0,1), (1,0), (1,-1)`, so the addressed record alphabet has 120
values and needs seven bits. A tape of `T` attempts therefore uses at most
`7T` bits. The retained preparation label and tracked-state index use seven
bits each, the counter uses `ceil(log2(T+1))`, and three status bits and two
distance bits suffice for the reported logical registers.

At `T=39` the retained-register bound is **298 bits**. An incoming packet
including its sequence number adds at most 13 transient bits, for **311
bits** altogether. Retaining every distance output as a second tape adds
78 bits. Immutable program/geometry storage and physical implementation
costs are separate. Overwriting a working register is not asserted to be
thermodynamically free. No failed message is silently accepted as a source
event or replaced by an inferred current.

```sh
python -m model.repair_sector_bridge.record_observer
python -m pytest -q -W error model/repair_sector_bridge
```

The command reexecutes the native branches and exact rational relaxation
calculation. [record_observer.py](record_observer.py) defines the bounded
controller and the logical register accounting.
