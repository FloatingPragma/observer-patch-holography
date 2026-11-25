# Exchanging two charged three-port cores

Two vacancy–excess–vacancy trimers on the registered Interval Echo material exhibit exchange phase **+1** on an eight-port ring. A second experiment releases both trimers together. At a delay of 65,536 native time units, **778,428 of 1,048,576** port-load records show both intact at a displaced pair of positions. Giving the two trimer centers negative exchange statistics, with identical local contact energies and hopping magnitudes, yields **46,005** such records. The instruments resolve the statistics of this moving, net-negative source. Its exchange disagrees with the electron assignment. [Simulated receipts](results.json) · [Instrument and native generator](measure.py) · [Independent labeled-center check](test_exchange.py)

## How the material carries the two charges

The eight physical addresses in owner zero are 0, 1, 5, 4, 2, 3, 6 and 7 in that cyclic order. Each trimer has loads (0, 2, 0) in a sea of ones, carries net load −1 and has one excess that marks its center. Two trimers fit at opposite positions with one sea port in each gap. Excesses outside the ring compensate the total board load. Each private contact on the ring receives a native (0, 0) silent write and every unused physical contact touching the ring receives a (1, 1) silent write. The receipt gives the neighboring swaps that present the correct loads and then restore the launch word. These heads annihilate the native generator for every later input on the ring. The underlying single-trimer construction also supports freely propagating bound cores on a fifteen-port route across five owners. [Single-core material and binding proof](../mobile_bound/README.md) · [Silent-head and blank-column verification](../conjugate_flight/measure.py)

## What separates exchange from circulation and internal rotation?

At the registered contrast coefficient eight, the blank native contact evolves for a time π to the ordinary swap on all nine possible load pairs. Advance one trimer one port by moving its leading vacancy into the gap, then its excess, then its trailing vacancy. Advancing both trimers four times exchanges their positions and preserves the orientation and order of all six constituents. The comparison branch holds the source. This is a controlled transfer, with intact trimer checks in the word calculation at each completed three-hop translation. The individual hops pass through intermediate noncompact configurations.

A separate experiment advances **one** trimer eight times. Both experiments use 24 physical swaps. On each of the eight directed bonds, each program carries a vacancy twice and an excess once. Consequently any fixed one-body phase per directed link and load species contributes exactly the same product to both sharp-source circuits. With a coherent path reference, both native circuits yield 1,048,576 bright outcomes in the real quadrature and zero dark outcomes. Reversing the reference sends all starts to the dark port. The calibrated exchange phase is +1. A comparison assigning half a flux quantum to the center changes both raw signs. A fermionic exchange changes only the two-center sign. The reference, address selection and clock belong to the instrument.

This circuit acts by an ordinary permutation on the tensor slots of the two three-port blocks. Its complete action on two aligned internal wavefunctions is the block flip. For two identical wavefunctions with amplitudes $u(a)$, the exchanged amplitude at $(a,b)$ is $u(b)u(a)=u(a)u(b)$. This holds for every coherent internal state carried rigidly by the slots, including states odd under local reflection. It gives a specific obstruction to obtaining fermionic statistics by adding internal structure to otherwise unchanged load transport. It leaves open charge carriers with a mobile history or field whose evolution supplies a path-dependent exchange phase. [Finite phase and block-flip lemmas](RigidExchange.lean)

A native control holds the same negative charge but resolves internal parity. Start each trimer with its excess at the left end. A path-controlled reflection, composed of three adjacent swaps, and a dark reference result prepare the normalized difference of the left-excess and right-excess words. Both preparations succeed in one quarter of starts on average. Exchanging these two identical odd states gives **262,237 bright and zero dark** responses, with 786,339 preparation refusals. Appending the same physical reflection on just one copy gives **zero bright and 261,778 dark**, with 786,798 refusals. Thus the detector resolves a genuine internal minus sign separately from the phase on exchanging whole charges. This control concerns addressed transport of the endpoint-excess superposition. Free binding was tested for the central-excess source.

## What happens on free release?

The generator sums the actual blank contact columns over all eight ring bonds, with all other incident terms silent. The retained bases have 168 sharp material words for one trimer and 420 for two. Every start is read at all eight ports by the ternary meter. The measured intact, displaced event requires a translation of the entire sharp launch word with different port addresses. The other addresses, including breakup configurations and returns, remain in the receipt. The table gives counts for delay 65,536 in a common seven-delay ladder. Each entry has 1,048,576 starts.

| Source | Native | Negative exchange of the centers | Extra half-flux | Both changes |
|---|---:|---:|---:|---:|
| One displaced intact trimer | 982,912 | 982,842 | 470,109 | 469,701 |
| Two displaced intact trimers | 778,428 | 46,005 | 46,014 | 777,387 |

At this delay the native two-trimer detector finds an intact separated pair at any address in **1,034,206** starts. The independent one-trimer exposure resolves the half-flux ambiguity of the two-trimer data. The two-center comparison keeps each contact energy and all local hopping magnitudes. Give each excess a private label, hop it with its load, forbid passage through another excess at the same port, and compress the resulting dynamics into symmetric or antisymmetric center-label amplitudes. In the basis ordered by center address, a single excess crossing the cut passes all other centers in that ordering. Its off-diagonal entry acquires the factor $(-1)^{c-1}$ for $c$ excesses in the antisymmetric sector. This definition also specifies the comparison dynamics during breakup. The independent tests build the full labeled cover and check that both projections intertwine the generators. Repeated squaring of the full propagator checks the discriminating free-flight probabilities. The result pertains to this finite ring, source and comparator. It supplies no three-dimensional exchange theorem.

From the repository root, generate an independent receipt at an unused path and check the committed measurement.

```sh
d=model/candidate_laws/interval_echo/openwave/row_11_spectrum/composite_exchange
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.composite_exchange.measure --output /tmp/oph-trimer-exchange.json
python -m pytest -q "$d/test_exchange.py"
(cd Lean && lake env lean "../$d/RigidExchange.lean")
```
