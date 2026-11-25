# Do two charges acquire a minus sign when they exchange and wind the written heads?

On a fifteen-port loop through five owners, exchanging two native negative charges and winding one charge around a parked partner produce the same final loads **and** the same fifteen reversed seam heads. The relative exchange phase is **+1**. A quantum superposition of the two winding records survives the comparison. Its odd winding state acquires a minus sign on either route. A separate autonomous release distinguishes native and negative exchange statistics through the two material port readings alone. [Native sources and measurement](measure.py) · [Full address and interferometer counts](results.json) · [Independent labeled-charge tests](test_exchange.py)

## Two histories ending in the same record

The route uses physical ports 3, 2, 10, 20, 15, 14, 24, 25, 32, 46, 36, 37, 50, 58 and 48, returning to port 3. A real blank contact joins port 25 to parking port 29. Start vacancies at ports 3 and 25. Two excesses on distant, separately isolated dimers compensate them. Each of the fifteen ring links has a writable even, unit-contrast head oriented to allow a forward crossing. An alternate preparation changes the parity of one such head. Private references on the route are zero-written and each unused incident seam is closed at unit filling. The stock replay checks the bypass, sharp writes, shutter parking and final loads. The native echo supplies the parity and direction preparations. [Native 99-state blocks](../live_trimer/measure.py) · [Even and odd writes](../two_winding_loops/measure.py)

On the exchange arm, advance the two vacancies in alternation around complementary arcs until they have traded starting ports. Each hop has one vacancy at its source and one unit of material at its destination. On the other arm, move the vacancy at 25 to port 29, wind its partner once around the ring and return the parked vacancy to 25. Retaining the order of two bookkeeping labels verifies that only the first arm exchanges the charges. Both arms move a vacancy forward over every ring contact exactly once. The extra two blank contacts are ordinary full swaps with unit amplitude. The arms thus match every one-body circulation phase as well as the complete final direction record.

For a directed crossing of seam $e$, let $A_e$ be the charge-transfer matrix element of its native contact pulse, acting on its head space. A forward-current herald keeps the prescribed destination and refuses the other output. On the accepted branch the two programs contain exactly the same factors $A_e$. Factors on different seam heads commute. This proves equality of the two complete accepted operators even on entangled input heads and an untouched quantum reference. A comparison with antisymmetric identical charges changes the sign of the exchanged-label arm. A one-body half-flux reverses both amplitudes. The Ramsey instrument supplies a coherent path control, addressed pulses, an exposure clock and the common current herald. The physical phase comparison assumes the route selection is coherently erased and that accepted detector records reveal no arm label. The free-release measurement below uses direct native evolution.

One sharp signed-contrast echo input, followed by a direction herald, prepares an equal superposition of even and odd contrast on the link between ports 2 and 10. The other fourteen links have even contrast. With an exposure of π/8 on each written link, the two parity components have different transfer amplitudes. In the native real-reference bank the counts are **514,771 bright, zero dark, 9,952 transfer refusals and 523,853 source refusals** in 1,048,576 starts. The comparably instrumented negative-exchange model gives zero bright and 513,749 dark. The resulting contrast information is retained when the charge interferes. Both source and current refusals count against the original exposure bank.

## A winding qubit and its own minus sign

For two charges on their starting ports, the accessible all-even sector contains two winding records. They differ on **all fifteen heads**. Call them $I$ and $F$. A completed forward arm acts on winding as $|F\rangle\langle I|$ times a common phase. Its reverse sends $F$ to $I$. Coherently select these directions with an equal path control, herald completion, and accept the plus output on recombination. In the matched filled-host convention the operator is $-X/2$ for **both** the exchange and the parked reference. Its success probability is one quarter on every winding input, including a state entangled with a held reference. Such an input can be prepared on the all-even guide by coherently selecting a hold or the deterministic full winding. The complementary projection rejects the unsuccessful path flags.

The relative odd/even winding eigenphase is −1 on each arm. After matching that internal operation, the native real-reference winding bank has **262,128 bright, zero dark and 786,448 refusals**. Reversing the reference gives zero bright and 262,666 dark. The antisymmetric comparison instead gives zero bright and 262,450 dark with the original reference. This is an exchange measurement on a shared winding degree of freedom. On this circuit the charge positions do not each possess a separate localized winding qubit. The exact source-connected sector has 82 charge-address pairs with one winding record each and two pairs with two records. The two-record pairs differ on the whole cycle. At fixed material addresses, equal local charge divergence forces the difference of two directed flows to be constant around the ring. [Finite flow and gauge lemmas](MatchedWinding.lean)

## Free evolution with unread heads

Simultaneously enable the ring and the spur, release the sharp source and count both terminal vacancy addresses. The negative-exchange comparison gives each vacancy a tag and transports that tag with its physical hop. Projection into symmetric or antisymmetric tag order leaves all diagonal interactions and all hopping magnitudes equal. The tests construct the labeled cover, compare both projections and independently square the full propagator. The release bank retains all address bins at thirteen power-of-two waits from 1/4 to 1,024. The following are the returns to the source address or address pair at time 512, the largest native-versus-negative difference for the all-even pair over this ladder. Each column uses 1,048,576 starts per setting.

| Source and device | Native | Negative exchange |
|---|---:|---:|
| One vacancy, open spur | 236,763 | 235,961 |
| Two vacancies, open spur, even ring | **437,557** | **114,130** |
| Two vacancies, open spur, one odd link | 335,858 | 203,314 |
| Two vacancies, closed spur, even ring | 163,454 | 163,361 |

The active, two-charge generator spans 86 reachable load-and-head words. With the parking contact closed, changing the statistics on the ring is exactly a diagonal rephasing by the last head direction. It preserves every finite-time charge histogram for this sharp source. A single vacancy also has identical predictions under both choices of exchange sign. Adding a half-flux on the ring cut is that same head rephasing for either exchange model with **or** without the spur. Opening the spur allows two histories with opposite label permutations to reach identical physical records, making their exchange sign measurable independently of circulation. The evidence identifies positive exchange on this source, geometry and accessible component. Other carriers with independent mobile internal states require their own exchange paths and matched record comparison.

From the repository root, run with one numerical thread and an unused output path.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.winding_exchange.measure --output /tmp/oph-winding-exchange.json
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/winding_exchange
(cd Lean && lake env lean ../model/experiments/interval_echo_angular/charged_orientation/winding_exchange/MatchedWinding.lean)
```
