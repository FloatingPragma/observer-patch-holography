# When a single repair produces two copies of a record

Open one native repair between two cells and keep its coin and signed head unread. Each owner then measures its own endpoint with two local echo standards. On all three boards an input with an even total load gave both readers the same answer on every launch. For input (0,1) the two owners each saw zeros and ones with equal probability, but on **all 65,536 starts per board their answers differed**. The smallest board gave 32,707 events of one ordering and 32,829 of the other. An observer who checks only the single-reader frequencies misses this distinction. The source, the selected gateway and its shutters are supplied. The repair itself produces the two local records without a read-conditioned instruction. [Frozen experiment](PROTOCOL.md) · [Complete record bank](records.json.gz) · [Acquisition and count decoder](measure.py)

The two observers were owners 0 and 6 on the 20-cell board and owners 0 and 3 on the two refinements. A balanced-pair pointer read against a zero standard and then a two standard gives words 01, 00 and 10 for loads zero, one and two. Undoing each detector echo preserves the measured load. The code checks those instrument columns on the full 99-dimensional seam. The word 11 is an inconsistent detector report. A missing acquisition has a separate error. The 160 settings contain 10,485,760 synthetic starts drawn from finite Born probabilities. Every endpoint pair, disagreement and clock reading is included. The simultaneous 99% half-width for the enrolled probabilities is 0.010038. [Native instrument and free generator](apparatus.py)

## Which source questions both ends can answer

For a sharp input $(a,b)$, write $s=a+b$, $l=\lfloor s/2\rfloor$ and $h=\lceil s/2\rceil$. The write places the two branches $(l,h)$ and $(h,l)$ in distinct private coin states. On even sources both endpoints therefore carry the same sharp value $s/2$. The three inputs (0,2), (1,1) and (2,0) all announce one. The other even inputs (0,0) and (2,2) announce zero and two. At an odd sum the two owners see the adjacent values in complementary order. In every case the two readings together recover the original sum.

This argument covers coherent and mixed inputs with an arbitrary retained reference. If a matrix element survives tracing the private head and the opposite endpoint, its two inputs share the signed contrast, coin and opposite endpoint load. Those three quantities determine the incoming pair uniquely. Thus either single-site reduced channel carries the diagonal matrix with weight one half at $l$ and $h$ for an input matrix unit on its diagonal, and gives zero on each off-diagonal input unit. The endpoint channel acts linearly on the corresponding blocks of a spectator reference. For an input confined to the even sector, the joint classical report plus reference is exactly the result of measuring the incoming sum and writing its value at both ends. This copies the measured value while preserving the conditional reference statistics. All 81 input matrix units are checked independently against both native output channels. [Lean identities for the traced matrix elements and sharp questions](Endpoints.lean)

For any enrolled subset of incoming labels, join two labels when one reader could obtain the same load from either. A sharp source question is readable at both ends precisely when its answer is constant on each connected component. Necessity follows from the shared possible reading. Sufficiency follows by assigning that answer to each reachable local level. Even an arbitrary local quantum measurement reduces to its diagonal on these single-site outputs, so it cannot distinguish a shared level with certainty. On the full nine-label domain the odd sums connect all three levels. Every exactly recoverable sharp source question is then constant, although the two marginal output states are equal for every density matrix.

| Maximum distinct sharp source answers at either end | Nonempty input supports |
| ---: | ---: |
| 1 | 393 |
| 2 | 111 |
| 3 | 7 |

All seven ternary supports include (0,0), (2,2) and at least one source of total two. The count checks all 511 supports and independently compares every binary coarsening. It describes the specified blank gateway and its two disjoint readers. A board-wide characterization of which record families the native process makes public remains open. [Classifier](classify.py) · [PR-130 publication experiment](../README.md)

## What a later reader can still recover

Every just-written column lies in the exact kernel of the guarded native gateway generator at any contrast. Writing equals applying the involution $U$ to a blank input $J$; the local generator is $UG_\kappa U$ and $G_\kappa J=0$. All isolated holds therefore preserve the endpoint readings. Both a parked write and a modified repair that stores the same signed head and coin while preserving the two incoming loads fail to broadcast the (0,2) source. These two controls have the same immediate endpoint readings.

Opening one further blank contact from the first owner's port to a vacancy changes the result. That arm and the original gateway remain active for the entire hold, and the clocks are read on separate calibrated copies. At the half-turn appointment and native contrast, source (1,1) gave unequal readings on all 65,536 starts. For (0,2), which gave the same two ones directly after repair, 32 starts disagreed here. Removing the contrast from the two free contacts raised that last count to 43,087. The direct output record does not predict its own later availability without specifying the continuation.

In a separate reversible comparison, either half-exchange phase of the (0,2) source gave two endpoint readings of one. Unwriting the gateway and reversing the preparation restored the source on every trial for both phases. Reading the private signed contrast in addition reduced the returns to 32,920 and 32,568. Agreement on the public sum is compatible with recoverable coherence within a sum sector.

Run from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -m model.experiments.interval_echo_consensus.publication_tree.repair_outputs.measure
python -B -m pytest -q -W error model/experiments/interval_echo_consensus/publication_tree/repair_outputs
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_consensus/publication_tree/repair_outputs/Endpoints.lean)
~~~
