# How far apart are the weakly bound excesses?

The two excesses on the twelve-port cell are usually two ordinary seams apart. In the native joined-energy window, 246,492 of 252,983 accepted launches had nonadjacent occupied ports. Two cold single excesses on independently prepared copies of the same cell were even farther apart when their observed port addresses differed. The joined pair has positive binding energy and **fails the enrolled relative-extension test**. [Prospective test](PROTOCOL.md) · [All original-start records](counts.json) · [Ground comparison](ground.py)

## An observer's ruler on the same clocked pair

On uniform sharp two-excess stocks a finite-width energy filter accepts at 84.72475 exchange units. The accepted and missed counts are both recorded. A supplied comparison of signed port echoes computes the shortest seam distance between the two excesses. This reversible pointer uncomputes intermediate addresses before its one read. A separate apparatus measures both individual ports and then forgets them apart from the near/far flag. Its backaction can distinguish those preparations. The cell has 30 adjacent port pairs, 30 next-nearest pairs and six antipodal pairs. The measurand is a graph distance in seams. [Original positive-stock energy survey](../README.md) · [Acquisition and decoder](measure.py)

At wait zero, accepted native counts in these three shells were 6,491, 246,056 and 436, from 16,777,216 original starts. The nonadjacent fraction was 97.434% (simultaneous 99% interval 97.283–97.580%). The unfiltered uniform pair gave 54.562% (54.504–54.619%). On two independently clocked cold single-excess copies, the same signed port indices were compared and equal addresses were kept as a separate outcome. Among 77,872 distinct-address coincidences, 77,259 were nonadjacent, or 99.213% (99.053–99.354%). The corresponding mean squared distances were 3.932, 3.092 and 3.986 seams squared. Most of the excess over the uniform-stock reference therefore follows from how each cold single occupies this material. Only 0.172% of joined bright pairs were antipodal. Halving contrast gave a joined nonadjacent fraction of 85.282%. At zero contrast all three nonadjacent fractions were compatible with the uniform geometrical value of 6/11.

A separate certificate gives the same ordering for the true native ground levels. Exact fraction-free inertia of twice the registered matrices isolates one eigenvalue below 84.80 on the joined pair and one below 42.46 on the cold single. The deposited positive vector x for each stock gives exact component ratios (Hx)ᵢ/xᵢ in [l,h]. For normalized x, the norm of (H−l)x is at most h−l. On the orthogonal complement of the ground this bounds the norm of x's component by (h−l)/(threshold−l), the trace distance between its pure state and the true ground. Any port predicate changes in probability by at most that distance. Two independently prepared singles cost at most the sum of their errors; the distinct-address coincidence must then be renormalized. This encloses the genuine nonadjacent probabilities in **[0.97427431, 0.97427534]** joined and **[0.99186620, 0.99186641]** on distinct addresses of the two cold singles. The enclosure evaluates rational numbers before reporting rounded decimal limits. [Exact inertia and port weights](ground.py) · [Deposited comparison vectors](../two_counts.json)

## First observed separation depends on the read schedule

In fresh trials the clock was followed by the adjacent/far pointer at zero. Of initially adjacent, clock-bright copies monitored every unit through time 16, 623 of 6,648 first reported a nonadjacent pair. Monitoring every four units instead gave 156 of 6,388. The respective simultaneous 99% intervals were 7.77–11.16% and 1.64–3.48%. A single final shell read at time 16 on another initially adjacent cohort found 0.672% nonadjacent. On initially far native pairs, 4,545 trials returned adjacent at that same appointment. These measurements distinguish monitoring schedules at a fixed total duration.

Without the first shell read, uniformly sampled sharp words have density I/66. A successful clock makes its unnormalized density proportional to a function of the same Hamiltonian, so the entire untreated shell distribution is stationary at every free wait. The five native waits gave 97.37–97.44% nonadjacent at each terminal. A first detection under monitoring is a recorded event. This closed finite graph has no outgoing fragment counter, and its recorded returns distinguish shell crossings from irreversible decay. The measurements establish neither a nuclear halo beyond the matched cold reference nor an intrinsic escape lifetime.

From the repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_17_nuclei.like_stock.outer_pair
python -Bm "$p.measure" --acquire /tmp/oph-outer-pair.json
cmp /tmp/oph-outer-pair.json model/candidate_laws/interval_echo/openwave/row_17_nuclei/like_stock/outer_pair/counts.json
python -Bm "$p.measure"
python -Bm "$p.ground"
python -Bm pytest -q --pyargs "$p.test_pair"
~~~
