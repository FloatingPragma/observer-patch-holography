# Two free orbital currents on a native Interval Echo cell

A source shares one unit of load among the five ports bordering a marked vertex. Its phase turns once or twice around that pentagon. A second instrument asks which way charge crosses a bond. It mixes the two endpoints for a quarter period, then uses separate occupancy standards to count **left, right and elsewhere**. A reader averages right minus left over the five consistently oriented bonds around each of twelve marks. Each bond uses fresh starts, including those with no load on the observed endpoints. No B14 projector, fitted angle or doublet gate enters the apparatus. The results concern the simulated observer's terminal records.

| Isolated cell, source 0, one unit of free exchange time | One turn | Two turns |
|---|---:|---:|
| Mean signed reading around launch mark 0 | 0.19085 | 0.11778 |
| Mean signed reading around adjacent mark 1 | +0.08561 | −0.05193 |
| Simultaneous 99% interval for mark 1 divided by mark 0 | [+0.403, +0.497] | [−0.520, −0.370] |
| Ratios across all 60 adjacent source/receiver pairs | +0.4388 to +0.4543 | −0.4611 to −0.4379 |
| Frequency from an independent opposite-ring arrival measurement | 1.6185 [1.6048, 1.6322] | 0.6177 [0.5937, 0.6409] |

## What the instruments measure

For amplitudes $a,b$ on an oriented native seam the local Hamiltonian in exchange units has hopping $-1/2$. The instantaneous continuity flux into $b$ is $\mathrm{Im}(\bar a b)$. A positive half-exchange rotation followed by the endpoint read gives twice this flux. Both physical signs were calibrated on every seam using opposed quarter-phase preparations and a separate short free flight. The largest signed calibration error was 0.00519, inside the simultaneous radius 0.01364. Each sharp occupancy is checked against the reversible native echo with a double-load reference. [Native apparatus](apparatus.py) · [Lean endpoint identity](PortCurrent.lean)

The first source occupies just one ring. The receiving apparatus interferes each endpoint with its antipode on routed native splitters, then reads the remaining difference amplitudes with the seam meter. The independent experiment launches a coherent difference across the two rings, lets it evolve freely, and uses the bare seam meter. It puts the sign selection at the source. All 360 separate endpoint proportions and all 48 circulation means agree within their simultaneous comparison widths on the isolated cell, for the two windings at exposures one and three. Mixing the two rings incoherently, with exactly the same initial port occupations, excludes 20 of 90 proportions and 10 of 12 circulations for each winding at exposure one. Reading the original one-ring launch directly at exposure zero excludes all 24 single-ring forecasts intended for the parity comparison. [Primary protocol](PROTOCOL.md) · [Independent source protocol](SOURCE_CHECK.md) · [Source and clock counts](source_counts.json)

## Why the signed patterns appear

Write $J$ for the native antipodal permutation. An oriented ring mode $f_m$ has equal modulus and phase increment $2\pi m/5$. The isolated adjacency commutes with $J$. On $f_m+Jf_m$ it acts by $-1$. On $f_m-Jf_m$ it acts by $2\cos(2\pi m/5)-2\cos(4\pi m/5)$, which gives $+\sqrt5$ for $m=1$ and $-\sqrt5$ for $m=2$. Free motion thus commutes with forming the antisymmetric source. The source/check equivalence includes all individual terminal probabilities on the two retained rails. The original single-ring source contains both parity components. Its total probability on the initially empty opposite ring is $\sin^2(\omega_m t/2)$, where $\omega_m=|1\pm\sqrt5|/2$ for $H=48I+L/2$. This population read uses ordinary port occupation with no parity sorter or circulation analysis. Both frequencies inferred from exposure one forecast the five held appointments at times two through six. Conjugate launches reverse the circulation and preserve the arrival frequency. Halving the kinetic exchange halves both frequencies against common external appointments. The two materials agree in units of their own separately calibrated exchange hands. Across the clocks' 16 appointments the largest phase error is 0.00676 rad against a simultaneous bound of 0.02047 rad. The clock supplies no conversion to physical rest energy.

The exact integer cyclotomic check applies native adjacency to every oriented launch. A second calculation adds the seam currents in $\mathbb Q(\sqrt5)$ and compares every marked pair to the independently constructed B14 projector. The deterministic adjacency isomorphism interchanges labels 9 and 10 between the two published chart conventions and fixes the other ten labels, including 0 and 1. Relative to the source diagonal, all twelve ring meters give the corresponding projector column $P_{qp}/P_{pp}$. The diagonal signed reading itself is $\sin(2\pi m/5)/5$. In the primary acquisition all 348 isolated absolute forecasts agree within 0.00610. The scalar sign here is the adjacent cosine. A spinor Ramsey experiment with two positive half-turns reads its negative. [Exact check and replay](test_ring_current.py) · [B14 marked-port comparison](../MARKED_PORTS.md)

## What transfers to an open federation

Opening every seam of the 20-cell and 80-cell boards preserves the two signs at the marked pair, for owners 0 and 13 at both appointments. It changes the currents' magnitudes. All 96 absolute forecasts at time one and 76 of 96 at time three fall outside the isolated-cell interval. Allowing each source its own measured amplitude still excludes 24 of 96 normalized profile comparisons at time one and 61 of 96 at time three. For the independent coherent source and bare detector, each of the sixteen open comparisons differs from late parity sorting on at least one ring and 34 to 53 of 90 terminal proportions. The crossed-cell comparison thus resolves a physical difference even when the local signs agree. The two locally accessible patterns provide orbital calibration for a future experiment that couples a free charged excitation to these currents. Neither a spinor coupling nor a unique bracket family is measured here.

## Replay

From the repository root run the following commands. Two independent streams contain 2,294 terminal banks of 131,072 simulated starts each. The decoder includes every empty or elsewhere outcome. Radii use a union bound over the enrolled readouts with a conservative integer-grid allowance. Overlapping reconstructed circulations do not require independent errors.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -m model.experiments.interval_echo_gauge.matter_class.ring_current.measure
python -B -m model.experiments.interval_echo_gauge.matter_class.ring_current.source_check
python -B -m pytest -q model/experiments/interval_echo_gauge/matter_class/ring_current/test_ring_current.py
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_gauge/matter_class/ring_current/PortCurrent.lean)
~~~

The first two commands compare every simulated count with its deposited receipt. [Primary counts](records.json) · [Terminal-only statistical readback](readback.py)
