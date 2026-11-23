# Does an arriving source slow both receiving clocks?

Two preparations of the same written hand at the same gateway give the same phase in vacuum. A sharp load released in the next cell causes a measurable lag only when the clock has a resident metering load. On each screen the incoming source carried about 51 native energy units in independent preparation measurements. Both clocks include every incoming and quiet outcome in their phase reading. [Design](DESIGN.md) · [Program](measure.py) · [Complete counts](results.json)

| Screen | Clock lag for a release one bridge away, in microradians | Simultaneous 99% interval | Clock lag for a release two bridges away, in microradians |
| ---: | ---: | ---: | ---: |
| 20 cells | −284.2 | [−334.1, −234.4] | +1.5 [−48.4, +51.3] |
| 80 cells | −284.9 | [−334.8, −235.0] | +4.1 [−45.8, +53.9] |
| 320 cells | −281.5 | [−331.4, −231.7] | −2.8 [−52.7, +47.0] |

These are occupied-minus-empty phase differences at the last appointment, when the source-free hand reads two radians. A one-bridge release resolved a lag at all three appointments and sizes. A direct launch onto the clock endpoint produced a lag of about 1.997 radians at the last appointment. The changed law retained mobile exchange and restored the clock energy during double occupancy. Every changed-law comparison with the empty meter was consistent with zero. The intervals come from a single 99% bound on all clock quadratures in all preparations and sizes. The two-bridge lag is unresolved at this precision.

## What the receivers measure

The sender prepares a single load at a specified gateway. Separate ensembles in that sharp state measure exchange parity on its six incident blank seams. Each odd result adds one unit to the six contrast costs of eight. The detector in cell zero knows the occupation of its parked endpoint. It makes a terminal number read with three alternatives for the mobile load, the clock gateway, another port of its cell or quiet. With motion parked after that read, metering loads restore the hand endpoints to one each. Inverse echo and a final sharp load read give the first quadrature. An exchange pulse before the sharp read gives the other. The hand angle is the argument of the two marginal signed means on independent ensembles. All source bins contribute to both means. The program checks the analyzer against the native 99-state operators. The enrolments, appointment pairing, terminal metering and gates are supplied services. [Earlier implementation of the written hand](../interval_echo_einstein_station/volume.py) · [Independent sector and propagation checks](test_response.py)

## Which interaction moves the hand?

Let $b$ be the clock gateway, and let $P_b$ project onto a mobile load at that port. The clock's sign parities have empty-contact energies 32 and 33. On each marked contact the native generator acts on the symmetric-coin word of energy $g$ as $g$ on all eight endpoint load words other than $(1,1)$, where it acts as zero. Zero-energy guards on the other seams of the resident endpoint keep its load in place. All contacts evolve at once during the exposure. The empty meter sees energy $g$ in either arm regardless of the source location. The material flight therefore factors out of its Ramsey quadratures exactly.

For the occupied meter, after removing the uniform phase, a mobile load instead evolves with $K-gP_b$ in arm $g$. Here $K$ is the carrier generator with the specified written seams. If $u_g$ is the source wave in this arm, the Ramsey contrast divided by its vacuum value is $C=\langle u_{33},u_{32}\rangle$. Its derivative in the audit exposure is $-i\langle u_{33},P_bu_{32}\rangle$, which gives

$$
|C(t)-1|\leq\int_0^t \sqrt{\Pr_{33}(b,s)\Pr_{32}(b,s)}\,ds.
$$

Thus every phase or visibility response of this meter requires an opportunity for the mobile and resident loads to coincide at the receiving seam. A traversal into the receiver really occurs. In the 20-cell occupied preparation at the last one-bridge appointment, the chance of finding the mobile load anywhere in the receiving cell was 0.000260, and the chance of finding it at the gateway was 0.000211. The alternate local rule restores both written energies at $(1,1)$ at every seam. This removes the clock shift while retaining ordinary propagation. [Native clock algebra and changed rule](measure.py)

This finite record identifies a material-contact contribution to clock delay. It excludes a single source-dependent change of elapsed time as the explanation of the entire shift in these two preparations. The apparatus supplies no independent curvature or proper-volume measurement. Repeating the local placement at three sizes tests the contact response across carriers and leaves a smooth metric limit open. A gravitational interpretation of other clocks or operating regimes requires their own observations.

From the repository root, the following commands replay the full acquisition and independently check the resident-load reduction and three-screen controls. The optional `--output PATH` writes a fresh receipt.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_metric_response.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_metric_response
```
