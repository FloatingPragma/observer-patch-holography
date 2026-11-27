# Where the two conjugate preparations can still differ

Before conversion, a charged vacancy arriving at port 160 can read the winding phase. With a material source and material detector, the positive preparation gave 41,065 phase hits in the positive analyzer and 1,551 in the negative analyzer. The negative preparation gave 1,573 and 41,061, respectively. Each setting used a separate bank of 2³⁸ original simulated starts. Misses and refusals are in the same records. [Earlier charged-vacancy acquisition](../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge/results.json)

After conversion and neutral release, every single terminal measurement on the active ring alone gives the same distribution for the two conjugate recipes. The independent finite audit uses the actual source pulses and the complete charged source generator. Its six arm permutations preserve all matrix entries. For each of eight sharp words at the tip, swapping two equal bits holds its measurement fixed and interchanges the winding signs. Their populations thus agree before conversion at any charged-flight time. Once the neutral pair is released, the bit on either written seam XOR the indicator of a nonunit load on its leaf is conserved. This holds on every one of the 3,472 directed exchanges in the 960-word release chart. The two leaves and their heads distinguish all four launch sectors. Ignoring these records removes interference between sectors at any later free-flight time. Each sector launches the same localized pair for both winding signs. These statements assume the original spectator and source records remain unavailable. They include pair misses and partial arrivals. [Finite support and source audit](audit.py)

A different instrument could compare a remote pair click with a coherent two-head read at the conversion site. We calculated the full distribution for such an explicitly supplied instrument. Each remote click includes every ordered pair at the chosen owner. At the source the two quadratures compare (01 ± 10)/√2 and (01 ± i10)/√2 on separate copies. The X and Y entries below are the differences between the two source preparations in positive-minus-negative analyzer probability, conditional on the remote click. The blank and doubly written local outcomes remain in the denominator. [Protocol and coherence assumptions](ERASER_PROTOCOL.md) · [Calculation](eraser.py)

| Remote owner | Time | Remote pair probability per original start | X difference | Y difference | Magnitude |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 12 | 256 | 3.8150 × 10⁻¹² | −0.03445 | +0.26356 | 0.26580 |
| 12 | 512 | 6.3526 × 10⁻¹² | +0.13233 | +0.19902 | 0.23900 |
| 12 | 2048 | 1.8490 × 10⁻¹¹ | −0.03774 | +0.20215 | 0.20564 |
| 17 | 256 | 8.5387 × 10⁻¹² | +0.04229 | −0.01861 | 0.04620 |
| 17 | 512 | 6.0474 × 10⁻¹² | −0.08968 | −0.00002 | 0.08968 |
| 17 | 2048 | 8.0838 × 10⁻¹² | −0.07883 | +0.15867 | 0.17717 |

The marginal pair predictions agree with all 54 held native receiver banks across both sites, three appointments, three character preparations and three terminal pulse settings, using a simultaneous 99% binomial check. The table's joint columns are forecasts for an additional measurement. A fixed phase difference through parking would rotate the two columns; decoherence between the retained sectors could erase both. The earlier count banks verify the *charged* local analyzer and the *neutral* remote pair arrivals separately. A physical two-head analyzer on the released material would be a new experiment. [Count comparison and sharp standards](test_eraser.py)

The subsequent two-head experiment tested this forecast with a borrowed vacancy and two material paths. The comparison at receiver 12, wait 256 resolved a 6.21 percentage point difference on the positive reference panel. A separate single-bypass construction supplied a finite bound for arbitrary stable parking phases. [Physical joint read](joint_receivers/README.md) · [One-bypass bound](theory/source_access/README.md)

Run the checks with one numerical thread from the repository root.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations
python -Bm "$p.audit"
python -Bm "$p.eraser"
python -Bm pytest -q --pyargs "$p.test_audit" "$p.test_eraser"
~~~
