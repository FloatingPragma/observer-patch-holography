# A receiving patch times its matter-phase observations

Cell 6 holds two material hands and the matter detector. At three scheduled appointments, its two local fuel readings reconstructed angular times of **0.24945, 0.49913 and 0.99806 radians** in the positive-mixer ensembles. In those same trials the patch also recorded a bit for matter on port 78. Comparing fresh ensembles under conjugate settings of the (82,78) mixer gives the matter contrasts below. Each simulated preparation yields one of eight joint three-bit outcomes. [Protocol](CLOCK_PROTOCOL.md) · [All counts and held checks](clock.json)

| Supplied flight appointment | Phase from patch records, positive mixer | Phase from patch records, conjugate mixer | Sampled matter contrast |
| ---: | ---: | ---: | ---: |
| 0.25 | 0.24945 | 0.24741 | +0.03866 |
| 0.50 | 0.49913 | 0.49895 | −0.00891 |
| 1.00 | 0.99806 | 0.99706 | −0.02839 |

The free matter sector on this preparation contains 2,730 reachable configurations. Neither endpoint of seam (76,77) or (80,81) occurs in any of them. Each clock seam can therefore hold a different sign preparation at the same fixed coin. Its symmetric and antisymmetric sign states have native energies 8 and 9 and remain confined to the unloaded seam for every matter configuration. All other terms act on the matter registers and their own collars. The free evolution factorizes exactly between the two hands and this matter sector. The clocks continue to evolve in a cell that the unit reaches. [Reachability check](results.json) · [Native generator and product check](clock.py)

The hands start at a sharp positive sign with one hand advanced by a supplied quarter turn during setup. A bare echo at empty endpoints leaves either sign unchanged and cannot reveal it through the loads. The receiver instead uses its local fuel and two admitted converters. The first couples the even sign record to blank with fuel present; the second couples the odd sign record to the same fuel. Their coherent paths interfere. The calculation includes the free generator during the finite pulses and checks their restriction on the full 198-state fuel and seam space. Four independent calibration preparations give a sinusoidal response amplitude of 0.49919. The two terminal fuel bits, together with the reference phase offset, give two quadratures across the sampled preparations. [Fuel converter](../../../contrast_memory_spherical_v1/collar.py) · [Local joint-phase detector](phase.py)

A negative-sign control shifts the inferred phase by π. Over both signs, all appointments and both mixer settings, the largest circular error was 0.002942 radians. Preparing parity eigenstates instead kept the two fuel-click rates stationary within a sampled range of 0.001163. The largest difference between a clock run's matter contrast and the independent prediction was 0.000539. Every calibration and terminal setting used 1,048,576 simulated preparations.

These are readings of an ensemble phase modulo 2π. Each individual preparation supplies only two hand bits. The controller selects release and terminal windows; no patch counted repeated ticks or triggered the appointments in this experiment. The readings put these local matter measurements on a physical phase scale without furnishing a matter frequency or an independently prepared wave number.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_30_arrows.clock
```
