# A released charge probed in freely exchanging material

On the 320-cell Interval Echo screen, a freely released vacancy grows from an rms radius of 0.9484 to 1.5892 links between exchange times one and three. A static probe then sees a strong scale dependence. At time three its signed interaction with the captured source is −7.56550, −4.12929 and −0.80899 exchange units in the first three shells about the launch site. A single independently calibrated contact strength and the separately sampled source profile predict −7.56414, −4.13253 and −0.80856. In the registered core-radius window the sampler intervals exclude an intrinsic 10% change. [Frozen protocol](PROTOCOL.md) · [Complete native probabilities and sampled records](records.json.gz)

**Full OW-23 remains Partial.** This test measures a material signed interaction and its spatial extent. Every captured configuration, every arbitrary blank-bath input common to the four sign settings and every static hold have the same four-sign mean coefficient on adjacent guard sites, and zero beyond them. A spatially unbiased bath with any fixed or mixed mobile inventory also has exactly zero terminal polarization response to either source sign. A differently prepared coherent bath does respond. Formation of a physical gauge charge and a change in its calibrated coupling near an independently identified intrinsic particle radius remain open. The released vacancy has a measured, age-dependent packet width and a recorded capture address. [Class-wide result and scope](PROOF.md) · [Published criterion and comparison](https://github.com/openwave-labs/openwave/blob/c531265eb6085725c6268e156804499c810ba5e4/MODELS.md#L127)

## Comparing an energy probe with an independently located source

The source starts as a sharp vacancy at port 0 or 2 on the 20- or 320-cell board, with one load at every other site. All material contacts act during release. Independent cohorts use eleven restored signed reads per owner to locate the vacancy. Its squared graph distance from the launch site gives the size measurement. A separate native dimer calibrates the supplied appointment times. In a fresh probe cohort the measured vacancy address becomes a static guard. A material standard sets its two possible signs, and a second guard is addressed uniformly over each enrolled shell. Coincident addresses return a named refusal scored as zero per original start.

A 256-position coherent clock on the guarded boundary contacts reads the four-sign energy difference. At native strength eight, the direct-contact calibration gives bin −65 and the separated calibration gives zero, corresponding to signed mean energies −16.25 and zero. Changing only the probed contrast to zero leaves a direct exchange contribution of −0.25. The bulk is free during the static hold. The clock isolates diagonal boundary terms during its interrogation. The all-open energy mean obeys the same signed identity. The source stock, sign guards, coherent clock, capture and shell addresses are supplied. [Native propagation and clocks](native.py)

| Cells | Age | Source rms radius, links | Shell | Energy from separate probe counts | Prediction from the radius cohort's address counts |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 1 | 0.94845 | 1 | −10.59989 | −10.60023 |
| 20 | 3 | 1.58962 | 1 | −7.56079 | −7.56387 |
| 20 | 3 | 1.58962 | 2 | −4.13456 | −4.13286 |
| 20 | 3 | 1.58962 | 3 | −0.80691 | −0.80585 |
| 320 | 1 | 0.94840 | 1 | −10.59988 | −10.60057 |
| 320 | 3 | 1.58916 | 1 | −7.56550 | −7.56414 |
| 320 | 3 | 1.58916 | 2 | −4.12929 | −4.13253 |
| 320 | 3 | 1.58916 | 3 | −0.80899 | −0.80856 |

The table includes every radius-window comparison. At age three the 320-cell intervals for the measured-to-predicted contact coefficient lie within [0.993,1.008] across these shells. The fourth shell, outside the window, has only 2,521 direct readings in 16,777,216 starts on that board, and its interval [0.863,1.137] leaves a 10% deviation unresolved. All four shells, the initial point source, the zero-contrast control, zero-energy readings and coincident refusals are archived. This comparison uses the full measured address distribution. An rms radius alone would not determine the shell energies.

## Does the surrounding material respond to the sign?

With the captured site guarded, a separate dilute bath holds one mobile vacancy on the remaining ordinary ports. Compare the uniformly mixed initial address with a coherently prepared uniform amplitude. Their initial sharp material-word distributions are identical. For either source sign, release all unguarded contacts, allow the guarded source's diagonal couplings to act, and count the mobile vacancy in each shell. The first-shell susceptibility is half the difference of the two sign-conditioned probabilities.

At time one with the source captured at its original port, the 320-cell incoherent bath gives +0.00000208 ± 0.00000792 and the coherent control gives −0.00024514 ± 0.00000779. Their calculated responses are zero and −0.00024538. On 20 cells the coherent response is −0.00398764 ± 0.00003089. The altered flight keeps exchange active and gives −0.00026276 ± 0.00000864 for the 320-cell coherent control. The coherent input is a prepared susceptibility control. Its change under a source-sign switch establishes that the detector responds to motion in material with the same initial occupation probabilities.

The archive also includes complete joint source-capture and polarization cohorts for every enrolled age, sign, law and appointment. In the native age-three 320-cell coherent cohort, 61,319 positive-source starts and 90,482 negative-source starts gave both a capture at the launch port and a subsequent first-shell bath vacancy. Each cohort had 536,870,912 original starts. With the positive source, 477,258,970 starts captured neither of the two enrolled source addresses and contribute to the refusal bin. The receiver therefore can recover both the all-start comparison and the conditional response at a given capture. [Acquisition and estimators](measure.py)

## Why the observed energy width runs differently

Even when the static polarization and the signed mean interaction vanish, the spectral variance can depend on the separation of two guards. On 20 cells, with half the active sites vacant, the signed variance at distance two is 251.17 ± 0.81 exchange units squared for a fixed total inventory, and 263.96 ± 0.83 for independent stocking at the same density. At distance three the respective readings are −16.14 ± 0.89 and +0.20 ± 0.91. A signed combination of four variances can be negative.

The two boundary heads read their overlapping vacancy counts on the same material word. Their independently calibrated signed energy scores give twice their covariance. A separate full 252-state microscopic calculation checks its equality with the four-sign spectral variance. For independent stocking, shared neighbours supply the entire signal. Fixing the total vacancy number adds a long-range anticorrelation even when the two neighbourhoods are disjoint. The observed scale dependence in this channel comes from fluctuations and the prepared inventory. [Exact comparison and tests](test_free_bath.py) · [Noncommuting moment identities](SignedMoments.lean)

All counts are PCG64-seeded multinomial draws from native-law probabilities. Each source, static probe and fluctuation cohort has 2²⁴ simulated starts. Each polarization cohort has 2²⁹. The reported 99% simultaneous allowances quantify Monte Carlo sampling resolution across at most 2,000 scalar means, using the bounded-score Bernstein inequality with variances from the calculated law and a finite sampling-grid allowance. Replay checks all deposited dyadic probabilities before reconstructing the seeded counts. They do not represent uncertainty in the chosen coupling, source preparation or timing model.

From the repository root, replay with one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_23_coupling.free_bath
nice -n 15 .venv/bin/python -B -O -W error -m "$p.measure"
nice -n 15 .venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling/free_bath
(cd Lean && nice -n 15 lake env lean -j1 ../model/candidate_laws/interval_echo/openwave/row_23_coupling/free_bath/SignedMoments.lean)
```
