# Two vacancies follow moving wells through a relative exchange

Two identical vacancies hop on an eight-port circuit while private-reference exposures move two wells. At the slowest enrolled rate on the native guide, **4,156,124 of 4,194,304** plus-arm starts end with both vacancies at the two intended destinations. The two arms exchange the destination assignments. Interference gives **27,793 dark outcomes in 4,194,304 unfiltered starts**, consistent with the ordinary positive sign. The negative-statistics comparison gives **4,184,719 dark outcomes** with the same pulse calendar and energy costs. [Frozen protocol](PROTOCOL.md) · [All address and interferometer records](results.json.gz)

| Guide law | Pair probability at destination | Pair dark probability | Single-winding dark probability |
| --- | ---: | ---: | ---: |
| Native | 0.990850 | 0.006572 | 0.009528 |
| Antisymmetric exchange | 0.990410 | 0.997674 | 0.009528 |
| Ordinary exchange, one-body $\pi$ flux | 0.990410 | 0.997674 | 0.998655 |
| Zero contrast | 0.003412 | 0 | 0 |
| Parked exchange | 0 | 0 | 0 |

The table gives calculated finite-law probabilities at duration 4,096 per moving step and split interval 1/4. Every recorded count is a PCG64 multinomial sample of the calculated port law. Each bank has 4,194,304 starts; a conservative simultaneous 99% radius for the reported probability marginals is 0.001523, including quantization of the sampling grid. In the native pair run the two final exit probabilities sum over all 28 vacancy positions. The interpolated dwell supplies a pair of moving wells. Closing exchange leaves the holes at their launch sites. Suppressing contrast makes both coherent arms propagate identically. A single vacancy then happens to reach the final port with probability 0.96799 at the longest appointment, but its probability at an earlier scheduled destination is below 0.000011. A terminal coincidence alone cannot certify the path.

## How the comparison works

The guide is the first eight-port cycle in the first owner for which all eight ports have a native private reference. It is embedded in the 20-owner board. Swap two sharp stock vacancies onto opposite guide ports of the unit-filled host. In each moving step the contact dwell shifts from one reference to its neighbor while a second deeper well remains available to the parked vacancy. Every guide bond opens together between successive half private pulses. Within the guide the native exchange and contrast act on all one- or two-vacancy address words, with no directional shutter that carries a hole to its target. The stationary-marker positions, shutter calendar, coherent direction arm and final balanced analyzer are supplied. The native Hamiltonian on a blank bond, the sharp terminal echo, and the separate calibration of the timing grid are verified against the 99-state seam. [Pulse sequence and arbitrary-position generator](exchange.py)

The two arms move from wells at indices 0 and 4 to wells at 2 and 6 by reflection-related routes. The source sees a single moving minimum on each leg and a parked partner well. Independent cohorts read the whole address word at every midpoint and endpoint, including all other addresses. On the native slow pair run the two vacancies follow the enrolled regions with calculated probability at least 0.99084 at all nine observed frames. In a separate cohort recombination proceeds with no midflight reads. Timing each step at 256 or 1,024 instead gives native final pair probabilities 0.03214 and 0.42026; the corresponding dark probabilities are 0.57132 and 0.33583. The two faster schedules fail the prospectively required 95% final localization. At 4,096 the half-sized split interval 1/8 gives final pair probability 0.99276 and dark probability 0.00399. The incoherent control gives half on each exit. A free blank bond calibrated the two interval settings as 0.25046 and 0.12444 exchange units.

The two alternative pair laws coincide exactly on a no-passing cycle. Reordering two vacancies gives the same boundary hopping sign as a one-body $\pi$ flux. The separate one-vacancy winding returns to the bright exit for the statistical-sign control and to the dark exit for the flux control. The two pair histories hit the same sharp destination word. Reflection fixes their interference sign on that word for *any* nonzero return amplitude. The held all-address probabilities and the independent material positions establish transport. [Complete run and readers](measure.py) · [Configuration-space proof](PROOF.md)

## Which native sources can acquire an adiabatic minus sign?

For arbitrary ternary loads and arbitrary written seam sectors the native Hamiltonians share one real basis in which every allowed nonnegative contact term has nonpositive off-diagonal entries. The geometric parallel transport of a separated lowest band along a closed parameter path is consequently nonnegative. In the basis of positive ground states on connected components it can permute those states. A simple ground has geometric phase +1, and an entire nonzero lowest band cannot acquire a central minus sign. Incoming states antisymmetric under a component permutation, excited spaces, finite-speed loops, and sign-changing interactions have distinct hypotheses. This places a limit on slowly exchanging ground-state cores across the full native positive-contact family. The measured source is two mobile vacancies in a protected guide. Its identification with free electrons and the full OW-07 minus-sign criterion remain **Open**. [Local integer certificate and finite-dimensional proof](PROOF.md)

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.experiments.interval_echo_angular.charged_orientation.adiabatic_exchange
nice -n 15 python -B -O -W error -m "$p.measure"
nice -n 15 python -B -W error -m pytest -q model/experiments/interval_echo_angular/charged_orientation/adiabatic_exchange
```
