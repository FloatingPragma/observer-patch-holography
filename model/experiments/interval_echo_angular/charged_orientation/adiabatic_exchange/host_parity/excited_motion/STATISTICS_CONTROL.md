# Changing the vacancy exchange rule on the same four-port pump

The sharp-source experiment starts with two vacancies on leaves 2 and 5 of a commissioned Y centered on physical port 4; port 9 is the third leaf. Both arms start with the same sharp stock, retain all six possible final addresses and meet at the original balanced detector. Introducing antisymmetric vacancy hopping reverses the once-around dark pattern for both central detunings. Changing the detuning also reverses the pattern with ordinary native hopping. A dark exit at a single detuning therefore cannot identify exchange statistics.

| Microstep | Central bias | Native dark, one turn | Antisymmetric dark, one turn | Native dark, two turns | Antisymmetric dark, two turns |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1/16 | +16 | 0.996812 | 0.003874 | 0.000293 | 0.000958 |
| 1/16 | −16 | 0.003581 | 0.997113 | 0.010144 | 0.000791 |
| 1/32 | +16 | 0.997006 | 0.003800 | 0.000894 | 0.001344 |
| 1/32 | −16 | 0.001810 | 0.998021 | 0.000043 | 0.001636 |

Each turn has duration 393,216. Its 49,152 or 98,304 appointments each use 128 split microsteps. The independently measured probabilities on the planned pair addresses at the four pinned frames are all at least 0.985899 for the changed rule. The three crossover readings retain at least 0.993348 on the two planned pair addresses. These position records use separate cohorts from the interference measurement. The numbers are exact-law probabilities evaluated numerically, with every sharp launch included. [All eight results and address records](statistics_control.json)

For a hop between ordered guide sites a and b, the changed generator multiplies the native off-diagonal matrix element by (−1) to the power of the number of other vacancies strictly between the two sites. The original diagonal contrast, private-bias pulses, pair-projector apertures, clock, arm paths and detector are unchanged. The phase convention for each sharp vacancy-pair basis vector can change where the negative link appears. The product around its six-state configuration cycle reverses for every such convention. A one-particle gauge phase on a spoke of the physical tree cancels on the exchange-versus-retracing loop. For one vacancy the spectator count is zero and each pair aperture vanishes. The two arms then have identical propagators. On the original sharp six-contact direct exchange, keeping the comparison arm's retracing order gives dark probabilities zero with native hopping and one with antisymmetric hopping. In the slow two-vacancy pump the programmed band geometry contributes its own sign. The introduced antisymmetry describes a comparison law. Spin, a response to physical rotation and an identification with electrons require separate measurements. [Sharp-source experiment](README.md) · [Native and comparison matrix elements](statistics_control.py)

The replay builds each split microstep at the original aperture angle, raises it to its actual integer service count and multiplies in temporal order. The reference arm retraces its first half. Transposition reverses the symmetric split-pulse sequence while retaining its dynamical phases. Independent tests project two labeled-particle generators onto symmetric and antisymmetric amplitudes, explicitly multiply both arm schedules and match the four native outputs to the original deposited acquisition. The changed-statistics run makes no phase-selected source or new clock-bin choice. [Independent comparisons](test_statistics_control.py)

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
p=model.experiments.interval_echo_angular.charged_orientation.adiabatic_exchange.host_parity.excited_motion
python -m "$p.statistics_control" --output /tmp/y-statistics.json
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/adiabatic_exchange/host_parity/excited_motion
```
