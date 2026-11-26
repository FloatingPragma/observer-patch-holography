# Two ordinary vacancies exchange in an excited band

Two identical vacancies move through three successive pairs of leaf ports on an unthreaded Y. Every other physical load is one. A balanced detector reads a negative relative phase in the slow, excited pair band. On **4,180,732 of 4,194,304 unfiltered simulated launches**, it returned dark for one exchange. The predicted dark fraction falls to **0.0293% after two exchanges**. A native direct exchange of the same two vacancies on the same four ports returns bright on every sharp launch. [All launch, position and exit records](records.json.gz) · [Derivation](PROOF.md)

| 393,216 native flight units per turn | Sharp dark, once | Sharp dark, twice | Least sharp probability at a planned pair address |
| --- | ---: | ---: | ---: |
| Positive central bias, excited pair band | 99.6812% | 0.0293% | 99.3264% |
| Negative central bias, lowest band | 0.3581% | 1.0144% | 99.7694% |
| Positive, finer clock | 99.7006% | 0.0894% | 98.8260% |
| Negative, finer clock | 0.1810% | 0.0043% | 99.5787% |

The probabilities use 49,152 aperture appointments with increment 1/16, and independently repeat them at 98,304 with increment 1/32. A source starts with two sharp vacancies on leaves A and B. On one arm the preferred pair moves from AB to BC to CA and back to AB. The other arm advances halfway and retraces. Reflection matches their instantaneous line energies and aperture service counts. Separate cohorts read both vacancies at thirteen appointments. At the three crossovers the sharp pair occupies one of its two allowed adjacent configurations with probability at least 99.23% in the main excited run. The phase-selected source measures a finite clock before the exchange, with every other bin counted as a refusal. It gives 99.9525% dark among accepted starts at the coarse setting. The raw sharp result above uses every launch. A single vacancy on this material sees a constant pair-aperture contribution and identical evolution on both arms. [Pulse and clock model](pump.py) · [Full detector acquisition](measure.py)

The original ±8 experiment gave a negative interferometer reading, while its sharp stock missed the 95% position requirement. Doubling the private bias was a deliberate change to improve localization, frozen before the held runs. The released records include both strengths, all three original running times, both follow-up times and both clock increments. The sign reverses on changing the static detuning with the same moving aperture. With +16, freezing the pair aperture while native swaps follow the same six-position exchange route gives zero sharp dark probability and 0.127% on the admitted incoming source. It is a faster, nonadiabatic control. [Original protocol](PROTOCOL.md) · [Comparison settings](CONTROLS.md) · [Held sharper localization](SHARP_PROTOCOL.md) · [Fixed-angle exchange](ROUTE_PROTOCOL.md)

### Identifying the transported cores

In a separate held run, one vacancy stays on its physical leaf throughout each third of the turn. Every routing operation, private exposure and active spoke on that leg can be executed without touching the parked port's load. During AB→BC the B vacancy stays put and the A vacancy goes through the centre to C. Parking C and then A on the next two legs interchanges the final assignments. **4,154,675 of 4,194,304** sharp launches passed the sequential address checks at +16. Independent copies bypassed the checks and gave dark probability **99.5916%**. At −16 the same checks passed on **4,158,347** launches with dark probability **0.5856%** on the interference copies. This provides a trajectory witness on each sign setting even if individual core labels are unobservable. [Continuous parking protocol](PARKED_PROTOCOL.md)

The prepared physical stock has ordinary unit hosts. A bare host permutation cannot supply the negative phase. The symmetric native clock instead gives the excited real line different reflection parities at the start and halfway through the drive. Their product is −1. The lowest comparison line gives +1. These signs characterize their gapped material processes. The pump uses supplied positive private exposures, timed native contact gates, coherent arm switching and a slowly rotated pair trap. The vacancy exchange is an operational finite success. Portable electron statistics and an intrinsic response to a physical 2π rotation remain unresolved.

### Reproduce

From the RER root, set numerical libraries to one thread. Replaying the full schedule takes longer than testing the checked-in data.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1
p=model.experiments.interval_echo_angular.charged_orientation.adiabatic_exchange.host_parity.excited_motion
python3 -B -m "$p.measure"
python3 -B -m pytest -q model/experiments/interval_echo_angular/charged_orientation/adiabatic_exchange/host_parity/excited_motion
```

Each counted cohort has 4,194,304 original starts. Terminal data order the clock refusal, bright and dark exits. Frame data order the refusal and six vacancy words. The continuously parked path records the source refusal, the first failed positional check and the final success. The published counts are simulated draws from the stated finite laws.
