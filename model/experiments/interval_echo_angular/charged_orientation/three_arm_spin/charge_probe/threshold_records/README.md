# A signed arrival and its enclosing boundary

The lowest mobile winding branch gives a negative material reading at port 160 while the complete boundary of its receiving owner has the same instantaneous readout as idle stock. The two observations can be made on the same arrival. Three disjoint observers along the guide also distinguish this single transported unit from independent deposits with the same individual detector statistics. OpenWave lists row 1 as untested. Its electrical quantization criterion remains on **HOLD**. [Official row](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L435) · [Emitter and threshold branch](../../recoil_search/THRESHOLD.md) · [Complete records](results.json)

## What the observers record

A terminal observer reads the written signed contrast $d$ and the sum $s$ of the two output loads after a native echo on a blank seam. The complete nine-input calibration gives $d=a-b$ and $s=a+b$ for sharp input loads $(a,b)$ in ascending address order. These outcomes identify the pre-echo pair through $a=(s+d)/2$ and $b=(s-d)/2$. The same detector gives $(s-2,d)=(-2,0)$ for two local vacancies and $(2,0)$ for two excess units. The coin in the output is immaterial to both readouts. A reference pair $(1,1)$ gives $(0,0)$. [Calibration and measurement code](measure.py)

One observer at owner 13 reads its interior spoke 158–160. A second reads the three enclosing inter-owner collars and six gateway loads. For the full boundary comparison, an active blank seam has its blank flag and its two reference loads. A held zero-write seam also allows the two loads and prepared zero-contrast collar to be checked by undoing the write. The raw local operators preserve these reference words. Negative and positive gateway-stock controls give different load and contrast readings on every boundary seam. The boundary observer compares the local readings with independently specified idle stock. Both observers retain the source and clock refusals in the original-launch bank. [Protocol](PROTOCOL.md) · [Native complement and contact calibration](../README.md)

The following counts are seeded draws from the finite-model outcome law. The printed upper endpoints show the sensitivity of a sampled boundary test, using a simultaneous 99% interval across the enrolled banks. The model equality comes from the exact support argument below.

| Physical steps | Negative tip arrivals per $2^{46}$ original starts | Boundary changes on those arrivals | Upper 99% change probability conditional on arrival |
| --- | ---: | ---: | ---: |
| 19 | 1,838 | 0 | 0.514% |
| 20 | 1,321 | 0 | 0.714% |
| 21 | 1,030 | 0 | 0.914% |

The computed signed tip rates are $2.551$, $1.908$ and $1.445$ in units of $10^{-11}$ per original start. The positive conjugate gives the same computed arrival rates and the same conditional idle boundary law. The receiver identifies each tip arrival from the two load records on the spoke.

## Finite implication for boundary comparisons

The native flight has one moving unit and a parked compensator. Every active magnitude-one seam is internal to a terminal owner. Exterior collars are either blank exchanges or held zero writes. The receiver tip is an interior port. The other observer's boundary support shares no port with the signed spoke or the parked donor.

The record $(d,s)$ identifying the marked charge at the receiver tip restricts the moving unit to that tip. Every microscopic history producing that record has unit host loads at all six boundary endpoints and the same three boundary collars as idle stock. In a state calculation, the normalized reduced boundary density operator conditioned on the tip record therefore equals the idle density operator. The same equality holds for arbitrary superpositions of the retained winding heads. Every measurement acting only on those instantaneous boundary registers has **exactly the idle outcome law**, including a coherent comparison of all three collars. The implication holds on the three deposited lengths whenever the specified signed tip event has positive probability under the uninterrupted native flight. The first-line reconstruction checks the geometry and reproduces the deposited tip probabilities. [Native flight reduction](../../apparatus.py) · [Inherited spectral preparation](../../recoil_search/threshold.py)

An instantaneous Gauss reading carried by this boundary would have to distinguish the enclosed signed arrival from idle stock after calibration. These measurements and the finite implication rule out that interpretation **on this aperture**. A boundary constructed from time-extended transport records or a field law with independent propagating boundary records requires its own attachment to the same arrival.

## Three observers on the arms

Three observers read the two ends of three disjoint blank contacts from the same middle cross section. The negative unit gives one signed deposit at most. A regional mean below one also arises if each observer receives its own independent mixture of occupied and idle inputs. Calibrating those mixtures to the actual one-observer laws preserves every local outcome probability. The joint distribution gives the difference.

| Physical steps | Multiple deposits per admitted clock read | Independent deposits, same local marginals | Total variation of the two joint laws |
| --- | ---: | ---: | ---: |
| 19 | 0 | 1.2772% | 2.5819% |
| 20 | 0 | 1.1445% | 2.3123% |
| 21 | 0 | 1.0548% | 2.1301% |

The first extended line fixes the table. The same full joint comparison at the three inherited two-tone appointments tracks interferometric transport with one occupied contact at most. A contrast-four flight on the unmodified 19-step geometry changes its location records at half beat by $7.20\times10^{-10}$ in total variation per original start. It preserves the exclusivity of the one-vacancy sector. The correlations identify a shared material deposit under the supplied preparation. Dynamical selection of an elementary electric unit requires physically distinguishable competing preparations. [Archived joint probabilities and refusal counts](results.json)

## Replay

Use the repository environment from the OPH root. The first command reconstructs the outcome law and performs the seeded draws. Use `--output` with a fresh path to retain a replay. The second command checks the local echo calibration, archived detector marginals, rejection event and the replay against the deposited probability laws.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
pkg=model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.charge_probe.threshold_records
.venv/bin/python -B -W error -m "$pkg.measure"
.venv/bin/python -B -m pytest -q --pyargs "$pkg.test_records"
```
