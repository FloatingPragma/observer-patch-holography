# Reading the phase of a split load at two stations

A source splits one unit between Alice and Bob and applies one of two opposite relative phases. Alice measures her own rails. Bob receives her recorded bits, selects his local analyzer and reports which phase he thinks the source used. They keep every simulated start. Each additional reference in the table is one independently split material unit delivered to the same two observers.

| Receiver setting | Correct on phase 0 | Correct on phase π | Correct overall | Predicted |
| --- | ---: | ---: | ---: | ---: |
| No reference | 65,563 | 65,934 | 50.162% | 50% |
| One split unit | 98,468 | 98,470 | 75.126% | 75% |
| Two split units | 111,651 | 111,945 | 85.295% | 85.355% |
| Two, with the signal path read early | 65,373 | 65,425 | 49.895% | 50% |
| Two, with exchange alone throughout | 81,960 | 49,345 | 50.089% | 50% |

Each phase column covers 131,072 starts per setting. The simultaneous 99% half-width is 0.604 percentage points for individual bank accuracies. Alice's words alone have equal predicted frequencies for either sign. Reading the signal path early erases the distinction. At zero contrast the phase writer also loses its action. In that control the receiver favors one guess, as predicted by the separate sign accuracies of 5/8 and 3/8. [Frozen protocol](PROTOCOL.md) · [Complete simulated records](records.json)

For a fixed independent reference of total load L, number-preserving local instruments and classical communication can distinguish the signs with success at most (1 + cos(π/(L+2)))/2. The native preparations and receivers reach the bound in probability for L=0,1,2. The proof covers arbitrary reference multiplicities. Its one- and two-unit inequalities also have a Lean check. This gives an operational extension of the fixed-inventory results for PR-03, whose general additivity premise still covers other effects, written records and free histories. The simulated draws use the registered quantum instrument. [Bound and attainment](PROOF.md) · [Lean inequalities](Bound.lean)

From the repository root, replay the terminal count archive and run the independent regression checks with a single numerical thread.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_consensus.inventory_effects.shared_reference.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_consensus/inventory_effects/shared_reference
```
