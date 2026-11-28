# Enclosing flux and signed transfer on the written four-face gateway

IT'S NEVER ABOUT THE STATE, IT'S ALWAYS ABOUT WHAT OBSERVERS MEASURE, COMPARE AND RECORD AS SHARED AGREED-ON CONSENSUS FACTS.

On the archived 14→24 four-face gateway, the receiver's load counter and the boundary-head observer agree on the transferred signed material unit. The exact half-contact control gives both observers the same marginal records as native transfer and yields disagreement on half the joint records. The same boundary observer can keep its signed read while the destination observer checks for the archived local doublet. These are simulated records on the written-ring extension and the charge is measured in its host's load units. [Design](DESIGN.md) · [Receipt](results.json) · [Original binding analysis](../binding_audit/README.md)

## Which records the identity ties together

For a region $R$, let $Q_R$ be the sum of the material observers' load changes relative to one occupant per port. Let $F_R$ count written heads that point into $R$ across its boundary. Take a reference copy of the same preparation before release and record $Q_R^{(0)}$ and $F_R^{(0)}$. Heads on seams internal to $R$ contribute a constant of one apiece to its total inward stock. On an allowed material step the signed load increment at each port equals the increment in the number of heads naming that port. An allowed ring step changes the head pattern and keeps those port counts fixed. Consequently the archived generator preserves the records

$$Q_R-Q_R^{(0)}=F_R-F_R^{(0)}.$$

The executable independently reconstructs the entire archived generator: for each sign it checks 352 directed material moves, including 32 directed bridge moves, and 608 directed ring moves among the 160 reachable records. All 4,096 cuts of the twelve active ports satisfy the identity exactly. The receiver checks the whole destination fan, then its gateway and its remaining five ports on separate fresh preparations when the inner seam reads can disturb coherent detection. The changes on the two disjoint parts add to the change on the whole fan.

The whole fixture contains 19 active heads and one supplied signed load. Its outer cut never changes. A putative neutral reference with the same local background would require 18 heads for a positive launch and 20 for a vacancy. The sealed nineteen-head fixture cannot provide that reference. The result tracks a transfer between observers. The original one-unit reservoir preparation and any exterior compensator are outside the measured record.

## Same marginals and different joint records

The calibrated rate-4 bridge contact applied alone for exposure $\pi/16$ gives the following exact outcomes for either archived input. For the other sign replace $q$ by $-q$.

| Receiver material read | Receiver head change | Native coupled contact | Independent half-pulses |
| ---: | ---: | ---: | ---: |
| 0 | 0 | ½ | ¼ |
| 0 | $q$ | 0 | ¼ |
| $q$ | 0 | 0 | ¼ |
| $q$ | $q$ | ½ | ¼ |

Each individual observer records a distribution with mean $q/2$ in both tests. Coincidence depends on the coupling. The changed address-local law keeps the same binary head, receiver, material unit, and pulse strength, allowing the two local flips to occur separately. The probabilities use exact rational arithmetic from the contact's two-rail pulse. Readout support includes both mismatch outcomes. A sharp load of $q/2$ requires an extended carrier and would be a different test.

## Original free release

The joint probabilities below include every launch of each declared input. Conjugate charges give the same probabilities within $3\times10^{-12}$ on replay. The exposure values identify the original simulation appointments. They have no conversion to a laboratory clock.

| Exposure | Signed load and flux in the right fan, inputs 1 and 2 | Whole receiver flux, gateway load, and receiver doublet, inputs 1 and 2 |
| ---: | ---: | ---: |
| 27.0625 | 53.6001%, 69.7690% | 10.6190%, 29.6613% |
| 30.8125 | 66.0695%, 47.5632% | 15.4251%, 15.0071% |

Each destination-doublet acceptance lies entirely within a single outer-flux result, so it can be compared jointly with that enclosure. Resolving every interior head on the same copy can disturb the coherent doublet. The parked contact gives zero probability of the right-fan transfer for either input at every free duration.

## Simultaneous suffix cuts on longer guides

The two-, three-, and four-fan guides attach to actual repeated contacts on the same 20-cell board. The independent replay carries each bridge head explicitly as part of the local step record. Starting from the original input columns and reference orientations, it recovers all 160, 800, and 3,840 words per sign, verifies every generator entry, and checks exact enclosed load and head increments. The signed increment at each port determines the identity for every subset. In the four-fan guide, the three joint bridge records comprise only `000`, `100`, `110`, and `111`, with `1` meaning that the recorded head points to the opposite endpoint from the source reference. The number of flipped bridge heads equals the net displacement in inter-fan steps from the source. Each suffix reads the signed unit exactly when the charge is inside it. [Path model](../fan_path/fixture.py) · [Path design](PATH_DESIGN.md) · [Path receipt](path_results.json)

At the preselected exposure 30.8125, observers on the final fan record the following joint probabilities for terminal load and the full bridge-change pattern. The two signs agree within numerical tolerance. The receipt includes each alternative occupied fan, with probabilities on all original starts.

| Fans | Physical bridges toward the receiver | Joint terminal probability, inputs 1 and 2 |
| ---: | --- | ---: |
| 2 | 14→24 | 66.0695%, 47.5632% |
| 3 | 14→24, 34→116 | 10.7175%, 15.3182% |
| 4 | 68→22, 14→24, 34→116 | 8.41762%, 8.41762% |

A unit of inward head change is inherited from binary heads and the declared one-load preparation. The experiment establishes when observers can compare this transfer to a flux record. It leaves the law selection of an elementary unit, the comparison with an independently identified electromagnetic field, and any measurement in physical electric units open. OpenWave row 1 stays on HOLD. [Official criterion](https://github.com/openwave-labs/openwave/blob/c9b5707ac03285eb4c885115245ac2cd5e90273d/MODELS.md#L101)

From the repository root, run

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.measure
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -m model.candidate_laws.interval_echo.openwave.electron_field.field_charge_records.paths
~~~
