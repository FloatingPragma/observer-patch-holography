# Can local gateway geometry produce a propagating matter band?

On the 20-cell federation, turning through the gateway farthest from entry produces the same three-port segment for **every starting port**. A detector three cells away receives a freely released load. On the 80- and 320-cell boards the farthest rule often introduces a fourth port. A second rule based on the 120° opening of a symmetric three-arm junction recovers long, repeating routes on all three boards, provided the preparer screens starts by their observed port degrees. These tests establish a geometrically selectable material sector. How the repair dynamics would select its preparation remains open.

## What the observer prepares and reads

Each cell has three gateways. A preparer queries the two possible exit gateways in the current cell, selects by local angular separation, marks the first shortest internal path in native port order, and crosses the actual external seam. The route ends before its first repeated cell. This uses access to cell adjacency, gateway identities, port order and a visited-cell record. No full cell itinerary, target mode or prescribed number of internal sites enters the growth rule. The material must still be prepared by setting one oriented even-sign arrow on each path link, empty writes off the path, and the ordinary degree references. The route-growth rule, collar writes, independent stationary reference and detector pulses are supplied operations. The simulation computes the list of marks from native adjacency and tests the resulting released state. It does not simulate spontaneous formation of the marked state.

A sharp launch occupies one addressed port. After release the occupation is read from its calibrated incident arrow-end records. Summation over ports of a cell gives that cell's detector outcome. On the longest small-board farthest route the third cell received 19,876 of 262,144 simulated launches at native exposure 4, from zero at launch (exact probability 0.07541475). The first cell after launch received 23,481. All four 99-state column checks underlying this reduction give zero residual. Direct time evolution and an independent spectral sum agree within 1e−10 before sampling. Equal-length symmetry cases use common seeds to control sampling variation.

## All starts and controls

| Board | Directed starts per tie rule | Longest true-antipode route | Longest farthest route, low / high tie | Lines inside the two trimer gaps on longest farthest route, low tie |
|---:|---:|---:|---:|---:|
| 20 | 60 | 0 | 14 / 14 | 0 / 0 |
| 80 | 240 | 2 | 36 / 38 | 4 / 10 |
| 320 | 960 | 4 | 69 / 68 | 2 / 10 |

All 60 starts on the small board produce degree word (6,5,6) under either tie rule, for 8 through 14 cells. A gateway at distance two has exactly one intervening marked site on the chosen geodesic. The other gateway cannot occupy this interior on the small board. Each larger board also has cells with two antipodal gateways. A turn through some of these cells contributes four sites. Nearest-gateway results, all route-length and segment histograms, chosen physical addresses and exact finite-operator audits appear in the [baseline receipt](results.json). Genuine straight-through propagation on an antipodal local gateway fails already at route construction on the small board.

For the exploratory 120° rule, the candidate angles in the inscribed icosahedron have cosines +1/√5, −1/√5 or −1. The preferred exit minimizes its cosine discrepancy from −1/2, including when the ideal choice is unavailable. This rule was proposed after the farthest-rule census. OPH supplies the gateway geometry. Equal three-arm tension is a hypothesis about how a preparer might use it.

| Board | Whole three-port routes, low / high tie | Longest whole three-port route, low / high tie | Longest unrestricted route, low / high tie |
|---:|---:|---:|---:|
| 20 | 60 / 60 | 14 / 14 | 14 / 14 |
| 80 | 37 / 62 | 29 / 21 | 44 / 34 |
| 320 | 82 / 117 | 53 / 27 | 84 / 66 |

Every directed entry contributes to each count. Routes can overlap. The counts measure preparation availability, not probabilities for a naturally formed state. On the longest unrestricted 320-cell route, the two putative gaps contain 4 and 4 interior levels with the low-port tie rule, or 7 and 7 with the high-port tie rule. [Exploratory census and complete address traces](junction.json)

## What the frequencies establish

On any complete three-site word, the bulk site diagonal is (−8,0,−8), the hopping is −4, and the bond between cells carries the Bloch phase. The two bulk gaps run from −12 to −8 and from −4 to −6+2√17. Repeated-cell momentum is measured in inverse cell spacings. Reflection at the prepared ends quantizes interior spatial indices at πj/M for j=1,…,M−1. Direct diagonalization on the finite physical routes verifies all 3(M−1) corresponding bulk lines. This is an operator cross-check using known cell count and degrees. On the low-tie screened routes it verifies 39, 84 and 156 distinct lines. Removing the degree references from the same 14-cell route puts 6 and 14 levels inside the two gaps.

A second experiment uses only one sharply localized central launch, an independent frequency reference on three ordinary ports in an untouched cell, and the earlier [frequency-selective arrival instrument](../row_30_route/receiver.py). Coherent holds filter the signal and the observer reads position. Each setting retains all 4,194,304 simulated trials, including rejections. On the longest screened low-port routes the detector three cells downstream recorded 2,786, 3,110 and 3,037 arrivals at −6 on successively larger boards. Every shot was rejected at −10 and −1 on these guides. A null result at one setting samples a narrow frequency window. On an unscreened 80-cell high-port route the −10 filter accepted 171 events including one at the third cell; finite eigenvalue counts supply the stronger full-interval control. [All receiver and rejection counts](filter.json)

The positive result is a particular propagating, frequency-selective excitation that can be prepared through local geometric choices. The three-site pattern differs from the earlier four-site material. These experiments do not measure a local rest frequency and independent characteristic speed for this three-site excitation, or show that it obeys Klein–Gordon dispersion in observer records. Neither a gap alone nor a deterministic route-growth service closes OW-30.

## Reproduce

Run from the RER checkout root. Each module verifies the checked-in receipt. The byte-independent comparison uses numerical values rounded to ten decimal places for operator and probability diagnostics. Integer detector counts are compared exactly.

```sh
.venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_30_local_routes.experiment
.venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_30_local_routes.junction
.venv/bin/python -m model.candidate_laws.interval_echo.openwave.row_30_local_routes.filter_read
.venv/bin/python -m pytest -q model/candidate_laws/interval_echo/openwave/row_30_local_routes
```

The [protocol](PROTOCOL.md) distinguishes the original census and the exploratory angular follow-up. The geometry test independently checks every angular separation against the coordinate icosahedron.
