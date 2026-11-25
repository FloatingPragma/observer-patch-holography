# A small positive formation gap that survives repetition

On the 144-port material route, a charged vacancy and its conjugate support the same lowest band. With one written barrier every 48 ports its infinite periodic minimum is **0.003786395** above matched filling. The actual 144-port charge-and-head component has minimum **0.003790409**. A finite clock prepares a packet with 99.9762% of its accepted weight in the first band. A detector one period away records **2** arrivals at release and **40,504** after 1,048,576 native time units, in 1,048,576 simulated starts per bank. The periodic mark spacing sets the material energy scale. None of the eight spacings reproduces both charged-lepton mass ratios. [Full records](results.json) · [Native sector and apparatus](measure.py) · [Exact spectral brackets](certify.py) · [Independent checks](test_periodic.py)

## Physical preparation and measurement

The [144-port route](../../marked_bands/material.py) has 48 eligible side contacts. Leave every second backbone contact blank and prepare each intervening one with an even magnitude-one head. Erase every ordinary private reference on the route. Write an even magnitude-two head on every selected side contact, then silence all remaining incident contacts. A compensator occupies an isolated dimer. The filled environment outside the route and its isolated stock factor from this sector.

The receipt identifies the side arms, neutral stock pairs and phase reference. The common bypass obtained by removing all 72 writing contacts on the backbone and all 48 eligible side contacts is biconnected. This supplies native paths for stock transpositions and for phase control while avoiding one specified endpoint. A transposition along a blank path and its reverse restores every intermediate load. Each side deposit spends one vacancy and one excess and leaves both endpoints filled. Each backbone deposit temporarily borrows one vacancy. A coherent write and load herald select either direction of the even head. A full native directed pulse corrects the other herald, and the vacancy returns to stock. The private zero writes follow the deposits. For a shorter transmitting component, reverse the direction of the active head at the stated right boundary. The closure routine checks that the ring and compensator factor. The writable magnitude-two head adds 32 to a charged port while leaving its load in place. [Stock and full native side-write column](../../marked_bands/stock.py) · [Coherent backbone write](../../../../../../experiments/interval_echo_angular/charged_orientation/winding_memory/measure.py)

The preparation assumes addressed contacts, a private phase rail, a coherent reference clock and the declared initial balanced stocks. Exchanging vacancy and excess reverses the required active direction. On the full component both uniform direction patterns occur at the launch end after a winding. The native column calculation for both charge signs returns the same relative generator and the same filled reference energy.

## Why a periodic estimate applies to a physical open component

A blank contact contributes diagonal 8.5 and hopping −0.5. An enabled magnitude-one contact contributes diagonal −4 and hopping −4, flipping its direction whenever the charge crosses. A blocked endpoint contributes zero from that contact. The interior therefore forms a weighted positive graph Laplacian with alternating bond weights 4 and 1/2. Each selected mark adds 32 on one site. A full winding reverses all 72 heads. Following all nonzero columns from the source gives positions −1 through 144, with their head histories, on a **146-state open path**. For a shorter length $L$ the corresponding path has $L+2$ states. Its interior agrees with the infinite periodic operator. Its two endpoint diagonals exceed the principal truncation of that operator by four each. Extending a finite vector by zeros bounds its energy below by the periodic minimum.

On an infinite repetition at least one positive onsite mark per period precludes a zero-energy mode. A zero of the quadratic form would require a constant modulus on all links and zero amplitude at every mark. The lowest periodic eigenvector has one sign. The scalar Jacobi recurrence has a real two-by-two transfer matrix of determinant one over a period. Its trace at a Bloch eigenvalue equals twice the cosine of the Bloch angle. Thus successive bands are intervals between ordered periodic and antiperiodic eigenvalues. The package brackets every band edge by exact integer inertia counts at rational endpoints. The leading path determinants obey the tridiagonal recurrence. Adding the closing bond introduces its two directed full cycles and the contribution of the endpoint transposition. Sign changes of the nonzero leading determinants count eigenvalues below the bracket endpoint. The same test brackets each true finite minimum. [Integer certificates](certify.py)

## All uniform spacings on the declared route

The spacing $s$ counts eligible middle ports. The combined period is the least common multiple of 2 and $3s$. Include every divisor of 48 with at least three combined periods on the physical route. The short length is the largest integral number of periods that fits within half the route. Compare that length with twice its value and the full route at fixed spacing. For the family test, compare all ordered triples of distinct band bottoms on one fixed material with the laboratory ratios 206.7682827 and 1776.93/0.51099895069. The last column minimizes the larger of the two absolute fractional errors, allowing any band to serve as light member. The finite spectra and separate comparisons anchored strictly to the lowest band are in the receipt.

| Side spacing $s$ | Physical period | Infinite minimum | Short minimum | Doubled minimum | Full minimum | Best joint miss |
|---:|---:|---:|---:|---:|---:|---:|
| No side mark | 2 | 0 | 0.001601401 | 0.000411523 | 0.000411523 | — |
| 1 | 6 | 0.492050086 | 0.492053765 | 0.492051074 | 0.492051074 | 97.84% |
| 2 | 6 | 0.221266138 | 0.221376012 | 0.221294078 | 0.221294078 | 95.19% |
| 3 | 18 | 0.088387633 | 0.088391704 | 0.088388769 | 0.088388769 | 87.96% |
| 4 | 12 | 0.058899930 | 0.058965731 | 0.058917052 | 0.058917052 | 81.93% |
| 6 | 18 | 0.026561653 | 0.026605348 | 0.026573357 | 0.026573357 | 59.94% |
| 8 | 24 | 0.015031815 | 0.015063357 | 0.015040556 | 0.015040556 | 29.21% |
| 12 | 36 | 0.006715885 | 0.006734626 | 0.006721512 | 0.006721512 | 58.43% |
| 16 | 48 | 0.003786395 | 0.003817151 | 0.003794340 | 0.003790409 | 29.54% |

Spacing 16 compares 48, 96 and 144 backbone ports. All other settings compare 72 and 144. The infinite minimum is positive and below 0.05694 at spacings 6, 8, 12 and 16. The unmarked control has a positive finite energy and zero infinite minimum. At fixed spacing 8 the best simultaneous band-bottom ratios are 1, 267.125962 and 2461.461419. The maximum relative miss is 29.215%. Even admitting any ordered pair of excited levels with the actual 144-port ground as light member in each preparation gives at least 29.25% maximum error. These systems contain a full orbital spectrum. The count of bands follows the written period.

## Does the low band carry a charge between periods?

For each setting, prepare a sharp charge inside the first period by addressed swaps. Count signed arrivals anywhere in the next complete period, tracing over every head history. The unfiltered bank includes higher energies. An independent bank uses a reference in an equal superposition of $M$ clock delays spaced by 1/8, then projects the reference back onto its starting state after compensation at the first-band midpoint. This finite filter has complex response equal to the mean of the $M$ compensated phases. Choose the smallest power of two $M$ at least 32 divided by the product of the delay step and the gap from the top of the first band to the bottom of the second. All settings use the same rule and wait ladder. Finite states below the midpoint of that gap define the first-band weight reported for the open component. The filter can spread a packet before release, so the zero-delay bank uses the same filter.

At spacing 16 the clock has 32,768 delays. The lowest periodic band spans 0.003786395 to 0.003803892. The filter admits 4.1628% of starts, with 99.9762% first-band weight. Immediately after preparation, 2 starts arrive in the neighboring period, 43,191 arrive elsewhere, and 1,005,383 are refused. At delay 1,048,576 those counts are 40,504, 3,103 and 1,004,969. Independent repeated squaring of the controlled binary clock gates and the full free propagator reproduces the joint arrival probability. The delay unit, mark separation, attached reference and spatial detector are material and apparatus inputs. A physical rest mass requires the dispersion and independent unit conversion of the same identified species.

From the repository root, run with one numerical thread and use an unused output path.

```sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -m model.candidate_laws.interval_echo.openwave.row_11_spectrum.light_charge_gap.periodic_marks.measure --output /tmp/oph-periodic-charge.json
python -m pytest -q model/candidate_laws/interval_echo/openwave/row_11_spectrum/light_charge_gap/periodic_marks/
```
