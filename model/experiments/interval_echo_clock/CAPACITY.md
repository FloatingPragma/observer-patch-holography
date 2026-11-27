# Clock accuracy at long handoff intervals

On the reciprocal unit-link fixture, a five-bit interval field suffices to bring the **expected** traveller/home rate within 0.0005 of the optical target for every positive integer handoff interval. At b bits the first overflow alone gives error 1−√(1−4⁻ᵇ). This is 0.00195504 for four bits and 0.000488401 for five. The finite-state calculation checks the representable intervals, and monotonicity bounds all later overflows. The five-bit margin is a design tolerance. The moving observer occupies 14 cells at rest and 15 across a transfer. The nine-bit elapsed word, prepared hands, guide, readout windows and local primary source have their original roles. [Frozen tolerance and experiment](CAPACITY_PROTOCOL.md) · [Measurements](capacity.json)

| Bits in each of age, length and address | First overflow | Intervals outside a 0.0005 rate tolerance | Largest absolute rate error |
| ---: | ---: | --- | ---: |
| 3 | 8 | 8–31 | 0.0078432584 |
| 4 | 16 | 16–31 | 0.0019550361 |
| 5 | 32 | none | 0.0004884005 |

## Why the walk approximates the optical rate

Let q count the service opportunities between two confirmed handoffs. On a closing positive the patch takes one compulsory signed step. With unlimited rejection sampling, a further step occurs with probability 1/q. For a symmetric walk the generating function of the first-return probabilities is F(z) = 1 − √(1−z²). This follows from the return generating function U(z) = Σ binom(2n,n) z²ⁿ / 4ⁿ = 1/√(1−z²) and the renewal identity U = 1 + FU. The compulsory first step makes the ideal conditional suppression probability qF(1/q). There is one closing opportunity in q. With the same fair primary rates for home and traveller, the ideal ratio of the expected elapsed counts is therefore 1 − F(1/q) = √(1−1/q²). The independent optical interpretation takes one reciprocal relay hop per q opportunities. The closing probabilities depend only on the patch’s interval field and native fair draws. The physical speed, unit length and equality of opportunities supplied by this fixture set the domain of the comparison. [Local clock program](../../candidate_laws/interval_echo/clock.py) · [Optical protocol](../../evaluation/ONE_LINK_CLOCK.md)

The actual instrument caps the walk at twelve steps and the address selection at twelve trials. For each representable q, the certificate propagates exact rational mass on the signed positions and removes first returns. It compares the resulting rate to the target by squaring rational interval bounds. On overflow a positive closing mark is accepted, giving mean rate one. The optical target increases monotonically toward one, so the maximum tail error is at the first overflow. For any prescribed completed itinerary the expected elapsed ratio is the duration-weighted average of its interval rates. The same 0.0005 bound applies under those fixture assumptions. Individual histories fluctuate. A material reunion also requires both elapsed words to remain within their nine-bit capacities. [Certificate and negative controls](capacity.py)

The material replay traverses intervals 8, 16, 31, 32, 2 and 1 on the 20-, 80- and 320-cell screens at two fixed seeds. It compares both endpoint words with the mark histories, checks every age reset and exercises the first-return control at both sides of each width boundary. The five-bit register represents 31 exactly. At 32 it intentionally saturates, posts a positive mark and resets locally on a handoff. The tolerance bounds the timing consequence of saturation. It does not extend the range of the physical age register.

Run a five-bit material journey from the repository root with the first command, or reproduce the all-interval certificate and material controls with the second.

```sh
python -m model run --law interval_echo --clock-trip --level 0 8 16
python -m model.experiments.interval_echo_clock.capacity
```
