# Can fixed detection losses explain the three winding records?

The two observers recorded which exit the source vacancy reached when the remote pair arrived. At owner 17 and reference setting +, every mixture of the three winding recipes gave the away-from-hub exit on at most 86.31% of intercepted pairs at wait 320, and on at least 91.97% at wait 640 (simultaneous 99% intervals). These disjoint ranges exclude a complete orthogonal triple travelling under common three-state unitary evolution and read by the same terminal instrument. The detector in that hypothesis may have arbitrary fixed unequal efficiencies and arbitrary missed outcomes. [Frozen test](PROTOCOL.md) · [Original-start count banks](counts.json) · [Decoder](measure.py)

| Source winding | Away at 320 / intercepted | Away at 640 / intercepted |
| --- | ---: | ---: |
| Zero | 23,092 / 27,070 | 29,431 / 30,050 |
| Positive | 20,995 / 24,738 | 23,244 / 25,054 |
| Negative | 20,903 / 24,942 | 24,110 / 25,276 |

Each cell is a separate cohort of 2⁵² simulated original launches. All seven trial outcomes are kept. The native pair yield across the 48 native cohorts ranged from 5.09 to 18.51 per trillion original starts. The comparison takes the recorded exit directly, before assigning a recipe to an event.

## Why detection efficiency need not be known

Write the three orthogonal prepared inputs as projectors Pⱼ that sum to the identity on the travelling space. At each time let U be any common unitary. For an arbitrary accepted exit effect E, summing its unconditional probability over the three inputs gives tr(E Σⱼ U Pⱼ U†) = tr(E). Let A be the sum of all accepted effects. In the equal delivered mixture the fraction of interceptions at exit E is therefore tr(E)/tr(A) whenever acceptance is positive. Common scalar transmission at an appointment cancels. A common random mixture of unitaries gives the same identity. The calculation also permits any three density operators summing to the identity.

At each time that fraction is a convex mixture of the measured single-recipe fractions, weighted by their intercept rates. Unknown positive production efficiencies do not affect a recipe's conditional fraction. Allowing *any* mixture at each time gives the hypothesis every choice of weights. Two disjoint ranges for the same exit exclude it. The finite proof checks the operator identity, the bounds on a weighted score, and separation for independent weights. [Lean certificate](StationaryDetector.lean)

Of 24 pairs of waits (six per owner and reference setting), the registered contrasts separated the ranges in 11 native comparisons, 19 after reading the first source direction before recombination, and one when only the flight contrast was halved. All are covered by the same simultaneous 99% counting bounds. Overlapping ranges leave their comparison undecided. The setup samples an orbital flight and consults material records stored at the conversion site. These additional records and time-dependent spatial acceptance remain possible physical explanations. Recovering a weak-flavour law demands independently measured charged partners and acceptance on the same arrival.

Use one numerical thread in the repository's scientific Python environment.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.fixed_detector
python -Bm "$p.measure" --acquire /tmp/oph-stationary-detector.json
cmp /tmp/oph-stationary-detector.json model/candidate_laws/interval_echo/openwave/row_10_oscillations/fixed_detector/counts.json
python -Bm "$p.measure"
python -Bm pytest -q --pyargs "$p.test_stationary"
~~~
