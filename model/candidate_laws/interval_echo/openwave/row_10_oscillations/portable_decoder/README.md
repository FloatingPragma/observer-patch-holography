# Does the same winding decoder work for another observer?

A material observer calibrated six click labels at one receiver and one time. The two observers then applied exactly those labels to fresh simulated trials at two other times, including a receiver on another cell. They compared a successful interception of the neutral pair with a sharp reading at the conversion site. The calibration and trial cohorts have different seeds. [Frozen design](PROTOCOL.md) · [Calibration counts](../joint_receivers/blind_counts.json) · [Held count banks](counts.json) · [Replay and decoder](measure.py)

### Held flight appointments

The first score gives each of the three winding preparations equal weight among intercepted pairs. The second uses an independently frozen binary decision on the two signs from the same calibration. Intervals have simultaneous 99% coverage across the enrolled comparisons. The last column averages the native pair yield over equally frequent prepared recipes and analyzer panels, counting every original simulated launch.

| Receiving owner | Flight time | Three recipes correct (%) | Two signs correct (%) | Pairs per trillion starts |
| --- | ---: | ---: | ---: | ---: |
| 12 | 384 | 33.97 [33.59, 34.37] | 50.71 [50.10, 51.32] | 5.86 |
| 12 | 768 | 33.80 [33.35, 34.27] | 49.91 [49.18, 50.64] | 9.45 |
| 17 | 384 | 33.48 [32.55, 34.43] | 50.18 [48.79, 51.58] | 4.62 |
| 17 | 768 | 34.67 [34.09, 35.26] | 51.37 [50.45, 52.29] | 6.68 |

Three of four joint three-recipe scores exceed uniform guessing. Two of the four binary scores resolve a sign preference in the enrolled direction. A prior sharp read of the first source direction erases the resolved binary preference at all four settings. The three-way decoder frequently labels both signs as the zero recipe. Applied to only the signed preparations, with a zero guess scored as wrong, its correct fraction ranges from 2.51% to 11.94%. On a fair random fallback for each absent pair, the three-way decoder gains just 0.014 to 0.273 expected correct labels per trillion original starts over guessing on every launch. The code computes this difference directly from integer counts. Conditioning on receipt with equal recipe weights and averaging over all launches answer different questions because pair yield varies with preparation.

Halving only the flight contrast at fixed preparation and detector gave binary scores of 53.31% to 57.38%. This control demonstrates dependence on the material during flight. The physical read uses the emitter's parked loan vacancy and conversion-site memories after classical comparison of records; its operation requires the loan's external contacts to be parked during the wait. Nothing in the label rule tags a weak interaction or measures a propagation mass. Interpreting the three scores as neutrino mixing angles would require an independently identified charged partner in both production and detection.

From the repository root, with its scientific Python environment:

~~~sh
export VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
p=model.candidate_laws.interval_echo.openwave.row_10_oscillations.portable_decoder
python -Bm "$p.measure" --acquire /tmp/oph-held-winding.json
cmp /tmp/oph-held-winding.json model/candidate_laws/interval_echo/openwave/row_10_oscillations/portable_decoder/counts.json
python -Bm "$p.measure"
~~~
