# Free propagation of the bound pair

The two-load pair from the stability experiment can reach another echosahedron with every seam on the sphere active. Its lowest 20 energy states form an isolated band. Under two cell-based spatial measurements the tested spatial profiles do not give a single resolvable Klein–Gordon frequency at each wave number. The arrangement of the three gateways in each cell accounts for most of the frequency spread. [Full-band protocol](PROTOCOL.md) · [Spectral weights and patch counts](results.json)

## What neighboring observers detect

Starting on seam (0,1) in cell zero, the pair evolved freely on the 20-cell sphere. At exposure 256 a terminal detector in cell 4 registered an adjacent pair in **387 of 262,144** simulated trials. Parking the pair gave zero distant clicks. Exchange alone gave 1,456 clicks at that receiver and much weaker retention in the birth patch. A distant detection establishes transfer; the reader accepts any adjacent pair in its patch. Each setting uses an independent ensemble, and the simulated records include quiet.

Repeating the release from the first cell of each gateway class gave transfer both within and between classes. At exposure 256, cell 2 detected 417 pairs from cell 1, of the same class, and cell 0 detected 345 from the same source, across classes. Each detector received 262,144 trials. The patch address identifies its gateway geometry. An adjacent-pair click carries no spectral-band label. [Every cell, both sources and all four appointments](classes.json)

A joint phase can also be read after free flight. Two disjoint seam words have identical single-load states for either source sign. In counts simulated from the predicted outcome distributions with 67,108,864 launches per sign, four of six registered neighboring routes pass the test for the source sign among admitted arrivals. All six routes have different predicted flag fractions for the two signs, and the two that miss the cutoff have the smallest differences. The full pair has an all-time adjacent probability of at least 47/64 after source acceptance. Crossing into a neighboring cell requires amplitudes on separated configurations. [Preparation, range of efficiencies and all six tests](PAIR_PHASE_PROTOCOL.md#what-the-receivers-record)

## Why a free mass shell fails this test

The full 28,680-state native Hamiltonian retains every separated configuration. The lowest 20 eigenstates each carry at least 98.62% adjacent-pair weight, with bandwidth 0.03585 and a gap of 0.37563 to the next state. A spatial amplitude uniform over each cell's 30 internal seams puts more than 91.6% of its weight into this band in every cell-graph momentum class. Its zero-momentum RMS spectral width is **0.01738** in native frequency units. The first nonzero momentum shifts its centroid by only **0.00002705**. Using each isolated cell's ground pair orbital in the spectral projection gives the same ratio of width to shift, about 642. A good fit to the centers of these broad distributions would provide little information about a propagating matter wave. [Full generator and analysis](experiment.py)

The eight cells with adjacent gateway ports have more than 99.9998% overlap with the lower frequency group when their internal seam orbital is projected and orthonormalized within the band. For the twelve cells whose three gateways all lie two internal steps apart, the corresponding upper-group overlap exceeds 99.9997%. The split between the groups is 0.03563. The largest off-diagonal coupling in these coordinates is 0.00004207. This spectral association motivated an additional test on each geometry. [Gateway diagnostic](GATEWAYS.md) · [Effective Hamiltonian and all overlaps](gateways.json)

## Holding the gateway geometry fixed

Both geometries were tested using their direct intercell seams and a second spatial metric that includes *all* routes through the complementary cells. The latter is the Schur complement of the original cell Laplacian, obtained by minimizing the graph energy over the complementary amplitudes. The free-flight generator is identical in all four comparisons. [Class-specific protocol](CLASS_PROTOCOL.md) · [Four comparisons and finite patch records](classes.json)

| Selected cells | Spatial metric | Spectral result |
| --- | --- | --- |
| 8 with adjacent gateways | Direct seams | Four separate dimers. Both frequency classes have RMS width below 0.12% of their separation. There is one nonzero squared wave number, enough to calibrate one shift and none to test it. |
| 8 with adjacent gateways | All routes | The shift is 0.000003844 at λ = 1.1315 and falls to 0.0000000616 at λ = 1.2. Any positive-mass Klein–Gordon curve anchored to the first shift predicts at least 64 times the second measured centroid shift. |
| 12 with separated gateways | Direct seams | The half-gap between the two class rest centroids is 0.017846. With speed calibrated at the first nonzero momentum, held shifts miss by 12.0% to 59.2% of the forecast. The rest linewidth exceeds the first shift. |
| 12 with separated gateways | All routes | The first shift has negative sign, giving negative inferred speed squared for the proposed positive-frequency branch. The zero-momentum width is three times the magnitude of this shift. |

Here λ is squared wave number in the stated graph metric. The half-gap and all bandwidths come from the generator in native exposure units. An independent rest-frequency reading, physical ranging and a branch-selective patch preparation are additional requirements for an observer-measured mass shell. The finite tests address these orbitals and these metrics on one screen. The one-dimensional driven-guide experiment has a different propagation mechanism. [Guided measurements](../row_30/README.md)

## A 56-cell connected test on the larger sphere

The 80-cell screen has **56 connected cells** whose gateway ports all lie two internal steps apart. Their geometry offers 56 distinct squared wave numbers, including zero. On the 320-cell screen the largest component of that class has 106 cells; its spectrum was outside the experiment. On the 80-cell sphere, the full two-load sector has 460,320 configurations and 2,520 active seams. More than 93.96% of each tested orbital lies in the upper cluster of 56 states. [Protocol and numerical bounds](REFINEMENT.md) · [Connectivity, band isolation and measured widths](refinement.json)

| 56-cell metric | First shift from zero momentum | Second shift | RMS width at zero momentum |
| --- | ---: | ---: | ---: |
| Direct intercell seams | −0.0000001574 | +0.0000116190 | 0.000026270 |
| All routes through the sphere | −0.0000079013 | +0.0000075927 | 0.000026270 |

The first shift has the wrong sign for positive speed squared on this upper branch in both metrics. The first two shifts have opposite signs, so changing a common energy zero or reversing the orientation of phase cannot make both positive. The spectral widths also fail the enrolled resolution gate. A radius-ten solve retained 93,548 configurations, with a maximum residual of 6.5 × 10⁻¹⁰ against the full operator, under 0.42% of the smaller first shift. The looser radius-six and radius-seven solves give the same signs but exceed the precision threshold on one or both metrics. This larger connected class resolves the shortage of wave numbers in the four dimers and leaves the free mass shell unestablished. The refinement measures a spectral prediction for the candidate patch orbitals. The observed propagation counts above use the 20-cell screen. [All 56 momentum classes in each metric](precision.json)

## Reproduction

From the RER root, or the extracted archive root, run:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -B -m model.candidate_laws.interval_echo.openwave.row_30_branch.experiment
python -B -m model.candidate_laws.interval_echo.openwave.row_30_branch.gateways
python -B -m model.candidate_laws.interval_echo.openwave.row_30_branch.classes
python -B -m model.candidate_laws.interval_echo.openwave.row_30_branch.refinement
python -B -m model.candidate_laws.interval_echo.openwave.row_30_branch.precision
```

The commands compare fresh computations with their receipts. Passing `--output PATH` creates a separate receipt. The tightest numerical solve used 655 MB of peak resident memory. The register paths and counts are local measurements on independently prepared copies. The branch calculations diagonalize the generator and do not count as a local frequency measurement.
