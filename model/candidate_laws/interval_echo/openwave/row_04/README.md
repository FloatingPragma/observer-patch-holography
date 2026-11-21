# Row 4 · Free stability of a localized Interval Echo pair

An initially adjacent pair of unit loads was released on the complete spherical federation with every seam active. After 150 native time units on 20 cells, a one-cell observer could still detect the pair on one of its 30 internal seams with probability **98.3881%**. Only **2.8306%** lay on its starting seam. The material distribution moved within the patch while the constituents overwhelmingly stayed together. Removing the contrast interaction lowered the same local detection probability to **14.9559%**. [Frozen test](PROTOCOL.md) · [Measurements, detector counts and addresses](results.json)

| Sphere and last exposure | Preparation or interaction | Adjacent anywhere, lowest over tested times | On a seam in the original cell, at last time | Both in original cell, at last time |
| --- | --- | ---: | ---: | ---: |
| 20 cells, 150 | Sharp adjacent release | 99.0688% | 98.3881% | 98.7525% |
| 20 cells, 150 | Native partial exchange before release | 95.0280% | 94.4304% | 95.8774% |
| 20 cells, 150 | Exchange-only control | 19.1675% | 14.9559% | 25.3091% |
| 80 cells, 32 | Sharp adjacent release | 99.0545% | 98.3465% | 98.6924% |
| 80 cells, 32 | Native partial exchange before release | 95.0696% | 94.5248% | 96.4974% |
| 80 cells, 32 | Exchange-only control | 24.2755% | 21.9452% | 34.4417% |

The table gives finite-unitary probabilities over every preparation, including those recorded as quiet. A separate sample had **16,093/16,384** native patch detections at time 150, against **2,474/16,384** with exchange alone. Each detector uses a fresh ensemble. It supplies no stabilizing intermediate measurements. At that exposure, the patch radius enclosing 90% of preparations was zero cell steps for the native sharp and perturbed pairs and four steps for exchange alone. These radii use the simulated position distribution. The one-cell detector uses the registered native echo on the 30 internal seams and can directly register either a seam address or quiet. [All appointments and detector outcomes](results.json) · [Source and terminal instrument](../../../../experiments/interval_echo_encounter/README.md)

## What stabilizes the pair

The registered native interaction on a blank collar contributes an exchange and a penalty proportional to the squared difference in endpoint load. Interval Echo includes the common degree-reference terms. With two distinct units, their total contrast is 10 on an occupied seam and 12 when separated. After subtraction of 12κ times the identity, the full binary-sector Hamiltonian is H = K − 2κP, where P detects an adjacent pair. The measurement uses the registered κ = 8. All internal and external exchanges operate throughout flight. The only supplied interventions are initialization, the optional initial partial exchange and the terminal read. Evolving with κ = 0 preserves the same load count and admits the large losses shown above. [Native law](../../../contrast_memory_spherical_v1/law.py) · [Exact sector restriction](../../../contrast_memory_spherical_v1/sectors.py)

The isolated low-energy band is invariant under free evolution. The registered all-size bound gives Pr(adjacent at time t) ≥ 1 − 40/(2κ−6)² = **0.60** for every real t and every normalized release supported on occupied seams. Its proof uses bounded degree, the native cell incidence and the external matching. This is a statement about constituent separation. Spatial containment of the pair in its original cell has the finite exposure and board limits in the table. The two screen sizes have different gateway placements and observation horizons. Repeating the local result on them gives a size check at their common appointments. [Binding proof](../../../contrast_memory_spherical_v1/PROOF.md)

Before release the sharp source places one unit at port 0 and one at port 1. A separate check used the 36 uniformly weighted, initially separated two-unit preparations on cell zero of the 20-cell sphere. After one unit of unassisted evolution, a sharp read at the specified seam prepares precisely the released state with probability **0.01891418%**. The remaining **99.98108582%** of those starts are quiet at that seam. The high persistence probabilities in the table apply after a successful preparation. The perturbed arm moves 4% of the preparation probability into a nonadjacent configuration by a single native partial exchange. [Independent encounter experiment](../../../../experiments/interval_echo_encounter/README.md) · [Present source and branch calculation](./__main__.py)

## Scope for an OpenWave submission

**The dynamic free-persistence test passes for this prepared two-load excitation.** The submission concerns a charge-unassigned composite on a discrete quantum federation. Its electromagnetic coupling, spin and identification with the electron have yet to be established. Derrick's continuum scaling theorem concerns a different class of field functionals. Here a spectral gap explains the measured relative binding. The simulation runs in the RER executable. Upstream validation and placement are separate decisions. [Public criterion](https://github.com/openwave-labs/openwave/blob/26fbb5fdd07049625d6bb7930af2f6bf545fb6fc/MODELS.md) · [Candidate specification](../../SPEC.md)

## Reproduction

Run from the RER repository root with its model dependencies installed. The first command rebuilds both boards, evolves all three arms and verifies the archived receipt. The second compares the sparse evolution against a direct dense cell exponential and validates the native preparation and readout columns.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04
OPENBLAS_NUM_THREADS=1 python -B -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_04
```

The binary basis contains 28,680 configurations on 20 cells and 460,320 on 80 cells. A zero of energy has been removed. Uniform Chebyshev propagation uses the enclosure [−2κ,12] and an analytic Bessel-series tail below 10⁻¹² at each appointment. For any s > 0, the contour representation gives |Jₖ(z)| ≤ exp(z sinh s − ks). Summing this bound geometrically controls the omitted tail without access to a desired trajectory. The exact single-cell exponential and a direct small-board sparse exponential provide separate numerical comparisons. The receipt records model probabilities, simulated finite counts and the supplied native exposure for each setting. [Implementation](./__main__.py) · [Independent checks](test_release.py)
