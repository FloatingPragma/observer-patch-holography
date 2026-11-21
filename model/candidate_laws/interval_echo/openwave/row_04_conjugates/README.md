# Can a port relabeling supply the partner of the persistent pair?

All five nonidentity permutations of the three port loads change the terminal records of the row 4 pair on the 20-cell screen. Exchange alone makes all six assignments indistinguishable in the calculated records. This tests an inexpensive candidate for a particle–antiparticle map. The permutations act on the ports and leave the zero-load reference registers fixed. Electric charge requires an independently calibrated interaction. [Protocol](PROTOCOL.md) · [All simulated counts](results.json)

The source places the exceptional load on ports 0 and 1 in a uniform background. After free evolution one patch reads the seam carrying the two exceptional loads, or quiet. It uses the same echo with its output labels mapped for each preparation. The largest total variation distances in the table occurred at native exposure 4. Each comparison includes the 30 addresses and quiet.

| Background | Exceptional load | Range of the generator difference, native energy units | Largest record distance from the row 4 source |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 0 | 0 |
| 0 | 2 | 48 | 0.08667 |
| 1 | 0 | 32 | 0.64079 |
| 1 | 2 | 32 | 0.46099 |
| 2 | 0 | 176 | 0.64422 |
| 2 | 1 | 64 | 0.70327 |

The 65,536 simulated preparations per setting independently resolve each alternative against the unpermuted source. On an isolated echosahedron every port has degree five, and all four assignments with a one-unit difference between background and exceptional load have identical free record laws.

For background load $b$, let $s$ be the exceptional load minus $b$, and let $r_i=6-d_i$ count the degree references at address $i$. At exceptional ports $i,j$ the contrast energy above the uniform background is

$$
\kappa\left[s^2(12-2\,\mathbf 1_{i\sim j})+2bs(r_i+r_j)\right].
$$

Exchange acts identically in the six preparations. On the full sphere, flipping the deviation changes the reference contribution at gateways relative to internal ports. Changing its magnitude changes the binding energy. This explains both the local null control and the observable difference on the complete screen. An excitation with a different internal state, transformed references or extended seam dressing falls outside this test. [Native sector and reference law](../../../contrast_memory_spherical_v1/SPEC.md) · [Implementation](experiment.py)

Run from the RER root or extracted evidence bundle.

```sh
OPENBLAS_NUM_THREADS=1 python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_conjugates.experiment
```
