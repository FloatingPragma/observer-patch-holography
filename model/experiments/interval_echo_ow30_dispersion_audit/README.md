# Resolving the dispersion of a traveling dimer reversal

On three exposed Interval Echo strands, observers can distinguish two frequency laws that agree on all 22 lines of the published validation. Sending shorter waves through the same material separates their forecasts. The 22 new frequencies agree with the lattice square root to within 0.000739 dimer-clock units. Every line differs from the alternative by at least 0.01102, against a tolerance of 0.0048. [Complete readings](results.json) · [Test and exposure history](PROTOCOL.md)

The alternative has hopping out to two neighbouring sites. With the squared lattice wave number $u=4\sin^2(k/2)$, its frequency is $m+c^2u/(2m)-c^4u^2/(8m^3)$. This is the expansion of $\sqrt{m^2+c^2u}$ through second order in $u$. Both use the same mass and speed readings from the independent coupling meter. Each already predicts the original 22 frequencies within 0.000575. At short wavelengths their higher-order terms separate. This compares two explicit frequency forecasts for the same counted responses. The range-two formula is a rival curve, with no claim to predict the full instrument. [Original validation](../../candidate_laws/interval_echo/openwave/submission/row_30_confirmation/validation.json) · [Coupling meter](../../candidate_laws/interval_echo/openwave/row_30_collective/coupling_meter.py)

| Cells | Measured wave number (radian/dimer) | Frequencies | Largest lattice error | Smallest range-two error | Odd-hold spatial residual / ceiling | Residual using launch labels |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 2.386–3.016 | 6 | 0.000734 | 0.01150 | 416.63 / 554.07 | 2,307,004 |
| 80 | 2.479–3.059 | 8 | 0.000739 | 0.01275 | 886.67 / 1137.74 | 641,366 |
| 320 | 2.355–3.043 | 8 | 0.000693 | 0.01102 | 820.31 / 963.37 | 1,487,012 |

This is a simulation in the exact reduced sector of a commissioned open strand. Native exchange on each two-port dimer moves its one load with the same amplitude. Each connecting seam retains its prepared coherent memory and contributes the same conditional neighbour energy. Exterior seams are idle and the two ends have fixed boundary energies. The resulting orientation Hamiltonian depends on strand length and these couplings. The execution repeats the full-board commissioning and native-column checks, checks the reduction against full orientation matrices of five and six dimers, and evolves its reduced two-time transfer matrix on each measured strand. It reuses the three validation geometries with new sign counts and new high-frequency launches. [Reduction and boundary conditions](../../candidate_laws/interval_echo/openwave/row_30_collective/PROTOCOL.md#local-reduction-and-exact-dispersion) · [Preparation and reduction receipts](results.json)

Short-wavelength sine launches reverse the sign on alternating receiver addresses relative to their long-wavelength partners. The intensity weights in the original spatial decoder therefore apply to both. Dimer-clock counts set the temporal step, pair-transition counts give the frequencies, and the fitted spatial phases from even holds predict all receiver addresses on odd holds. Substituting the high launch labels fails on every board. With the native wave turned off by zero contrast, the largest amplitude at a candidate line is 0.000210 and the calibration amplitude at frequency 2 remains above 0.478. The largest held pair-sum error is 0.000176. [Executable experiment](measure.py)

The short-wavelength result concerns the finite lattice dispersion. The continuum Klein–Gordon comparison uses ordinary $k$ in the registered window $0<k\leq0.50$. All the new lines lie outside it. They distinguish the two microscopic frequency curves on this commissioned one-dimensional strand. The earlier long-wavelength data, their tolerances and the separate upstream decision on the state and decoder retain their original scope. This is a same-group supplementary measurement on exposed boards. [Candidate and terms](../../candidate_laws/interval_echo/openwave/row_30_collective/confirmation/TASK.md)

At OpenWave revision `0df73dc6`, the canonical M10 pin was `4fe10f9a`. The frozen v3 source at `92370140` returned exit 5 from its read-only upstream admission check. The official history had no landing for the v3 terms. The row was untested and no public selectors were drawn. [Canonical pin](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#L12-L17) · [Official criterion and cell](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L467)

From the RER root run the following command with the model dependencies installed. It regenerates the acquisitions and compares the report to its saved values. Add `--retain /path/to/new-directory` to keep every raw binary count array (about 285 MiB).

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.experiments.interval_echo_ow30_dispersion_audit.measure
```
