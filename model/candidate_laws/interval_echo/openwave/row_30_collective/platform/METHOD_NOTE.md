# M10.5 · Measuring a matter wave on paired material

$$
H_{\rm strand}=J\sum_{i=1}^{n-1} Z_i Z_{i+1}-h\sum_{i=1}^{n}X_i+J(Z_1-Z_n)+C,
\qquad \omega^2=4(J-h)^2+4Jh\,[2\sin(k/2)]^2 .
$$

The first equation is the restriction of the free seam generator to a supplied material sector. The second is a forecast for a mobile defect in the alternating orientations of loaded dimers. The sign of the boundary fields follows the strand orientation. The additive constant cancels from transition frequencies. For the spatial phase measured along the strand, the difference between squared lattice momentum and $k^2$ begins at order $k^4$. The ordinary small-$k$ continuum test and the finite-lattice test have separate gates. [Proposed task](../tasks/m10_5_task_details.md)

A port-pair observer reads which endpoint holds the load. Two seeded matchings commission the material by physical seam operations. The observer enrolls the longest open strand before measuring a frequency. Fresh local clock and phase records from outside the enrolled band determine $h$ and $J$. Changes in wall orientations reach detectors in other cells. Pair lines $s_{ij}$ determine individual frequencies through $\omega_1=(s_{12}+s_{13}-s_{23})/2$ and $\omega_i=s_{1i}-\omega_1$. Other pair sums test the reconstruction. Spatial scans independently measure $k$. On the boards of this experiment all enrolled momenta are finite. The intercept $2|J-h|$ is a prediction from the two local meters.

## Equation-to-code map

| Quantity | Executed source in the platform tree |
| --- | --- |
| Native echo, exchange and squared-contrast columns used in the reduction | [Native seam operators](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/contrast_memory_spherical_v1/law.py#L67) |
| Ternary stock, physical matched collars and full-board load conservation | [Commissioning procedure](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/interval_echo/openwave/row_30_collective/formation.py#L87) |
| Finite defect generator and direct physical spin-space checks | [Reduced generator and exact small chains](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/interval_echo/openwave/row_30_collective/experiment.py#L97) |
| $h,J$, bulk intercept and characteristic squared speed from local records | [Local coupling instruments](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/interval_echo/openwave/row_30_collective/coupling_meter.py#L120) |
| Pair spectroscopy, spatial scans and free wall transport | [Count generation and decoding](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/interval_echo/openwave/row_30_collective/experiment.py#L237) |
| Wave numbers, line residuals, rival tangent and changed-contrast response | [Measured dispersion tests](https://github.com/openwave-labs/openwave/blob/main/openwave/xperiments/m10_oph/research/scripts/m10_5/vendor/model/candidate_laws/interval_echo/openwave/row_30_collective/holdout.py#L15) |

## Reproduction

The OpenWave port carries the frozen author implementation, native sphere fixtures and measured local calibration records. It re-executes the dynamics and the observer apparatus using those platform-local files. The author's Git tree serves as a byte reference. The same finite-sector compression is used in both implementations, so matching them checks portability and custody. Independent evidence for the reduction comes from diagonalizing the physical spin generator and testing its two-time detector responses on small chains. The platform entry point records the origin of every imported model module. [Platform runner](../scripts/m10_5/run.py) · [File map](../scripts/m10_5/layout.json)

From an OpenWave checkout with the port installed, select the Python environment specified by the shipped requirements. In the commands, the source checkout supplies the release named in the task's companion source pointer. Use a fresh output directory for each acquisition.

```sh
runner=openwave/xperiments/m10_oph/research/scripts/m10_5/run.py
python -I -B "$runner" pilot --source /path/to/oph --release /path/to/pilot
python -I -B "$runner" acquire --source /path/to/oph --release /path/to/confirmation
python -I -B "$runner" replay --source /path/to/oph --release /path/to/confirmation
```

The pilot always uses the exposed seed 927016. It writes the measured records and checks their count-stream hashes against the saved holdout. On that source the author reports 19 lines across three board sizes with a largest finite-lattice discrepancy of 0.000585 clock units. The absolute gate is 0.0048. These pilot results were known before the proposed registration. After the exact task lands, the owner can execute the port on an OpenWave branch descending from the registration. The formal path retains every count array, including failed settings, and evaluates every registered gate. The execution command prints the source and platform revisions alongside its import inventory.

**Scope.** This test concerns an effective one-dimensional wave in a commissioned material, with length in dimer spacings and frequency read by a material hand. Its coverage excludes an isotropic three-dimensional vacuum dispersion, an absolute particle mass in laboratory units, and a zero-momentum transition measured directly on these finite open strands.
