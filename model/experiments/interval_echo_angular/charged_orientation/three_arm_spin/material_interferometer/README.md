# Two ways to emit and detect the charged winding state

The three-path pointer and a separately constructed local circuit using the vacancy itself prepare the same two winding states. Either detector reads the phase of a vacancy emitted by either source. After 32 units of free evolution, the pointer's archived counts predict all 16 crossed measurements. A second set of 16 measured intervals agrees at the independently specified formation clock centered on **16.43120649** native units. Each set includes both source types, both detectors, and both winding senses at each end. [All original-start count banks](results.json) · [Count prediction and full free evolution](measure.py)

## How the vacancy serves as its own interferometer

At the emitter three blank spokes join hub 4 to three marked triangles. A sequence of three blank contacts splits the vacancy into equal-amplitude paths, one at each triangle corner. The dwell times are $2\arcsin(1/\sqrt3)$, $\pi/2$ and $\pi$. Once the vacancy is at the corner, eight pulses alternating between the two blank edges of that triangle return it to the same corner with phase $\omega=\exp(2\pi i/3)$. A triangle containing only unit loads stays filled and contributes phase one. The eight dwell times were calibrated against the isolated three-port gate alone. Their inverse writes the conjugate phase using positive dwell times. The emitter winds each of the three triangles once with a full written hop and two blank hops. It then reverses the blank splitter and requires a vacancy at the original hub. Other source outcomes count as refusals. [Native pulse sequence, stored dwell times and full local blocks](apparatus.py)

Let $\ket{j}$ mean that only direction head $j$ is reversed. Conditional on the hub read, the two preparations are $(\ket{0}+\omega^{s}\ket{1}+\omega^{-s}\ket{2})/\sqrt3$ for $s=+1,-1$. They have identical charge address and head populations. Their inner product is zero. Each source accepts one third of its starts. Hub 160 has its own three spokes and triangles. When a vacancy arrives there, the detector reverses its physical pulse sequence, reads the hub and tests all three heads for the unwound orientation. Each head test routes the vacancy to the near contact, records whether a directed crossing occurs, and reverses the test. A failed flag counts as a refusal. No arm label is recorded on an accepted coherent route. All source heads are discarded at the remote read.

The comparator instead uses a separate supplied coherent path pointer, pointer erasure, and a pointer measurement. The charge-only instruments use addressed native contact exposures and sharp local load heralds. They require the same prepared guide and remote compensating stock, and use no additional probe pair. The local calculations include all 80 linearly independent vacancy-position and head words at each terminal. The resulting success amplitudes at the hub agree with the pointer states to $2.6\times10^{-15}$ per entry up to an overall phase. These are independently constructed instruments for the **same** input states. [Pointer source and detector](../RECIPROCAL.md)

## Three detector readings on fresh launches

The table gives the phase-matched counts. Each row uses independent banks for the two analyzer settings and includes all source and detector refusals. The charge prepares itself and reads itself at the opposite owner through the local circuits above.

| Propagation and original starts per bank | Emitted winding sense | Analyzer $+$ | Analyzer $-$ |
| --- | ---: | ---: | ---: |
| Free, $2^{38}$ | $+$ | 41,107 | 1,514 |
| Free, $2^{38}$ | $-$ | 1,559 | 41,092 |
| Formation clock, $2^{42}$ | $+$ | 21,244 | 85 |
| Formation clock, $2^{42}$ | $-$ | 92 | 21,161 |

At free time 32 the bare charge counter accepts $7.25848\times10^{-7}$ per original start for either sense. Matched and opposite phase probabilities are $1.49650\times10^{-7}$ and $5.57103\times10^{-9}$. Reading and forgetting the source winding word before free flight erases this contrast. The material analyzers then give 14,659 and 14,740 successes in the two senses, with overlapping simultaneous intervals.

The supplied formation clock coherently combines $2^{25}$ delays spaced by $1/32$. With its original center, both senses have clock admission $6.03650\times10^{-5}$ per original start. Across the four material settings the simultaneous 99% intervals have endpoints between $6.0345\times10^{-5}$ and $6.0381\times10^{-5}$. The bare receiving charge appears with probability $1.68103\times10^{-8}$ per start. For the calibrated pulses, projecting the prepared amplitude onto the two symmetry coordinates used for the clock changes its norm by less than $3.1\times10^{-15}$. Subsequent conditional operations are contractions, which bound the omitted amplitude by the same amount. The matching winding read yields $4.84771\times10^{-9}$, and the opposite read yields $2.16610\times10^{-11}$. A clock-selected state can extend to the receiver during preparation. The ordinary free run starts with the vacancy at the emitter.

For the forecasts, each local winding sense corresponds to a known $Y$ eigenstate of the calibrated pointer. The archived pointer counts give a Clopper--Pearson interval for each receiving outcome. In the narrow-clock archive the retained reference is conjugated, so opposite $Y$ signs on the reference and receiver correspond to a matched physical input and read. The source curves obey arm-reflection symmetry, interchanging both senses and leaving the matched probabilities equal. The count comparison uses a union bound for all training and test intervals. The material data recover the guided excitation's **internal winding phase**. Its transformation under a spatial rotation, helicity, and three-dimensional rest mass require different experiments.

From the repository root, use one numerical thread.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.material_interferometer.measure --write
python -m pytest -q model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/material_interferometer
```
