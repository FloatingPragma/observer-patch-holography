# Do clocks share a change of elapsed time beside a freely moving pair?

Two internal observers start with independently calibrated clock lines. They compare a binary Ramsey outcome on identical source preparations, using a shared start appointment and analyzer phases. An excess and a vacancy, launched separately on the unit-filled board, continue to move while the hands run. At the first receiver the two line records exclude *every common distribution of elapsed times* as the explanation of the entire change. The clock hands have different local material interactions with the source. These observations measure how the hands respond; a common proper curvature would require an independently calibrated optical and stress measurement. [Counts from all launches](results.json) · [Active-law operations and replay](measure.py)

## What happens before a hand is read

Two launch detectors report an excess at 10 and a vacancy at 1. Their loads have no common contact initially. After 1/8 exchange-time unit of free evolution on the unmarked 20-cell host, a pair of source detectors can announce the ordered encounter (0,1). Its probability is 0.00092907. All surviving and failed launches appear in the heralded record. A separate acquisition makes no intermediate measurement and continues **every** source directly to the clock. Neither receiver uses a terminal from the initial or heralded pair. Their two seams are (2,11) and (12,13). These identifiers specify operations and give no physical distance or velocity.

Sharp load and swap-parity measurements on separate copies give initial energy $[101.99997,102.00003]$ above the uniform background in exchange units. Free flight conserves that energy distribution. On the sharp stock *conditional on the encounter reading*, the same tests give $[132.99997,133.00003]$. The extra conditional energy includes the selection and preparation by the encounter measurement. The unmarked sharp encounter stock has a previously bounded probability of at least 85.47% of remaining on enrolled adjacent high-energy addresses at every free dwell. In the clock experiments both constituents remain present, interact with the clock and then undergo an independent contact-coincidence read. Over 98% of the heralded stock gives a final contact at each of the three near appointments and all three clock settings. Contact here means a joint material detector event. Assigning a physical rest mass requires a dispersion measurement in calibrated observer units. [Free-pair spectrum, retention bound and material readouts](../../../../../experiments/material_inertia/README.md)

## Exact joint clock and matter dynamics

On the one written seam, use a contrast of magnitude two. Choose reversal parity $\pi=\pm1$ and a symmetric or antisymmetric coin. Direct application of the full 99-state native operator to all nine endpoint words gives a scalar $g=32$ or $33$ on the antisymmetric coin for the two parities. On the symmetric coin the scalar is $g$ when at least one terminal differs from unit load, and zero on the word $(1,1)$. These columns are invariant exactly. The other seams are unmarked and continue evolving through the exposure. The marked seam changes the local matter Hamiltonian in every enrolled arm; its interaction is explicitly included.

Denote by $K$ the ordered excess--vacancy Hamiltonian with the blank term of that one seam removed and by $Q$ the projector onto the clock terminals having an excess or vacancy. The two coin branches have generators $K+gI$ and $K+gQ$. The unit standard compares opposite parities with the antisymmetric coin in both arms. Its two generators are $K+32I$ and $K+33I$. Hence its gap is exactly one for **every** material state, including a moving, entangled pair. In the filled source-absent preparation the other two standards have gaps 32 and 33. The relative analyzer phases are set to these integer multiples of the same blank calibration phase. Source-absent null outcomes follow from the joint generator. The published acquisition records one quadrature, the right or wrong port, for each hand.

For every incoming wave $\psi$, the probability of the wrong port for line $g=32,33$ at the nominal blank phase is
\[
 p_g=\tfrac14\Vert (e^{-it(K+gQ)}-e^{-itK})\psi\Vert^2 .
\]
The protected line has $p_1=0$. The real part of an interference amplitude determines this probability, so a wrong port can reflect a phase or a visibility change. For the heralded source $Q\psi=0$ and $\|QK\psi\|^2=1/4$ on the first receiver. Taylor expansion gives $p_g=g^2t^4/64+O(t^5)$ at small $t$. The leading response already requires a source path to a clock terminal. It is a change in a measurable hand under source contact.

## Comparing the records without choosing a time coordinate

Let the shared elapsed time be any random variable $T$ and the common calibration phase be $\phi$. If elapsed time alone accounts for the changes then $p_g=\mathbb E\sin^2[g(T-\phi)/2]$. From $|\sin(gu)|\le g|\sin u|$ for positive integers $g$, every distribution and every phase convention satisfy
\[
 p_{32}\le1024p_1,\qquad p_{33}\le1089p_1.
\]
The bound permits arbitrary winding, stochastic clock offsets common to the hands and arbitrary correlation between the time and the source. It compares like-for-like preparations. Different clock-specific material interactions are an explicit way to violate its hypothesis.

At the first receiver after a quarter-unit exposure, $p_{32}$ is in $[0.0028452,0.0028582]$ on the heralded stock. The unit line implies $1024p_1<0.0000101$. With **no herald measurement**, the line-32 wrong-port probability lies in $[0.0104641,0.0104650]$, against $1024p_1<0.000000122$. All bounds hold at 99% simultaneous coverage. Thus selection of a rare encounter is unnecessary for this effect. On that unselected acquisition the final pair-contact coincidence probability is about 5.96% for line 32; it cannot be interpreted as a read exclusively on bound pairs. The more distant receiver in the wiring has much smaller responses over these appointments. A universal lapse would have to account for co-located standards at the first receiver alone.

Each setting uses $2^{40}$ simulated starts. Four terminal coincidence bins retain both source outcomes and both hand outcomes, with the nonheralded count recorded when a herald is required. The interval calculation uses Clopper–Pearson bounds conditional on the actual number of successes for each bank, and a union bound across all reported comparisons. The independent sharp-stock energy measurements use bounded swap-parity outcomes. The sparse, Hermitian joint generator has 57,360 ordered source configurations. A scaled Taylor series bounds the norm error in each evolved vector; the replay independently checks the formation amplitude against matrix-exponential propagation. All measurements use internal detector labels. These finite observations alone assign neither an affine ray coordinate nor a universal gravitational coupling.

~~~sh
OPENBLAS_NUM_THREADS=1 python -m model.candidate_laws.interval_echo.openwave.row_27_metric.formed_pair_clocks.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_27_metric/formed_pair_clocks
~~~
