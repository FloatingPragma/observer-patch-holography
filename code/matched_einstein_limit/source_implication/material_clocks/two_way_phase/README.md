# Can two stations register a boost on the entropy equilibrium?

On the 80-cell equilibrium two material standards measure the phase between arriving pulses on the gateway $1\leftrightarrow38$. The signals travel on the free board, and each pulse has its own first-hit record. These are separate A-to-B and B-to-A first-receipt cohorts. The forward and reverse records jointly give a fractional radial velocity within **0.000524 of zero** at both shutters for both probe charges, under a collinear Doppler interpretation. Giving the second pointer a known extra wait after its recorded arrival yields **0.799926 to 0.800072** by the same inference. Both experiments leave the detectors parked. A change in the pointer appointment alone can imitate the $4/5$ axial speed in the ten-observer construction. [Native flight and clock reader](measure.py) · [Simulated sign counts for every start](results.json)

## Which clocks and events are compared?

Each trial prepares the same normalized two-head equilibrium on four registered ports $(0,1,38,39)$. Choose the direction and sign of a fresh orbital probe before launch. Record its first hit at the opposite gateway at half-tick appointments, or a miss. In each direction one cohort launches at zero and a separate one at $\Delta$. The two settings use $(\Delta,H)=(1/2,3),(1,11/2)$ in exchange units, with $H$ the deadline after each launch. The parked head generator is scalar on each reversal block, so waiting before launch leaves this block-diagonal source invariant. Changing direction or charge can change the entire transit-time distribution. [Equilibrium and its preparation](../../full_algebra/PROOF.md) · [Native return and relay kernels](../../full_algebra/reuse/edge_clock.md)

One pointer is an isolated blank dimer whose two orbitals have unit angular-frequency separation. The second uses the sign of a parked marked head with even and odd reversal energies separated by one. Assign the dimer to port 1 and the sign reference to port 38 in both directions. Two ensembles per epoch measure cosine and sine with analysis phases separated by $\pi/2$. The native local operators on all 99 seam states were checked for both standards. Independent auxiliary copies carry the hands, avoiding the measured shift of source energy that an initial sharp sign measurement on the flight's original left head would cause. Appointments and phase gates are supplied. Operationally the script replays the actual native transit-time distribution and reads an independent local-clock copy at the recorded appointment. This uses no geometry-dependent proper-time coupling of that clock to an in-flight body. [Preparation, read disturbance and native clock operators](../PROOF.md) · [Independent propagation and gap checks](test_measure.py)

A first hit at tick $j$ contributes a signed Ramsey reading at the actual first arrival $j/2$ after launch. A miss contributes zero on each quadrature. In the prompt acquisition, integer-evolution first-hit probabilities $p_j$ give the unconditional complex mean for launch epoch $t$,
\[
 Z(t)=e^{it}\sum_{j=1}^{2H}p_j e^{ij/2}.
\]
The product $Z(\Delta)\overline{Z(0)}$ has argument $\Delta$ whenever the sum is nonzero. For the smaller deadline its absolute visibility $|Z(0)|$ is 0.1231–0.1727 across the four directed probe settings (rounded outwards). At the larger deadline it is 0.0548–0.0745. A separate cohort reads the sending clock at both launch epochs. Every endpoint count includes all starts. Phase offsets and stable detection losses cancel in each elapsed phase. The recorded schedule selects the unwrapped branch $0<\Delta\phi<\pi$. An unlogged extra winding requires a different clock calibration. The same prescription applies to all four directed probe settings.

## From the two frequency ratios to a radial velocity

For each direction form $D_{AB}=\Delta\phi_{A,\rm send}/\Delta\phi_{B,\rm receive}$ from its own two pairs of clock readings. In a collinear inertial optical comparison, with local oscillator gains $g_A,g_B$ relative to proper time and recession rapidity $\xi$, the reciprocal pulse records obey
\[
 D_{AB}=(g_A/g_B)e^{-\xi},\qquad D_{BA}=(g_B/g_A)e^{-\xi}.
\]
Thus the product $P=D_{AB}D_{BA}$ gives $\xi=-\log(P)/2$ and $\beta=(1-P)/(1+P)$, independent of the clock units at either end. Interpreting $\beta$ as a physical velocity requires those optical paths and clocks. Even for counterpropagating rays at one event the product measures only the longitudinal component if the velocity has an unknown transverse part.

| Native probe | Launch interval / shutter | Prompt radial velocity | Extra reception wait |
| --- | ---: | ---: | ---: |
| Vacancy | 0.5 / 3 | [−0.0003931, 0.0003449] | [0.7999448, 0.8000594] |
| Vacancy | 1 / 5.5 | [−0.0003836, 0.0004021] | [0.7999464, 0.8000538] |
| Excess | 0.5 / 3 | [−0.0004803, 0.0005069] | [0.7999264, 0.8000709] |
| Excess | 1 / 5.5 | [−0.0005232, 0.0005082] | [0.7999352, 0.8000720] |

The extra wait is $2\Delta$ on the later cohort in each direction; both cohorts use identical first-hit laws and apertures. Its population value is exactly $P=1/9$. If the reader kept only these frequency ratios and discarded its arrival and shutter records, the apparent radial Doppler shift would match an ideal collinear recession at speed $4/5$. The retained hit and read appointments distinguish the explanations. Completely suppressing flight makes the Ramsey first-hit signal zero and leaves velocity undefined.

The raw archive preserves the plus and minus counts at each successful first-hit tick and the separate miss count for each analyzer. Joint sampling of these bins uses their native probabilities. A returning pulse and a common-event moving observer are outside this acquisition. Each tabulated phase comes from independent batches of $2^{38}$ simulated starts per quadrature. With per-start score in $[-1,1]$, the simultaneous error for each real component is at most $\sqrt{2\log(2\cdot256/0.01)/2^{38}}+7\cdot10^{-7}$ at 99% coverage. The union includes the local and remote quadratures for all eight panels. The calibrated equilibrium's Hilbert–Schmidt uncertainty is below $2\cdot10^{-7}$; since each component effect has operator norm at most one, its state-induced error is at most $\sqrt{8}\,2\cdot10^{-7}$. The remaining allowance covers the integer-propagator recurrence through eleven first hits, conversion of its four-dimensional blocks to double precision and the $2^{-40}$ probability grid. The inherited propagator bounds each half-tick block error below $10^{-19}$; the first-hit recurrence sums fewer than $2^{11}$ error terms times a factor of eleven. The angular error of a complex mean in a disk of radius $r<|Z_{\rm obs}|$ is at most $\arcsin(r/|Z_{\rm obs}|)$. [Native propagator](../../full_algebra/native_precision.py) · [First-hit error recurrence](../../full_algebra/reuse/clock_history.py)

The preparation is the same finite entropy equilibrium used for the gravity comparison. Both input devices occupy fixed addresses and there is no registered meeting of independent moving clocks here. A proper timelike-frame measurement would have to carry a hand along a trajectory, establish its coincidence with another observer, verify the reception appointments locally, and calibrate relative motion there before comparing regional stress or area. A separate travelling-clock experiment achieves material handoff and reunion on a different routed source. [Travelling clocks and reciprocal pulse timings](../../../../../model/experiments/interval_echo_clock/README.md)

## Does the arriving probe itself carry a clock hand?

In this flight a single vacancy or excess of fixed charge $q$ occupies one of 960 ports. At a specified arrival port $j$ there is exactly one load word for the mobile: $1+q$ at $j$ and one everywhere else. The microscopic basis in each conserved reversal sector is $|j\rangle\otimes|c_L,c_R\rangle$. The two coin bits belong to the **parked source seams** $(0,1)$ and $(38,39)$. Blank hops act on position times the identity on both heads. In particular the actual gateway block is $2H_{38,1}=-I_4$ for both charges and both left reversals. Each marked contact acts on the coin at its own seam. [Native sector embedding and generator](../../../../../model/experiments/interval_echo_einstein_station/mobile_collar_source/measure.py)

For any mobile-only effect $E$ and sharp arrival projector $P_j=|j\rangle\langle j|$, one has $P_jEP_j=\langle j|E|j\rangle P_j$. Thus after the hit its conditional expectation is independent of which path brought the load, its dwell time and its correlations with the heads. There is no second internal mobile state that could accumulate a readable relative phase during this flight. The hit *probability* and the conditional head state can depend on all three. A local detector acting jointly on the arrived load and a parked head reads that larger system.

On an isolated blank dimer the two position eigenmodes of exchange have energies 8 and 9, with a phase detectable by retaining both orbitals. Compression of either energy projector to either endpoint gives $P_j/2$. Reading one sharp endpoint therefore erases that internal-phase interpretation. The two possible probe charges travel in separate total-load sectors; the source and terminal read admit no coherence measurement between them. A clock encoded in a coherently collected set of ports or in a multi-load object would use a different flight and arrival measurement. [Independent native fiber checks](test_measure.py) · [Clock phase carried by a moving pair on a different source](../../../../../model/candidate_laws/interval_echo/openwave/row_03_clock/moving_pair/README.md)

Run from the RER root with one numerical thread.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 PYTHONPATH=code:.
python -B -W error -m matched_einstein_limit.source_implication.material_clocks.two_way_phase.measure
python -m pytest -q code/matched_einstein_limit/source_implication/material_clocks/two_way_phase
~~~
