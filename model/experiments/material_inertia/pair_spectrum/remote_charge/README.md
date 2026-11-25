# A different owner reads the charge of the moving pair

A three-dimer material guide carries the parity and sign of a regional charge to a receiver on the far side of a cell boundary. The receiver reads just the two loads on its last dimer. On the registered law two measurement settings separate an excess from a vacancy and assign the same neutral response to the ordinary host and a compensated pair. With the first bridge cut, all four sources give a fair coin at the receiver. This adds an operational long-range read of the existing material charge to OW-19. The experiment uses a scheduled polarizable guide and native material standards. [Frozen protocol](PROTOCOL.md) · [Source and records](results.json) · [OpenWave criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#simplest-test-per-criterion)

## The measured response

OpenWave lists M10 dark matter as untested. Its stated test asks for a stable neutral massive state. On the same 20-cell sector already studied in OPH, a freely moving excess–vacancy pair occupies a separated 24-mode cluster on a unit-filled host. The previously measured charge of a region is the sum of load departures from this host. Three additional guarded dipoles now mediate the response between this source and a separate observer. Only the head of the instrument sees the sampled ports. Every write route avoids the receiving owner. Once writing finishes, the only way for the receiver to learn the answer is through the guide interactions. [Earlier stability and transport](../README.md) · [Regional calibration](../native_source_probe/regional_census/README.md) · [Polarizable guide construction](../../../../candidate_laws/interval_echo/openwave/row_20_electric/README.md)

Each calibration stratum has 65,536 new simulated launches. The table gives the count of far-end word 20 at the second tick. Every other launch gave word 02. The two settings use separate preparations. The port ordering is fixed before the simulation.

| Material at the source | Regional charge | Parity setting | Sign setting |
|---|---:|---:|---:|
| Vacancy | −1 | 65,536 | 65,536 |
| Ordinary host | 0 | 0 | 32,774 |
| Excess | +1 | 65,536 | 0 |
| Internal excess and vacancy | 0 | 0 | 32,687 |

At ticks zero and one the probability is exactly one half for each setting and source. After the full second tick a signed source changes the far reading. Cutting the first bridge produced 32,515 and 32,886 word-20 outcomes on host and excess, respectively, in the parity setting. Halving both bridge dwells gave 16,346 and 49,387, against probabilities one quarter and three quarters. Removing contrast from every blank seam, while retaining the private reference operations, preparations and schedule, gave 32,886 and 32,917. The complete receipt also retains the separate control in which only the guide flight loses contrast. All 120 calibration banks satisfy the simultaneous 99% Hoeffding allowance. [Native word and detector programs](measure.py) · [Independent checks](test_remote.py)

## Why a distant read works

Use the two basis states 20 and 02 of each protected dimer, with $Z=+1$ on 20. Each internal contact gives the matrix $(32+1/2)I-X/2$. On a link from the second rail of one dimer to the first rail of the next, crossing exchange leaves the guarded sector and its projection is zero. Its retained diagonal is $J(I+ZZ)$ with $J=65/4$. The two ports within each dimer have equal degree so the degree-reference cost is a scalar. Each dimer starts in its even free mode. The deposited six-port layouts provide three mutually disjoint candidate guides for each of the 60 owner regions. The first all-unit guide can be used without a guide-failure cut: two material defects can spoil no more than two of the three candidates. Reading the initial guide occupations commutes with the regional charge. Replacing the six units by three neutral dipoles conserves total load; the supply and replacement are part of the measuring apparatus. [Board and route construction](geometry.py)

The source writes via refocused diagonal contacts. One target of load $n=1+q$ imparts a relative $Z$ rotation $\exp[-i(65t/4)qZ]$ on the head. Its other phase depends on the target alone. The complete three-qutrit pulse verifies this action and returns every extremal and ordinary operand. Consecutive target pulses therefore write $Q=\sum q$ without looking up a classical charge label. At the parity setting, $t=2\pi/65$, followed by a fixed private phase on the head's first rail. At the sign setting, $t=\pi/65$. A final head half-contact converts the resulting phase into a signed orientation. The charge-write, tomography and travel clocks are supplied meter operations. [Source refocusing and its full ternary action](../native_source_probe/PROOF.md)

One subsequent tick opens both bridges for $\pi/65$, then each dimer for $\pi/2$. Apart from scalar phases write this as $U=\prod_j\exp(i\pi X_j/4)\prod_i\exp(-i\pi Z_i Z_{i+1}/4)$. The receiving orientation satisfies $U^\dagger Z_3 U=-Z_2 X_3$, then $(U^2)^\dagger Z_3 U^2=Z_1 X_2 X_3$. The first of these expressions has no head support. At the second tick, with dimers two and three initialized in their even modes, the far read measures the head orientation. Direct native matrices give the same identities. Thus the probabilities of word 20 at that appointment are $|Q|$ in the parity setting and $(1-Q)/2$ in the sign setting for the candidate alphabet $Q\in\{-1,0,1\}$. These are statements about marginals. The terminal operations can still disturb coherent states within a charge fibre.

## Transfer to the free band

For sharp inputs the exact far parity effect is the projector on $Q\ne0$. The charge write is diagonal in the target loads, so the same marginal effect applies to all coherent inputs, mixtures and untouched quantum references. Summing over guide eligibility returns this projector on the initial ensemble. Consequently the original measured parity counts and their uncertainty predict this remote parity probability on the same freely arrived preparation. At one million exchange units, each of 180 native panels had 67,108,864 original launches including the initial spectral-filter refusals. The largest simultaneous 99% upper limit from their charge counts is 0.00099103. Including the numerical source allowance gives a worst-case limit below **0.000996**. The previous operator bound on different owners gives **0.001070534** for every state of the 24-mode band at every free exposure and every region with the stated external guide. At the three launch-owner regions the matched-marginal controls instead gave at least 27.28% charged. The receipt transfers the original 540 native and control histograms and retains the signed source counts and refusals. The new guide calibration used fresh draws; the pair estimates reuse the archived data and have no additional statistical power. [Original readout counts and numerical error bounds](../native_source_probe/regional_census/results.json)

The polarization read establishes transmission of a signed response calibrated against independent material excesses and vacancies. An identification of those standards with electric charge in a Maxwell field requires a further physical measurement. The original independent kinetic test on this same band also distinguishes different inferred masses at different spatial profiles. A single physical rest mass and a full OW-19 submission remain open. [Kinetic measurement and failed cross-profile prediction](../sum_rule_inertia/README.md)

From the repository root, run both commands with single-threaded numerical libraries.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
nice -n 15 .venv/bin/python -B -O -W error -m model.experiments.material_inertia.pair_spectrum.remote_charge.measure
nice -n 15 .venv/bin/python -m pytest -q model/experiments/material_inertia/pair_spectrum/remote_charge
~~~
