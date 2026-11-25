# Native complex instruments with retained state updates

One common native control program implements the eight committed sharp phase-table slots and an unsharp version of each. The weak effects have eigenvalues 3/4 and 1/4. Their conditional maps are the positive-root Lüders maps, rather than only matching first-read probabilities. Reading an extra address bit gives the same grouped effects with different disturbance. The source, neutral pointers and final read architecture are shared across contexts and charge signs. [Native compilation](measure.py) · [Complete acquisitions](results.json) · [Operator and readback tests](test_instruments.py)

## The table and its scope

The [committed phase receipt](../../../../../code/phase_operation_producer/PHASE_OPERATION_RECEIPT.v1.json) has eight slots but four distinct projectors:

| Slots | Projector |
| --- | --- |
| Diagonal, conjugated 0 and 1 | \(|0\rangle\langle0|\) |
| Conjugated 2 and 3 | Projector on \((1/2,-\sqrt3/2)\) |
| Conjugated 4 and 5 | Projector on \((1/2,\sqrt3/2)\) |
| Phase | \((I+Y)/2\), projector on \((1,i)/\sqrt2\) |

For each sharp projector \(P_c\), the unsharp yes effect is \(E_c=I/4+P_c/2\), with complement \(I-E_c\). The table spans all four Hermitian qubit coordinates, including the imaginary quadrature. The literal projectors are comparison targets. The chosen contexts and weak strength are apparatus settings; this construction does not derive their selection, PR-04 or the probability rule.

The input is one vacancy or excess on logical rails 0 and 1 of a unit-filled data cell on the registered 20-cell board. Its sign is unknown to the controller. Initially unit-filled rails 2 and 3 provide the dilation. Four fixed neutral material pointer services, their hosts and the balancing-reference port 239 are exactly those of the [full-cell address instrument](../full_address/README.md). Source data port 7 is an additional unit-filled host during the output phase correction. All other contacts are parked. Local pulse timing, protected isolation, feedback/request implementation and the final sharp occupation semantics are supplied.

An arbitrary coherent charge-balanced input is allowed. Source sign is an audit-block index only; no gate word is selected by it. The same phases and full branch maps on the two signs are checked, including cross-sign Choi entries and untouched references.

## Explicit native compilation

An isolated blank source contact has block \(17I/2-X/2\) for either sign. A reference at degree-five data port 0 for \(\pi/16\) gives the same relative phase \(i\) on a marked rail for both signs; \(3\pi/16\) gives its inverse. These produce real rotations through \(\pm\pi/3\) with a contact and two reference pulses. The phase frame is a contact for \(3\pi/2\), whose first row resolves +Y. The receipt contains every frame word. No separate sign control or abstract unitary-reachability step is used.

After the frame \(U_c\), couple rails 0 and 2 for \(\pi/3\), and rails 1 and 3 for \(2\pi/3\). The logical/ancilla routing is \(j\mapsto j+2y\). The nonincident pairs are implemented by actual native swaps to an adjacent contact and their exact reverses. All transit ports are restored.

Writing \(\theta_0=\pi/6\), \(\theta_1=\pi/3\), the two output blocks before phase correction are

\[
D_0=\operatorname{diag}(e^{-17i\theta_j}\cos\theta_j),\qquad
D_1=i\operatorname{diag}(e^{-17i\theta_j}\sin\theta_j).
\]

Their relative phases matter for the instrument. A fully refocused native contact to unit host 7, with total dwell \(5\pi/51\), applies relative phase \(e^{-5\pi i/6}\) on source rail 0 for either sign. This removes the relative phase of both blocks. Conditional native permutations return outcome-1 rails 2 and 3 to logical rails 0 and 1; all other rails are restored. The inverse native frame gives, up to retained outcome scalars,

\[
K_{c,0}=\sqrt{E_c},\qquad K_{c,1}=\sqrt{I-E_c},
\]

\[
\sqrt{E_c}=\tfrac{\sqrt3}{2}P_c+\tfrac12(I-P_c),\qquad
\sqrt{I-E_c}=\tfrac12P_c+\tfrac{\sqrt3}{2}(I-P_c).
\]

The two physical scalars are \(h^{14}e^{i\pi/3}\) and \(ih^{14}e^{i\pi/3}\), with \(h=e^{-133\pi i/260}\) from the material pointer. They are common across charge signs. The receipt retains them; they cancel only after forming the labeled outcome maps \(K\rho K^*\).

For the sharp table use dwell zero on the first dilation pair and \(\pi\) on the second. The same native record and recovery steps give \(P_c\rho P_c\) and \((I-P_c)\rho(I-P_c)\). Thus all eight literal sharp slots also have explicit finite instruments on this apparatus.

## The pointer read determines the disturbance

Encode the source address in the four material bits. For the coarse instrument, reverse the six calls that wrote logical bit zero before reading anything. This restores that pointer coherently and erases its information about the logical rail. Keep the outcome bit and the two high address bits; the latter retain possible escape from the four dilation rails as explicit rejection records. The operation is checked on all 384 physical signed address/register inputs.

Read all four material bits and hosts with the same complete seventeen-word alphabet as the address protocol. Ideal coarse outcomes are raw words 0 and 2. Fine readout keeps bit zero and therefore exposes raw words \(r=j+2y\) in 0 through 3. All other words, including the malformed-pointer/host complement, remain in the acquisition. Nothing is discarded or accepted conditionally.

The four fine effects are

| Raw word | Effect |
| --- | --- |
| 0 | \(3P_c/4\) |
| 1 | \((I-P_c)/4\) |
| 2 | \(P_c/4\) |
| 3 | \(3(I-P_c)/4\) |

Fine words {0,1} and {2,3} have the same yes/no effects as the coarse unsharp instrument. Their conditional maps are sums of separate Kraus maps and dephase in the \(P_c\) basis. Summing amplitudes would give the wrong instrument.

The alternative grouping {0,3} versus {1,2} gives effects \(3I/4\) and \(I/4\). A direct native coin uses equal \(\pi/3\) couplings on the two dilation pairs and leaves the logical input unchanged in either outcome. It has exactly those same effects. The grouped fine map instead dephases the source. These are actual different programs and different poststates for the same effect; neither effect data nor a context-normalized table selects the map.

## A second native read tests the poststate

Prepare the logical +X state with a native half contact and reference pulse. In the phase context, first apply the coarse weak, fine or direct-coin instrument. Retain the raw first record. Public-record-conditioned \(\pi\) swaps reset the neutral pointer orientations. A common native X analyzer and the full address encoder then produce the second record. The same physical pointer architecture is reused. All \(17\times17\) joint record bins are retained.

| First instrument | Ideal probability that the X analyzer returns plus |
| --- | ---: |
| Coarse weak Lüders | \((2+\sqrt3)/4=0.9330127\ldots\) |
| Fine weak, grouped afterward | \(1/2\) |
| Direct coin | 1 |

For example, the vacancy banks gave 30,624, 16,559 and 32,768 plus records out of 32,768 starts. Excess banks gave 30,611, 16,437 and 32,768. The grouped fine event {0,3} has joint X-plus probability 3/8, while the direct coin's same-effect yes event gives 3/4. These readings distinguish the conditional maps rather than only their effects.

An orientation control replaces the phase frame's \(3\pi/2\) contact by \(\pi/2\). On the same +Y input its unsharp yes probability changes from 3/4 to 1/4. The alphabet, source, routes and terminal reader are unchanged; this is an instrument-setting control, without an emergence claim for a new physical law.

The primary bank has 250 cohorts of 4,096 starts, totaling 1,024,000. Six poststate cohorts add 196,608 starts. Every source preparation has an explicit common native word. Sampling uses exact ideal-law trees after comparison with unnormalized native propagation centers. No clipping or renormalization hides drift, and no count represents literal native-event calls. This is a construction replay, without a prospective physical-data claim.

## One common calibrated input, not a historical source transfer

The [tagged material preparation](../tagged_preparation/README.md) supplies one common reduced input for all eight contexts. A calibrated native source contact gives address amplitudes \(\sqrt{111/179}\) and \(i\sqrt{68/179}\). A fixed neutral tag at ports 76 and 77 records the source position coherently through an explicit routed native word. Tracing the unread tag gives \(\operatorname{diag}(111/179,68/179)\) for both signs. The measurement word avoids the tag, its host 79 and the charge reference.

The historical fibre census supplies that calibration ratio. The material preparation is a fresh software construction; it is not identified with the old fibre trajectory. Its 64 joint cohorts reproduce every sharp benchmark weight and the corresponding weak weight on the same prepared input. Skipping the tag leaves the real-context readings fixed but changes sharp +Y from 1/2 to 0.98535881 and weak +Y from 1/2 to 0.74267940. The source/readback semantics, calibrated fresh input and clock remain supplied. [Preparation and joint source/record acquisition](../tagged_preparation/measure.py)

## Verification and limits

The largest native-versus-analytic branch operator discrepancy is \(2.46\times10^{-14}\); the same-word sign blocks differ by at most \(1.29\times10^{-14}\). Exact symbolic checks establish the table, positive square roots and complex rank. Full branch and reference checks test trace preservation and all output states. Their numerical residuals are regression diagnostics, not physical measurement precision.

Each accepted first-measurement history publishes its complete native wait list and norm budget. The largest weak history has 972 waits and generator-norm budget 31,624. Preparation half-trace error \(\epsilon_p\) and absolute wait error \(\delta\) give complete first-measurement classical/quantum output error at most \(\min(1,\epsilon_p+31624\delta)\). The readback experiment appends its published reset/analyzer/read word and has a separate longer budget. The bound assumes the final sharp measurement rule; physical pulse or pointer errors require their own calibration.

This realizes finite complex sharp and nonprojective instruments, alternative record coarse grainings and their operational state updates on the specified source. It supplies neither a law-selected measurement menu nor a derivation of Born probabilities, full-grammar additivity, physical clocks, laboratory gate fidelity or external producer/run custody. PR-03, PR-64 and PR-65 retain those broader obligations.

```sh
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m model.experiments.interval_echo_consensus.signed_effects.phase_table.acquire
OPENBLAS_NUM_THREADS=1 .venv/bin/python -m pytest -q model/experiments/interval_echo_consensus/signed_effects/phase_table model/experiments/interval_echo_consensus/signed_effects/tagged_preparation
```

The operator API is `branch_matrices(context, mode="sharp"|"weak"|"coin", sign=...)`; sign chooses a block for audit and is absent from the physical word. `raw_branch_maps` retains all seventeen source-output maps. The acquisition's `--write` exclusively creates a missing receipt; replay preserves saved evidence and compares discrete fields exactly with bounded floating serialization tolerance.
