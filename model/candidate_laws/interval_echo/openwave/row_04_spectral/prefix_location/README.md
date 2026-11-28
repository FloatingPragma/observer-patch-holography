# Where can an observer read the two spectroscopy preparations?

An observer with access to all material ports and every seam memory except three has identical outcome probabilities on the two ideal persistent preparations. With equal prior probabilities the best success in telling them apart on a single launch is exactly **50%**. Reading the direction of the single seam 0–1 distinguishes them exactly. Both statements hold after every undisturbed free wait on the twenty-cell spectroscopy material. The apparatus supplies the same initial arrow chart and withholds the source program label. [Original selection and signed detector](../provenance/README.md) · [Executable comparison](verify.py)

| Original pulse column | Earlier path to port 5 | The same two terminal alternatives |
| --- | --- | --- |
| 14 | 0 → 5 | 5 → 9 and 5 → 11, with opposite amplitudes |
| 15 | 0 → 1 → 5 | 5 → 9 and 5 → 11, with opposite amplitudes |

The physical material records are the occupied port and the seam heads reversed along its path. On either tip, the two preparations differ on precisely 0–1, 0–5 and 1–5. Denote those three head registers by $T$. After factoring the fixed collar records, the two ideal preparations have the form $|t_{14}\rangle_T|\chi\rangle$ and $|t_{15}\rangle_T|\chi\rangle$. The head words are orthogonal. The same $|\chi\rangle$ contains the load and the two final hop alternatives, including their relative minus sign. For any terminal effect $I_T\otimes E$, the compression to the two-column space is $\langle\chi|E|\chi\rangle I_2$. Thus measurements supported on the complement also give identical probabilities for arbitrary coherent combinations of the columns. The two exact integer eigen-equations give the same free phase, so the identity holds at any terminal exposure. A detector that later interacts with $T$ measures a different effect.

The verifier derives the physical partition from the reachable words, checks the reduced diagonal and off-diagonal operators in integer arithmetic, and checks the distinguishing head on the same source. Each tip has occupation probability one half and the total load is one. Every free hopping-current operator between reachable words has zero compression to this ideal two-column space. This statement concerns material transport by those hops. An intrinsic magnetic response would require its own coupling and meter.

The release inputs are the board, pulse words and clock banks bound in OPH revision `e2810645c67ca634a05ea005f84bf95dc922d66d`. A second calculation used the release replayer’s independently constructed history graph, diagonalized its free generator and traced the physical heads by explicitly summing pairs of words. Restoring the head at 0–1 gave unit ideal discrimination. With all three heads withheld its largest success on the 80 finite-pulse comparisons agreed within $3\times10^{-15}$. [Release replayer](../../../../../../openwave/interval_echo/row_04/spectroscopy_release/headline.py)

On the 32 finite-accuracy pulse preparations the maximum calculated discrimination success outside $T$ is **50.01496%** rounded upward, conditional on commissioning. This compares the two programs at each of sixteen resolved strengths and five waits, including launch. The common commissioning refusal conveys no information. A conservative bound at every free wait follows directly from the signed launch measurement. If both preparations have original-start ideal signed weight at least $F$, their trace distances from their respective ideal inputs are at most $\sqrt{1-F}$. After discarding $T$ the ideal outputs coincide, so the equal-prior success is at most $\tfrac12+\sqrt{1-F}$. The archived calculated weights give an upper bound of **50.071%**. The simulated witness counts give **57.621%** at simultaneous 99% coverage. The zero-strength clock remains inconclusive. These are probability calculations and an inference from the archived counts. [Input bounds and acquisition](../DONOR_RELEASE.md)

For electron identification, the preparation distinction is demonstrably stored on $T$ in this experiment. The two-tip signed excitation is a separate observable common to both inputs. An experiment could transfer a record from $T$ to a remote charged arrival. Calibrating the remote detector on a fresh, initially independent receiving material and reading it without information from the launch apparatus would test that transfer. The present experiment establishes free persistence at its original site. [Electron identity and the competing hypotheses](../../../../../experiments/interval_echo_angular/charged_orientation/three_arm_spin/ELECTRON_IDENTITY.md)

From the repository root run

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python -B -W error -m model.candidate_laws.interval_echo.openwave.row_04_spectral.prefix_location.verify
~~~
