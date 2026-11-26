# A measured noise floor for two overlapping repair triangles

A finite Interval Echo observer process contracts expected distance to four exact normal forms. Each block consists of eight offers on two triangles. For every starting state and every order fixed at block start, its expected terminal distance is
\[
\mathbb E[D'\mid\mathcal F_{\rm start}]
   =\frac{27}{4096}D+\frac{939}{2048}.
\]
Repeated blocks approach a mean of \(6/13\) unresolved positions. Retaining a single fault setting throughout a block changes the terminal mean from \(255/512\) to \(3/2\) when all six positions initially disagree. Both implementations have the same isolated one-offer law. [Frozen predictions](PROTOCOL.md) · [Full material branch banks](records.json.gz)

## What the observers read

Each triangle carries a protected boundary bit and an active three-bit word. Its one-hot bank address is the number of positions that disagree with the boundary. The service exposes two independently routed reads of the boundary and the active word, including a direct mismatch probe at each offered position. After the command, a different observer reads the published bank and checks its complete word against its address. In disagreement coordinates the state consists of two boundary bits and six disagreement flags. Hamming distance counts the coordinates on which two such states differ. The nearest normal state has the original boundaries and six zero flags, so its distance is precisely the number of disagreements. The [exact repair construction](../README.md) and its [termination, local diamond and normal-form proof](../Repair.lean) give this target for all 256 states. The noisy run retains the same quotient and isolated appointments.

A spare extremal dimer initially holds 20. On the native seam its two orientations differ only by exchange: up to a common phase, the amplitude crossing after a dwell \(t\) has magnitude \(\sin(t/2)\). At \(\pi/3\) the stay and crossing probabilities are \(3/4\) and \(1/4\). The material pointer resolves the two outcomes. A full exchange after a crossing restores the dimer, including for the second draw. A first stay commands correction. A first crossing and a second stay preserve the old covered bits. Two crossings command the opposite of the boundary at every covered bit. These three commands have weights 12, 3 and 1 out of 16. Their outcomes are written into two disjoint balanced dimers and read back before use. The controller performs a new draw for every nominal offer. A projected, reset clock has the same next-outcome law conditional on earlier readings, including the choice of block. A shared clock command correlates the two or three positions touched by an offer. [Native hardware, projected clock and published rank](instrument.py) · [Independent scalar operator replay](scalar.py)

The simulation retains every positive joint command branch, with its Born weight, and every queried pointer result throughout each block. There is no finite-shot sampling. The four primary banks each contain \(3^8=6{,}561\) conditional launches. Each has total weight \(16^8=2^{32}\) in integer units. The additional scalar runs use the original routed single-launch meter and actual projected 99-state native seam pulses. They check mixed disagreement words under all four protected boundary choices and selected branches of the packed observations. [Scalar pointer tapes](scalar.json.gz) · [Acquisition and independent kernel composition](measure.py)

## How a fair block contracts

Within a triangle the four labels are its three pairs and its triple. Each block contains all eight labels exactly once, exposing and acting on every class of mismatch. The supplied controller may select any of the 40,320 orders using all earlier public records before the first draw of that block. Time counts completed appointments; the receiving station reads both published words at each such time. Block endpoints are deterministic multiples of eight appointments. In the ideal setting every offer applies the exact repair. An active mismatch can survive a noisy appointment as a measured failure.

Let \(p\) be the probability of one address disagreeing. An offer covering it changes its marginal to \((1+3p)/16\). A correct address is corrupted with probability \(1/16\); a wrong address survives with probability \(1/4\). Other addresses keep their marginal. After \(k\) exposures its signed deviation from \(1/13\) is multiplied by \((3/16)^k\). Every site receives three exposures in a fair block. Summing the six marginals gives the displayed equality without assuming independent site outcomes. The additive intercept is attained at a normal input; with that intercept fixed the slope is also attained at every nonnormal input.

For a strict prefix, an initially correct site contributes at most \(313/4096\) after three exposures, or \(304/4096\) after two. At least two sites belong to any remaining offer and hence have received at most two exposures. An initially wrong site's contribution is at most its initial value plus the same allowance. This gives \(\mathbb E[D_t\mid\mathcal F_{\rm start}]\le D+465/1024\). Equality occurs on the reverse schedule, seven offers after a normal start. The constants for the fair-block theorem are therefore \(L=8,\ A=1,\ \beta=465/1024,\ \lambda=27/4096,\ \varepsilon=939/2048\). At block \(m\) the mean is exactly \(6/13+(27/4096)^m(D_0-6/13)\) for a fixed initial state. The [Lean kernel proof](Contraction.lean) composes the full three-outcome update and proves the coordinate law, fair endpoint, strict prefix bound and repeated-block recursion. The archived table checks every one of 256 offer subsets at every one of 64 flag configurations. Arbitrary independent fault models or schedule choices made using future clock outcomes require their own kernels.

## What changes under the controls

| Initial disagreement and instrument | Expected final count at the receiver |
| --- | ---: |
| None, fresh fault per offer, forward or reverse | 939/2048 |
| All six, fresh fault per offer, forward or reverse | 255/512 |
| All six, one retained fault word across eight offers | 3/2 |
| All six, only the first triangle serviced | 13551333/4194304 |
| None, preserve and corrupt commands interchanged | 2457/2048 |

From a normal start the nominal forward and reverse schedules give probabilities \(3565164681/2^{32}\) and \(3130738209/2^{32}\) for *zero* terminal mismatches. Their expected distances agree. The measured noise correlates sites within an offer, and the order changes the complete terminal law. The retained-fault comparator reuses the same physical latch at all eight appointments. With probability \(3/4\) it corrects everything. Otherwise all six initial mismatches survive. Holding one triangle out keeps its three mismatches, and the direct comparison of its published word and count address confirms the failure of coverage. [Regression checks and source verification](test_blocks.py)

PR-121 is Demonstrated for this finite simulated repair process, its stated metric, and its registered conditional fault law. The result allows any block order chosen from the public history at its start. The native error source, isolated appointments, and calibrated readers specify the implementation to which these constants apply.

Replay from the repository root with single-threaded numerical libraries.

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -O -W error -m model.experiments.interval_echo_consensus.transaction_component.confluence.noisy_blocks.measure
python -B -O -W error -m model.experiments.interval_echo_consensus.transaction_component.confluence.noisy_blocks.scalar
python -B -m pytest -q model/experiments/interval_echo_consensus/transaction_component/confluence/noisy_blocks
(cd Lean && lake env lean -j1 -DwarningAsError=true ../model/experiments/interval_echo_consensus/transaction_component/confluence/noisy_blocks/Contraction.lean)
~~~
