# One marked commutator probes the five-dimensional sector

Current measurements around one marked pentagon at a time resolve both triplet shapes. To separate the four viable B14 families, an observer can instead compare two **orders of perturbing an independently prepared probe**. There is an exact marked pair for which the predicted result distinguishes the two families with an active five-dimensional sector irrespective of the unknown coupling strength.

Let $J$ exchange antipodes and $P_5=(I+J)/2-\mathbf1\mathbf1^T/12$. Drive $u=P_5e_0$ and $v=P_5e_1$. In each drive the marked port and its antipode have weight $5/12$, and the ten others have weight $-1/12$. For $[u,v]$, the fixed B14 Reynolds tensors give the following coordinates in the registered port chart. Sites 4 and 8 have the same labels in both charts.

| Bracket family | Output at 4 | Output at 8 | Ratio 8/4 |
|---|---:|---:|---:|
| P1 or P2 | 0 | 0 | undefined |
| F2 | $-b$ | $b(-1+\sqrt5)/2$ | $(1-\sqrt5)/2$ |
| G2 | $-b$ | $b(-1-\sqrt5)/2$ | $(1+\sqrt5)/2$ |

Here $b$ is the actual five-to-triplet channel coefficient in the pinned convention. Its sign and magnitude are free subject to the compact stratum. The two outputs are signed coordinates of the same bracket. Two independently calibrated channels at 4 and 8 suffice if their relative gain and polarity are known. On F2 the entire output is in the plus triplet and on G2 it is in the minus triplet. The code contracts the published exact tensors with the two physical source weights. On the compact families the five-five channel is linear in $b$ alone. Division by the nonzero sample coefficient gives its unit vector for the entire family. [Exact contraction, native commutators and all twelve coordinates](five_commutator.py) · [Frozen family equations](../../../../../code/b14_jacobi/classify.py)

## What the present instruments do

The two separately measured sets of pentagon currents admit exact operator checks on the one-load Hilbert space. Give each directed edge $a\to b$ the real skew matrix $E_{ba}-E_{ab}$. Write $K_q$ for the sum on the ring at $q$. The bare signed terminal mean has observable $iK_q/5$. The late odd sorter has observable $i(I-J)K_q(I-J)/10$ on the original source. These formulas follow from the calibrated endpoint read. Their commutators were calculated exactly. The existing count banks contain expectations of single currents.

* **Bare bond meters.** The twelve skew matrices span 11 directions. Even adjoining $[K_0,K_1]$ increases this to 12. All first commutators span 44. For the specified five-sector pair and the $1/5$ mean normalization, the skew-generator commutator has entry $(4,8)=-1/100$. Every bare ring operator has zero there because 4 and 8 share no seam. Thus identifying the bare ring operators directly with a linear representation of a closed B14 port bracket fails on this pair.
* **Odd-sorted meters.** Their span has dimension six, reproduces itself under brackets, kills the even subspace and commutes with native adjacency. On the odd subspace the adjacency has two three-dimensional eigenspaces, with eigenvalues $\pm\sqrt5$. The currents therefore fill the entire six-dimensional block algebra $\mathfrak{so}(3)\oplus\mathfrak{so}(3)$. This is the P2 **derived algebra type for these observable operators**. Both five drives are identically zero on this instrument. A doublet restricted to one block sees only one of the two factors.
* **Local occupation phases.** Driving $u$ and $v$ through private port phase/degree controls gives two diagonal matrices, which commute identically for any weights. A null with those controls alone has no power to distinguish P1/P2 from F2/G2. Free propagation inserted between the controls can itself create an order effect and needs a matched clock control.

## Independent probe protocol

1. First bind each signed port drive to a repeatable response of a separate freely prepared carrier. Calibrate polarity and relative gain using individual drives, known positive and negative controls, and signed standard probes. Use equivalent independently mounted readers at 4 and 8. Check sensitivity to both triplet sectors and recover held linear combinations of the port drives. An unresponsive carrier or a blind receiver is an inconclusive run.
2. Apply the two weighted drives in both orders to fresh copies using equal timings. Flip either sign, repeat at smaller pulse strengths, and run both one-drive controls and an equal-time parked pair. Use an independent reference arm and both analyzer phases to read the signed difference, with right, left, empty and refusal retained per launch. Swapping the order reverses the second-order term. Isolate its bilinear coefficient by alternating the two independent signs and checking the scaling with their product.
3. Tomograph the resulting infinitesimal action against the calibrated twelve single-drive responses. Check residuals on independent preparations, including nonadjacent ports. Only a closed, independently observable port response supports the table. Use the signed 4 and 8 output ratio once the denominator resolves and check the other ten coordinates on held preparations. If the probe observes only a quotient, label the result for that quotient. The bare current's 4--8 entry supplies a direct out-of-span control. Prepare port 8 and recombine 8 with 4 in a calibrated reference analyzer on both pulse orders.
4. A null constrains the product of coupling scale and detector sensitivity. Without a measured lower bound on that product there is no finite-sensitivity exclusion of F2/G2. If the five response is excluded over a declared sensitivity range, test self-brackets of noncollinear marks 0 and 1 in each triplet. One active triple is the P1 prediction and two are the P2 prediction, subject to the same sensitivity and kernel checks.

No independent-probe commutator counts have been taken here. The exact calculations can be reproduced from the repository root with the following command.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python -B -O -Werror -m model.experiments.interval_echo_gauge.matter_class.ring_current.five_commutator
~~~
