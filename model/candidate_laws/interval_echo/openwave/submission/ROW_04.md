# A free pair on the Interval Echo sphere

Two prepared unit loads remained detectable as an adjacent pair in their original echosahedron during free evolution on a 20-cell sphere. The local receiver detected them in 4,028 of 4,096 releases at native exposure 37. Removing the contrast interaction gave 824 detections with the same receiver. The experiment concerns an excitation of unspecified charge and spin. [Complete terminal records](row_04_receipt.json)

## The law and the instrument

In the basis of two distinct occupied ports, the native free generator, up to a scalar, is

\[
H_\kappa=\tfrac12 L_{\rm moves}-2\kappa P_{\rm adjacent}.
\]

Each occupied endpoint can exchange with an adjacent empty endpoint. These moves give the configuration-graph Laplacian in the first term. The reference-compensated contrast of the two units is 10 if their endpoints touch and 12 otherwise. Subtracting the common \(12\kappa I\) gives the second term. Endpoint swaps preserve both units and admit every separated configuration. The 630 contacts supply 28,680 two-load configurations on 240 ports. The inherited interaction and the compensated generator are part of the candidate's specified blank-memory evolution. [Interval Echo source](https://github.com/FloatingPragma/observer-patch-holography/blob/7a11a741dcc8de00a2ca826b3800909566467d62/model/candidate_laws/interval_echo/SPEC.md) · [Local generator](https://github.com/FloatingPragma/observer-patch-holography/blob/7a11a741dcc8de00a2ca826b3800909566467d62/model/candidate_laws/contrast_memory_spherical_v1/law.py) · [Sector reduction](https://github.com/FloatingPragma/observer-patch-holography/blob/7a11a741dcc8de00a2ca826b3800909566467d62/model/candidate_laws/contrast_memory_spherical_v1/sectors.py)

The designer supplies the adjacent preparation, blank collars, uniform native exposure, initialization of the local receiver, and terminal read windows. During exposure every seam stays active. At each terminal appointment a fresh release is measured once. The native echo maps each seam's joint occupancy to a readable pair flag. Copying the flag and undoing the echo realizes the terminal projector. In the binary two-load sector, at most one seam can report the pair. The home patch returns one of 30 occupied seam addresses or quiet. Quiet includes separated and out-of-patch pairs. The independently written checker builds the allowed moves from the fixture, computes the coherent evolution, and samples these 31 outcomes. The receipt contains every sampled address in draw order. [Preregistered design](ROW_04_DESIGN.md) · [Independent execution](row_04_check.py)

The same release can also be prepared by a native encounter. Start the two individual units uniformly on the 36 **separated** port pairs in cell 13. After one free exposure, measure the joint flag at ports 156 and 157. Its projector has rank one in the two-unit sector, so every detection prepares exactly the state used in the release above. The calculated probability is 0.01891573%, with 99.98108427% quiet. The user supplies the two initial units and chooses when and where to read. This low preparation yield applies before the conditional free-persistence experiment. A supplemental computation gives the probabilities for all 36 inputs and independently recovers the earlier cell-0 encounter yield. [All encounter probabilities](row_04_formation.json) · [Preparation check](row_04_formation.py)

## What the patch recorded

Both arms start at ports 156 and 157 in cell 13. The reference arm uses the candidate value \(\kappa=8\). The control uses \(\kappa=0\). Counts are per 4,096 releases. The probability columns give the corresponding finite model values.

| Exposure | Native detections | Native probability | Control detections | Control probability |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 4,048 | 0.990983 | 1,780 | 0.450002 |
| 7 | 4,052 | 0.989653 | 1,182 | 0.275769 |
| 37 | 4,028 | 0.980272 | 824 | 0.202906 |

All three native frequencies have simultaneous one-sided lower bounds above 95.48%. At the last appointment the control's upper bound is 22.98%. The design fixed a common 2.857-percentage-point Hoeffding radius for these four decisions. It required native lower bounds above 80% and the final control upper bound below 65%. Cell 13 tests a new address. The source receipt discloses that binding and the earlier cell-0 result were already known. For that historical release, the independent implementation reproduces probabilities 0.983880525 and 0.149558605 at exposure 150. [Timing and hashes](ROW_04_FREEZE.json) · [Results including quiet outcomes](row_04_receipt.json)

### Between appointments

A separate calculation after the registered test bounds the probability of a home-seam report **at every real exposure from 0 through 37 by at least 94%** for the cell-13 release. It rules out an interval of wholesale dispersal followed by reformation. Different terminal measurements use fresh preparations and do not measure a sequence of positions on a single member of the ensemble. [Window calculation](row_04_window.py) · [Window result](row_04_window.json)

For a state \(\psi(t)\), the distance to the home-seam subspace in projective Hilbert space is \(\arccos\sqrt{p(t)}\). The triangle inequality bounds its change between two times by the length of the unitary path, at most \(\Delta H\,|\Delta t|\) in native units. Here \(\Delta H=1.5\), obtained directly from the off-diagonal column of the launch state. A mesh of spacing 1/8 finds a minimum probability 0.978883326 at 34. Every time in the interval lies within 1/16 of a mesh point.

The mesh uses degree-24 Taylor steps for \(H+2I\). Its entrywise absolute matrix is symmetric with maximum row sum 14, so its spectral norm is at most 14. With step 1/8, the truncation error per step is at most \(e^{1.75}1.75^{25}/25!\). For IEEE double arithmetic the code uses \(\gamma_{128}=128u/(1-128u)\) per sparse matrix action, scaling and accumulation, with \(u=2^{-53}\). Each row contains at most 13 terms. Propagating the errors of the 24 Taylor terms and 296 steps gives a probability allowance below the reserved \(10^{-6}\). The calculation subtracts that allowance before applying the angle bound and rounds the final lower bound down to a whole percent. It also checks the three common appointments against the independently computed terminal probabilities.

The native all-time result on any admissible sphere bounds adjacency anywhere by 60% for every initially adjacent release at this coupling. The stronger 94% result uses this particular home cell, seed and finite interval. Neither result determines electromagnetic quantum numbers. [All-size binding theorem](https://github.com/FloatingPragma/observer-patch-holography/blob/7a11a741dcc8de00a2ca826b3800909566467d62/model/candidate_laws/contrast_memory_spherical_v1/PROOF.md)

## Reproduction

With Python, NumPy and SciPy installed, run these commands in this directory. They use the carrier from the pinned RER checkout. A standalone export can supply the same file with its path after `--board`.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
board=../../../../executable_model/carriers/FEDERATION_LEVEL0.npz
python -B -W error row_04_check.py --board "$board"
python -B -W error row_04_window.py --board "$board"
python -B -W error row_04_formation.py --board "$board"
python -B -W error row_04_check.py --board "$board" --remove-contrast
```

The first three commands verify the archived receipts. The last intentionally exits 1 when the mutated candidate loses local persistence. The dated review records the clean-copy timing and the failure message. The public “Read, compare, repair, record” scene illustrates local reading in a related teaching implementation. The numbers above come from the finite native generator in the submitted scripts. [Review](ROW_04_REVIEW.md) · [Teaching simulator](https://simulation.floatingpragma.io/patches?scene=repair)
