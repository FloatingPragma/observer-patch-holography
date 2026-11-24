# Can a stationary load medium change the rate of both seam clocks?

An independently prepared blank seam has an initial exchange beat of one in every surrounding state. In a uniform stationary binary medium at filling $f$, its second initial phase rate on a degree-five pair is $50-68f$. A thermal medium with occupation coherence can read both the exchange and contrast energies through its two-port modular operator. The two responses per unit of native energy differ. [Native law](../../../candidate_laws/contrast_memory_spherical_v1/law.py) · [Calculation](measure.py) · [Probability records](results.json)

## The instantaneous clock comparison

Let the two ports be A and B, with the measured seam memory blank and independent of the environment at launch. The initial derivative of the reduced density is the commutator with the environmental expectation of the Hamiltonian. Every other contact touches at most one of A and B. The effective first-order generator therefore has the form

~~~math
h_{AB}+L_A\otimes I+I\otimes L_B+c I .
~~~

The outside state may contain correlations among loads and memories. Between the symmetric and antisymmetric words $s,a=(|10\rangle\pm|01\rangle)/\sqrt2$, all one-port contributions to the difference of diagonal matrix elements cancel. The exchange term gives one. Operationally, launch on A and read the quadrature $2\,\mathrm{Im}(\psi_{10}^*\psi_{01})$, averaged over every exterior outcome. Its initial slope is one because the two-port hopping element is (-1/2). A half contact followed by the endpoint occupancy difference reads this quadrature. This argument concerns the initial rate of a fresh local probe. Later evolution couples its state to the environment.

The reference medium has a unit load on each port independently with probability $f$, with all memories blank. The probability of each global word depends only on its conserved load number, so the density is stationary. For a degree-five pair, each port has four other contacts and one reference. On each outside contact, averaging the native exchange and contrast gives a one-port energy difference $(\kappa+1/2)(1-2f)$. The reference adds $\kappa$. The second probe interferes the empty pair and the symmetric one-load word using a parked isolated reference rail to keep total load fixed. With zero reference energy its initial quadrature slope is $2\kappa+4(\kappa+1/2)(1-2f)$. A reference offset common to every preparation subtracts the same constant from these frequencies.

The population calculation uses seam (0,5) on the complete 12-port cell, all 4,096 binary words, and the full native degree correction of the standalone cell. Both selected ports also have degree five on the 20-cell federation, where the initial-rate derivation applies. The two finite-exposure readouts use a supplied launch, release time and terminal quadrature. For the second probe the analyzer measures $i(|v\rangle\langle s|-|s\rangle\langle v|)$, where $v$ places the load on the reference. Orthogonal outcomes score zero. The calculation gives the probabilities' signed differences, without sampling or conditioning on successful return.

| Filling | Exchange read at $t=1/32$ | Parked-rail-to-symmetric read | Neighbour's occupancy response |
| ---: | ---: | ---: | ---: |
| 1/4 | 0.030149 | 0.775045 | 0.00023976 |
| 1/2 | 0.029802 | 0.414912 | 0.00023839 |
| 3/4 | 0.030149 | −0.036775 | 0.00023976 |

The last column compares two fresh releases with the launch port A occupied or empty and all other ports drawn at the given filling. Its leading term at B is $t^2/4$: a neighbouring occupied and empty port exchange with amplitude $it/2$, and all other diagonal-population contributions cancel in the comparison. This gives a local propagating density response. The exposure uses native exchange units. These data supply neither a macroscopic dispersion relation nor a material proper-time calibration. The exchange-only control keeps the alphabet and contacts and moves the second initial slope to $2-4f$.

## What a regional modular measurement reads

For any faithful regional density diagonal in sharp binary occupations, $K_B=-\log\rho_B$ is diagonal. Moving probability from $s$ to $a$ changes only off-diagonal entries. Its modular slope vanishes while the native bond energy increases by one per unit transfer. Both single-port marginals are fixed. Choosing the tangent tensored with an independent exterior density keeps every other bond energy fixed, and disjoint interior transfers can be summed. Positivity holds at sufficiently small amplitude for a faithful global preparation. For a fixed-number family use a spectator word of the requisite load in both arms and restrict to the occupied support. On a connected blank federation every stationary sharp-diagonal binary load state has constant probabilities within each fixed-number sector: the commutator with the Hamiltonian equates the probabilities of any two words one neighbouring swap apart, and those swaps connect the sector. The diagonal-reference obstruction also applies to their regional densities on the common support.

For a reference with regional coherence, take the native stationary grand Gibbs density on the standalone cell with chemical potential $\mu=\kappa=8$. Load reflection on this degree-five graph leaves $H-\mu N$ invariant and gives half filling. The code diagonalizes all 13 load-number blocks and traces out ten ports. A second transfer compares an equal mixture of 00 and 11 to $s$. Both endpoint marginals are fixed and the native bond energy rises by eight. A sufficiently small signed perturbation of the full-rank thermal density realizes each tangent, with a normalized spectator density on the other ten ports.

| Inverse temperature | Smallest regional probability | Modular slope per unit exchange energy | Modular slope per unit energy, aligned mixture to symmetric |
| ---: | ---: | ---: | ---: |
| 0 | 0.25 | 0 | 0 |
| 1/32 | 0.20156 | 0.030944 | 0.043854 |
| 1/16 | 0.11099 | 0.059356 | 0.148252 |
| 1/8 | 0.00581 | 0.125661 | 0.539646 |

At inverse temperature $1/16$ the two-port modular observable detects the exchange transfer. Its conversion from energy differs between the two preparations by a factor of 2.50. These are thermal reduced-state probabilities evaluated from the declared equilibrium. The two-port calculation has no proper radius. A conformal vacuum interpretation requires distinct evidence on a refining observer algebra and calibrated balls.

The discrepancy also has an analytic high-temperature test. At half filling $H-\mu N$ is, up to a constant, a sum of bonds $h'_{ij}=\sum_{\alpha}J_\alpha\sigma_i^\alpha\sigma_j^\alpha$ with $J_x=J_y=-1/4$ and $J_z=-(1+2\kappa)/4$. Expand the normalized partial trace of the exponential to second order and take its negative logarithm. The square of the interior term cancels in the logarithm. Single external bonds have scalar traced squares. Two different exterior bonds contribute a nonconstant term only if they meet a common neighbour of A and B, since a Pauli matrix on an unpaired traced site has trace zero. With $m=2$ common neighbours this gives

~~~math
K_{AB}=cI+\beta h'_{AB}
-m\beta^2\sum_{\alpha} J_\alpha^2\sigma_A^\alpha\sigma_B^\alpha+O(\beta^3).
~~~

The swap transfer changes the expectation of $XX+YY$ by $-4$ and of $ZZ$ by zero. The aligned-to-symmetric transfer changes them by $2$ and $-2$. After division by the respective native energy changes, the second modular slope minus the first is $m(1+2\kappa)\beta^2/4+O(\beta^3)$. Its leading coefficient is $8.5$. Analyticity of the finite full-rank reduction excludes a common conversion on this region for every sufficiently small positive inverse temperature.

### Moving the same seam into the interior

At $\beta=1/16$ trace the same 12-port density onto the edge, the edge with both common neighbours, and the union of the edge and all its nearest neighbours. Use the same full-cell tangent for each comparison, the signed two-port transfer tensored with the normalized identity on the other ten ports. Its restriction to any of the regions is the same transfer tensored with the normalized identity on the included spectators. The native energy changes are one and eight at every region size. Each tangent leaves the mean load number fixed. The logarithm is taken after the partial trace. This computes equilibrium population responses for specified graph regions.

| Ports retained | Swap per unit energy | Aligned to symmetric per unit energy | Second minus first |
| --- | ---: | ---: | ---: |
| 2, endpoints | 0.05935639 | 0.14825170 | +0.08889531 |
| 4, both common neighbours included | 0.05755803 | 0.07015837 | +0.01260033 |
| 8, every neighbour of either endpoint included | 0.06253325 | 0.06250043 | −0.00003281 |
| 12, full cell (identity) | 0.06250000 | 0.06250000 | 0 |

The discrepancy shrinks and reverses sign. With all contacts of the excited edge inside the region, each response is within 0.054% of the full-cell coefficient. The full-cell equality follows directly from $K=\beta(H-\mu N)+cI$. It is a thermal control. These three traces have fixed mesh size and use two variations on one bond. Establishing a small-ball modular stress law requires a vacuum reference, proper balls, its radius-dependent kernel and uniform estimates on the admitted variations.

There is also a uniform check for this tangent embedding. Let $M$ be the normalized spectator trace of $K_R$ minus $\beta h_{AB}$. Subtract its identity and single-port parts to obtain the two-port interaction

~~~math
M_{\rm int}=M-\tfrac12 I_A\otimes{\rm Tr}_A M
 -\tfrac12({\rm Tr}_B M)\otimes I_B+\tfrac14({\rm Tr}\,M)I_{AB}.
~~~

For every Hermitian endpoint transfer $D$ with both one-port partial traces zero, its embedded modular-response error is $\operatorname{Tr}(D M_{\rm int})$, bounded in absolute value by $\|D\|_1\|M_{\rm int}\|_\infty$. Direct diagonalization gives norms 0.343793, 0.0318690 and $2.66704\times 10^{-5}$ on the two-, four- and eight-port regions. This bound includes imaginary coherences and changes of the endpoint number distribution at fixed marginals. A relative error for a tangent with vanishing energy is undefined. The norms are double-precision estimates. The test independently reconstructs them from the two signed transfers using the symmetries of this thermal cell.

Run from the repository root.

~~~sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m model.experiments.seam_program.collective_rate.measure
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python -m pytest -q model/experiments/seam_program/collective_rate
~~~
