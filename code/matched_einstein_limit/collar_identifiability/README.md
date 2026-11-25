# Can a positive collar fluctuation supply the signed gravitational source?

A material observer can read a positive collar amplitude whose mean grows as the square root of the number of independent sources. The same preparations give identical complete parked seam-energy histograms. On a directly measured change of history persistence the amplitude increases while the signed modular energy of the observed records decreases. A positive conversion factor from the amplitude to that modular variation fails on this family. [Native preparation, reset and word detector](../../../model/experiments/interval_echo_einstein_station/backreaction.py) · [Finite acquisition and controls](measure.py)

## Retaining the history of a word

Each service writes two three-level collars. A fair material tag selects one of the two cyclic word laws \(P_\pm=P_0\pm D\). The [native source experiment](../../../model/experiments/interval_echo_einstein_station/backreaction.json) supplies preparation, feed-forward, marked memory and clocked word readout. Conditional on the retained tag, each reset uses a fresh independent draw. The control selects a fresh tag for every write. The three matrices have uniform one-collar marginals. In numerator units of \(1/36\), \(P_0\) has entries 4 on the diagonal, 6 on the forward cycle and 2 on the reverse cycle. \(D\) has entries \(0,-1,+1\) respectively. The two parked collar levels have single-mark energies 0, 32 and 33. Summing those energies gives the identical histogram in each arm with multiplicities \(4,8,8,4,8,4\) at totals \(0,32,33,64,65,66\). This equality holds conditional on the tag. Consequently every sequence of those energy readings on these fresh writes has the same law in both arms.

The native echo, fixed port phase and terminal seam pulse identify all three word labels. Let \(Y=1\) for a forward pair and zero for its complement. The conditional click probabilities are \(7/12,5/12\). The centered word score \(-\log(9P_0)\) has the same *conditional mean* as \(-\log(3)(Y-1/2)\) after subtraction of its unconditional mean. This follows by transferring probability between the two directed cycles. Suppress reversal exchange during the terminal seam pulse with the preparation, ports, read schedule and score fixed. The word channel then merges the two nonzero-magnitude labels, giving forward probability \(1/3\) in both arms. Each detector's common offset is obtained on its own reference law.

## Square-root scaling and scatter

Prepare \(N\) banks independently. Read \(K\) fresh words per bank, retaining its underlying tag, and evaluate
\[
 B_{N,K}=\log(3)\left|\sum_{i=1}^N \left(\frac{1}{K}\sum_{j=1}^K Y_{ij}-\frac12\right)\right|.
\]
The observer uses the word records alone. Write \(S_N\) for a sum of \(N\) independent fair signs and \(a=\log(3)/12\). Conditional expectation of the unrectified score is \(aS_N\) up to an overall sign. Convexity, the conditional variance \(35/144\) of each click, and the Lipschitz bound for absolute value give
\[
 0\leq {\mathbb E}B_{N,K}-a{\mathbb E}|S_N|
 \leq \frac{\log(3)}{12}\sqrt{35N/K}.
\]
The binomial identity \({\mathbb E}|S_N|=N2^{1-N}\binom{N-1}{\lfloor(N-1)/2\rfloor}\) follows by pairing the two tails and telescoping with \(j\binom Nj=N\binom{N-1}{j-1}\). Dividing by \(\sqrt N\) gives a limit of \(\sqrt{2/\pi}\). Fresh tags give \({\mathbb E}B_{N,K}=\log(3){\mathbb E}|S_{NK}|/(2K)\); their finite reading tends to zero at fixed \(N\). Without reversal exchange the retained and fresh mean readings coincide at every \(N,K\).

| Sources | Retained tags, 4096 reads each | Fresh tags, 4096 reads each | Retained tags in the infinite-read limit |
|---:|---:|---:|---:|
| 1 | 0.091551 | 0.006848 | 0.091551 |
| 4 | 0.142391 | 0.013696 | 0.137327 |
| 16 | 0.292966 | 0.027393 | 0.287662 |
| 64 | 0.587493 | 0.054785 | 0.582099 |

Entries are expectations in nats, calculated by convolving the complete finite count distributions. A single bank collection has intrinsic scatter even at infinite reading precision. Since \({\mathbb E}S_N^2=N\), its relative standard deviation tends to \(\sqrt{\pi/2-1}\), about 0.756. Averaging \(m\) independent collections reduces this by \(1/\sqrt m\). Reusing the identical tags across \(m\) radii preserves the relative scatter in the sum. Relating \(N\) to baryonic mass, identifying a series of physical cuts and selecting their joint history are separate requirements.

## The sign and support of the source

Use two successive emissions and hide the tag. Their fresh joint density is \(P_0\otimes P_0\). Mix the persistent preparation with the fresh one in proportion \(s\), \(0\leq s\leq1\). The resulting diagonal four-register density is \(R_s=P_0\otimes P_0+sD\otimes D\). Both partial traces of \(D\) vanish. Every proper marginal of these four word registers therefore agrees at all \(s\). Every sum of observables supported on at most three of these registers has zero response. A four-register read resolves the difference: the trace distance between \(R_1\) and \(R_0\) is \(1/72\), and the forward-click covariance increases by \(1/144\). This support bound applies to word records with the tag invisible. A retained tag jointly read with both collars can distinguish a single conditional write. For any number of independent banks, expanding the tagged product after averaging the sign leaves only even powers of \(D\) in each bank. At least two whole emitted pairs from the same bank must occur in any distinguishing untagged marginal.

For one source and two writes the positive score \(B_{1,2}\) equals \(\log(3)/2\) when the click bits agree and zero otherwise. Its expectation has derivative \(+\log(3)/144\) along \(R_s\). The modular operator of these same four records at their unperturbed reference \(R_s\) is \(-\log R_s\). On the same transfer its expectation has first variation
\[
 -{\rm Tr}\big[(D\otimes D)\log R_s\big]
 =-\frac{1}{144}\log\frac{(36+s)(4+s)}{(12-s)^2}.
\]
The nine forward-forward, nine reverse-reverse and eighteen mixed directed pairs give this expression. The signed charge is strictly negative for \(s>0\). The same apparatus can reverse the tag before the second write. This realizes the anticorrelated endpoint $R_{-1}$, and mixing with fresh tags gives negative $s$. The modular slope then changes sign, while the positive score has the same derivative. At the independent reference $s=0$ the modular slope is zero. The reference and source direction therefore carry observable information. Every additive proper-marginal subtraction contracts to zero. At an interior \(s\), both tangent signs are admitted. For independent replicas and an additive regional modular generator the two response differences add and keep their opposite signs.

The finite construction thus supplies an observable positive fluctuation and distinguishes it from a signed cap-energy transfer. Ordinary calibration of a positive coupling on the parked rest energy cannot transfer this amplitude by that energy meter: all its input laws coincide. The native four-register reference also gives the opposite sign for the proposed identification. A gravitational explanation through a different physical ball algebra or pressure response requires its independently measured source and area variation on the same states. This finite experiment measures neither proper aperture nor lens deflection. [Einstein source conditions](../PROOF.md) · [Exterior profile and optical closure](../../cosmology/rar_deep_regime/finite_halo/README.md)

## Restricting the modular operator to material records

On the protected two-seam family each marked-word Hamiltonian acts as a scalar on every endpoint input. The same property holds on the four modes storing a quarter-unit physical energy contrast in a single marked memory. Free evolution and classically scheduled material probes give one common history law in both arms, including all adaptive choices from earlier material outcomes. This follows by canceling the label-dependent phase before each instrument and conditioning on its result. [Exact source columns and scope of the instruments](../seam_overlap/README.md#an-exact-invisible-correlation-on-the-native-law)

Let \(\Phi\) produce the actual material output record and let \(\sigma\) be a reference faithful on its output support. Restrict the states first. The observed modular response on either protected transfer is exactly \(-{\rm Tr}\,\Phi(\rho_1-\rho_0)\log\Phi(\sigma)=0\), since the output difference itself vanishes. This holds for each admissible protocol on every number of copies. Direct source-word observations or coherent control of alternative service histories change the output algebra. The restriction removes the apparent contradiction between an invisible metric variation and a nonzero hidden modular source. It also removes the proposed dark source on those transfers, including the quarter-unit Hamiltonian difference. Universal coupling to that full Hamiltonian energy through the restricted readouts would still fail.

A [dynamically formed two-head collar](../../../model/experiments/interval_echo_einstein_station/collar_pairing/README.md#energy-and-round-trips-on-the-same-collar) is visible by this prescription. Use the 80-cell excess formation, an excess probe, and two units of native exposure per leg. The material register reports an outward miss, a return-leg miss or a completed return with probabilities \((0.90690362,0.05425742,0.03883896)\) on the paired state and \((0.88569202,0.06109865,0.05320932)\) on the product of its marginals. With these probability vectors \(p,q\), the signed transfer \(D=p-q\) gives \(-D\cdot\log p=-0.06454239\) nats at the paired reference. Subtracting the independently measured product-source output generator gives \(-D\cdot\log(p/q)=-0.00583834\) nats. Both signs reverse on weakening the pairing. The latter expression equals minus the sum of the two classical relative entropies. It is a contrast between specified outcome distributions, with no native area coefficient or vacuum ball weight attached. Suppressing contrast during the return probe with the formed source fixed gives zero for both responses to numerical precision. [Full outcome replay](material_algebra.py)

~~~sh
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.collar_identifiability.measure
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m matched_einstein_limit.collar_identifiability.material_algebra
OPENBLAS_NUM_THREADS=1 PYTHONPATH=code python -m pytest -q code/matched_einstein_limit/collar_identifiability
~~~
