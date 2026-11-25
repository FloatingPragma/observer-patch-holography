# How two arms of a bound trimer carry conjugate colour

An observer on one arm reads one of three orbital labels. Observers on the other two arms each read a label from their own detector. Their antisymmetric pair carries the **conjugate** three-dimensional rotation. One can predict all three outcomes from an independently launched single load that traverses the inverse of the pair's control word. Four held tables agree with that prediction. Using a fundamental action for the pair fails at three settings. A second experiment reads opposite centre phases on the constituent and the pair with the same material clock. [Frozen transport protocol](PROTOCOL.md) · [All address counts](counts.json) · [Centre protocol](CENTRE_PROTOCOL.md) · [All clock counts](centre.json)

## Source, routes and records

The source is the [energy-selected bound three-load face](../../../candidate_laws/interval_echo/openwave/row_13_baryon/arm_symmetry/README.md). One sharp face and a four-stage energy filter on the guarded native Hamiltonian select a simple line at 100.8894585 in exchange units. Its sign under spatial arm permutations is determined by the energy line. On the one-load-per-arm three-label detector the ideal admitted state is the normalized determinant. The main decoder omits the fourth column of each arm. Each source run has a filter acceptance of 0.08364493 and conditional orbital acceptance of 0.00040717, about 34 admissions per million original launches. The calculation uses the actual finite filter in each new source measurement. The earlier all-open experiment bounds binding and fidelity after release. The present transport and clock experiments use protected receiver contacts.

A native π contact on a blank seam exchanges its endpoint loads. Transposing two addresses along a path and restoring all intermediate addresses transports the three labels of each arm into owners 1, 4 and 6 on the twenty-cell board. The three outbound words contain 23, 21 and 27 native swaps. Each word returns after the receiver operation. The receiver triangles have three private degree-five ports. The controller times a half exchange on the first two addresses (A), or a half exchange on the last two followed by a phase on the third (B). A degree-five private phase on each of the two interacting rails removes the common energy shift. All three rotations, the return routes and the terminal projective reads are specified apparatus. Qutrit number is conserved on each receiver, so excluded occupations cannot enter the counted triplet through the controls. Every launch ends in one of 29 recorded bins, including the two sources of refusal. [Native gates and measured words](measure.py)

Let V be the word applied to the first constituent, and U the word applied to each of the other two. Word letters below read in time order. The one-load calibration sends each of three sharp labels k through the inverse U word and then V. It counts which label a emerges. In the trimer read, when the two other arms return distinct labels b and c, k is the missing label. The joint probability for the ordered triple a,b,c, conditional on detector admission, is one sixth of the calibration probability of a at k. The callee receives independent sharp launches and the prediction uses its observed counts. Transposing U in that calibration instead implements the ordinary-representation foil. It is experimentally distinguishable on three of the four enrolled settings.

| First arm | Other two arms | Admitted out of 2³² starts | Repeated labels | Largest discrepancy from inverse-word prediction | Largest discrepancy using transpose |
|---|---|---:|---:|---:|---:|
| AB | AB | 146,426 | 0 | 0.0017 | 0.1674 |
| AB | BA | 146,940 | 30,916 | 0.0013 | 0.1182 |
| BAB | AB | 146,313 | 68,473 | 0.0015 | 0.1183 |
| A | B | 145,512 | 84,936 | 0.0012 | 0.0012 |

The simultaneous 99% allowances are below 0.0078 in each table. They include counting uncertainty from every held and calibration bin. Each sharp calibration uses 262,144 starts. The finite-filter vector differs from the ideal line by at most 1.09×10⁻⁸. Normalizing after triplet admission increases the possible conditional probability difference by at most 2ε/(√p−ε), below 0.0000011. With the original labels measured and forgotten before the common AB, 109,762 admitted starts had repeated labels. Parking one of the remote routes gave 85,380. Both controls fail the calibrated prediction at every setting. With zero contrast throughout the source and receiver law the declared source filter accepts approximately 1.75×10⁻⁴⁹ of starts. The corresponding empty admission banks provide no conditional transport measurement.

## Why the conjugate appears and what the centre measures

Identify the missing-label space with the antisymmetric pair by the isometry

\[
 J|k\rangle=2^{-1/2}\sum_{bc}\varepsilon_{kbc}|b,c\rangle .
\]

Then \((U\otimes U)J=(\det U)J\bar U\) for every U(3) word, by the cofactor identity. Hence the sharp joint probabilities above involve \(|(VU^\dagger)_{ak}|^2/6\). The determinant factor disappears in these sharp reads. On SU(3) the transported pair is exactly the antifundamental. This source attachment follows because the measured energy line selects the determinant before the independent rotations are applied. The output test asks only for three ordinary address reads. It uses no final singlet projection.

The complete local-unitary stabilizer also follows. On any pair of arms the reduced density matrix is one third of the projector onto the antisymmetric subspace. If X⊗Y preserves this projector it commutes with the swap, so X⊗Y=Y⊗X and X is a scalar multiple of Y. Applying this to both pairs forces all three local frames in a ray-preserving word to have the form Uᵢ=zᵢU. Conversely such a word multiplies the determinant by z₁z₂z₃ det U. The complete vector stabilizer requires this phase to equal one. For three independently timed SU(3) centre words, the neutral subgroup thus consists exactly of the nine triples of cube roots with product one. This describes the centre action on this trimer. It gives no electroweak doublet, hypercharge assignment or identification of the six-axis gluing loops in PR-35. A physical colour connection in PR-46 also requires a source for the transported frames themselves. [Finite theorem on the different electroweak quotient](../../../../Lean/Screen/TraceBalancedKernel.lean)

For a separate phase reading, a single two-rail material clock visits the selected receivers by native transpositions. A contact on two binary ports, refocused by a private half-turn and corrected on each individual port, realizes exp(iφ n₁n₂) exactly. Apply these clock-controlled phases across three labels in the ratios (1,1,−2). A sharp single load independently calibrates each label. At θ=2π/3 the word is the centre of this traceless orbital rotation. The simulated readings give the following conditional clock quadratures. Each entry has simultaneous 99% uncertainty below 0.0119.

| Clock controls | X | Y |
|---|---:|---:|
| One constituent | −0.4997 | +0.8650 |
| Complementary pair | −0.4989 | −0.8676 |
| All three | +1.0000 | −0.0052 |

At θ=π/4 the constituent gives (0.4739, 0.1387) and the pair gives (0.4700, −0.1406). The three-arm trace stays at one. With uniform material-number phases (1,1,1), the central exposure gives the same reading, while the π/4 and π settings distinguish the two controls. All 27 clock settings agree with the predictions for their respective words. The same clock visits each selected arm. An inverse half contact and, for X, a private phase turn each quadrature into an ordinary terminal port read on an independent cohort. Their sampling uses the declared Born instrument. [Native clock coupling](centre.py) · [Independent checks](test_transport.py)

## A singlet certificate from four sets of local readings

Can three observers certify the state without an apparatus that projects onto the proposed singlet? Each observer switches among four local qutrit bases. They compare their final address labels in the same basis. Three settings use the discrete Fourier transform following a quadratic phase, and the fourth uses the original addresses. Each change of basis compiles into private phase waits and native pair contacts on the receiver triangle. Four independent cohorts of 2³² original starts are simulated per case. The protocol specifies a common source and basis-independent orbital admission. Counts use the Born probabilities of the registered generator with the stated selection and detector. [Frozen local protocol](LOCAL_PROTOCOL.md) · [Scheduled basis changes and decoder](local_bases.py) · [Complete counts and intervals](local.json)

For each pair of arms and each basis, read the frequency of equal labels. Let D be the sum of the twelve equal-label probabilities. The four qutrit bases form a projective two-design. On each pair the four equal-label operators therefore sum to I + swap. This identity also follows by expanding the quadratic Fourier phases and summing the cube roots. On all three arms the resulting operator is 3I + swap₁₂ + swap₁₃ + swap₂₃. It has eigenvalues 0 on the single determinant line, 3 on the mixed-symmetry space of dimension 16 and 6 on the symmetric space of dimension 10. Thus for **any mixed state** on the admitted three qutrits its singlet fidelity F obeys 1 − D/3 ≤ F ≤ 1 − D/6. The tests verify the two-design identity using the compiled native matrices and independently check all 27 eigenvalues. The measurement uses local position detectors throughout.

| Source and receiver words | Admissions across all four bases | Conditional singlet fidelity, simultaneous 99% | Singlet weight per million original starts, simultaneous 99% |
|---|---:|---:|---:|
| Idle | 585,142 | [0.9997491, 1] | [33.70, 34.38] |
| Common AB on all three | 585,045 | [0.9997491, 1] | [33.82, 34.51] |
| First AB, pair BA | 586,469 | [0.7811826, 0.8974710] | [26.40, 30.94] |
| Measure and forget, then common AB | 585,094 | [0, 0.5087639] | [0, 17.53] |

The conditional intervals use simultaneous two-sided exact binomial coverage over the 48 pair frequencies and a 10⁻⁷ tolerance for numerical gate synthesis and rounded probabilities. The final column assigns zero weight to either refusal. It combines the same bound with an exact admission interval on the first cohort of each case, allowing for all 52 intervals simultaneously. The product uses no independence assumption between the estimated admission rate and fidelity. For common AB, 15,742,824,242 starts fail the filter and 1,436,459,897 fall outside the orbital sector, out of 17,179,869,184 original starts. The observers can identify the prepared colour singlet to high fidelity once admitted. Its preparation yield is about 34 per million launches.

## Reproduction

From the RER root, run the following commands.

~~~sh
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.experiments.interval_echo_gauge.trimer_transport.measure
.venv/bin/python -B -W error -m model.experiments.interval_echo_gauge.trimer_transport.centre
.venv/bin/python -B -W error -m model.experiments.interval_echo_gauge.trimer_transport.local_bases
.venv/bin/python -m pytest -q model/experiments/interval_echo_gauge/trimer_transport
~~~

Each receipt retains every launch-probability vector in sampling order, with the 14-decimal sampling grid. Numerical libraries differ in the last bits of these probabilities, and one such bit can change every later multinomial draw. Replay recomputes each vector from the native evolution and requires agreement to 2 × 10⁻¹¹ relative plus 10⁻¹⁴ per bin. It then draws the archived counts from the recorded vector with the original seed.
