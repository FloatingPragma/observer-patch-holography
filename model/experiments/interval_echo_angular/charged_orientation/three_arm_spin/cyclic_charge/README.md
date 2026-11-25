# Three charged winding characters reach the other observer

A vacancy emitted from port 4 reaches port 160 with **three coherently recoverable winding states** under the registered Interval Echo law. Each owner can operate as emitter or receiver. Two complementary measurements certify transmission of a qutrit after all three emitting head records are discarded. Across six choices of direction and independent apparatus, the simultaneous 99% lower bound on entanglement fidelity ranges from **0.73246 to 0.73641**. A system carrying at most two entangled levels has fidelity at most $2/3$. The two opposite winding characters have exactly the same energy spectrum. These records identify three states of one guided charge and give only two distinct source spectra. [Native records and controls](results.json)

## From a finite symmetry to two local measurements

The three equivalent arms of the [physical fixture](../README.md) can be cycled together. The induced permutation $C$ of the 3,584 charge and head words obeys $C^3=1$ and commutes exactly with the free generator. Eight words are fixed. The other words form triples, so the three character spaces have dimensions 1,200, 1,192 and 1,192. The three locally preparable characters are the equal-phase sum of the single-wound arms and the two sums with phases advancing or retreating by $2\pi/3$. All three start with exactly one vacancy on the common physical port.

On the material emitter, blank spokes split the actual vacancy over the three marked triangles. Each visited triangle writes its own direction and the inverse splitter returns the vacancy to the common port with probability $1/3$. Omitting the relative phases prepares the trivial character. Opposite native phase loops prepare the other two. The receiving owner uses its own inverse circuit and three local unwound checks. A separate controlled-path pointer supplies the same preparations and measurements. An independent sharp occupation read tests charge arrival before the winding analyzer. For the complementary sharp-arm basis, an isolated native blank-spoke pulse attenuates a single winding preparation to the same $1/3$ acceptance. The finite group and matter program select these bases without a mass target. [Emission and physical detector](instruments.py) · [Independent pointer and local phase calibration](../material_interferometer/README.md)

The common coherent extractor first splits the arriving vacancy. Each branch approaches its local head in reverse. On a matched crossing its vacancy parks at one of three distinct material addresses, 164, 162 or 156 at the far end. Three unit-occupancy checks at the unused exits accept all successful paths together. A held auxiliary vacancy from a compensated pair checks the three local heads for zero using touches disjoint from the storage ports. A stock transposition restores every intermediate load even when the stored vacancy occupies the route. Replaying all eight input head words gives the three-by-eight matrix $D[I_{\rm one\text{-}hot}]/\sqrt3$, where $D$ contains the calibrated phases of the three output ports. Thus every accepted one-hot superposition and its untouched reference reach three material charge locations coherently. The reverse observer uses ports 6, 8 and 0. The apparatus supplies its held stock, addressed contacts and appointments. Free flight uses the unmodified native generator.

## What the count records certify

Each preparation waits 32 native units without a clock filter. For the maximally mixed three-state input the average decoded probability is $1.61898\times10^{-7}$ per original start including both one-third apertures. Each input and analyzer setting has $2^{38}$ original simulated starts. Every bank retains source refusal, missing charge, other charged arrivals and accepted phase. The phase comparison uses independent banks for each binary projection. Numerical calibration of the registered pulses agrees with both intended bases to better than $10^{-12}$ in state norm. Each original-start probability interval includes an additional $2\times10^{-11}$ allowance. [Crossed instruments and sampling](measure.py)

For a normalized successful channel Choi state, let $p_Z,p_X$ be the respective same-label fractions in the sharp-arm and Fourier bases. The original-reference basis is conjugated for preparation. The two match projectors have product $|\Phi_3\rangle\langle\Phi_3|$, where $\Phi_3$ is the maximally entangled qutrit pair. Their orthogonal complements outside this common line give $F=\langle\Phi_3|\rho|\Phi_3\rangle\geq p_Z+p_X-1$. A normalized pure Schmidt-rank-$r$ state overlaps the maximally entangled state with squared magnitude at most $r/3$ by Cauchy–Schwarz on its $r$ Schmidt coefficients. Convexity preserves that bound for a mixture.

| Experiment | Sharp-arm match | Three-character match | Lower fidelity |
| --- | ---: | ---: | ---: |
| Exact free-channel probabilities | 0.90046 | 0.92204 | 0.82250 from the two bases |
| Material emission and detection, 4 to 160 | 0.90004 | 0.92084 | 0.73246 at simultaneous 99% |
| Material emission and detection, 160 to 4 | 0.89923 | 0.92227 | 0.73323 at simultaneous 99% |
| Source with the sharp-arm label dephased | 0.89970 | 0.33515 | 0.16765 by the same 99% bound |

The directly evaluated Choi fidelity is 0.88841. The two-basis certificate uses only the count intervals. Its confidence bound sums lower intervals for the matched outcomes and upper intervals for all outcomes before dividing. The physical charge source with a pointer at the destination and a pointer source with a physical detector also clear the two-level ceiling independently, in both directions. The analytic group selection and fidelity threshold are checked in [Lean](CharacterRead.lean).

A second probe allows the emitting observer to measure its own residual head word after the first vacancy has arrived at the far port. When both local one-hot Fourier reads succeed, $k_{\rm donor}+k_{\rm receiver}=k_{\rm input}\pmod3$. When the donor is all unwound its contribution is zero. All nine allowed one-hot splits appear. The independent scans give 122,971 accepted joint records forward and 122,539 reverse, all in allowed bins. Ignoring the donor produces a nonzero wrong-character rate at the receiver. Since the sharp hub is fixed by $C$, the measured sum rule follows by applying $C$ to the input ket and the two local output bras. Outcomes with different characters have zero amplitude.

## What if the observers use different phase conventions?

The emitter shifts the origin of its three arm labels and phases the one-hot addresses by $1,\omega,\omega^2$. The other observer independently reverses its arm ordering with a shifted origin and phases $1,\omega^{-1},\omega^{-2}$. The product action extends to all local head masks, including the residual source memory. We propagated the full Hamiltonian again after transforming its 3,584 basis words, the two source preparations and both local detection maps. The same random draws produce exactly the same **90 original-start count banks**, with a largest probability difference of $1.3\times10^{-18}$. Both descriptions have 245,233 accepted character outcomes, all respecting the sum rule. Keeping the receiver's old labels in the new frame yields 204,720 forbidden-label outcomes. A physical arm intervention changes what the source does. These compensated coordinate changes express the same experiment in different local conventions. A Standard Model gauge representation requires its own physical current and coupling measurements. [Frame transformation](frames.py) · [Held comparison](frames.json)

## Why these states do not fix three rest masses

Reflecting the entire arm apparatus exchanges the two winding characters and leaves the real native Hamiltonian invariant. Every spectral projector therefore gives exactly the same weight on the two conjugate preparations. More generally a real free generator gives equal spectral measures to a state and its complex conjugate, even without the arm reflection. The symmetric character can differ. At moment ten of $2H$, the sharp-source diagonal is 121592396787244320 and each off-diagonal is 128. All earlier off-diagonal moments vanish. The symmetric character has moment $d+2c$, while either winding has $d-c$. The exact difference is 384. [Exact integer recurrence and spectral readout](spectrum.py)

| Measured source response | Symmetric character | Either winding character |
| --- | ---: | ---: |
| Mean formation energy | 25.5 | 25.5 |
| Energy variance | 0.75 | 0.75 |
| Lowest observed line with source weight above $10^{-8}$ | 8.9307090985 | 8.9381807248 |
| Its source weight | $6.10\times10^{-8}$ | $9.09\times10^{-4}$ |

The lower lines depend on the specified guide and threshold. An independent exact certificate bounds every energy of this entire preparation between 8 and 27, so neither charged-lepton heavy/light ratio fits under one scale. The cyclic symmetry is a symmetry of the three complete arms. Measuring an individual intervening vacancy address generally disturbs its character. These phase-resolved charged states give an operational test for a possible family basis elsewhere in the material. A lepton-family assignment still needs three distinct rest poles, spatial spin, the appropriate gauge currents and flavor-changing responses on the same carriers. [Full spectrum bound](../flavor_search/README.md) · [Finite three-line geometry](../../../../interval_echo_lepton_register/README.md)

From the repository root run the following commands.

~~~sh
p=model.experiments.interval_echo_angular.charged_orientation.three_arm_spin.cyclic_charge
d=model/experiments/interval_echo_angular/charged_orientation/three_arm_spin/cyclic_charge
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m "$p.measure"
python -m "$p.frames"
python -m pytest -q "$d"
(cd Lean && lake env lean -j1 "../$d/CharacterRead.lean")
~~~
