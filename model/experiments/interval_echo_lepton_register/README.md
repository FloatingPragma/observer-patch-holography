# What does the three-line detector measure?

A single vacancy on a filled twelve-port cell produces three interference intensities in the ratio **1 ∶ 206.770316 ∶ 3477.472837** under registered Interval Echo contacts. A free-release test on the same preparation distinguishes these intensities from rest energies. On the marked edge, every nonzero combination of the three detector modes couples to other modes when the guards open. The exact predictions, complete simulated port counts and independent checks are in the [experiment](measure.py), [records](results.json) and [tests](test_measure.py).

The apparatus declares an internal edge and labels its two incident face tips. A path of native contacts and private load phases spreads one vacancy uniformly over twelve ports. Eight are idle, the two endpoints supply the central arm, and the two tips supply the oriented arms. Its four position probabilities are 2/3, 1/6, 1/12 and 1/12. Independently repeated preparations give probability 2/9 that exactly one of “first tip on C” and “endpoint on E” occurs. Copying one address instead gives 1/4. The endpoint and tip names describe a geometric read. Identifying them with the repair events of the distinct 432-state lepton round requires an instrument map. [Protocol](PROTOCOL.md) · [Register and field constructions](../../candidate_laws/interval_echo/openwave/electron_field/lepton_family/theorem_interfaces/README.md)

On a fresh vacancy the apparatus sets opposite tip phases in radians, first from the exact exclusive probability and independently from the observed exclusive fraction. This conversion sets the exposure. It applies a three-port Fourier reader to the symmetric endpoint pair and the two tips. Native exchanges route each pulse. A degree-five port supplies the private phase; swaps return the rail to its address. Two independent external owners can read each terminal vacancy through the echo. Keeping every output, the three character probabilities are

$$
P_k=\frac{1}{36}\left[\sqrt{2}+2\cos(\theta+2\pi k/3)\right]^2,\qquad k=0,1,2 .
$$

The three probabilities sum to 1/3 on the coherent uniform input. At the count setting the three bracketed amplitudes are positive and their intensities have Koide quotient 2/3. All thirty choices of internal edge give the same three intensities. In 4,194,304 simulated starts on the marked edge the three channels received 370, 78,119 and 1,318,931 clicks in ascending intensity order. The six-digit ratios above are predictions of the pulse calculation. The finite calibration observed 931,464 exclusive events in 4,194,304 starts. Feeding that frequency back as the actual angle gives intensity predictions 1 ∶ 205.109227 ∶ 3451.978496. The small channel is sensitive to this calibration. Dephasing the vacancy before the reader gives three equal probabilities of 1/12. Reversing the angle exchanges the two smaller responses. Setting the angle to zero makes them equal. Tripling the exchange strength at the frozen exposures gives quotient 0.460044. [Complete ports and controls](results.json)

For an independent free-return test, reverse the reader from each sharp character output, release it into the twenty-cell federation, then reapply the same reader. The terminal bank keeps each of the twelve home addresses and a thirteenth outcome for an absent vacancy. At one native time unit the three summed character return probabilities, in ascending initial intensity order, are 0.370760, 0.371419 and 0.456450. The constant host energy cancels in these reads. A vector in the marked three-mode space has equal amplitudes $x$ at endpoints 0 and 1, and tip amplitudes $y,z$ at ports 5 and 7. Outside this space the free-generator rows at ports 4, 6 and 8 give $-y/2,-z/2,-(x+z)/2$. Their common kernel is zero. The test checks an integer rank witness on all thirty edges. Thus the exact Koide intensity pattern can coexist with the absence of even one free eigenstate supported in its three-mode space. [Lean coefficient and leakage lemmas](../../../Lean/ObserverPatchHolography/LeptonRegisterCharacterisation.lean)

From the repository root run the following commands.

```sh
OPENBLAS_NUM_THREADS=1 python -m model.experiments.interval_echo_lepton_register.measure
OPENBLAS_NUM_THREADS=1 python -m pytest -q model/experiments/interval_echo_lepton_register
```
