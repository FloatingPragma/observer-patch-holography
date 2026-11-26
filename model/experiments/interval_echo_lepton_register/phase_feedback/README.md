# Turning three geometric click probabilities into three measured phase rates

A single vacancy in a guarded twelve-port Interval Echo cell carries three cyclic addresses. A fixed sequence of native pulses gives three distinct, independently read low-frequency phases on this material. The sequence is unitary and requires only a final measurement. All three characters receive the same applied phase. Their frequency ratios approach **1 ∶ 206.770316 ∶ 3477.472837** at weak drive. A second implementation uses a herald at each tick and resolves its individual phase rates. Releasing the same preparation into the freely evolving federation leaks every vector of the three-address space. The measured frequencies characterize the periodic drive.

## Encoding and preparation

An internal edge and its two face tips determine a four-port tile. A proper order-three rotation makes three disjoint copies covering the twelve ports. The cyclic addresses are the three translated tiles, the input ports are 0, 2 and 9, and each tile has center, positive-tip, negative-tip and idle paths. The same three paths carry a single vacancy with charge −1 in all preparations. The apparatus supplies a compensating stock, the addressed rotation and the phase reference.

A uniform native vacancy read on the marked edge gives path weights \(1/6,1/12,1/12,2/3\). Independent tip and edge records have exclusive-event probability \(2/9\). The phase program converts that exact probability to radians. This is the exposure convention of the earlier vacancy experiment. The historical 432-state register has its own event proof, and its free Dirac compiler uses a specified family mass coupling. Equal event fractions across these constructions do not identify their native clocks or interactions. [Marked vacancy preparation](../README.md) · [Theorem and compiler interfaces](../../../candidate_laws/interval_echo/openwave/electron_field/lepton_family/theorem_interfaces/README.md)

On each tile, routed native rotations prepare the four amplitudes from the square roots of the path weights. The two tip paths receive opposite phases and are shifted by one tile in opposite directions. On each output tile the coherent equal sum of its three active branches is brought back to its central address. The program checks the whole twelve-column unitary \(U\), its physically compiled inverse, and the return block on the three central addresses. With \(R\) the cyclic shift, this block is

$$
A=PUP=(\sqrt2 I+e^{i\theta} R+e^{-i\theta}R^\dagger)/6.
$$

Its intensity effect \(E=A^\dagger A\) is diagonal in the three characters. The historical 3,600-address coherent analyzer has active mask counts 500, 250 and 250 with amplitude denominator \(60\sqrt3\). Comparing these coefficients gives exactly \(E=(6/5)Y_{\mathrm{compiler}}\) after matching the cyclic orientation. This identifies the normalized effect on the three-character module. The two full preparations, kinetic generators and clocks remain distinct. [Original compiler and mask counts](https://github.com/FloatingPragma/observer-patch-holography/blob/3bcacf1276ddb3ead35be2b9adcd63b6643e2d35/model/lepton_mass_compiler/PROOF.md)

 The measured inputs to this construction are the spatial flag weights. Equal erasure, shift directions and one radian per exclusive-event fraction are apparatus settings. The archived finite flag frequency would give ratios **1 ∶ 205.1092 ∶ 3451.9785**; its 99% interval has endpoints mapping to middle/light ratios 199.23 and 211.24. The precise clock results below use the prescribed exact \(2/9\) setting. [Native pulses and free release](circuit.py)

## The same pulse on every character

Start on the central addresses. Apply \(U\), give every central output a relative phase \(-\epsilon\), apply \(U^\dagger\), then check the nine complementary addresses for the vacancy. All-clear preserves their unobserved complement coherently. A third external gateway carries the reference arm. The other two external owners can implement the vacancy pointer via the established echo and reverse echo. The relative negative phase is attainable by a positive periodic private dwell. All three characters take the same time on a successful tick.

The exact successful operator is
$$
K=P U^\dagger[I+(e^{-i\epsilon}-1)P]U P=I+(e^{-i\epsilon}-1)E.
$$
For an eigenvalue \(p\) its phase rate is \(\omega=-\arg(1-p+p e^{-i\epsilon})\) per tick and its retention probability is \(1-2p(1-p)(1-\cos\epsilon)\). A calibrated nonreal kick determines \(p\) from its measured quadrature. For \(0<\epsilon<\pi\), phase alone gives \(p=\sin\omega/[\sin\omega+\sin(\epsilon-\omega)]\). The local pulse has an exact weak-drive limit. In one tick the difference from \(e^{-ip\epsilon}\) is bounded by \(p(1-p)\epsilon^2/2\): expand the two phases about their weighted mean, cancel the linear terms, and bound each remainder by half its squared angle. Both amplitudes have modulus at most one, so \(N\) ticks multiply the error bound by at most \(N\). Holding \(N\epsilon\) fixed gives an effective evolution generated by \(E\) in accumulated kick angle. Each tick also requires routing and nine pointer checks, so this limit does not fix an energy in seconds. [Checked retention and identification lemmas](MeasuredDwell.lean)

## Coherent drive with a final read only

A fixed alternative cycle replaces every intermediate check with a native π phase at the nine complementary addresses. Its unitary is $W=(2P-I)[I+(e^{-i\epsilon}-1)U^\dagger P U]$. Each character spans a two-dimensional invariant plane with its component outside the central addresses. The restricted determinant is $-e^{-i\epsilon}$ and its trace is $(2p-1)(e^{-i\epsilon}-1)$. The low phase rate is therefore
$$
\omega=\epsilon/2-\arcsin[(1-2p)\sin(\epsilon/2)].
$$
Its ratio to the exposure tends to $p$. This is a Floquet phase per repetition of the complete pulse schedule. On the full twelve rails there are three distinct low modes, three partner modes near a half turn and six modes exactly at a half turn. The code checks that spectrum directly on the compiled physical unitary.

At common exposure 1/128 the low rates in ascending order are about 0.000000706651, 0.000146114 and 0.002457350 radians per tick. All intermediate amplitudes are retained. The original character input has weight above 0.9999967 on its low mode. A final interferometer compares its return against the external reference. The exact unwanted branch weight is bounded before clock inference: $(W-I)$ on that input has squared norm $4p\sin^2(\epsilon/2)$, and the partner eigenvalue has distance at least $2\cos(\epsilon/2)$ from one. Its contribution to the return relative to a single unit-modulus oscillation is bounded by $2\tan^2(\epsilon/2)$ for every delay. Including this nuisance and the counting errors, the final-only records give middle/light **206.769702–206.770883** and heavy/light **3477.45035–3477.47012** at 99% simultaneous confidence for this experiment. All final output failures are recorded. The exact finite-drive ratios are 206.770258 and 3477.459652. [Full native cycle](coherent.py) · [Final-only counts](coherent.json)

## Reading the rates independently

Prepare a Fourier character in one arm of a balanced interferometer and park the reference vacancy outside the cell. The matching global charge in the two alternatives preserves the phase comparison. After each power-of-two tick count up to \(2^{24}\), two analyzer phases record both outputs and every failed return. The simulation uses \(2^{42}\) virtual starts per setting. The estimator uses only these counts, the common drive angle, and \(0\leq\omega\leq\epsilon\). It resolves the phase windings across the delays and returns simultaneous 99% intervals. The native circuit is used to generate the outcomes and to check them separately. [Complete count banks and calibrations](readings.json) · [Reader and replay](measure.py)

At common exposure 1/128 the measured middle/light phase ratio is **206.770362–206.770491** and the heavy/light ratio is **3477.4981–3477.5002**. Inverting the measured phase map gives effect-ratio intervals **206.770246–206.770375** and **3477.47170–3477.47380**. Increasing the common exposure to 1/16 changes the measured phase ratios to approximately 206.7777 and 3479.161 while recovering the same two effects. The middle/light weak-limit ratio exceeds the repository's frozen laboratory ratio 206.7682827 by 9.83 parts per million. That difference is the prediction of the supplied Koide–Brannen circulant at 2/9 radians. An electron anchor converts its ratios to energy units. [Frozen laboratory comparison and calibration](../../../candidate_laws/interval_echo/openwave/electron_field/lepton_family/hierarchy/README.md)

Two controls isolate what the measurements identify. At a common phase of π all three returns have real, positive amplitude and zero phase, although the static geometric probabilities and their Koide quotient stay identical. Omitting both the return checks and the complementary π phase gives a free-to-repeat twelve-port unitary conjugate to the bare common phase gate, whose only eigenphases are zero and the supplied phase. The static records and three character labels alone therefore cannot distinguish these timed experiments.

With all addressed controls removed the independent native free-release calculation gives \(4PH(I-P)HP=3I+2\mathbf1\mathbf1^T\), whose smallest eigenvalue is three. At one native time the three-character preparations retain only 22.5%, 56.4% and 56.4% probability on the central addresses. The theorem predicts the ratios of the returned, weakly driven record. A charged-lepton assignment would require a free, independently tagged energy–momentum pole for each character, spatial spin and exchange on that excitation, and an independent energy calibration. [Existing free-pole readout](../pole_meter/README.md)

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -m model.experiments.interval_echo_lepton_register.phase_feedback.measure
python -m model.experiments.interval_echo_lepton_register.phase_feedback.coherent
python -m pytest -q -W error model/experiments/interval_echo_lepton_register/phase_feedback
(cd Lean && lake env lean -j 1 ../model/experiments/interval_echo_lepton_register/phase_feedback/MeasuredDwell.lean)
```
