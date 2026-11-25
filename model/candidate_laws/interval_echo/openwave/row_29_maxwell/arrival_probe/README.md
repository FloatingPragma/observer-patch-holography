# What does an independent probe feel when the polarization arrives?

An owner three cell boundaries from the emitter can read a circulating material current and detect the arriving load with a second, independently prepared load. At injection, reversing the probe velocity leaves every initial address acceleration exactly equal. The measured circulation would give a nonzero velocity-odd acceleration under a magnetic interpretation with nonzero scalar gain. This rules out that identification for this probe, on the specified arrival and preparation. [Original polarization flight](../polarization_channel/README.md) · [Protocol](PROTOCOL.md)

## Two measurements on the same arriving preparation

The emitter launches one unit on the intact 80-cell screen. Use the same route \(0,3,1,18\), exposure 64 and port frame as the six-state polarization experiment. Run all six transverse inputs, the axial input and a blank sham. The owner at cell 18 measures clockwise and anticlockwise material probability current on each of its twenty triangular faces, each on fresh preparations. Weight the sign by its current eigenvalue \(\sqrt{3}/2\) and the outward unit face normal, then average the twenty vectors to obtain the axial reading \(M\). Other local outcomes and empty cells also count.

On separate launches the receiving owner shutters its external exchanges and reads whether the arriving load occupies one of probe ports 0 and 1, another port in the cell, or the exterior. The first outcome blocks injection. In either other case a supplied source launches an amount-two probe with amplitudes \((|0\rangle\pm i|1\rangle)/\sqrt2\). These have opposite initial velocity along the same unit direction \(u\). A third source uses a real plus sign. The observer records the final probe port. A matched control takes the same admission reading and stores the arrived single outside the interaction before injecting the probe. The probe estimator subtracts that control and divides by all source launches, including blocked and exterior outcomes. A separate matched circulation bank imposes the same admission gate and retains its failures. Complete face and address histograms use \(2^{38}\) simulated multinomial draws per setting. [Counts, inputs and routing](records.json) · [Acquisition](measure.py)

## Velocity reversal at injection

For a local double-address projector \(P_j\), the first two derivatives of its probability are the expectations of \(i[H,P_j]\) and \(-[H,[H,P_j]]\). The two-load register has 132 ordered positions for a double and a single. Up to a scalar, twice the native Hamiltonian in this register is an integer matrix. For each probe seam \(a,b\), restrict the single to the other ten ports. The first commutator has the same probe blocks as the empty receiver. In the second derivative, the \(a,b\) block is a real symmetric ten-by-ten matrix. Reversing the probe's relative phase replaces its expectation difference by a multiple of that block minus its transpose, which vanishes. The calculation covers every admitted single density matrix, all twelve probe addresses, all thirty injection seams, and the native and zero-contrast laws. It also holds for the solo stored-probe control. A detector check supplies a phase on another side of the injection triangle. The third port then has initial probability accelerations 1 and 0 for the two probe signs. [Integer operator check and detector control](onset.py)

The matched face records give \(\|u\times M\|=0.00028983\) and \(0.00024359\) for the two circular inputs. The simultaneous 99% norm error bound is 0.00001371, from bounded outcomes over eight sources, twenty faces and three components. The probe has exactly zero velocity-odd initial acceleration on both inputs. Hence no nonzero scalar gain maps \(u\times M\) to this acceleration. These are two separately executable instruments on identically prepared arrivals. Their units use the supplied port chart and exchange exposure.

## Finite contact and the scope of the field interpretation

Let \(D\) be half the difference between the control-subtracted terminal displacement vectors of the two probe signs. In a half-unit contact window the two circular inputs give the following readings per original launch, in cell circumradii. The simultaneous 99% displacement half-width is 0.00001865.

| Incoming circular mode | Fraction admitted locally | Native \(D\cdot u\) | Zero-contrast \(D\cdot u\) |
| --- | ---: | ---: | ---: |
| Positive | 0.015059 | −0.00168023 | −0.00009878 |
| Negative | 0.012446 | −0.00222796 | −0.00012919 |

Both windows, \(1/8\) and \(1/2\), resolve a negative longitudinal displacement for each of the seven occupied native inputs. The blank sham gives zero. The finite-window rule \(D=g(u\times M)\) predicts zero longitudinal displacement for every gain. For comparison, a charged particle launched perpendicular to a uniform magnetic field at speed \(v\) has longitudinal displacement \(v[\sin(\omega t)/\omega-t]\) relative to free flight. That quantity reverses with the incoming velocity even under the exact Lorentz force. The finite-window result tests a displacement formula. The double-commutator test above tests acceleration at injection.

Velocity-even and real-source displacements are recorded for every input. The probe's electric charge, an independent electric-field calibration, a magnetic-field calibration with the correct force response, and an independently clocked wave speed are unresolved. The experiment identifies a local interaction with arrived material on a shuttered cell. The OpenWave requirement is recovery of Maxwell dynamics and transverse waves traveling at \(c\). [Official criterion](https://github.com/openwave-labs/openwave/blob/0df73dc6be3bde9f0b8ced6044effc79f0d92bd6/MODELS.md#L135)

From the RER repository root:

~~~sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_29_maxwell.arrival_probe.measure
.venv/bin/python -B -W error -m model.candidate_laws.interval_echo.openwave.row_29_maxwell.arrival_probe.onset
.venv/bin/python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_29_maxwell/arrival_probe/test_probe.py
~~~

The test rebuilds the full local Hamiltonians from native seam columns, propagates the incoming load independently on the board, reconstructs the face projectors, and checks all outcome banks and normalizations. It shares the native law, board and port chart with the producer. [Reconstruction](test_probe.py)
