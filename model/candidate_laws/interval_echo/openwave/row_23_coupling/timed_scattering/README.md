# Does scattering begin to resolve the interaction at the measured core radius?

Two freely evolving material doubles give a measurable collision signal when the observer interferes an early contact with a late one. How much momentum transfer is needed to see that signal depends strongly on how long the material flies between the contacts. On the narrower core, the sixteen-period signal already exceeds the registered 10% threshold below half the inverse core radius. This protocol leaves OW-23 partial. [Declared test](PROTOCOL.md) · [Decoded measurements](results.json) · [Terminal counts](counts.npz)

## What the observer measures

The source is the native two-double line on the commissioned 214-port guide of the 80-cell Interval Echo federation. Its independent endpoint counts give rms radii of 6.588 guide-port spacings at winding 0 and 2.376 at winding 13. In each run a coherent switch puts the same addressed contact before or after a period of free motion. Both paths reach the detector at the same time. The observer reads the address pointer's Fourier bin and the time switch's recombination port. The preparation filter, matching schedule, coherent switches, off-guide vacancy test, Fourier analyzer and timing reference are supplied instruments. All reported starts are simulated. [Source and independent radius read](../mobile_form_factor/README.md) · [Native gate used in the guide](../../row_30/experiment.py)

Let $D_j$ be the no-spill contact followed by Fourier read $j$. With $L=214$, $q=2\pi j/L$, no-spill occupied-contact amplitude $a$, and occupied guide addresses $x$, its action on a material word is diagonal,

$$
D_j=\delta_{j0}+\frac{a-1}{L}\sum_x e^{-iqx}.
$$

For each accepted initial orientation the early and late recoil states are $E=U^t D_j\psi$ and $B=D_j U^t\psi$. The acceptance probabilities are retained in the unnormalized source states. The two recombination ports measure $\|E\pm e^{i\phi}B\|^2/4$ for phases zero and $\pi/2$. Their differences give the two components of $\langle E,B\rangle$, averaged over the two source orientations before normalization. Separate early-only and late-only banks give both arm intensities. The observable $F=\langle E,B\rangle/\sqrt{\langle E,E\rangle\langle B,B\rangle}$ divides out the single-contact strength and the individual recoil intensities. The detector keeps every recoil bin, both time ports, spills and rejected sources. [Executable acquisition and decoder](measure.py)

Both flight laws start from the same native prepared state and use the same native contact. Switching the contrast from 8 to 0 *during flight* changes an isolated double by a global phase. On the occupied–vacant gate the contrast term is $4\kappa I$. Two doubles on separate matched gates acquire the product of those phases. A doubly occupied gate with one double at each end has zero generator. Consequently only an encounter on a common gate changes the relative flight phase. The first matching lasts $\pi$, giving the relative factor $e^{8i\kappa\pi}=1$ at the two enrolled contrasts. The second matching's collision phase is diagonal in the output occupation words. Thus for one period $U_8$ equals $U_0$ followed by a diagonal word phase, up to a common phase. That final phase commutes with every $D_j$ and cancels in the interference. Over multiple periods intervening motion mixes those words and makes an earlier collision observable. [Independent complete unordered-pair check](test_timed_scattering.py)

## Records

The comparison uses $|F_8-F_0|$ at fixed source, delay and Fourier bin. A bin is resolved if its simultaneous interval excludes zero. The last column includes every possible first 10% crossing on the **sampled** positive bins. If an earlier bin is statistically unresolved at 10%, the conservative bracket starts at the end of the initial contiguous run certified below 10%. Interpolation between discrete modes is not measured.

| Source | Flight periods | Resolved bins / 107 | First bin certified above 10% | Possible first-crossing $qR$ |
| --- | ---: | ---: | ---: | ---: |
| Winding 0 | 0, 1 | 0 in each | none | all bins below 10% |
| Winding 13 | 0, 1 | 0 in each | none | all bins below 10% |
| Winding 0 | 4 | 103 | 100 | 11.73–19.46 |
| Winding 13 | 4 | 99 | 65 | 0.78–4.81 |
| Winding 0 | 16 | 107 | 3 | 0.384–0.584 |
| Winding 13 | 16 | 106 | 4 | 0.196–0.296 |
| Isolated double | 0, 1, 4, 16 | 0 in each | none | all bins below 10% |

At four periods the closest certified-below bins preceding the first certified-above bins are 98 and 62, respectively. Earlier unresolved bins explain the wider brackets. At sixteen periods the preceding bins 2 and 3 are certified below and all earlier bins are also below. The largest binwise lower bounds on the effect reach 0.622 and 0.324 for the two sources after sixteen periods. The isolated-double upper bounds stay below 0.0038 at every delay. The forward read is a constant on each fixed-load sector and its exact timing overlap is one at either law. [Complete curves and bin classifications](results.json)

Each of the 96 acquisition settings uses $2^{40}$ multinomial simulated starts. The simultaneous detector intervals use at most 300,000 one-sided count tails, each budgeted $0.01/300000$, and include probability quantization in the sampler. Signed component intervals and both independent arm intensities are propagated through the normalization and the contrast subtraction. Detector-only family coverage is at least 99%, and the union with the earlier radius family has coverage at least 98%. This comparison measures a two-time collision response, conditional on the stated apparatus and source. A gauge-coupling interpretation and a common core-controlled onset remain open. [Forward-tail checks](../../../../../experiments/count_intervals.py)

## Reproduction

From the repository root with NumPy, SciPy and pytest installed, write fresh receipts to a separate directory and run the independent word-evolution, detector-control and archived-count checks.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
python -B -W error -m model.candidate_laws.interval_echo.openwave.row_23_coupling.timed_scattering.measure --output /tmp/oph-timed-scattering
python -B -W error -m pytest -q model/candidate_laws/interval_echo/openwave/row_23_coupling/timed_scattering
```

The NPZ keys give source, flight contrast, delay and analyzer. Single-arm arrays contain the 214 recoil bins, spill, then source rejection. Each recombination array contains 214 plus-port bins, 214 minus-port bins, spill, then source rejection. The quarter-turn analyzer uses $+\pi/2$ in the late-arm amplitude, so its port difference reads the negative imaginary part of the overlap.
