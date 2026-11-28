# One outside head reads the current of the arriving charge

An observer at 1045 can infer the material current on 1053→1052 from the direction of the outside head on 1045–1052. On the late arrival with the source report unread, the exterior estimate is **[−0.021844, −0.021226]** in native current units. An independent counter on the material side gives **[−0.021806, −0.021326]**. The inference holds on every one of the 30 wheel configurations and their coherent combinations. It needs neither an angular label nor access to the internal material counter during the outside experiment. [Prospective procedure](DESIGN.md) · [Count banks](records.json)

## What observers do

The same outside detector has also been run on the negative-charge arrivals that first passed a calibrated guide rest-line measurement. The origin read preserves the two source words. A native four-face handoff restores the required wheel and exterior heads. [Sequential measurements and calibrated collector](COLLECTOR.md#what-the-returned-charge-and-the-outside-observer-establish-together)

On separate launches the source at 98 writes one member of a fresh reference pair. At 1045 the observer heralds an arrival in either of two face words without resolving which one arrived. Five native pulses load that word onto the wheel. A calibrated energy filter can select its stationary angular doublet. The two adjoining outside arrows retain their arrival directions. An observer flips a fair setting bit and opens the single rim contact for $\pi/16$ or $3\pi/16$. An idle wait equalizes the elapsed times. The observer then opens just the adjoining face for $\pi$ and reads the original outside head. An agreement is a changed head after the short pulse or a preserved head after the long pulse. Each accepted event gives one bit. [Written preparation and energy filter](../README.md) · [Physical contact and flip](instrument.py)

In the positive chart the two initial outside arrows follow 1052→1045→1053. The face can flip exactly when the shared head points to 1052. Complementing the three arrows and using vacancy transfer gives the same change-of-head detector for the negative arrival. For an enabled material pair write $a$ for the source amplitude at 1053 and $b$ for the destination amplitude at 1052. A positive quarter-cycle gives destination amplitude $(b+ia)/\sqrt2$ up to overall phase. Three quarters gives $(b-ia)/\sqrt2$. Their probability difference times the calibrated contact rate $t=4$ is $2t\,\mathrm{Im}(a^* b)$, precisely the oriented material-number current in the continuity equation. Disabled words contribute equally in both arms. The exact identity for the complete arrival sector is

\[
  \langle J_{1053\to1052}\rangle
    =4(p_{\rm short,flip}-p_{\rm long,flip})
    =4(2p_{\rm agree}-1).
\]

The Lean proof sums this identity over unresolved words or reference outcomes. The microscopic replay independently constructs both pulse arms, checks every entry of the two detector effects, and reads the current through a separate source/destination counter. Halving only the outside-face rate at the same wait halves the contrast. Closing the face makes the agreement probability exactly one half. Stationarity and material conservation equate the currents on all five rim sides in the isolated doublet. [Finite proof](../../../../../../../../../../Lean/SourceLaw/ExteriorHeadCurrent.lean) · [Operator checks and controls](measure.py)

| Late source record at 98 | Current from the outside head | Independent material counter |
| --- | --- | --- |
| Positive reference report | [−0.048887, −0.047991] | [−0.048856, −0.048147] |
| Negative reference report | [+0.010848, +0.011867] | [+0.010962, +0.011776] |
| Report unread | [−0.021844, −0.021226] | [−0.021806, −0.021326] |

The late gate admits 0.1305953% of original launches. The table uses simultaneous 99% simulated-count intervals rounded outward. Counts before the energy selection are included in the receipt. The native current sign refers to motion of the material excitation. Its electric current would also carry its electric charge.

## Comparing the exterior detector with the original reference

Within the energy doublet, three settings of this same outside head, together with normalization, span the four real Hermitian observables. In the first setting the observer opens only the adjoining face. In the second, the observer first applies full material exchanges on 1048–1053 and 1053–1052. The third setting uses the randomized signed-current read. The least singular value of the calibrated four-row matrix exceeds **0.01179773**. Independent source axes at 98 give nine joint settings. Full rank discharges the span premise of the finite observer theorem. Every linear comparison on the doublet then has a decoder in these outside records. The aperture-closed control has rank one. [Finite observability theorem](../../../../../../../../../../Lean/Geometry/SourceTemporalObservation.lean) · [Detector calibration and interval calculation](record_banks.py)

The chosen decoder measures a partial-transpose witness. Its expectation for every separable preparation is nonnegative. The direction was fixed by the previously published incoming face/reference witness and the calibrated transfer. Count banks on the late arrival give **[−0.01018938, −0.00780364]**. The same comparison at the early appointment gives **[0.09564351, 0.09801931]**. Randomizing the late source reference gives **[0.20649213, 0.20887788]**. The intervals are simultaneous at 99% across the three witness comparisons and carry the arrival, angular-frame and numerical allowances. Both the energy gate and its refusals precede the analyzer settings. These trusted measurements certify entanglement between the selected material-and-rim record before interrogation and the reference kept at 98. The final exterior read is destructive and returns a classical bit. [Source witness](../../../WITNESS.md) · [Nine-setting receipts](records.json)

These are finite records of the written-field extension. This experiment closes the observability premise for the stationary current and its source correlation. The electron spin magnetic moment requires further attachment. In the existing Pauli square theorem, three kinetic momenta, a magnetic field, their calibrated commutators and the spin factor occur on one construction. A single circulating material side gives only one of the required responses. [Pauli square](../../../../../../../../../../Lean/EventAlgebra/PauliMinimalCoupling.lean)

## Replay

From the repository root run

```sh
python -B -W error -m model.candidate_laws.interval_echo.openwave.electron_field.fan_path.local_launch.material_loop.five_faces.exterior_meter.measure
```

The current detector uses $2^{42}$ original starts per comparison and the witness uses $2^{45}$ per setting. The programs draw these banks from the finite model probabilities. The source at 98, face at 1045 and the guide are the same throughout this test. The probe pulses, energy appointment and detector choices are supplied.
