# When two local modular clocks give a joint clock

Alice and Bob can each assign the same modular time to an accepted pulse of their local hand. If they combine records, correlations in the state they used for calibration can invalidate that assignment. On this 20-cell experiment the separate clocks agree in both preparations. The joint clock prediction succeeds when the reference hands are prepared independently. One unread fair exchange of two different references changes the joint prediction. The physical clock pulses are the same in both tests. [Protocol](PROTOCOL.md) · [All simulated cohorts and forecasts](records.json)

## Reading the reference and the hand

The source and detectors are the native protected two-rail hand and third-rail contact apparatus of the [two-owner exchange experiment](../README.md). Each source has one double per admitted hand. A 32-phase average makes its reference diagonal in the hand's swap-parity basis. In the shared preparation a native fair bit selects either a hold or an exchange of differently prepared hands. In the independent control each owner obtains its hand from a separate realization of this service. Enrolling both attempted heralds gives acceptance probabilities 3/4 and 9/16, respectively. The protected routes, independent preparation of the inputs, fair service and coherent phases are supplied. The calibration and clock trials are separate.

For the two swap parities the third-rail plus probabilities are 5/8 and 1/8. Their two-by-two response is invertible. Joint terminal counts therefore determine the four reference populations within this prepared diagonal family. In the order ++, +−, −+, −−, the shared reference has weights (5/12, 1/4, 1/4, 1/12). Both local marginals are (2/3, 1/3). Independent references with these marginals have weights (4/9, 2/9, 2/9, 1/9).

On fresh probes Alice starts with a localized hand and Bob with an independently phase-averaged localized hand. Sixteen fair native services publish a count. Each accepted service rotates both hands by π/16. The two final readouts are Alice's hand position, with either of two analysis phases, and Bob's third-rail bit. Both bits, all 17 count groups and all misses are archived. A separate control retains the contrast term during the exposures and removes exchange. The calibrated clock then stops in every count group. [Native-contact acquisition and count decoder](measure.py)

## What changes when the records are joined

Let a,b,c,d be strictly positive normalized joint populations in the two-hand energy basis, ordered by rows. Under the joint modular unitary, Alice's two possible transition rates are log(a/c) and log(b/d), depending on Bob's parity. They agree exactly when ad = bc. For normalized binary tables this is exactly the condition that the joint state factor into the product of its marginals. In that case each common conditional rate is its own marginal rate. This finite criterion concerns continuous modular flow on a diagonal reference family. A finite set of detector readings may also agree at isolated phase recurrences. [Lean proof](../../../../../Lean/Time/JointClockReference.lean)

At eight accepted pulses the observed covariance of Alice's unshifted hand bit and Bob's detector bit was −0.0000175 ± 0.001203. The shared reference predicted +0.07465 ± 0.00304. Its two inferred rates relative to the marginal rate were 0.7374 and 1.5835. Across counts 4, 8 and 12 and both phase settings, 22 of 24 joint probability intervals excluded that reference's prediction. The independently prepared reference gave 24 of 24 compatible predictions. Those entries share cohorts; the bounds apply simultaneously. The programs use signed componentwise propagation from binomial Hoeffding bounds, including sampling quantization, with at least 99% family coverage. The calibration uses only the simulated reference counts and the specified analyzer. [Independent full-matrix unitary and control checks](test_joint_time.py)

This determines when local modular timing extends to the combined two-hand algebra for the references tested. A physically selected reference and an attachment to proper time are premises for using that timing beyond the protected hands.

From the repository root run:

```sh
python -m model.experiments.interval_echo_clock.reference_exchange.joint_time.measure
python -m pytest -q model/experiments/interval_echo_clock/reference_exchange/joint_time
```
