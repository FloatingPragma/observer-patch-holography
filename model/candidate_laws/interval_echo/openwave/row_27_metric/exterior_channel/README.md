# A freely moving source changes a uniform exterior probe population

On the 20-cell Interval Echo federation an observer uniformly illuminates all ten sites of an exterior radar belt. Four exchange units after release the marked receiving gate collects 13.92475% of unit probes for an occupied pole and 12.38075% for an empty pole. The same free mediator, uniform probe preparation, detector geometry and release appointment are used in both cases. Any passive model in which the source selects a lossless propagation on these ten orthogonal positions, including any positive mixture of such propagations, predicts 10% at every site. [Preparation menu](PROTOCOL.md) · [All end-site counts and certified populations](results.json)

### Which records the source changes

The preparation uses the two protected belts and occupied or empty pole of the native exterior timing experiment. A sharp mediator on its own belt is released from position 0 or 5, facing the receiver or the pole respectively. It evolves for twelve exchange units before the exterior probe launches. Both particles subsequently move. The observer chooses a sharp start uniformly among the ten exterior sites independently of the mediator. The final read resolves every site of the closed probe sector; the two particles cannot exchange belts. The pole lies three to six links from the surveyed exterior sites. The direct material coupling acts through the mobile mediator and its shared contact with receiver site 2. [Native preparation, contact and timing calibration](../EXTERIOR_PROTOCOL.md)

Let the columns of the joint wave matrix be the ten sharp probe launches. The marked-gate probability is one tenth of the sum of squared column entries with that receiving address, summing over **every** final mediator position. For the boundary-facing launch and wait four, integer Taylor propagation bounds the exact probabilities by

| Pole | Certified gate probability |
| --- | ---: |
| Empty | 0.1238073531 to 0.1238076048 |
| Occupied | 0.1392474282 to 0.1392476800 |

Thus the occupied-minus-empty increase exceeds 1.5439 percentage points before sampling. There is a further signed measurement at wait two. For the two primary launch addresses the source reduces the gate fraction by 2.1662 and 0.2429 percentage points in separate simulated samples. For a uniform choice among the other eight addresses it *raises* the fraction by 0.3023 points. The simultaneous 99% uncertainty for each occupied-minus-empty difference is ±0.1249 points. [Counts, source controls and uncertainties](results.json)

### A test independent of the clock calibration

Write the input density on the exterior belt as I/10. A source-dependent, possibly time-dependent Hamiltonian generates a unitary U. For every such propagation, U(I/10)U†=I/10. Integrating over random Hamiltonians, durations, bends or detector-independent classical source trajectories preserves this identity. Any fixed subsequent detector calibration observes the same source contrast, zero. This argument requires a closed ten-state receiver, sharp normalized site preparations and an input-independent distribution of the proposed propagators. The experiment excludes every such explanation of the **complete material signal** here. A shared geometric contribution could still coexist with probe-dependent interactions. The spatial distribution of a selected beam and the curvature of a common null congruence require additional measurements.

The joint two-belt evolution gives the changes directly. In a one-unit probe sector the protected contact contributes −32 to the joint diagonal energy when both particles are at their facing gates. When the mediator is parked its ten possible positions are conserved. Tracing over it gives a positive mixture of conditional probe unitaries. The control then preserves uniform illumination exactly, even if the stored mediator was coherent. A double probe has zero contact shift. Removing the contrast also gives zero. Both controls factorize. All recorded site frequencies in the three controls are within 0.0251 percentage points of 10%.

The complete mixture over **mediator** launch addresses provides a positive prediction for the mobile experiment. Independent uniform choices of both starts produce the joint I/100. For either pole state this is invariant throughout the joint run, including the preparatory twelve-unit mediator evolution. Let p₀(x) and p₅(x) be the observed terminal radar probabilities for the two primary sharp mediator starts. The predicted probability at site x when the source launch is uniform on the other eight addresses is

$$p_{\rm eight}(x)=[1-p_0(x)-p_5(x)]/8.$$

The source preparation and native unitarity fix this prediction without learning any of those eight responses. The largest absolute error across 800 held-site probabilities was 0.0230 percentage points. The simultaneous 99% radius is 0.0781 points. Each setting uses 2²⁴ independent simulation-generated preparations. The samples are multinomial aggregate draws from the finite dynamics. [Implementation](measure.py)

### Error control

The source and probe belt hopping amplitudes are −1/2 in the sharp-address basis. Irrelevant scalar energies cancel. Twice each Hamiltonian in the two certified gate calculations is an integer matrix. Both reconstructed matrices are also checked against the native contact routines. The interval calculation uses a 44-bit dyadic mesh, 128 substeps per exchange unit and sixteen Taylor terms per substep. The absolute row-sum bound gives a scaled operator norm at most one. At every Taylor term the complex vector has rounding error at most √(2d) times the mesh width in dimension d. Subsequent terms amplify this contribution by at most e. For an input of norm at most two the Taylor remainder is at most 6/17!. The sum over all substeps stays below 1/4, so unitarity inductively maintains the norm bound. The twelve-unit preparation error and four-unit probe error add. Multiplying their sum by three bounds the error in each projected terminal probability. The regression independently applies the sparse matrix exponential to the sharp starts. [Independent checks](test_channel.py)

Run the two commands at the RER root with the model environment active.

~~~sh
python -m model.candidate_laws.interval_echo.openwave.row_27_metric.exterior_channel.measure --write
python -m pytest -q -W error model/candidate_laws/interval_echo/openwave/row_27_metric/exterior_channel
~~~
