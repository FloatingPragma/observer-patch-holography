# Independent review of the proposed collective matter-wave measurement

**Superseded release assessment.** A subsequent audit tested the ring source, memory commitment, spatial ruler and verifier error handling. The earlier release recommendation is withdrawn. [Subsequent audit](../../submission/row_30_cautions/README.md)

The local records support a positive long-wavelength dispersion measurement for the traveling orientation defect on a commissioned one-dimensional strand. The public comparison cell is untested. A registered run on fresh public settings and a maintainer decision are outstanding. [OpenWave row and its explanation](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L467) · [Criterion](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/MODELS.md#L136)

## Independent numerical check

I read the 19 measured Bloch wave numbers and pair-inferred frequencies from the local holdout receipt. I combined the measured wave numbers with the two local meter coefficients without using the stored frequency predictions or scores. The ordinary spatial phase gives the continuum predictor √(m² + c²k²). A nearest-neighbour difference gives the lattice predictor √(m² + 4c²sin²(k/2)). The table gives the largest absolute frequency errors in local clock radians per unit of clock time. The meter intercept is 15.0005108939 and the squared propagation coefficient is 16.0015591352 in dimer-spacing units. [Raw local calibration](../coupling_meter.json) · [Measured frequencies and wave numbers](../holdout.json)

| Sphere cells | Distinct measured k | Distinct k at or below 0.50 | Lattice, all | Ordinary k, at or below 0.50 | Ordinary k, all |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 20 | 5 | 3 | 0.000526 | 0.001959 | 0.011567 |
| 80 | 7 | 4 | 0.000584 | 0.002261 | 0.016665 |
| 320 | 7 | 4 | 0.000567 | 0.001608 | 0.011260 |

The exact native coefficients m=15 and c²=16 also give a uniform continuum error bound. Since 0 ≤ k² − 2(1 − cos k) ≤ k⁴/12, rationalizing the two square roots gives 0 ≤ ω_cont − ω_lattice ≤ c²k⁴/(24m). This is at most 1/360 clock units for |k| ≤ 0.50. The bound concerns the dispersion functions on this strand and excludes sampling error.

At low k, all eleven qualifying lines pass the proposed 0.0048 tolerance. A tangent parabola also passes that absolute low-k tolerance. The distinction comes from measuring curvature across all wavelengths on the same lattice ruler. Centering each forecast and the readings on the first enrolled frequency gives tangent-to-square-root RMS error ratios of 43.94, 31.43 and 20.85. These are ratios of frequency errors, not confidence levels. The open strands have positive minimum wave number. The meter-predicted bulk intercept is an extrapolation to zero wave number. A separate count-only check on the archived development bank used an eightfold zero-padded Fourier peak followed by continuous frequency refinement. Three pair sums gave individual frequencies of 15.020314, 15.080336 and 15.177368. [Representative count bank](../example_counts.npz) The routing, material preparation, parity-allowed pair readout and supplementary packet measurements set the scope. [Local holdout](../HOLDOUT.md) · [Packet records](../packets.json)

## Source custody

A clean checkout at the published author base passed the frozen-source comparison. Changing the candidate wrapper was sufficient to reject the checkout. In a separate fixture with two local Git remotes, the verifier accepted the first exact landing with the correct source pointer. Withdrawing and restoring the identical bytes barred acquisition. The verifier also rejected a pointer to public source whose task text differed from the upstream task. Those fixture tests exercised custody without requesting any draw or running an experiment. The checks against the actual public repositories returned no eligible registration. [Guard](verify.py) · [Prospective terms](TASK.md)

The frozen material and receiver use the adopted Interval Echo generator at 4fe10f9a. The proposed state and decoder have their own task and leave the row-4 proposal under review separately. The maintainers adjudicate this scope. [Canonical pin](https://github.com/openwave-labs/openwave/blob/d4fc3e4046c472f4ad38468df8dddfa1f97bac97/openwave/xperiments/m10_oph/research/m10_theory_canonical.md#1-specification-of-record-and-the-pin) · [Row-4 proposal](https://github.com/openwave-labs/openwave/pull/749)

## Verdict on the staged registration

PASS for independent review of the proposed M10.5 registration on 2026-11-24. The audited task SHA-256 is `748eecfa31fbd3c1d08851896bec4dece665683cf38a8eb893987c38df720b9f` and the bound source-manifest SHA-256 is `3789a7be528ab9d20471486816df3f3ee39b003eb54c09eab3cd50d65629d4c9`. The proposal assigns the numbered task in the M10 roadmap. The platform implementation is staged separately under that task ID. Its entry point checks the source bytes, registration ancestry and working task before acquisition, runs the headless experiment from the platform's vendored modules, and checks the module origins. The formal acquisition and official cell decision are outstanding. The owner controls release of this version. [Platform entry point](../platform/run.py) · [Local registration instructions](README.md)

The separate ring study selected two multicell rings after inspecting all twelve local matching patterns. For these distinct closed preparations, three parity-allowed pair lines reconstruct a constituent frequency of 15.000003 and 14.999997 at zero wave number. Halving the free-flight contrast gives 7.000012 and 6.999983 in the same detector. These observations extend the local evidence, with the geometry selection and pair reconstruction disclosed. They contribute no prospective open-strand settings. [Ring study, source and limitations](../rest/README.md)
