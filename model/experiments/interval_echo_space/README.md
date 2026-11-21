# Three spatial response components read by neighbouring patches

On the intact 20-, 80- and 320-cell screens, three neighbouring one-cell observers can jointly distinguish three independent variations in the sender’s local load distribution. Their recorded counts identify all 144 held patterns across nine sender sites and two exposure settings. At three additional sender sites, a ten-cell patch predicts the neighbours’ measurements at three untrained exposures. Its largest error in an individual detection probability is 6.32 × 10⁻⁷. [Response measurements](result.json) · [Held station forecasts](forecast.json)

## Preparation and readout

The free generator is the registered Interval Echo one-load generator on all ports and seams. A sender prepares one sharp load per trial. Twelve anonymous source ports define the possible preparations. A recipient reads one of its twelve ports or its quiet outcome. A distinct ensemble is prepared for each recipient. The joint reader pools the three local inboxes by reported preparation and exposure. The three recipient readings on a row come from independently prepared ensembles. Each source, neighbourhood and encounter has its own finite aperture. Cohort pooling and estimation are supplied analysis. [Candidate implementation](../../candidate_laws/interval_echo/candidate.py) · [Acquisition and inverse](instrument.py)

At a short calibration dwell the source measures the probability table for its own ports. The gap between two sets of off-diagonal probabilities gives the internal incidence. Its first nonconstant Laplacian eigenspace has dimension three. The source prepares small positive and negative variations of the uniform port probabilities along an orthonormal frame of that measured subspace. The three signs in a held word select a direction between the calibration axes. None of the preparations requires a coherent multiport source. The simulation uses complete Born probabilities for the short calibration and separately draws finite counts for the indicated held tests. No embedding coordinates enter the inference.

The lab supplies patch membership, addressed ports, the mixtures, the native-time exposure and local detection. Between preparation and detection the excitation evolves freely. The elapsed-word clock is idle. In the one-load sector the generator is the federation Laplacian divided by two, plus 6κ times the identity. The constant contributes a common phase. The same measurements therefore apply to the contrast-memory reference and are insensitive to κ in this sector. [One-load reduction and band](../../evaluation/repair_band.py) · [Specification](../../candidate_laws/contrast_memory_spherical_v1/SPEC.md)

## Direct test of three independent responses

At each of the nine sender sites, each recipient's calibration includes the uniform preparation and the six signed single-axis perturbations. Eight different source mixtures are held out for each of two exposure settings. The rank-two control keeps the two singular components of the calibrated receiver response that minimize its squared error. A separate control supplies an unmodulated, uniform source for all intended words. All counts, including quiet, remain in the receiver distributions. [Frozen protocol](PROTOCOL.md) · [Exact probabilities and finite counts](records.npz)

| Native exposure | Sender sites | Smallest / largest singular value of joint response | Relative held-signal RMS error after retaining two components | Decoded patterns |
| --- | ---: | ---: | ---: | ---: |
| 2 | 9 | 0.248–0.499 | 17.7–40.3% | 72/72 |
| 4 | 9 | 0.205–0.286 | 16.3–24.7% | 72/72 |

Each local setting uses 4,194,304 simulated preparations. The sham decodes 16 of 144 words, with 18 expected for uniform guessing. The exact three-component interpolation agrees within 3.4 × 10⁻¹⁶, as expected from the linearity of mixtures. The physical tests are the observed third singular component, its effect on held predictions, and decoding with independent counts. Pooling matters. Only one of the 54 individual recipient panels has a smallest-to-largest singular-value ratio of at least 0.05. All 18 three-recipient panels meet that threshold. A receiver's probability of detecting a load ranges from 0.0372 to 0.0937 over the held settings. The operator fit and the held predictions use the absolute, unconditional counts. [Per-site results](result.json)

## Forecast at three additional sender sites

A ten-cell observer patch consists of the sender and two layers of neighbours. One very short survey at external dwell 1/32 gives the conditional port readings in that patch, including the quiet complement. The reader infers contacts from the off-diagonal readings and total port degrees from the return probabilities. It builds the finite Hamiltonian from these local observations and the single-cell hop coefficient 0.5 inferred at exposure 1/8. The simulator keeps all external links during the actual survey and the subsequent flight. The finite forecast uses the projected Hamiltonian on the surveyed patch. Propagation out of the patch and back can produce a discrepancy. [Held forecast protocol](FORECAST_PROTOCOL.md) · [Forecast calculation](forecast.py) · [Independent calibration](../../evaluation/patch_spatial_readout.json)

An incidence-isomorphism check preserves total port degrees and the sender, inner and outer patch roles. The three sender cells have the same marked geometry. Adding immediate neighbours distinguishes the 20-cell site, while the 80- and 320-cell stars agree. The complete ten-cell surveys distinguish all three sites. Their similar responses reflect shared local structure. [Independent geometry and forecast checks](test_forecast.py)

The three enrolled senders are cell 10 on the 20-cell screen, cell 40 on the 80-cell screen and cell 160 on the 320-cell screen. Each prediction uses only the short survey, preparation, declared form of the walk and exposure. All three receiver cells at each site, every held exposure and all eight words enter the comparison. The exact conditional probabilities used by the survey are a supplied input to this test. Its finite-count check draws 1,048,576 preparations per held setting per recipient.

| Held dwell | Largest absolute error in any port or quiet probability across boards | Largest relative RMS error in the three-component signal | Relative RMS error with two components | Patterns decoded |
| --- | ---: | ---: | ---: | ---: |
| 0.75 | <3 × 10⁻¹² | <2.0 × 10⁻¹⁰ | 18.76% | 24/24 |
| 1.5 | <4.3 × 10⁻⁹ | <2.4 × 10⁻⁷ | 18.50% | 24/24 |
| 3 | <6.32 × 10⁻⁷ | <2.22 × 10⁻⁵ | 18.04% | 24/24 |

[Held station measurements](forecast.json) · [Surveys, forecasts and counts](forecast_records.npz)

The measured directions describe three components of a prepared occupation distribution. Their transmission is detected across one intercell step. Establishing an effective spatial metric requires independent ranging and elapsed times on material patch clocks. A continuum comparison also requires matching physical regions across refinements. Neither experiment measures transport of an arbitrary coherent quantum state in a three-component wave band.

## Reproduce

From the repository root, both commands regenerate the source, probabilities, counts and analysis and compare them to the retained files.

```sh
python -m model.experiments.interval_echo_space
python -m model.experiments.interval_echo_space.forecast
python -m pytest -q model/experiments/interval_echo_space
```

Each command accepts `--output DIRECTORY` to write an independent copy. In `records.npz` keys identify the board size, sender, dwell and record type. The last receiver outcome is quiet. `forecast_records.npz` adds the local survey, model forecast and uniform-source records. The JSON files include the relabellings, detector cells, data selection and protocol hashes.
