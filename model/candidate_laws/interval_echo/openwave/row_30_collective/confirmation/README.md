# OW-30 confirmation on the coherent open strand

The candidate is a reversal in the alternating orientations of loaded dimers after controlled commissioning. Separate local measurements of exchange and conditional link phase predict its frequency as a function of the spatial phase. The observer can challenge a reported wavelength with fresh receiver readings on the same board. All counts, including the complementary signs, enter the calculation. [Candidate declaration](candidate.json) · [Acquisition terms and prior exposure](TASK.md) · [Spatial reader](spatial.py)

At each stimulated mode pair, the detector uses spatial records at even-numbered holds to fit the two wave numbers. The odd holds test the predicted spatial shape at every receiver. These even-trained wave numbers enter every frequency comparison. The same odd records score this prediction and the wavelengths from the sine templates used to launch the waves. The receipt also gives the earlier interior-recurrence reading on the same records. Passing requires the receiver profiles to predict the shape and the source-only substitution to miss. The fixed end of the strand and the two free amplitudes per scan are disclosed apparatus assumptions. Counts from the same local clock and separate temporal counters set the demodulation frequency in both halves.

The formal acquisition requires OpenWave's affirmative decision on the exact state and decoder, the same canonical RER pin throughout the registration's first-parent history, and an executable verified against both repositories. The sources of the retired comparison remain available in the parent package. The prospective procedure freezes the entire experiment and preserves each release's draw, numerical verdict, measurements and count hashes. [Binding checks](binding.py) · [Acquisition and replay](acquire.py) · [Manifest](files.json)

Run with the pinned NumPy and SciPy versions in `model/requirements.txt` and a clean checkout of the source commit. From that checkout, the command prefix is

```sh
V=model/candidate_laws/interval_echo/openwave/row_30_collective/confirmation
python -I -B "$V/run.py" check-sources
python -I -B "$V/run.py" status --openwave /path/to/openwave
python -I -B "$V/run.py" acquire --mode pilot --release /path/to/pilot
python -I -B "$V/run.py" replay --mode pilot --release /path/to/pilot
```

The public port uses the same launcher bytes in `openwave/xperiments/m10_oph/research/scripts/m10_5/run.py`. Pass `--source /path/to/the/frozen-rer-checkout`. Inspect a saved release with `inspect` and the corresponding mode. To run the exposed v2 geometries with the v3 reader, select `--mode repeat`. For the three fixed v3 local validation selectors, replace `--mode pilot` with `--mode validation` and use a separate directory. The public mode takes its selectors exclusively from authenticated registration history.

| Exit | Meaning |
| ---: | --- |
| 0 | All physical gates passed, or the requested provenance check succeeded |
| 2 | Complete measurements or enrolled geometry rejected the claim |
| 3 | Acquisition or required records incomplete |
| 4 | Invalid binding, damaged data or disagreement on replay |
| 5 | The exact candidate and terms await an upstream decision |

The exposed source comparisons and supplementary ring results have different preparations. The prospective public claim concerns the specified open-strand source and its decoder. [Counts and full instrument audit](../../submission/row_30_cautions/README.md)
