# Finite amplified source interval certificate

All 256 signed-source, resolution, orientation, turn and law cases pass the
independent outward calculation. The actual finite mixed state and all four
signed native-clock weights are retained. Native selected-clock probability
is at least 0.03305905697 before the declared service perturbation, and at
least 0.03305142757 with it. The finest pair-energy envelope is contained in
[1020.9461578383,1022.4156829702] keV under the single limiting calibration.

The [proof](PROOF.md) states the full source, joint operator, calibration and
error hypotheses. [SOURCE_INPUTS.json](SOURCE_INPUTS.json) binds the actual
calibrated source columns, finite phase words, complete producer index and
two independently admitted sources. No producer module is imported by the
interval calculation. Its counter/arithmetic ancestry is explicit in
[ANCESTRY.json](ANCESTRY.json). The source-blind inverse reader receives no
target frequency; that value is checked only after decoding.

Run the complete 32-case offline headline from this directory:

```sh
python3 -B certificate.py
```

It uses Python 3 and mpmath 1.3.0, reads only local bounded inputs, writes
nothing and compares the saved [HEADLINE.json](HEADLINE.json). The first
execution took 12.295 seconds. The full 256-case confirmation took 71.107
seconds and is separately replayed with:

```sh
python3 -B certificate.py --k 0 1 2 3 4 5 6 --verify CERTIFICATE.json
```

Only explicit `--record-to NEW_FILE` creates evidence, exclusively. Failed
scientific predicates remain active under optimized Python. The full source
and independent reviews retain all implementation failures and mutations;
the outward certificate does not replace those complete-instrument checks.

The 40-bit single-source and 52-bit pair service envelopes are conditional on
the declared complete-service and joint-input error bounds. They are not
claims of measured hardware accuracy. Rare control ports retain their costs
and nominal clocks without a precision guarantee larger than their
acceptance. The exact finite record laws are distinct from sampled data.
The calculation is same-group verification and does not pass external G6 or
authorize a row award or submission.
