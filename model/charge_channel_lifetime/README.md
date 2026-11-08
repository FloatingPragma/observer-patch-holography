# How long a created charge signal persists

A local two-outcome measurement creates two internal-charge branches on
the selected eight-site source. Their return probabilities first fall by
1% at source times about 0.172621 and 0.104457. Exact charge conservation
also supplies a local-detector contrast that persists after that overlap
is lost. Every herald, detector complement and source level is retained.
The [proof](PROOF.md) and [certificate](results.json) distinguish this
finite retention measurement from physical particle decay; the continuum
particle and resonance identification is not established.

Run `python -m model.charge_channel_lifetime --check` for the retained
certificate, or pass a receipt path to `--check PATH`. The producer uses
exact source matrices, complete spectral residual bounds and rational
Taylor remainders. `python -m pytest model/charge_channel_lifetime -q`
checks the source attachment, instrument, echo normalization, survival
intervals, conserved detector and fail-closed replay.
