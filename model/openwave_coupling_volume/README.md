# Colour coupling on a larger box

Two additional colour ensembles test the native loop-coupling trend on
side10, alongside the complete side6/8 records. The same bare action
and loop scheme are used throughout. Correlated block errors, hot/cold
starts, split histories and temporal-correlation diagnostics expose finite
volume and sampling effects. Physical scale and continuum running require
separate matching.

Analyze with `python3 -m model.openwave_coupling_volume --write` or `--check`.
Regenerate a raw chain with the pinned runtime from the parent package:
`python -m model.openwave_coupling_volume.production --side 10 --seed 4101 --start cold --warmup 400`.
The other retained new run uses side10/hot/4102/warmup500.
Both use400 records and beta6, adding800records to2400total.
The [proof](PROOF.md) states the statistical scope.
