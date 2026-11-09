# Momentum from a shared observation origin

A retained random origin selects where a complete history is observed on
the specified source. The same origin controls every step. Averaging these
complete records over an expanding bulk gives a translation-invariant
limiting observation state, with commuting momentum and nonnegative source
energy generators. The [proof](PROOF.md) keeps the original source,
reference and clock, and bounds the boundary loss of populated energy bands.
Nonzero momentum response and a particle dispersion relation require
further spectral information.

Run `python3 -m model.particle_primary_origin_momentum --check` to replay
the [exact receipt](results.json). Large source matrices and observation
registers are specified without allocating them.
