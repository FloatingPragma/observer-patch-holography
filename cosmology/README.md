# OPH Cosmology Papers

This directory holds OPH cosmology papers inside the public core repository.

## Released Papers

- `oph_dark_matter_paper.tex`
  - included in the normal paper release pipeline
  - owns the dark/anomaly stress branch, finite covariant parent contract, galaxy limit, cluster/cosmology contracts, and the dark-sector source interface consumed by later CMB/LSS work

## Staged Papers

These files are public working papers but are not wired into the normal release build yet.

- `oph_cosmology_finite_source_cmb_program.tex`
  - umbrella program paper for finite-source cosmology, physical CMB gates, scale bridge, simulator receipts, and claim boundaries
- `oph_inflation_without_inflaton_observer_screen_synchronization.tex`
  - theorem template for the inflation-free branch: flatness, horizon coherence, geometric screen spectrum, screen release amplitude, radial lift, and hot initial data
- `oph_cosmological_vacuum_and_structure_formation.tex`
  - theorem template for OPH-native vacuum gates, fluctuation ensembles, proto-object/worldline formation, and structure-seed receipts
- `oph_cosmology_data_likelihood_contracts.tex`
  - technical companion template for frozen source artifacts, no-data-use receipts, pooled reducers, Boltzmann transfer, and official likelihood comparisons
- simulator reference: https://github.com/muellerberndt/oph-physics-sim
- visualization companion: https://oph-universe-explorer.lovable.app

The staged CMB program distinguishes diagnostic proxies, conditional physical artifacts using a
frozen imported FLRW packet, and OPH-native physical artifacts derived from the quotient carrier.
Flat-sector work is likewise split: spatial Levi--Civita holonomy identifies the \(\kappa=0\)
branch, while direct theorem, conditional CMH, or explicit-assumption labels select it.

## Paper Split

Two cosmology papers would be too compressed. The split is five papers:

1. **Dark sector and structure**: the released dark-matter paper, including galaxy phenomenology, transported stress, cluster behavior, and the dark/anomaly source contract.
2. **Finite-source CMB prediction program**: staged program paper focused on source-only inputs, scale calibration, Boltzmann transfer, frozen likelihoods, and simulator receipts.
3. **Inflation without inflaton**: staged theorem paper for observer-screen synchronization, flatness, horizon coherence, geometric screen spectrum, screen release amplitude, radial lift, and hot initial data.
4. **Cosmological vacuum and structure formation**: staged paper for OPH-native vacuum gates, fluctuation ensembles, proto-object/worldline formation, and visualization-facing receipts.
5. **Data and likelihood contracts**: staged technical companion for CMB/LSS/BAO/SN/WL/growth comparison protocols.

Only item 1 is release-published from this directory. The staged program paper may cite and coordinate the later tracks, but it should not become a release-bundle paper until the theorem gates close.
