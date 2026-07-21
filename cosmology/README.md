# OPH Cosmology And Dark Gravity

This directory develops the cosmological continuation of Observer Patch Holography: repair-charge dark gravity, finite-source primordial structure, observer-screen synchronization, Boltzmann transport, and direct comparison with public cosmological data.

The main OPH papers establish the structural observer, geometry, gravity, and matter branches. These focused companions carry the additional source laws, transfer maps, calibration, and likelihood contracts needed to turn that structure into cosmological observables.

## Paper Map

| Paper | Main contribution |
| --- | --- |
| [Dark Gravity](oph_dark_matter_paper.pdf) ([source](oph_dark_matter_paper.tex)) | Repair-charge condensate action, dust-like normal phase, deep-galaxy radial-acceleration branch, and coherent-source coupling |
| [Finite-Source CMB Program](oph_cosmology_finite_source_cmb_program.pdf) ([source](oph_cosmology_finite_source_cmb_program.tex)) | Conditional source-screen spectrum, finite receipt contract, transfer requirements, and CMB promotion path |
| [Inflation Without an Inflaton](oph_inflation_without_inflaton_observer_screen_synchronization.pdf) ([source](oph_inflation_without_inflaton_observer_screen_synchronization.tex)) | Observer-screen synchronization, conditional edge-center tilt, one-shell radial obstruction, physical dilation and tomography routes, flatness, and horizon coherence |
| [Cosmological Vacuum And Structure Formation](oph_cosmological_vacuum_and_structure_formation.pdf) ([source](oph_cosmological_vacuum_and_structure_formation.tex)) | Vacuum boundary, fluctuation ensembles, proto-objects, worldlines, and structure seeds |
| [Cosmology Data And Likelihood Contracts](oph_cosmology_data_likelihood_contracts.pdf) ([source](oph_cosmology_data_likelihood_contracts.tex)) | Boundary between the conditional source theorem, finite source instantiation, transfer, nuisance treatment, and official likelihood comparison |
| [Boltzmann Transport Derivation](oph_boltzmann_transport_derivation.pdf) ([source](oph_boltzmann_transport_derivation.tex)) | Finite transport interface between OPH sources and observable distribution functions |

The [radial-lift theorem note](RADIAL_LIFT_SOLUTION.md) gives the one-shell
non-identifiability proof, the physical source-dilation theorem, radial
cross-covariance tomography, the exact amplitude conversion, and finite-window
bound. The [simulator contract](SIMULATOR_RADIAL_CONTRACT_330.md) specifies the
fail-closed evidence split used by finite runs. The
[coverage audit](RADIAL_LIFT_COVERAGE.md) maps each radial result to its formal
paper theorem and supporting artifact.

The [physical CMB theorem program](physical_cmb_theorem_program.md) collects the remaining source, lift, stress, abundance, transfer, and likelihood obligations in one place.

## Dark-Gravity Structure

The dark-sector action yields a pressureless dilute background and, through its cubic condensed link energy, a spherical deep-acceleration law with baryonic Tully–Fisher scaling. The same action supplies a repair current, stress channel, and coherent-source coupling weighted by $\chi_\nu^{\rm can}$.

Physical promotion requires the canonical repair pair, a relativistic constitutive completion, abundance and lensing maps, CMB and cluster transfer, Solar-System response, and calibrated laboratory receipts to close on the same branch. Existing analytic and numerical comparisons are evidence for developing those maps; the papers state their precise status.

## Reproducibility

Companion code lives in [`../code/dark_matter/`](../code/dark_matter/). The larger simulator and visualization surfaces are:

- [OPH physics simulator](https://github.com/muellerberndt/oph-physics-sim)
- [Interactive simulation](https://simulation.floatingpragma.io)

The cosmology papers are focused research companions and are not part of the core release bundle.
