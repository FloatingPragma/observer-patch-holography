# OPH Dark Matter Supplement Code

This directory contains reproducibility and diagnostic code for:

- [Observer-Patch Holography and the Dark Matter Phenomenon](../../cosmology/oph_dark_matter_paper.pdf)

The package implements public-SPARC comparisons, screen-capacity arithmetic,
repair-matrix diagnostics, and CAMB wiring checks. These calculations test code
paths and empirical parameterizations. They are not predictions or likelihood
scores of the repair-charge condensate proposal.

Neutrino masses are external cosmology inputs, not OPH predictions. The
default is the conventional minimal-normal reference sum. The
`rejected_weighted_cycle_compare_only` scenario is target-informed,
incompatible with the NuFIT 6.1 correlated profile, and blocked from promotion.
Unless `--mu-eq` is supplied, the calculator derives its diagnostic
`rho_A/rho_b` ratio from the homogeneous run with the declared external input.

## Layout

- `scripts/`: Python calculators and diagnostic likelihood scaffolds.
- `data/external/`: local copies of the public SPARC tables used by the galaxy tests.
- `data/observational_comparisons.json`: external cosmology comparison rows;
  none is an OPH prediction target.
- `outputs/`: generated local scorecards and parent-grid JSON.
- `requirements-boltzmann.txt`: optional CAMB dependency list.

## Main Command

From the repository root:

```bash
python3 code/dark_matter/scripts/dark_empirical_scorecard.py --quiet
```

This writes:

```text
code/dark_matter/outputs/dark_empirical_scorecard.json
code/dark_matter/outputs/dark_parent_collar_grid_diagnostic.json
code/dark_matter/dark_empirical_implementation_status.md
```

The generated Markdown file records implementation diagnostics. Its CAMB and
cluster outputs have no physical score or prediction status. A declared
external neutrino scenario is required.

## Optional CAMB Dependency

The Boltzmann wiring check uses CAMB. It does not define CMB, BAO, growth, or
`S8` predictions. If CAMB is missing, install it in the active environment:

```bash
python3 -m pip install -r code/dark_matter/requirements-boltzmann.txt
```

## Scope

Implemented diagnostics:

- SPARC RAR comparison.
- Fixed and nuisance-profiled SPARC rotation-curve scaffolds.
- Z6/Poisson collar-coefficient arithmetic.
- Homogeneous scalar-load and CAMB wiring checks.
- Finite repair-matrix and timescale arithmetic.

Audited in the maintained companion comparison suite:

- fixed unit and conditional-Z6 RAR branches;
- error-aware BTFR slope and pivoted normalization;
- Cassini external-field quadrupole with benchmark reproduction;
- analytic `n_s=1-P/48` CAMB/Planck TT diagnostic;
- rejection of a physical score from compressed-cosmology wiring output.

Physical dark-sector tests require the repair-charge action, its source
normalizations, a relativistic stress tensor, source-derived abundance and
initial perturbations, disk and lensing solvers, a Boltzmann module, cluster
maps, and joint likelihoods. These derivations are work in progress.

## License

This code surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE).
