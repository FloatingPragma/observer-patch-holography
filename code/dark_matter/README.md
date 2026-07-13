# OPH Dark Matter Supplement Code

This directory contains the reproducibility code for the supplemental paper:

- [Observer-Patch Holography and the Dark Matter Phenomenon](../../cosmology/oph_dark_matter_paper.pdf)

The package is a pre-likelihood empirical surface. It reproduces the scorecard
rows in the supplement from public SPARC tables, the OPH screen-capacity
calculator, the repair-matrix diagnostic, and a compressed CAMB comparison.

The corrected cross-repo RAR, BTFR, Cassini, and analytic scalar-tilt audit is
summarized in
[`ROUND_ONE_PUBLIC_DATA_AUDIT.md`](ROUND_ONE_PUBLIC_DATA_AUDIT.md). Its
maintained error-aware likelihood and external-field implementations live in
the sibling `oph-physics-sim` repository and are hash-bound by the generated
best-of-public-data bundle. This separation prevents the older generated
scorecard from silently inheriting corrected numbers that it does not compute.

Neutrino masses are external cosmology inputs on this surface. The default is
the conventional minimal-normal reference sum, not an OPH prediction. The old
`0.0900119296` electron-volt weighted-cycle sum remains available only as the
explicit `rejected_weighted_cycle_compare_only` CAMB scenario; it is
target-informed, rejected by the NuFIT 6.1 correlated profile, and blocked from
all promotion paths. Unless `--mu-eq` is supplied explicitly, the main
scorecard derives its diagnostic `rho_A/rho_b` ratio from the homogeneous run
using the declared external neutrino input; it no longer reuses the old ratio
indirectly.

## Layout

- `scripts/`: Python calculators and diagnostic likelihood scaffolds.
- `data/external/`: local copies of the public SPARC tables used by the galaxy tests.
- `data/observational_targets.json`: compressed cosmology comparison rows.
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

The generated Markdown claim-boundary file gives the compact empirical table used in
the supplement. A checked-in scorecard generated with the rejected weighted-cycle
input is invalidated until it is rerun with an explicitly declared external
neutrino scenario.

## Optional CAMB Dependency

The compressed CMB, BAO, growth, and `S8` rows use CAMB. If CAMB is missing,
install it in the active Python environment:

```bash
python3 -m pip install -r code/dark_matter/requirements-boltzmann.txt
```

## Scope

Implemented diagnostics:

- SPARC RAR comparison.
- Fixed and nuisance-profiled SPARC rotation-curve scaffolds.
- Z6/Poisson collar coefficient diagnostic.
- Flat capacity-saturated homogeneous anomaly diagnostic.
- Finite repair transition matrix and cluster timing gate.
- Compressed CAMB rows for `Omega_m`, `sigma8`, and `S8`.

Audited in the maintained companion comparison suite:

- fixed unit and conditional-Z6 RAR branches;
- error-aware BTFR slope and pivoted normalization;
- Cassini external-field quadrupole with benchmark reproduction;
- analytic `n_s=1-P/48` CAMB/Planck TT diagnostic;
- explicit invalidation of the legacy compressed-cosmology evidence row.

Paper-grade empirical tests require the modules listed in the supplement:
a conservative disk-potential solver, hierarchical SPARC likelihood, finite
collar-packet parent evaluator with stress and exchange-current closure,
physical-clock receipts for any promoted `Gamma_rec`, CAMB or CLASS anomaly
module, cluster map forward model, joint likelihood runner, and reproducibility
harness with fixed data hashes.

## License And Patent Policy

This code surface is part of the OPH public repository. See the main
[LICENSE](../../LICENSE) and
[OPH Open Use And Anti-Patent Covenant](../../PATENTS.md).
