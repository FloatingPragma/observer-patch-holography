# OPH Dark Matter Supplement Code

This directory contains the reproducibility code for the supplemental paper:

- [Observer-Patch Holography and the Dark Matter Phenomenon](../../cosmology/oph_dark_matter_paper.pdf)

The package is a pre-likelihood empirical surface. It reproduces the scorecard
rows in the supplement from public SPARC tables, the OPH screen-capacity
calculator, the repair-matrix diagnostic, and a compressed CAMB comparison.

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
the supplement.

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
