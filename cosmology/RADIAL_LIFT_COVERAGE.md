# Radial-lift coverage audit

The formal theorem surface is
`paper/tex_fragments/RADIAL_LIFT_THEOREMS_330.tex`, included by
`paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex`. The surrounding primordial
fragment supplies the physical scale, mode, source, ancestry, and forward
projection gates. The reader-facing derivation is
`cosmology/RADIAL_LIFT_SOLUTION.md`, and the executable contract is
`cosmology/SIMULATOR_RADIAL_CONTRACT_330.md`.

| Result | Formal paper location | Supporting material |
| --- | --- | --- |
| Exact screen projection | Radial theorem fragment, opening projection formula | Reference implementation and simulator contract |
| Correlation restriction | One-shell correlation restriction theorem | Reader-facing radial-lift note |
| Exact radial null space | Low-mode-removed kernel corollary | Null-space receipt and schema |
| Positive non-identifiability | One-shell exact no-go theorem | Positive counterexample tests |
| Radial tomography | Unrestricted tomographic uniqueness theorem | Tomography receipt contract |
| Source-refinement transfer | Source-refinement orbit producer theorem | Dilation-intertwiner receipt |
| Finite-to-continuum limit | Source-refinement orbit producer theorem and proof | Strong-convergence fields in the schema |
| Power-law forcing | Source-dilation uniqueness theorem | Dilation and planted-wiggle tests |
| Restricted-family uniqueness | Source-family injectivity corollary | Multipole amplitude consistency check |
| Exact amplitude conversion | Flat thin-shell lift theorem | Amplitude round-trip tests |
| Finite-window stability | Finite-window bound theorem | Quadrature and certified-bound tests |
| Multipole consistency | Source-family injectivity and multipole consistency corollary | Forward-residual contract |
| Finite singular-value and prior dependence | Finite singular-value and prior-continuation theorem | Prior formula, resolution, and null receipts |
| Branch-typed radial theorem | Branch-typed radial boundary theorem plus the primary forward operator | Curved-branch contract |
| Complete promotion rule | Promotion rule and primordial common gates | Fail-closed receipt schema |
| Temperature and polarization firewall | Transfer-firewall corollary | Boltzmann and likelihood papers |

The compact paper states the one-shell obstruction, both uniqueness routes,
flat thin-shell amplitude conversion, finite source boundary, and transfer
firewall. The main observer synthesis contains the complete formal fragment.
The screen-microphysics paper identifies the finite readouts needed by the
source producer.

The book uses the spherical-snapshot analogy and the tomography analogy. It
omits theorem labels, issue numbers, receipt identifiers, singular-value
formulas, and operator-limit details.
