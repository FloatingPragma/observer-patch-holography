**Paper release:** `r1549`
**Released:** July 18, 2026

# Comparison Classes

Every reported number belongs to one of four classes:

1.  **source-derived**: all physical inputs follow from the declared OPH construction without using the comparison data;

2.  **imported**: established background or microphysical inputs are declared explicitly;

3.  **fit**: one or more parameters are inferred from the same data;

4.  **diagnostic**: the calculation tests shape, scale, code behavior, or an interface without supporting a physical OPH claim.

The cosmology comparisons in scope are diagnostic or imported. No source-derived joint cosmology likelihood exists.

# Required Model Information

A physical comparison must state:

1.  the covariant action and source map;

2.  the background solution and species content;

3.  the gauge-invariant initial conditions;

4.  the perturbation, collision, and exchange equations;

5.  the observable projection and numerical tolerances;

6.  every external input and fitted parameter.

For the repair-charge dark-sector proposal, the dilute dust-like equation of state is a conditional action consequence. The abundance, relativistic stress, initial perturbations, lensing map, clusters, and nonlinear structure equations are work in progress.

# Required Data Information

A comparison must identify the data release, selection, masks, covariance, calibration model, nuisance parameters, priors, and combination rule. Dataset overlap and shared calibration must be included. A diagonal sum is not a joint likelihood when the covariance is non-diagonal.

The model and data layers must remain separate. Measured values that select a source formula, normalization, branch, or applicability domain are inputs to that comparison rather than OPH outputs.

# Observable Contracts

<div class="center">

| Observable family | Required theory output | Required comparison object |
|:---|:---|:---|
| CMB | $`C_\ell^{TT}`$, $`C_\ell^{TE}`$, $`C_\ell^{EE}`$, lensing, foreground model | map or band-power likelihood with covariance and nuisance model |
| BAO and supernovae | $`H(z)`$, $`D_A(z)`$, luminosity distance, sound horizon | survey likelihood with calibration covariance |
| Weak lensing | metric potentials, nonlinear matter power, intrinsic alignments | tomographic likelihood with masks and covariance |
| Growth and RSD | $`P_m(k,z)`$ and $`f\sigma_8(z)`$ | windowed survey likelihood with cross-bin covariance |
| Clusters | mass function, selection, observable–mass map, lensing calibration | count and calibration likelihood |
| Galaxies | disk dynamics, gas and stellar nuisances, external environment | galaxy-level likelihood with selection and covariance |

</div>

# Diagnostic Results

The repository contains CMB overlays, compressed background comparisons, and galaxy scaling checks. Their formulas or inputs are selected with knowledge of the comparison data, omit required covariance or nuisance structure, or lack a physical source map. They test numerical paths and conditional formulas. They do not establish a cosmological prediction.

# No-Data-Use Boundary

A source-derived calculation must expose its dependency graph. The source side may use mathematical constants, declared OPH finite artifacts, and explicitly imported physical inputs. It may not use the observed quantity being claimed as an output, a fit to that quantity, or a proxy calibrated from the same dataset.

This provenance rule is enforced through the declared dependency graph and data-use classification.

# Status

The screen-spectrum and radial-lift mathematics have a conditional theorem packet. One-shell inversion is nonunique; physical source dilation and radial cross-covariance tomography are the two uniqueness routes. Finite source instantiation, repair-charge cosmological source, Boltzmann bridge, and joint likelihood are work in progress. All cosmological data products remain diagnostic. Cosmology is outside the OPH falsification program.
