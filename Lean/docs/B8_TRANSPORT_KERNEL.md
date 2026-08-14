# B8 finite transport kernel

The B8 package is an exact finite theorem layer. It does not identify a
repair step with physical time, an edge with physical distance, or any typed
quantity with a laboratory unit.

## Reversible Green--Kubo layer

`Thermodynamics/GreenKubo.lean` fixes the positive generator convention
`L = I - K`. For a finite stochastic kernel with nonnegative weights and
detailed balance, it proves the exact Dirichlet identity and nonnegativity of
`<f,Lf>`. A linear Poisson solver on centered currents then gives:

- `greenKuboPair_symm_of_poisson`, the Onsager symmetry;
- `greenKuboPair_finite_matrix_psd`, nonnegativity of every finite quadratic
  form of the coefficient matrix;
- `greenKuboPair_eq_integratedCorrelation_add_remainder`, the exact cutoff
  formula
  `<j,Rk> = sum_(n=0)^N <j,K^n k> + <j,K^(N+1)Rk>`.

Literal tail-zero statements are finite-extinction corollaries. They are not
presented as generic decay or as an infinite-time limit.

Under a strictly positive reference weight, the conditional-resampling heat
bath is idempotent. Its positive-lag correlation is therefore constant. A
fibre-centered current is killed in one
step and retains only its equal-time variance; a current with nonzero
one-step correlation has linearly growing partial sums. Thus one full-fibre
projector cannot serve as a nontrivial long-memory transport dynamics. A
physical OPH model with a nonzero decaying positive-lag tail needs a
source-derived nonidempotent evolution, such as a justified local or
random-scan composition, and an independently calibrated clock.

## Graph constitutive layer

`Thermodynamics/GraphDiffusion.lean` uses separate types for concentration,
temperature, particle amount, energy, particle and heat flux, edge distance,
clock increment, cell volume, heat capacity, and conductance. On a finite
closed oriented graph it proves:

- cancellation of total divergence and exact summation by parts;
- Fick and Fourier flux-gradient identities with nonpositive pairing under
  nonnegative conductance;
- amount--concentration and energy--temperature bridges;
- canonical one-step particle and energy updates with exact source balance
  and source-free global conservation;
- explicit two-vertex Fick and Fourier steps, plus negative-conductance
  controls.

The wrappers prevent accidental interchange of quantities; they do not
implement dimensional-unit arithmetic. The file proves no positivity
preservation, CFL condition, stability, convergence, entropy production,
boundary-reservoir law, hydrodynamic limit, or measured coefficient.

## Physical promotion boundary

A physical transport claim must additionally bind, on one source-derived
observer system, the transition generator, equilibrium reference, conserved
quantity, graph or continuum geometry, distance, clock, constitutive
coefficient, boundary conditions, and readout. Source-side common-reference and
collar work remains under #739, while energy/clock calibration is the empirical
import PR-15; #730 owns the remaining quantum attachment layers and the closed
#731 milestone supplies only a conditional mechanics interface. B8 supplies no frozen
prediction and changes no prediction-ladder row.
