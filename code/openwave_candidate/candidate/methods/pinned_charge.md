# Finite pinned-charge component on the native uniform-fiber instrument

`oph_openwave.pinned_charge.run()` executes the canonical 1x1 ring (87 port/transit positions) plus a fully enumerated modular triangle check. `run(include_larger=True)` adds the existing 2x2 benchmark (348 positions). No private workspace/Fable input is loaded.

The supplied model is H=H_m+κ diag(U), m=.23, κ=.29, modulus6, source port0. For each port z use ρ_z=δ_0−δ_z; for transit seam(a,b), use ρ_z=δ_0−(δ_a+δ_b)/2. Set G=(L^+)^(1/2), U_z=ρ_z^T Gρ_z/2. Full finite self terms are retained and U_0=0. This is a charge-only energy, not raw residue field energy. The modular triangle explicitly checks an integer Gauss-compatible lift, generator intertwining and every local C/E Kraus map. Position-only energy preserves the uniform-fiber code.

Every eigenbranch is retained, including negative-energy and transit flat branches. Select a positive eigenstate with maximum source-port overlap; this is supplied spectral preparation, not a positive-band projection of the dynamics. Store its complete amplitudes, full spectrum, position weights, residual and a phase-error bound t times the residual. The small eigen-residual is a floating diagnostic, not an interval eigenvalue theorem. H_m+U does not inherit the free squared-generator identity; a supplied port-only H0=L clock gives different probabilities.

Coherent eigenstate stationarity does not imply monitoring stability. Both the complete Bernoulli channel and a fixed source-seam native commit are executed with every outcome retained. Monitoring creates negative-branch weight and destroys the selected pure-state stationarity. The same total Bernoulli p across sizes is not fixed per-seam rate.

This is a PARTIAL finite-model observable. The literal row31 outcome is NOTTESTABLE: no certified infinite-screen bound/scattering threshold, hydrogen series, physical source/mass assignment, spatial refinement theorem or monitoring stability is established. A source-enhanced high-energy lattice eigenvector is not identified as an atom. No official criterion is promoted.
