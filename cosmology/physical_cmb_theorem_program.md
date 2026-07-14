# Physical CMB Theorem Program

This contract keeps the physical-CMB promotion boundary machine-auditable. OPH-FPE CMB curves remain diagnostics until a bounded finite observer-patch parent emits local state, declared ports and boundaries, readback records, feedback or repair moves, and a public evidence bundle that closes every contract below.

## Required contracts

- **Finite source:** source-only receipts for `A_zeta`, scalar shape and tilt, `q_IR`, `ell_IR`, freezeout, and `N_CRC`, with no measurement-tuned constants.
- **Physical scale:** receipts for physical `k`, angular `ell`, scale factor or redshift, freezeout, and the common primordial/anomaly mode basis.
- **Covariant stress parent:** finite stress-energy closure, recipient stress for nonzero exchange, gauge-independent variables, causal response, and refinement convergence.
- **Dark/anomaly kernels:** physical `B_A(k,a)`, `rho_A(a)`, and `Gamma_rec(k,a)` from one declared parent; `Gamma_rec` also requires a physical clock, active-fiber receipt, conserved-sector decomposition, and common-parent response pole.
- **Transfer and likelihood:** frozen source, solver, and likelihood hashes; custom-parent and Standard-Model-off regressions; CDM-limit recovery; blinded comparison; full observable likelihood execution; and a global rather than shard-local reduction.

The executable schema and fail-closed validator live in the sibling simulator repository at `oph-physics-sim/oph_fpe/cosmology/physical_cmb_contract.py`; its regression tests live at `oph-physics-sim/tests/test_physical_cmb_contract.py`. This document is the research-repository theorem contract, not a duplicate implementation.

## Owning theorem fragments

- `paper/tex_fragments/SCREEN_SPECTRUM_THEOREMS.tex` — screen scalar, precision operator, MaxEnt covariance, tilt receipts. Status: screen-level receipts open; the `P/48` tilt is an analytic candidate with the `kappa_rep` certificate pending.
- `paper/tex_fragments/PRIMORDIAL_BRIDGE_THEOREMS.tex` — screen-to-primordial radial lift, amplitude conversion, null-space and forward-projection receipts. Status: open.
- `paper/tex_fragments/PHYSICAL_SCALE_BRIDGE_THEOREMS.tex` — physical `k`, angular `ell`, scale factor, freezeout certificates. Status: open.
- `paper/tex_fragments/FINITE_COVARIANT_PARENT_THEOREMS.tex` — covariant stress parent, stress/exchange closure, kernel receipts. Status: open.
- `paper/tex_fragments/FINITE_QUOTIENT_ENSEMBLE_THEOREMS.tex` — quotient ensembles and claim tiers E0–E5. Status: shared surface, in use.
- `cosmology/oph_boltzmann_transport_derivation.tex` — transport reduction with declared imports; transfer receipts `Boltz-R1/N1/X1/V1`. Status: R1/V1 open, N1/X1 declared imports.
- `cosmology/oph_cosmology_data_likelihood_contracts.tex` — frozen source/solver/likelihood hashes and no-data-use receipts. Status: zero frozen physical-prediction receipts.

## Issue ownership

- [#371](https://github.com/FloatingPragma/observer-patch-holography/issues/371): finite-source contract hardening.
- [#372](https://github.com/FloatingPragma/observer-patch-holography/issues/372): physical scale bridge and mode calibration.
- [#373](https://github.com/FloatingPragma/observer-patch-holography/issues/373): finite covariant collar-packet stress parent.
- [#374](https://github.com/FloatingPragma/observer-patch-holography/issues/374): physical dark/anomaly kernels.
- [#363](https://github.com/FloatingPragma/observer-patch-holography/issues/363): Boltzmann transfer and frozen likelihood closure.

GitHub issue state is canonical; regenerate the open-problem ledger before using this list as a schedule.

## Claim boundary

Until the complete contract passes, OPH-FPE CMB outputs are measurement-facing diagnostics, not physical CMB predictions. Good shape agreement identifies useful work; it does not replace finite source, scale, stress, transfer, and frozen-likelihood receipts.
