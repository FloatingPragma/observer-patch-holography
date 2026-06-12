# OPH Hierarchy Proof Artifacts Bundle

Bundle ID: `oph-hierarchy-proof-artifacts-r1473-2026-06-12`

Generated: `2026-06-12T04:17:17Z`

This bundle packages the proof-obligation artifacts for the OPH electroweak hierarchy lane:

```math
P_\star \to \alpha_U(P_\star) \to
\frac{v}{E_\star} =
P_\star^{-1/2}\exp\left[-\frac{2\pi}{4\alpha_U(P_\star)}\right].
```

## What is included

- `certificates/R_P_public_pixel_certificate.json`
- `certificates/R_P_source_audit_pixel_certificate.json`
- `certificates/Pi_U_frozen_source_packet.json`
- `certificates/DAG_U.json`
- `certificates/R_U_interval_certificate.json`
- `certificates/R_U_krawczyk_certificate.json`
- `computations/hierarchy_numeric_witness.json`
- `computations/hierarchy_recompute.py`
- `certificates/R_HT_declared_surface_certificate.json`
- `certificates/RG_Higgs_naturality_defect_certificate.json`
- `certificates/R_WZ_boundary_certificate.json`
- `certificates/R_gamma_noG_DAG_certificate.json`
- `certificates/R_N_global_repair_tick_certificate.json`
- `certificates/local_global_resonance_audit.json`
- `computations/derive_global_repair_tick_lemma.py`
- `validators/validate_bundle.py`
- `validators/validate_manifest.py`

## Closure status

Closed inside this bundle:

1. The declared D10 source packet ledger.
2. The declared `DAG_U` forbidden-path check.
3. The `R_U` interval/Krawczyk inclusion witness.
4. The public and source-audit hierarchy computations.
5. The claim-tier boundary artifacts for `W/Z`, `R_HT`, `R_gamma`, and the local/global resonance.
6. The `R_N` global repair-tick lemma: under the declared area-law counting
   model the readback fixed-point equation is equivalent to the closure
   transport `G_N(1) = rho_star`, giving the full-cycle multiplier
   `(N_CRC/pi)^(-1/2)`; the declared 24-round branch yields
   `|g_*'| = (N_CRC/pi)^(-1/48)`.

External/source gates outside this bundle:

1. Source-only public Thomson endpoint transport `A_T(P)`.
2. Formal outward-rounded interval log from a certified interval stack.
3. Raw D10/D11 interval box for Higgs/top internals.
4. Concrete RG/Higgs naturality defect bound `epsilon_H`.
5. Full `R_gamma` stack for SI gravity/clock hierarchy.
6. Theorem-grade `W/Z` promotion.
7. Finite-machinery verification that `nf_{r,N}` delivers a single well-defined
   effective readback resolution, representation-to-spectrum derivation of the
   24-round repair count, plus electroweak tick projection, joint `(P,N)`
   stability, and RG/coarse-graining naturality for theorem-grade local/global
   resonance promotion.

## Run validators

From this directory:

```bash
python3 validators/validate_bundle.py
python3 computations/hierarchy_recompute.py
```
