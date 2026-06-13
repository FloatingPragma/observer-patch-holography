# OPH Hierarchy Proof Artifacts Bundle

Bundle ID: `oph-hierarchy-proof-artifacts-r1473-2026-06-13`

Generated: `2026-06-13T00:00:00Z`

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
- `certificates/R_EW_tick_projection_certificate.json`
- `certificates/R_local_global_hierarchy_resonance_closeout_335.json`
- `issue_332_rg_naturality_certificate.json`
- `certificates/R_PN_joint_fixed_point_certificate_report.json`
- `certificates/local_global_resonance_audit.json`
- `computations/derive_global_repair_tick_lemma.py`
- `verify_issue_332_rg_naturality.py`
- `verify_issue_335_local_global_resonance.py`
- `verify_issue_337_electroweak_projection.py`
- `verify_joint_fixed_point_certificate.py`
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
7. The product-branch joint `(P,N_CRC)` fixed-point theorem: on the
   product-separated source map `J(P,x)=(Gamma(P),C_hat(x))`, component
   contraction certificates imply a unique joint branch and stability with
   `q=max(qP,qN)<1`. A genuinely coupled source map requires the recorded
   weighted-sup derivative condition.
8. The electroweak tick-projection bridge: the unique projection map is
   `Pi_EW(P,N)=24*pi/(alpha_U(P)*log(N/pi))`; the target projection
   `Pi_EW(P_star,N_CRC)=4P_star` is equivalent to
   `B_EW(P,N)=alpha_U(P)*log(N/pi)-6*pi/P=0`. The certificate records the
   exact bridge target `N_EW(P_star)` and keeps the rounded `3.31e122` capacity
   label as a diagnostic.
9. The RG/Higgs naturality square for the selected exact branch:
   `epsilon_H=max(epsilon_n,epsilon_h)=0`, with measured weak-scale, Higgs,
   W/Z, gravity, Planck-area, and Lambda inputs excluded.
10. The issue-#335 close-out certificate: the prerequisite records are
   accounted for, but the full local/global `N_CRC` resonance is not promoted.
   The closed result is the exact surviving conditional statement: if an exact
   global capacity source certificate satisfies `B_EW(P_star,N_CRC)=0`, and the
   finite-readback and 24-round gates close, the target relation follows.

External/source gates outside this bundle:

1. Source-only public Thomson endpoint transport `A_T(P)`.
2. Formal outward-rounded interval log from a certified interval stack.
3. Raw D10/D11 interval box for Higgs/top internals.
4. Full `R_gamma` stack for SI gravity/clock hierarchy.
5. Theorem-grade `W/Z` promotion.
6. Finite-machinery verification that `nf_{r,N}` delivers a single well-defined
   effective readback resolution and representation-to-spectrum derivation of the
   24-round repair count for theorem-grade local/global resonance promotion.
7. Coupled-map extension of the joint `(P,N_CRC)` theorem, if OPH source work
   introduces cross-feedback between the local pixel and global capacity maps.
8. Exact global capacity source certificate satisfying `B_EW(P_star,N_CRC)=0`
   for the full local/global resonance theorem.

## Run validators

From this directory:

```bash
python3 validators/validate_bundle.py
python3 computations/hierarchy_recompute.py
python3 verify_issue_332_rg_naturality.py --check --output issue_332_rg_naturality_certificate.json
python3 verify_issue_335_local_global_resonance.py --check --output certificates/R_local_global_hierarchy_resonance_closeout_335.json
python3 verify_issue_337_electroweak_projection.py --check --output certificates/R_EW_tick_projection_certificate.json
python3 verify_joint_fixed_point_certificate.py --output certificates/R_PN_joint_fixed_point_certificate_report.json
```
