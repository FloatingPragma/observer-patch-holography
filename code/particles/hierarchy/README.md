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
- `certificates/R_screen_sieve_icosahedral_certificate.json`
- `certificates/R_EW_global_capacity_certificate.json`
- `certificates/R_readback_resolution_certificate.json`
- `certificates/R_m_rep_24_certificate.json`
- `certificates/R_local_global_hierarchy_resonance_closeout_335.json`
- `certificates/R_pixel_screen_resonance_summary.json`
- `issue_332_rg_naturality_certificate.json`
- `certificates/R_PN_joint_fixed_point_certificate_report.json`
- `certificates/local_global_resonance_audit.json`
- `computations/derive_global_repair_tick_lemma.py`
- `verify_issue_332_rg_naturality.py`
- `verify_issue_335_local_global_resonance.py`
- `verify_issue_337_electroweak_projection.py`
- `verify_screen_sieve_theorem.py`
- `verify_issue_342_readback_resolution.py`
- `verify_issue_343_m_rep_24.py`
- `verify_issue_344_exact_capacity.py`
- `verify_pixel_screen_resonance_summary.py`
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
6. The `R_N` global repair-tick theorem: under the area-law counting model the
   readback fixed-point equation is equivalent to the closure transport
   `G_N(1) = rho_star`, giving the full-cycle multiplier
   `(N_CRC/pi)^(-1/2)`. The issue-#343 representation-to-spectrum theorem
   supplies `m_rep=2*(8+3+1)=24`, hence
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
9. The icosahedral screen-sieve theorem: on the declared triangulated `S^2`
   screen branch, `q_v=6-deg(v)` obeys `sum_v q_v=12`; convex defect cost
   selects twelve unit fivefold defects; edge-center collars expose them as
   central ports; `A5/C5` gives the 12-vertex orbit. This supplies the
   geometric origin of the `P/12` exponent.
10. The issue-#344 EW-refined exact-capacity certificate:
   `C_EW(P,x)=(1-lambda)*x+lambda*6*pi/(P*alpha_U(P))` is a contraction at
   `lambda=1/2`; its unique fixed point gives
   `N_CRC^EW=pi*exp[6*pi/(P_star*alpha_U(P_star))]` and
   `B_EW(P_star,N_CRC^EW)=0`.
11. The issue-#342 finite readback-resolution certificate:
   the fixed-cutoff pipeline `F_r(N)=Cap_read(Obs(nf_{r,N}(U_{r,N})))`
   has one selected positive central capacity atom, so
   `rho_read(r,N)=sqrt(pi/F_r(N))` is a singleton and
   `rho_read(r,N_CRC) -> (N_CRC/pi)^(-1/2)` in the positive-root
   refinement limit.
12. The issue-#343 representation-to-spectrum round-count theorem:
   the observer-visible product adjoint has dimensions `8+3+1=12`; reversible
   write/verify orientation doubles this to `m_rep=24`. The SU(5) adjoint has
   the same single-orientation integer for the wrong support because it
   includes X/Y mixed gauge channels excluded by the product branch.
13. The RG/Higgs naturality square for the selected exact branch:
   `epsilon_H=max(epsilon_n,epsilon_h)=0`, with measured weak-scale, Higgs,
   W/Z, gravity, Planck-area, and Lambda inputs excluded.
14. The issue-#335 close-out certificate: the prerequisite records are
   accounted for and the full local/global `N_CRC` hierarchy-resonance theorem
   closes on the selected branch.
15. The pixel-screen resonance summary receipt: the selected `(P_*,N_CRC^EW)`
   branch emits `K_cell=4*N_CRC^EW/P_*`, checks
   `K_cell*(P_*/4)=N_CRC^EW`, and records the dimensionless de Sitter
   coordinate pair `Lambda_CRC*l_star^2=3*pi/N_CRC^EW` and
   `Lambda_CRC*a_cell=3*pi*P_*/N_CRC^EW=12*pi/K_cell`. This is a receipt-level
   composition of existing certificates, not a new primitive-carrier or SI
   Lambda theorem.

External/source gates outside this bundle:

1. Source-only public Thomson endpoint transport `A_T(P)`.
2. Formal outward-rounded interval log from a certified interval stack.
3. Raw D10/D11 interval box for Higgs/top internals.
4. Full `R_gamma` stack for SI gravity/clock hierarchy.
5. An independently source-closed physical `E_star`; this bundle fixes
   `v/E_star`, not a weak scale in GeV.
6. The D10 QT1--QT5 quotient-path certificate: finite quotient
   canonicalization, explicit path enumeration, exact rational measures and
   central trace, fibre Gram/residual derivation, and a positive MAR gap.
7. A concrete frozen RG/matching/threshold/scheme receipt with truncation
   intervals.
8. W/Z/H complex-pole, residue, analytic-sheet, and uncertainty certificates.
9. A hash-bound prospective source DAG and branch-rigidity receipt.
10. Coupled-map extension of the joint `(P,N_CRC)` theorem, if OPH source work
   introduces cross-feedback between the local pixel and global capacity maps.

## Run validators

From this directory:

```bash
python3 validators/validate_bundle.py
python3 computations/hierarchy_recompute.py
python3 verify_issue_332_rg_naturality.py --check --output issue_332_rg_naturality_certificate.json
python3 verify_issue_335_local_global_resonance.py --check --output certificates/R_local_global_hierarchy_resonance_closeout_335.json
python3 verify_issue_337_electroweak_projection.py --check --output certificates/R_EW_tick_projection_certificate.json
python3 verify_screen_sieve_theorem.py --check --output certificates/R_screen_sieve_icosahedral_certificate.json
python3 verify_issue_342_readback_resolution.py --check --output certificates/R_readback_resolution_certificate.json
python3 verify_issue_343_m_rep_24.py --check --output certificates/R_m_rep_24_certificate.json
python3 verify_issue_344_exact_capacity.py --check --output certificates/R_EW_global_capacity_certificate.json
python3 verify_pixel_screen_resonance_summary.py --check --output certificates/R_pixel_screen_resonance_summary.json
python3 verify_joint_fixed_point_certificate.py --output certificates/R_PN_joint_fixed_point_certificate_report.json
```
