# Artifact Status Ledger

| Artifact | Status | Meaning |
|---|---|---|
| `R_P_public_pixel_certificate.json` | conditional public endpoint | Public endpoint branch recorded; source-only endpoint transport required for a no-measured-Thomson claim. |
| `R_P_source_audit_pixel_certificate.json` | source-audit witness | Avoids upstream measured Thomson endpoint; does not close public endpoint without same-scheme hadronic transport. |
| `Pi_U_frozen_source_packet.json` | supplied | Records frozen D10 packet. |
| `DAG_U.json` | supplied / passes declared graph check | Excludes forbidden measured-data paths on the declared graph. |
| `R_U_interval_certificate.json` | supplied witness | Endpoint signs and derivative interval enclosure are included. |
| `R_U_krawczyk_certificate.json` | supplied witness | Krawczyk image lies inside `I_U` under the derivative interval. |
| `hierarchy_numeric_witness.json` | supplied | Public and source-audit branch computations. |
| `R_HT_declared_surface_certificate.json` | partial | Formula/output ledger supplied; raw interval input box missing. |
| `RG_Higgs_naturality_defect_certificate.json` | closed exact selected branch | Defines `epsilon_H` and points to the issue-332 certificate, which supplies `epsilon_H=0` and `epsilon_H in [0,0]` on the selected exact source-to-Higgs branch. |
| `R_WZ_boundary_certificate.json` | compare-only | Prevents accidental promotion of `W/Z`. |
| `R_gamma_noG_DAG_certificate.json` | skeleton only | No-G rule and missing components are explicit. |
| `R_N_global_repair_tick_certificate.json` | closed lemma on declared rounds | Derives the closure transport from the corpus readback map under the declared area-law counting model (`F(N)=pi/rho_read^2`, so `N=F(N)` is equivalent to `G_N(1)=rho_star`), hence the full-cycle multiplier `(N_CRC/pi)^(-1/2)` and `|g_*'|=(N_CRC/pi)^(-1/48)` on the declared 24-round branch without EW inputs; the proof takes the counting model and round count as declared inputs. |
| `R_EW_tick_projection_certificate.json` | closed projection bridge | Defines `Pi_EW(P,N)=24*pi/(alpha_U(P)*log(N/pi))`; the target `Pi_EW=4P` is equivalent to `B_EW(P,N)=alpha_U(P)*log(N/pi)-6*pi/P=0`. The rounded `3.31e122` capacity label is diagnostic and fails the exact bridge residual. |
| `R_EW_global_capacity_certificate.json` | closed exact capacity fixed point | Defines `C_EW(P,x)=(1-lambda)*x+lambda*6*pi/(P*alpha_U(P))` with `lambda=1/2`; the unique fixed point gives `N_CRC^EW=pi*exp[6*pi/(P_star*alpha_U(P_star))]` and `B_EW(P_star,N_CRC^EW)=0`. |
| `R_PN_joint_fixed_point_certificate_report.json` | closed product-branch theorem with coupled boundary | Defines the joint space `I_P x log I_N` and the source map `J(P,x)=(Gamma(P),C_hat(x))`; component contractions imply a unique stable joint fixed point. A genuinely coupled source map requires derivative bounds with `max(a+b/r,d+r*c)<1`, otherwise the coupled branch remains residual freedom. |
| `issue_332_rg_naturality_certificate.json` | closed exact selected branch | Verifies the RG/Higgs naturality defect on the selected exact branch: `epsilon_n=epsilon_h=epsilon_H=0`, with measured weak-scale, Higgs, W/Z, gravity, Planck-area, and Lambda inputs excluded. |
| `R_local_global_hierarchy_resonance_closeout_335.json` | closed conditional close-out | Accounts for #336, #337, #344, #338, and #332; closes #335 as the exact bridge statement with finite-readback and 24-round derivation gates recorded. |
| `local_global_resonance_audit.json` | not used | Computes mismatch for rounded `N_CRC`, records the closed tick, projection bridge, exact-capacity fixed point, joint product-branch, and RG/Higgs naturality components, and lists the remaining finite-readback and round-count gates. |
