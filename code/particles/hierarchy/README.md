# OPH Hierarchy Proof Artifacts Bundle

Bundle ID: `oph-hierarchy-proof-artifacts-r1473-2026-06-13`

Generated: `2026-06-13T00:00:00Z`

This bundle packages the proof-obligation artifacts for the OPH electroweak hierarchy lane:

All capacity "closure" language in this bundle is internal to the selected
bridge packet. Physical cosmic
capacity remains conditional on F and CP-1 to CP-3, and comparison with the
Lambda-located value additionally requires the joint cosmological posterior.
The 6.6 percent central-value mismatch is not an unconditional contradiction.

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
- `certificates/R_U_outward_rounded_interval_log.json`
- `computations/hierarchy_numeric_witness.json`
- `computations/hierarchy_recompute.py`
- `computations/generate_ru_outward_rounded_log.py`
- `tools/outward_interval.py`
- `tools/ru_formula_stack.py`
- `validators/verify_ru_outward_rounded_log.py`
- `certificates/R_HT_declared_surface_certificate.json`
- `certificates/R_HT_interval_input_box_log.json`
- `tools/ht_formula_stack.py`
- `computations/generate_ht_interval_box_log.py`
- `validators/verify_ht_interval_box.py`
- `HT_INTERVAL_EXTENSION.md`
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
9. The conditional icosahedral screen-sieve theorem: `sum_v(6-deg(v))=12`;
   a strict unit-splitting cost selects twelve unit defects; inverse pairing
   and the D-optimal vector/quadrupole selector give the unique icosahedral
   six-axis frame; edge-center collars expose twelve central ports. This
   supplies the geometric origin of `P/12` on the declared branch.
10. The issue-#344 conditional EW-refined bridge-capacity certificate:
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
   conditional on active current algebra `su(3)+su(2)+u(1)`, the adjoint has
   dimension `8+3+1=12`; reversible orientation gives `m_rep=24`. Screen selection,
   the global quotient, and matter-menu completeness are not premises. The physical-current
   identification is upstream.
13. The declared-map RG/Higgs commuting square:
   `epsilon_H=max(epsilon_n,epsilon_h)=0` for the supplied comparison maps,
   with measured weak-scale, Higgs, W/Z, gravity, Planck-area, and Lambda
   inputs excluded. This is a packaging identity, not a source-derived
   naturality result. The antecedent-only receipt has same-fingerprint exact
   completions with defects zero and one; source-derived maps and independent
   residual evaluation remain open.
14. The issue-#335 certificate composes the exact conditional local/global
   resonance. The screen supplies `X/12`, the product adjoint independently
   supplies `m_rep=24`, and `HIERARCHY-SCREEN-READOUT` is required to identify
   `log(E_cell/v)` with `Gamma_screen`.
15. The issue-#331 outward-rounded interval log for the `R_U` witness:
   `certificates/R_U_outward_rounded_interval_log.json` re-evaluates the full
   `R_U` formula stack in directed-rounding IEEE-754 binary64 interval
   arithmetic (`tools/outward_interval.py`, one-ulp outward step after every
   correctly rounded primitive, in-module exp/log with explicit series
   remainder bounds, no libm exp/log assumption). The log records the
   interval image of every named formula-DAG node at both `I_U` endpoints,
   at the center, and over the full interval with the derivative enclosure,
   plus the Krawczyk inclusion `K(I_U)` strictly inside `int(I_U)` and the
   witness interval. `validators/verify_ru_outward_rounded_log.py`
   reproduces every serialized bound bit-exactly from the declared
   structural inputs (`P_fwd`, cutoffs, one-loop coefficients, seed) without
   importing measured endpoint constants, and checks that the
   high-precision root and Krawczyk image lie inside the outward-rounded
   witness.
16. The issue-#333 raw interval input box for the Higgs/top declared
   surface: `certificates/R_HT_interval_input_box_log.json` publishes the
   eleven declared branch inputs of the frozen R_HT formula stack (candidate
   D10 tuple plus declared D11 core/Jacobian constants) with centers, boxes,
   units or dimensionless normalization, and per-input provenance tags, and
   re-evaluates the full stack in the same directed-rounding backend as the
   `R_U` log (`tools/ht_formula_stack.py`, per-node interval extension
   documented in `HT_INTERVAL_EXTENSION.md`). The log records the interval
   image of every formula node at the centers and over the full box, the
   forward-mode Jacobian enclosure of `pi_y`, `pi_lambda`, `m_H`, and `m_t`
   in all eleven inputs, the diagonal chart-block non-singularity
   certificate with determinant bounded away from zero, exact rational
   inclusion of the recorded `(m_H, m_t)` decimals, and a uniqueness scope
   restricted to the declared surface (single-valued arithmetic and
   injective diagonal readout; no criticality-system existence or
   uniqueness claim). `validators/verify_ht_interval_box.py` reproduces
   every serialized bound bit-exactly, enforces the input allowlist and
   provenance classes, and scans the verify path for measured-constant
   markers. The inputs are declared branch inputs, so the enclosure closes
   the arithmetic grade of the declared surface while the source-root,
   scale, QT1--QT5, RG/scheme, rigidity, provenance, uncertainty, and
   complex-pole gates stay open upstream.
17. The pixel-screen resonance summary receipt: the selected `(P_*,N_CRC^EW)`
   branch emits `K_cell=4*N_CRC^EW/P_*`, checks
   `K_cell*(P_*/4)=N_CRC^EW`, and records the dimensionless de Sitter
   coordinate pair `Lambda_CRC*l_star^2=3*pi/N_CRC^EW` and
   `Lambda_CRC*a_cell=3*pi*P_*/N_CRC^EW=12*pi/K_cell`. This is a receipt-level
   composition of the listed certificates, not an independent primitive-carrier or SI
   Lambda theorem.

External/source gates outside this bundle:

1. Source-only public Thomson endpoint transport `A_T(P)`.
2. `HIERARCHY-SCREEN-READOUT`, the physical identification
   `log(E_cell/v)=Gamma_screen`, and source production of the unit cost,
   inverse pairing, and D-optimal screen selector.
3. Proof-assistant formalization of the `R_U` interval witness. The
   directed-rounding outward interval log itself is supplied inside this
   bundle at `certificates/R_U_outward_rounded_interval_log.json` with
   verifier `validators/verify_ru_outward_rounded_log.py`, under the stated
   IEEE-754 assumption set; the remaining external step is a mechanized
   proof-assistant replay of that log.
4. Source derivation of the eleven declared branch inputs of the Higgs/top
   declared surface. The raw interval input box, per-node interval
   extension, Jacobian enclosure, and non-singularity certificate are
   supplied inside this bundle at
   `certificates/R_HT_interval_input_box_log.json` with verifier
   `validators/verify_ht_interval_box.py`; the remaining external step is
   source emission of the D10 tuple and the D11 core/Jacobian constants.
5. Full `R_gamma` stack for SI gravity/clock hierarchy.
6. An independently source-closed physical `E_star`; this bundle fixes
   `v/E_star`, not a weak scale in GeV.
7. The D10 QT1--QT5 quotient-path certificate: finite quotient
   canonicalization, explicit path enumeration, exact rational measures and
   central trace, fibre Gram/residual derivation, and a unique target-free source rule.
8. A concrete frozen RG/matching/threshold/scheme receipt with truncation
   intervals.
9. W/Z/H complex-pole, residue, analytic-sheet, and uncertainty certificates.
10. A hash-bound source DAG and branch-rigidity receipt.
11. Coupled-map extension of the joint `(P,N_CRC)` theorem, if OPH source work
   introduces cross-feedback between the local pixel and global capacity maps.

## Run validators

From this directory:

```bash
python3 validators/validate_bundle.py
python3 computations/hierarchy_recompute.py
python3 computations/generate_ru_outward_rounded_log.py
python3 validators/verify_ru_outward_rounded_log.py
python3 computations/generate_ht_interval_box_log.py
python3 validators/verify_ht_interval_box.py
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
