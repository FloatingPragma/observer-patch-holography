# Hierarchy Bundle Integration

This directory carries the frozen electroweak-hierarchy proof bundle supplied
for the OPH hierarchy audit. It is separate from the rounded `1.63094`
calibration carrier used by the older D10 particle-code path.

The bundle’s public endpoint branch uses

```text
P_C = 1.630968209403959324879279847782648941
alpha_U(P_C) = 0.041124336195630495
v/E_star = 2.0199803239725553e-17
```

The source-audit branch keeps the public Thomson endpoint out of the upstream
path and records

```text
P_cand = 1.63097209569432901817967892561191884270169
alpha_U(P_cand) = 0.04112424744557487
v(P_cand)/E_star = 2.0198114150099223e-17
```

The claim closed here is local electroweak hierarchy closure on the declared
source graph:

```text
P_star -> alpha_U(P_star) -> v/E_star
```

The `R_U` Krawczyk certificate proves a unique source zero inside the supplied
interval for the declared formula stack. The bundle also records its boundary:
source-only public Thomson transport, raw outward-rounded interval logs,
Higgs/top interval boxes, theorem-grade W/Z promotion, the full no-G clock
stack, and the theorem-grade local/global resonance theorem are outside this
certificate.

The `R_N_global_repair_tick_certificate.json` artifact closes the global
repair-tick lemma on the declared 24-round resonance branch: with the corpus
readback map `F = Cap_read(Obs(nf))` and the declared area-law counting model
`F(N) = pi/rho_read^2` (a D6-consistent modeling identification of
`Cap_read`), the fixed-point equation `N_CRC = F(N_CRC)` is equivalent to the
closure transport `G_N(1) = rho_star`, which forces the full-cycle multiplier
`(N_CRC/pi)^(-1/2)`; the declared 24-tick decomposition then yields
`|g_*'| = (N_CRC/pi)^(-1/48)` with no electroweak inputs. The counting model
and the 24-round count are declared branch inputs, and the corpus marks the
finite readback map as schematic. The finite-machinery
verification of the readback resolution and the representation-to-spectrum
round-count derivation remain open before the full resonance relation is
promoted.

The issue-#337 artifact closes the electroweak tick-projection bridge. The
projection map is

```text
Pi_EW(P,N) = 24*pi/(alpha_U(P)*log(N/pi))
```

and the target projection `Pi_EW(P_star,N_CRC)=4P_star` is equivalent to

```text
B_EW(P,N) = alpha_U(P)*log(N/pi) - 6*pi/P = 0
```

The certificate records the exact bridge target `N_EW(P_star)` and keeps the
rounded `3.31e122` capacity label as a diagnostic. The full local/global
resonance theorem requires the exact global capacity source certificate that
satisfies `B_EW(P_star,N_CRC)=0`.

The issue-#335 close-out certificate accounts for the closed global tick,
projection bridge, joint product branch, and RG/Higgs naturality records. It
does not promote the full local/global `N_CRC` resonance. It closes the umbrella
issue as the exact surviving conditional statement: if the exact global
capacity source certificate supplies `B_EW(P_star,N_CRC)=0`, and the finite
readback-resolution and 24-round derivation gates close, then the target
local/global hierarchy relation follows.

The joint `(P,N_CRC)` artifact defines the product branch map
`J(P,x)=(Gamma(P),C_hat(x))` on `I_P x log I_N`. Component contraction
certificates imply a unique stable joint fixed point. If a later source branch
contains genuine cross-feedback, the package records the necessary
weighted-sup derivative condition `max(a+b/r,d+r*c)<1`. The coupled branch
has residual freedom unless that condition is supplied.

The issue-#332 artifact closes the RG/Higgs naturality square on the selected
exact branch. The verifier emits

```text
epsilon_n = epsilon_h = epsilon_H = 0
epsilon_H in [0, 0]
```

and forbids measured weak-scale, Higgs, W/Z, gravity, Planck-area, Lambda, and
tuned bare-Higgs/cutoff-counterterm inputs. The optional `N_CRC` and repair-tick
checks are diagnostic unless supplied by the upstream resonance records.

Run the local guard with:

```bash
python3 validators/validate_bundle.py
python3 computations/hierarchy_recompute.py
python3 -m pytest test_hierarchy_bundle.py
python3 verify_issue_332_rg_naturality.py --check --output issue_332_rg_naturality_certificate.json
python3 verify_issue_335_local_global_resonance.py --check --output certificates/R_local_global_hierarchy_resonance_closeout_335.json
python3 verify_issue_337_electroweak_projection.py --check --output certificates/R_EW_tick_projection_certificate.json
python3 verify_joint_fixed_point_certificate.py --output certificates/R_PN_joint_fixed_point_certificate_report.json
```
