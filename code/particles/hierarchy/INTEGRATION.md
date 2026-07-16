# Hierarchy Bundle Integration

This directory carries the frozen electroweak-hierarchy proof bundle supplied
for the OPH hierarchy audit. It is separate from the rounded `1.63094`
calibration carrier used by the older D10 particle-code path.

Status correction (2026-07-16): the bundle certifies internal algebra on a
selected bridge branch. Its `N_CRC^EW` value is conditional on the physical
readback map F and premises CP-1 to CP-3. It neither closes the physical cosmic
capacity nor makes the 6.6 percent central-value comparison with the
Lambda-located capacity a contradiction or significance test.

The bundle’s public endpoint branch uses

```text
P_C = 1.630968209403959324879279847782648941
alpha_U(P_C) = 0.041124336195630495
v/E_star = 2.0199803239725553e-17
```

The source-audit branch keeps the public Thomson endpoint out of the upstream
path and records

```text
P_cand = 1.63097209569432901817967892561191884270169 (predates the 2026-07-14 converged rerun; the certified P_fwd is 1.630972095858897..., ledger row CL-6)
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
repair-tick theorem on the selected resonance branch: with the corpus readback
map `F = Cap_read(Obs(nf))` and the area-law counting model
`F(N) = pi/rho_read^2` (a D6-consistent modeling identification of
`Cap_read`), the fixed-point equation `N_CRC = F(N_CRC)` is equivalent to the
closure transport `G_N(1) = rho_star`, which forces the full-cycle multiplier
`(N_CRC/pi)^(-1/2)`. The issue-#343 representation-to-spectrum theorem supplies
`m_rep=2*(8+3+1)=24`, so the one-tick multiplier is
`|g_*'| = (N_CRC/pi)^(-1/48)` with no electroweak inputs.

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
rounded `3.31e122` capacity label as a diagnostic.

The issue-#344 artifact supplies the conditional EW-refined bridge-capacity
identity. In log-capacity coordinates,

```text
C_EW(P,x) = (1-lambda)*x + lambda*6*pi/(P*alpha_U(P))
```

with `lambda=1/2` is a contraction. Its unique fixed point gives

```text
N_CRC^EW = pi*exp[6*pi/(P_star*alpha_U(P_star))]
B_EW(P_star,N_CRC^EW) = 0
```

using only `P_star`, source `alpha_U(P_star)`, `pi`, and `exp`.

The issue-#342 artifact supplies the finite readback-resolution certificate.
At fixed cutoff,

```text
F_r(N) = Cap_read(Obs(nf_{r,N}(U_{r,N})))
rho_read(r,N) = sqrt(pi/F_r(N))
```

The selected finite branch has one positive central capacity atom and zero
variance, so the delivery resolution is a singleton. In the positive-root
refinement limit,

```text
rho_read(r,N_CRC) -> (N_CRC/pi)^(-1/2)
```

with weak/Higgs/gravity/cosmology calibrations excluded.

The issue-#343 artifact supplies the representation-to-spectrum round count.
The observer-visible product adjoint has dimensions

```text
dim(su(3) + su(2) + u(1)) = 8 + 3 + 1 = 12
```

and reversible repair uses write/verify orientation, giving

```text
m_rep = 2 * 12 = 24.
```

The SU(5) adjoint has the same single-orientation integer for a different
support; its X/Y mixed gauge channels are excluded by the OPH product branch.

The issue-#335 close-out certificate accounts for the closed global tick,
projection bridge, exact capacity fixed point, finite readback resolution,
representation-to-spectrum round count, joint product branch, and RG/Higgs
naturality records. It closes the umbrella issue as the full local/global
`N_CRC` hierarchy-resonance theorem on the selected branch.

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
python3 verify_issue_342_readback_resolution.py --check --output certificates/R_readback_resolution_certificate.json
python3 verify_issue_343_m_rep_24.py --check --output certificates/R_m_rep_24_certificate.json
python3 verify_issue_344_exact_capacity.py --check --output certificates/R_EW_global_capacity_certificate.json
python3 verify_joint_fixed_point_certificate.py --output certificates/R_PN_joint_fixed_point_certificate_report.json
```
