# Preregistered Flavor-Orbit Selector Candidates

Frozen menu under SELECTOR_SPEC.md clause S4. Every closed form on this page was
written before any numeric value was computed. The count is frozen at **12**.
Evaluation order is the listed order. No candidate is added, removed, or edited after
the first evaluation run.

## Declared structure (SPEC clause S2)

```
P      = P_fwd = 1.630972095858897        computed fixed point (CL-6 forward value)
aU     = alpha_U(P_fwd) = 0.041124247441816685   certified source-audit companion of P_fwd
V      = 12    icosahedral vertices (C5 curvature ports)
F      = 20    icosahedral faces (C3 fibers of the carrier)
E      = 30    icosahedral edges
G      = 60    face-corner flags (regular A5 torsor)
N_wc   = 24    write/check ports (x2 orientation doubling of the 12 ports)
z6     = 1/6   Z6 rank-one trace coefficient
eps    = P/24  Z6 shared scalar-reserve density; reserve fraction e^(-eps)
b      = P/4   pixel (screen-cell) entropy budget
w5     = 8/5   W5 quadrupole normalization, Q*Q = (8/5) P_5
N_g    = 3     generations
```

## Sector attachment (SPEC clause S5, declared once)

The up sector attaches to the face structure (the C3 fibers of the charged-lepton
face-corner carrier theorem; the up ray carries the adjacent-gap ratio `rho_ord > 1`).
The down sector attaches to the vertex structure (the twelve C5 ports). Where a
candidate uses a total budget with one split, the split is `F : V = 20 : 12 = 5 : 3`,
up on the larger share. This assignment is fixed for all twelve candidates.

## The menu

| # | Name | sigma_u | sigma_d | Structural rationale |
|---|---|---|---|---|
| C-01 | Z6-traced face/vertex budget | `P*F/6 = 10P/3` | `P*V/6 = 2P` | Each span is the Z6 rank-one trace (1/6) of the P-weighted cardinality of its carrier orbit: faces for up, vertices for down. |
| C-02 | Z6-traced edge/face budget | `P*E/6 = 5P` | `P*F/6 = 10P/3` | The same trace functional on the other adjacent pair of the icosahedral f-vector (30 edges, 20 faces). |
| C-03 | Z6-traced write-check/port budget | `P*N_wc/6 = 4P` | `P*12/6 = 2P` | The x2 write/check orientation doubling (12 to 24) as the up/down asymmetry; the Z6 trace supplies 1/6. |
| C-04 | Pixel-budget orbit product | `(P/4)*F = 5P` | `(P/4)*V = 3P` | Each span is the pixel entropy budget P/4 accumulated once per carrier-orbit element. |
| C-05 | Port-total pixel budget, F:V split | `(P/4)*24*(5/8) = 15P/4` | `(P/4)*24*(3/8) = 9P/4` | The 24 write/check ports each carry one pixel budget; the total 6P splits by the declared F:V bipartition. |
| C-06 | Reserve-withheld traced budget | `(10P/3)*e^(-P/24)` | `2P*e^(-P/24)` | C-01 with the Z6 shared scalar reserve e^(-P/24) withheld from the readable span in each sector. |
| C-07 | Reserve-restored traced budget | `(10P/3)*e^(+P/24)` | `2P*e^(+P/24)` | C-01 with the reserve returned to the span: the traced budget is the post-reserve remainder and the physical span is the pre-reserve total. |
| C-08 | Orbit log-cardinality | `ln(G) = ln 60` | `ln(F) = ln 20` | Gibbs reading: the span is the log-cardinality of the carrier orbit, flags for up (the regular torsor), faces for down. |
| C-09 | P-scaled log-cardinality | `P*ln 60` | `P*ln 20` | C-08 with the fixed point as the per-nat conversion of orbit entropy into log-Yukawa span. |
| C-10 | Transmutation depth, F/V split | `(5/3)*ln(1/aU)` | `ln(1/aU)` | The down span is the log-depth of the unified coupling at the computed fixed point; the up span carries the F/V ratio of the sector attachment. |
| C-11 | Quadrupole-read traced budget | `sqrt(8/5)*(10P/3)` | `sqrt(8/5)*2P` | C-01 read through the W5 quadrupole map, which scales family-plane norms by sqrt(8/5) (Q*Q = (8/5) P_5). |
| C-12 | Flag-total traced budget, F:V split | `P*G/6*(5/8) = 25P/4` | `P*G/6*(3/8) = 15P/4` | The Z6 trace of the full 60-flag torsor budget 10P, split by the declared F:V bipartition. |

**Frozen count: 12. Menu closed before evaluation.**

Each candidate emits both moduli in closed form from the declared structure with zero
continuous freedom (SPEC S1, S2). Implied dimensionless Yukawa ratios are computed by
`evaluate_candidates.py` only after all twelve pairs are recorded, using the ray shape
datum `rho_ord = 1.2942849363777058` for the per-generation gap display.
