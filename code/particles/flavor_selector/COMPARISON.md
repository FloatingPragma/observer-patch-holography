# Flavor-Orbit Selector Menu: Comparison Record

Written after `runtime/candidates_evaluated.json` recorded all twelve candidate
outputs. This is the single comparison pass permitted by SELECTOR_SPEC.md clause S3.

**Label: compare-only and non-blind.** The two-modulus window below is assembled from
the documented audit surface, and every point in it sits downstream of PDG quark rows
(SELECTOR_SPEC.md Section 4). No blind window exists for these moduli; that absence is
exactly the non-identifiability gap this module addresses.

## Window

```
sigma_u in [5.407843949508826, 5.5905]
sigma_d in [3.300314452061615, 3.4210589139721543]
```

Reference point for offsets: the exact current-family sigma target
`(5.579692209267639, 3.300314452061615)`
(`../runs/flavor/quark_current_family_exact_sigma_target.json`).

## Standing of every candidate

Offsets are relative to the reference point. In-window flags are per modulus.

| # | Name | sigma_u | offset | in | sigma_d | offset | in | Standing |
|---|---|---|---|---|---|---|---|---|
| C-01 | Z6-traced face/vertex budget | 5.436573653 | -2.565% | yes | 3.261944192 | -1.163% | no | excluded |
| C-02 | Z6-traced edge/face budget | 8.154860479 | +46.153% | no | 5.436573653 | +64.729% | no | excluded |
| C-03 | Z6-traced write-check/port budget | 6.523888383 | +16.922% | no | 3.261944192 | -1.163% | no | excluded |
| C-04 | Pixel-budget orbit product | 8.154860479 | +46.153% | no | 4.892916288 | +48.256% | no | excluded |
| C-05 | Port-total pixel budget, F:V split | 6.116145359 | +9.614% | no | 3.669687216 | +11.192% | no | excluded |
| C-06 | Reserve-withheld traced budget | 5.079393418 | -8.966% | no | 3.047636051 | -7.656% | no | excluded |
| C-07 | Reserve-restored traced budget | 5.818870611 | +4.287% | no | 3.491322367 | +5.788% | no | excluded |
| C-08 | Orbit log-cardinality | 4.094344562 | -26.621% | no | 2.995732274 | -9.229% | no | excluded |
| C-09 | P-scaled log-cardinality | 6.677761732 | +19.680% | no | 4.885955745 | +48.045% | no | excluded |
| C-10 | Transmutation depth, F/V split | 5.318595616 | -4.679% | no | 3.191157369 | -3.307% | no | excluded |
| C-11 | Quadrupole-read traced budget | 6.876782164 | +23.247% | no | 4.126069298 | +25.020% | no | excluded |
| C-12 | Flag-total traced budget, F:V split | 10.193575599 | +82.691% | no | 6.116145359 | +85.320% | no | excluded |

## Outcome

No candidate places both moduli inside the documented window. The menu of twelve is
exhausted, twelve exclusions are recorded, and the flavor-orbit selector remains open.

The nearest exclusion is C-01, `(10P/3, 2P)`: its sigma_u lands inside the window and
its sigma_d falls 1.163% below the window floor. C-01 is recorded as excluded on the
frozen acceptance rule of SELECTOR_SPEC.md Section 3; the rule requires both moduli
inside the window, and no post-hoc tolerance is granted. A future menu that revisits
the Z6-traced budget family reopens the count under clause S4 and starts from
exploratory standing.

The window itself is PDG-derived, so even an in-window landing would have conferred
compare-only standing, pending the frozen registration and the lattice out-of-sample
test stated in SELECTOR_SPEC.md Section 3.

## Ledger effect

- The paper's open-object row "source-derived flavor-orbit selector" stays open.
- Twelve structural closed forms are permanently excluded for the documented window:
  the Z6-traced orbit budgets (C-01, C-02, C-03, C-12), the pixel-budget partitions
  (C-04, C-05), the reserve-shifted variants (C-06, C-07), the log-cardinality
  readings (C-08, C-09), the transmutation-depth split (C-10), and the
  quadrupole-read budget (C-11).
- Any future selector claim from these families cites this record first.
