# Specification of the Flavor-Orbit Selector S

Formal acceptance specification for any candidate construction of the source-derived
flavor-orbit selector for the quark mass sector, the fifth open object of the particle
lane. The restricted source-spread non-identifiability theorem
(`paper/deriving_the_particle_zoo_from_observer_consistency.tex`, paragraph "Restricted
source-spread non-identifiability theorem") proves that the light-quark source equations
fix two ordered Yukawa profile rays while leaving their two endpoint spans free; the
no-extra-axiom MAR theorem in the same paper proves that Axioms 1-5 plus fixed P admit a
free `(R_>0)^2` family of equal-score quark spectra. This file states what a candidate
selector must be and what it must pass. It constructs nothing and selects nothing.

**Status: specification plus a frozen preregistered menu. No candidate is promoted.
The selector remains open.**

## 1. Object

The selector is a map

```
S : (declared structure) -> (sigma_u, sigma_d) in (R_>0)^2
```

with zero continuous freedom, whose output fixes the two free spread moduli of the
ordered quark Yukawa rays.

### 1.1 Domain and codomain

- **Domain.** The declared structural inventory of Section 2b, and nothing else. The
  selector takes no runtime numeric input beyond the computed fixed point `P_fwd` and
  declared integers.
- **Codomain.** `(R_>0)^2`, the fiber of the non-identifiability theorem: the pair
  `(sigma_u, sigma_d)` of endpoint spans on the ordered rays `E_u = sigma_u * v_u`,
  `E_d = sigma_d * v_d`, where `v_u`, `v_d` are the unit-span zero-trace profiles with
  adjacent-gap ratios `rho_ord` and `1/rho_ord` (paper, restricted theorem; run artifact
  `../runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json`,
  `exact_ray_classification`).

### 1.2 Convention

`sigma_q` is the half endpoint span of the ordered log-Yukawa profile,
`sigma_q = (1/2) ln(y_q3 / y_q1)`, the convention of the audit surface
(`../runs/flavor/quark_current_family_exact_sigma_target.json`,
`../runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json`). The ray shape
datum `rho_ord` enters only the downstream display of implied per-generation gaps; it is
domain data of the rays, never a selector output and never a scale input.

## 2. Required properties

A candidate is admissible for evaluation only if all of S1-S5 hold.

- **S1. Both moduli fixed.** The candidate emits a single pair
  `(sigma_u, sigma_d) in (R_>0)^2` in closed form. A candidate that fixes one modulus, a
  ratio only, or a one-parameter family is not a selector candidate.
- **S2. Declared structure only (2b).** The closed form is a function of exactly:
  - the computed fixed point `P_fwd = 1.630972095858897`
    (`code/P_derivation/FULL_DERIVATION.md`, CL-6 converged forward value), and its
    certified source-audit companion `alpha_U(P_fwd) = 0.041124247441816685` where a
    candidate declares it;
  - A5 orbit geometry: 12 vertices (C5 ports), 20 faces (C3 fibers), 30 edges,
    60 face-corner flags (the regular A5 torsor);
  - the port counts 12 and 24 (write ports; write/check orientation doubling);
  - the Z6 data: trace coefficient 1/6, reserve density `P/24`, reserve fraction
    `e^(-P/24)`;
  - the pixel budget `P/4`;
  - W5 Weyl data: `dim W5 = 5`, quadrupole normalization `Q*Q = (8/5) P_5`;
  - the generation count 3.
  No other constant, count, or function may appear. Zero continuous freedom: the closed
  form contains no adjustable coefficient.
- **S3. Blindness.** The evaluation cone of a candidate contains no PDG mass, no
  measured ratio, no running-mass target, no fitted spread, and no artifact listed in
  the `forbidden_ancestors` block of
  `../runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json`. The
  comparison surface of Section 4 is read once, after every candidate output is
  recorded in `runtime/candidates_evaluated.json`. A candidate revised after that read
  is a retrodiction and is excluded permanently.
- **S4. Menu discipline.** The candidate list is written in `CANDIDATES.md`, with a
  structural rationale per candidate, before any candidate is evaluated. The list is
  frozen at its recorded count, at most 12. The 219,615-member quark denominator
  grammar (COMPRESSION_SCORECARD.md, discrete structural selections; self-rejected) is
  the controlling precedent: a selector that needs a large menu certifies nothing, and
  any landing from such a menu is worthless. Adding a candidate after evaluation
  reopens the count and demotes every prior standing to exploratory.
- **S5. Sector attachment declared once.** The assignment of the up sector to one A5
  structure and the down sector to another is declared once, in the menu header, before
  evaluation, and is shared by every candidate. Per-candidate reassignment is a hidden
  binary dial and is excluded.

## 3. Acceptance and promotion

- **Evaluation standing.** A candidate that satisfies S1-S5 is evaluated; its pair and
  implied dimensionless Yukawa ratios are recorded. Evaluation confers no standing
  beyond "recorded".
- **Selector-candidate standing.** A candidate whose recorded pair lands inside the
  documented two-modulus comparison window (Section 4) at the window's stated width is
  a selector-candidate. The window is PDG-derived, so this standing is compare-only and
  non-blind; it is a basin location, never a confirmation.
- **Promotion.** Promotion requires, in order: (1) a frozen registration of the single
  closed form with an external timestamp; (2) an out-of-sample test. The registered
  pair `(sigma_u, sigma_d)` together with the ray shapes fixes all six dimensionless
  light-quark Yukawa gap ratios, in particular `ln(y_s/y_d) = 2 sigma_d / (1+rho_ord)`
  and `ln(y_c/y_u) = 2 sigma_u rho_ord / (1+rho_ord)`. Future lattice determinations of
  `m_u/m_d`, `m_s/m_d`, and `m_c/m_s` at higher precision (FLAG-class averages) test
  those ratios without touching the registration. A registered candidate whose
  predicted ratios fall outside the tightened lattice bands is excluded permanently;
  relabeling does not reopen it.
- **Exhaustion.** If no candidate reaches selector-candidate standing, the menu is
  exhausted, the exclusions are recorded, and the selector remains open. That outcome
  is recorded in `COMPARISON.md` and is final for this menu.

## 4. Comparison surface

The documented two-modulus window is assembled from the audit surface, all PDG-derived:

| Point | sigma_u | sigma_d | Source |
|---|---|---|---|
| Exact current-family sigma target | 5.579692209267639 | 3.300314452061615 | `../runs/flavor/quark_current_family_exact_sigma_target.json` |
| Current theorem-grade pair | 5.5905 | 3.3049 | same artifact |
| Unselected formula point | 5.578418804072826 | 3.4210589139721543 | `../runs/flavor/quark_sigma_source_nonidentifiability_obstruction.json` |
| Zero-coefficient counterpoint | 5.407843949508826 | 3.346451912983049 | same artifact |

Window: `sigma_u in [5.407843949508826, 5.5905]`,
`sigma_d in [3.300314452061615, 3.4210589139721543]`. Every value in this table sits
downstream of PDG quark rows; the comparison is therefore compare-only and non-blind,
and it is labeled that way wherever it appears.

## 5. Tests

- **T1. Determinism.** Repeated evaluation emits byte-identical artifact content; no
  timestamps or run identifiers in `runtime/candidates_evaluated.json`.
- **T2. Schema.** The artifact carries the frozen count, one entry per candidate with
  closed forms, both moduli positive, and the implied ratio block.
- **T3. Blindness tokens.** The evaluator source contains no PDG, CODATA, or measured
  reference token; its numeric constants are `P_fwd`, `alpha_U(P_fwd)`, `rho_ord`
  (display only), and declared integers.

## 6. What this specification does not do

It does not derive any candidate from the axioms; it does not promote any candidate; it
does not close the flavor-selector row of the paper's open-object ledger. The menu in
`CANDIDATES.md` is an exclusion instrument: each evaluated candidate either becomes a
registered selector-candidate or a recorded exclusion.

**Status: specification plus frozen menu. The flavor-orbit selector is open.**
