# Capacity Readback (Generator G2)

Specification and executable schema for the OPH capacity readback map
`F = Cap_read ∘ Obs ∘ nf`, the object whose construction and certified execution moves
closure rows CL-3, CL-4, and CL-7 (generator G2 in
[CONSISTENCY_STACK.md](../../docs/CONSISTENCY_STACK.md); rows in
[CLOSURE_LEDGER.md](../../docs/CLOSURE_LEDGER.md)).

Standing of the lane in one paragraph: the equation that would force the
measured capacity is not known. The conditional bridge
`N = pi*exp(6*pi/(P*alpha_U))` is the only closed form on the table; it lands
6.6 percent above the Λ-located value, its balance premise CP-1 is a declared
input, and the corrected-balance candidates
([CP1_CORRECTED_BALANCE_CANDIDATES_2026-07-17.md](CP1_CORRECTED_BALANCE_CANDIDATES_2026-07-17.md))
are underived. The working capacity `N_Λ = 3π/(Λ·l_P²) ≈ 3.313e122` is a basin
location under SL-4, an input everywhere it appears. The propagated Planck
posterior for the comparison lives in
[planck_posterior/](planck_posterior/PROPAGATION_RECORD.md).

## What this directory is

- [`F_READBACK_SPEC.md`](F_READBACK_SPEC.md): the formal acceptance specification any
  candidate F must satisfy before CL-7 can close: domain/codomain, the three factors
  with their paper sources, required properties P1–P5 (well-definedness, monotonicity,
  growth bounds for a contraction interval, count-density coherence with
  `log|Ω^sc_N| − N`, non-triviality), the V-08 blindness requirement, and the
  acceptance-test checklist A1–A7 with the certificate schema.
- [`toy_readback.py`](toy_readback.py): an executable toy model (strings over a
  three-letter alphabet; sorting as nf; a marker-cell filter as Obs; a declared log
  normalization as Cap_read) running the full schema end to end: solve `N = F_toy(N)`,
  certify the fixed point with an interval Banach check (`F_toy(I) ⊆ I` and
  `|F_toy'| ≤ L < 1`, evaluated with `mpmath.iv`), and emit
  `runtime/toy_readback_certificate.json`.
- [`test_toy_readback.py`](test_toy_readback.py): tests where enumeration matches the
  closed-form sector count, the toy fixed point exists with a certified enclosure, the
  certificate schema keys are present, output is deterministic, and a deliberately
  non-contracting variant is rejected.
- [`F_CONSTRUCTION_2026-07-14.md`](F_CONSTRUCTION_2026-07-14.md): the recorded
  construction run: candidate counts `|Omega^sc_N|` derived from the declared structure only
  (equal-area cells `4N/P`, Z6 reserve `P/24` in poisson and presence readings, the
  twelve-port load `log(N/pi)/12`, oriented 24-slot register, A5/C3 face-corner
  factors), carried as a full branch lattice with menu sizes and certified per branch.
- [`F_candidate_capL.py`](F_candidate_capL.py),
  [`F_candidate_capP.py`](F_candidate_capP.py),
  [`F_candidate_capK.py`](F_candidate_capK.py),
  [`compare_references.py`](compare_references.py): the three executed candidate
  families (log-count, port-inversion, cell-count readback) and the stage-3
  comparison; certificates in `runtime/F_candidate_*_certificates.json` and
  `runtime/F_construction_comparison_2026-07-14.json`.
- [`CP2_INVERSION_FORM_REDUCTION_2026-07-17.md`](CP2_INVERSION_FORM_REDUCTION_2026-07-17.md):
  the GAP-A3 reduction record: CP-2 decomposed into load mediation plus read
  consistency, the surviving freedom classified exactly (the load gauge), and the
  conditional theorem that P4 coherence plus the count-side premise CP-4 forces the
  fixed-point location across the whole gauge freedom without pinning the gauge.
- [`P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md`](P4_COUPLED_COHERENCE_REDUCTION_2026-07-17.md):
  the GAP-A8 reduction record: the clause-4 reading dichotomy (aggregate reading closes
  P4 degenerately by support collapse; per-configuration reading reduces it to CP-4 by
  the coherence-transfer theorem), with positive and negative controls.
- [`cp2_p4_premise_reduction.py`](cp2_p4_premise_reduction.py),
  [`test_cp2_p4_premise_reduction.py`](test_cp2_p4_premise_reduction.py): the interval
  machine checks for both reduction records (gauge-family axiom verdicts, support
  collapse, CP-4 controls, recorded-lattice separation); artifact
  `runtime/cp2_p4_premise_reduction_check.json`, deterministic across reruns.

## Construction run: outcome class (b), CL-7 stays open

190 branch rows executed blind (P only as the certified forward enclosure; no
measured Lambda, no SL-4 estimate, no CL-3 bridge value). 102 rows carry certified
interval Banach fixed points, all in `[1.4686, 1452.33]` nats; zero rows are
P4-coherent; no row lands within an order of magnitude of either reference capacity
(shortfall at least `10^119.36`). The constructed candidates are excluded and the
exclusion is structural: every count assembled from the declared 12/24-port,
Z6-reserve, and P/4-budget combinatorics has `1 - rho` set by the reserve fraction
(order `0.17..0.5`), so both displays of the target close at O(1)-O(10^3) nats. The
named missing theorem is gap item G2-GAP-1 in
[`F_CONSTRUCTION_2026-07-14.md`](F_CONSTRUCTION_2026-07-14.md): a source-side
derivation tying the port load `log(N/pi)` to the D10 inner observation step, the
same object as the CL-3 resolution path. The ledger is untouched.

## What this directory is not

No physical readback map is constructed here. The toy model carries no physical
content: its state space, repair rule, observer sector, and readback normalization are
declared toy choices with no relation to the OPH grammar. The toy certificate does not
move CL-7, and every artifact in `runtime/` says so in its own keys
(`physical_content: false`, `moves_cl7: false`, `cl7_status: "open"`).

## Promotion path

1. Construct the physical `U_N`, `nf`, `Obs`, `Cap_read` on the OPH normal-form grammar
   satisfying spec properties P1–P5 (finite cutoff first, then the refinement limit).
2. Pass the V-08 dependency-cone audit: the evaluation cone contains neither measured Λ
   nor the SL-4 estimate `N = 3.31e122` nor the CL-3 electroweak-bridge value.
3. Emit the stage-2 interval contraction certificate (acceptance tests A3–A5) with a
   directed-rounding backend for theorem-grade promotion, mirroring
   [`code/P_derivation/`](../P_derivation/).
4. Compare the blind output once against the SL-4 basin (A7). Inside: CL-7 closes and
   CL-3/CL-4 become evaluable at the one certified N. Outside: the verdict is recorded
   permanently (STRANGE_LOOP_PRINCIPLES.md rule 7).

## Usage

```bash
python3 toy_readback.py                          # emit runtime/toy_readback_certificate.json
python3 toy_readback.py --variant non_contracting --output /tmp/neg.json
python3 -m pytest test_toy_readback.py -q
```
