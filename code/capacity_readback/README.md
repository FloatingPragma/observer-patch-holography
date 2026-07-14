# Capacity Readback (Generator G2)

Specification and executable schema for the OPH capacity readback map
`F = Cap_read ∘ Obs ∘ nf`, the object whose construction and certified execution moves
closure rows CL-3, CL-4, and CL-7 (generator G2 in
[CONSISTENCY_STACK.md](../../CONSISTENCY_STACK.md); rows in
[CLOSURE_LEDGER.md](../../CLOSURE_LEDGER.md)).

## What this directory is

- [`F_READBACK_SPEC.md`](F_READBACK_SPEC.md) — the formal acceptance specification any
  candidate F must satisfy before CL-7 can close: domain/codomain, the three factors
  with their paper sources, required properties P1–P5 (well-definedness, monotonicity,
  growth bounds for a contraction interval, count-density coherence with
  `log|Ω^sc_N| − N`, non-triviality), the V-08 blindness requirement, and the
  acceptance-test checklist A1–A7 with the certificate schema.
- [`toy_readback.py`](toy_readback.py) — an executable toy model (strings over a
  three-letter alphabet; sorting as nf; a marker-cell filter as Obs; a declared log
  normalization as Cap_read) running the full schema end to end: solve `N = F_toy(N)`,
  certify the fixed point with an interval Banach check (`F_toy(I) ⊆ I` and
  `|F_toy'| ≤ L < 1`, evaluated with `mpmath.iv`), and emit
  `runtime/toy_readback_certificate.json`.
- [`test_toy_readback.py`](test_toy_readback.py) — tests: enumeration matches the
  closed-form sector count, the toy fixed point exists with a certified enclosure, the
  certificate schema keys are present, output is deterministic, and a deliberately
  non-contracting variant is rejected.

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
