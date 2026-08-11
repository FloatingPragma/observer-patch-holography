# Repair-current orientation verifier

This directory is a hermetic independent audit packet for the simulator's
schema-v3 repair-current payload. It imports no simulator or OPH module. The
vendored payload is byte-identical to the regenerated sibling artifact and has
SHA-256
`7f8ea7ef9c92a50e23207c2fe85d09ed2bce1c1aa539ae9914a9b9edd0df26d6`.

The verifier fails closed on the raw payload hash, schema and post-hoc status,
pinned source hashes and metadata, complete 8-by-8 count table, observer and
transition totals, row sums, antisymmetric current, serialized designations,
the 64 Lean count literals, and unexpected fields. It independently recomputes:

- raw pair maximizer `(3,4)`, gap `1343`;
- row-normalized pair maximizer `(2,3)`, gap `167/172`;
- raw cycle maximizer `(3,4,5)`, gap `1239691068`;
- row-normalized cycle maximizer `(3,4,5)`, gap
  `9391599/57188378`;
- exact pair/cycle sign reversal under count-table transposition; and
- the synthetic symmetric control, including exact detailed balance and zero
  raw pair, raw cycle, and normalized cycle gaps.

The normalized pair statistic is deliberately not treated as a reversible
orientation test. Row normalization of the symmetric control produces a
nonzero pair difference when row sums differ, even though detailed balance and
every cycle-product equality hold exactly.

This packet authenticates the hashes recorded in the vendored payload; it does
not contain or rehash the simulator's large raw inputs. Its result remains a
post-hoc diagnostic and is explicitly ineligible as preregistered validation.

Run from the repository root:

```bash
python3 code/thermodynamics/repair_current_orientation/verify_repair_current_orientation.py
python3 -m pytest -q code/thermodynamics/repair_current_orientation/test_verify_repair_current_orientation.py
ruff check code/thermodynamics/repair_current_orientation
```
