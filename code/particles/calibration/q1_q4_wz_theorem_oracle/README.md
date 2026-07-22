# QFT-Q1--QFT-Q4 and W/Z theorem oracle

This directory is a public, non-promoting oracle for the algebraic implications
used by the papers. It is not a QFT implementation and does not manufacture an
OPH-native action, chiral measure, Hamiltonian, BV restoration transcript,
physical current amplitude, nonperturbative observable tower, resonance sheet,
or numerical prediction.

The machine namespace deliberately uses `SM_QFT_*` identifiers so these tiers
cannot be confused with the older particle-pipeline `Q1`--`Q3` receipt names.
The graph records that QFT-Q2 and QFT-Q3 are parallel descendants of QFT-Q1.
The strict fixed-parameter W/Z calculation belongs to QFT-Q3; it has no edge
to QFT-Q4. QFT-Q4 consumes a separately supplied nonperturbative tower and its
resonance continuation consumes a separately supplied analytic-sheet packet.

Run:

```bash
python3 verify_symbolic_q1_q4_wz.py
python3 -m unittest -v test_theorem_oracle.py
```

`upstream_import.json` content-addresses the audited handoff. The oracle was
rewritten rather than copied because the supplied promotion checker accepted
forged status fields; only its theorem algebra and source hashes are retained.

