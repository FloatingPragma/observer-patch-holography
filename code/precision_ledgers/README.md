# Precision Ledgers Code

Numerical verification for the two-ledger error accounting of the synthesis
paper, Section 5.14 ("State--observable error propagation and the dynamical
bridge", issue #537).

Canonical runner:

- `gw_ledger_separation.py`: verifies (A) the ledger-separation counterexample
  (zero trace distance, identical bounded-observable expectations, yet a
  nonzero fractional group-velocity error), (B) the conditional bridge lemma
  bounds (Davis--Kahan projector estimate, band-uniform speed estimate with
  `L_B = 1 + 4 M_B`, and the group-delay bound) on randomized C^1 fiber
  families, and (C) the dual-norm inequality behind the state--observable
  ledger.

Run from the repo root:

```bash
python3 code/precision_ledgers/gw_ledger_separation.py
python3 -m pytest code/precision_ledgers/test_gw_ledger_separation.py
```

The current emitted receipt is `runs/gw_ledger_separation_current.json`
(deterministic: fixed seeds, no timestamps).

Scope: this suite validates the mathematics of the ledger separation and of
the conditional acceptance gate only. It does not construct the OPH TT
generator and proves no estimate `eta_dyn <= C_B r_N`; that OPH-specific
bridge remains open, as stated in Section 5.14.
