# Fractional Quotient-Sector Receipts

This folder mirrors the fractional exciton/FQAH sandbox in the paper-stack
particle code tree. It emits receipt scaffolds for quotient-sector normal
forms, Hamiltonian promotion checks, topological ledgers, optical line-fan
identifiability, refinement compatibility, and no-target-leak guards.

It is diagnostic only. A real material claim remains blocked until a
material-specific Hamiltonian/source-law proof and frozen sample comparison are
available.

```bash
python3 fractional/build_fractional_quotient_receipts.py
```

## Lean verdict certificate

The four quotient-correctness gates — `CANONICALIZER_IDEMPOTENCE`,
`REPRESENTATIVE_INVARIANCE`, `QUOTIENT_LUMPABILITY`, and
`NO_ORBIT_SIZE_BIAS` — are read from
`fractional_quotient_certificate.json`, whose verdicts are computed by the
certified decision procedures in
`Lean/ObserverPatchHolography/QuotientLumpability.lean` on the exact sandbox
instance. The Lean module pins the rendered certificate with `#guard_msgs`,
so instance or verdict drift fails the Lean build; the builder verifies the
checked-in file against that pin before reading any verdict and fails closed
(all four gates `False`) on any mismatch.
`fractional_quotient_negative_control_certificate.json` carries the pinned
rejecting verdicts and their checker witnesses; pass it via `--certificate`
to see the gates go `False` and the claim drop to `DIAGNOSTIC_ONLY`. The
remaining sandbox gates are declared scaffold, labelled
`DECLARED_SANDBOX_SCAFFOLD` in `receipts.json`.
