# Finite collar-gap certificate

This directory contains the executable finite certificate layer for the
source-defined collar generator

\[
L_{r,b}=\sum_{v\in V_r}c_{v,r,b}(I-P_{v,r,b}).
\]

All numerical proof inputs are exact rationals.  The checked calibration table
has 244 active types, `c_floor = 1`, `eta_upper = 1/2`, and consequently
`gap_lower = 1/2`.  It is deliberately a finite Ising calibration, not a
physical compact-simple-gauge Yang--Mills receipt.

```bash
python3 code/yang_mills_gap/collar_gap_certificate.py certify \
  --manifest code/yang_mills_gap/manifests/atomic_4d_ising_calibration.json \
  --output code/yang_mills_gap/receipts/atomic_4d_ising_calibration.receipt.json
python3 code/yang_mills_gap/collar_gap_certificate.py verify \
  --manifest code/yang_mills_gap/manifests/atomic_4d_ising_calibration.json \
  --receipt code/yang_mills_gap/receipts/atomic_4d_ising_calibration.receipt.json
```

The `physical_compact_gauge_uninstantiated.json` manifest must fail closed.
It has no source type table, physical provenance, or continuum receipts.
