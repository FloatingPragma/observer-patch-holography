# Yang--Mills collar-gap certificate

`verify_collar_gap_certificate.py` checks the finite source-type contract used by the
issue-306 theorem. It uses exact rational arithmetic. A valid certificate must include every
field in the collar source signature, positive rate lower bounds, a closed refinement-transition
table, and exact conditional-probability rows whose total-variation influences have common upper
bound strictly below one. The reported analytic
floor is

```text
gap_lower = c_floor * (1 - eta_upper).
```

Run the bundled contract witness with:

```bash
python3 code/yang_mills/verify_collar_gap_certificate.py \
  code/yang_mills/certificates/issue_306_theorem_contract_witness.json
python3 -m pytest code/yang_mills/test_collar_gap_certificate.py
```

The bundled JSON is deliberately marked `theorem_contract_witness`. It verifies the checker and
the explicit constant `3/8`; it is not evidence that the OPH compact-gauge regulator tower has
those kernels or influence bounds. A Clay-facing receipt must be generated from the actual finite
transfer matrices, use `scope: physical_source_receipt`, and independently pass the continuum,
OS/noncollapse, and transfer/intertwiner gates named in the papers.
