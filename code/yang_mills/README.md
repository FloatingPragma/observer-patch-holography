# Yang--Mills collar-gap certificates

## Source-contract verifier

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

## Finite calibration fixture

`finite_collar_gap_certificate.py` independently checks a finite Ising
calibration family. Its exact rational table has 244 active types,
`c_floor = 1`, `eta_upper = 1/2`, and `gap_lower = 1/2`. It is deliberately
not a physical compact-simple-gauge Yang--Mills receipt: the physical
placeholder manifest must fail closed for missing source, continuum, and
transfer evidence.

```bash
python3 code/yang_mills/finite_collar_gap_certificate.py certify \
  --manifest code/yang_mills/manifests/atomic_4d_ising_calibration.json \
  --output code/yang_mills/receipts/atomic_4d_ising_calibration.receipt.json
python3 code/yang_mills/finite_collar_gap_certificate.py verify \
  --manifest code/yang_mills/manifests/atomic_4d_ising_calibration.json \
  --receipt code/yang_mills/receipts/atomic_4d_ising_calibration.receipt.json
python3 -m pytest code/yang_mills/tests/test_finite_collar_gap_certificate.py
```

## Finite Z2 transfer-receipt diagnostic

`z2_finite_transfer_receipt.py` evaluates the paper's finite
ground-state-transform and cross-fiber receipt on Z2 lattice gauge theory
(L x L periodic spatial torus, Gauss-law sector, one heat-bath collar per
spatial link). The ground-state transform is the Doob transform by the Perron
vector. Two transfer objects are tested: the Wilson transfer matrix with
`H = -log(T / lambda_max)` and the Kogut-Susskind Hamiltonian.

Result on the committed receipt (`receipts/z2_finite_transfer_receipt.json`):

* the receipt is exact at `beta_s = 0` (every rate equals `log coth beta_t`, the dual coupling);
* for every tested interacting Wilson point
  `beta_s = beta_t in {0.1,0.3,0.5,0.7,1}` it fails: the best constant-rate
  fit leaves a 4% to 24% relative residual at `L = 3`;
* for the Kogut-Susskind Hamiltonian the single-flip form is exact with
  fiber-dependent rates `c_l(o) = lam (r + 1/r)`, `r = Omega(o)/Omega(X_l o)`,
  so the scalar cross-fiber equality fails (rate spread 1.04 to 2.05 at
  `L = 3`) while `c_l(o) >= 2 lam` remains an analytic positive floor;
* among the interacting tested points, the Dobrushin sum `eta_*` is below one
  at the weakest Wilson and Kogut-Susskind couplings.

This is a finite-grid diagnostic on a toy gauge system (`physical_clay_receipt:
false`), not a universal no-go or a compact-simple-gauge receipt. Anisotropic
couplings, the anisotropic Wilson-to-Hamiltonian limit, alternate transfer
objects, and fiber-dependent positive rates remain viable. The next proof
target is the quotient-space variable-rate approximate-tensorization lemma.

```bash
python3 code/yang_mills/z2_finite_transfer_receipt.py --L 2 3 \
  --output code/yang_mills/receipts/z2_finite_transfer_receipt.json
python3 -m pytest code/yang_mills/tests/test_z2_finite_transfer_receipt.py
```
