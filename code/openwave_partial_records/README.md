# Partial OpenWave results on native records

Native records of the observer construction carry part of the test for
these OpenWave rows. Each item links to its section of the [proof](PROOF.md).

- **[Row 3, de Broglie clock](PROOF.md#row-3-de-broglie-clock).** A native
  internal reference on the compact charged state: the walk W = S(2AA*-I)
  rotates span{u,w} by theta = acos(209/210), and a native triangle cycle
  makes the theta harmonic visible in ordinary coordinate records, with all
  source and coordinate outcomes retained.
- **[Row 4, particle stability](PROOF.md#row-4-particle-stability).** The
  exact invariant span of that walk stays inside one twelve-port cell at every
  integer time.
- **[Row 6, spin](PROOF.md#row-6-spin).** The exact CAR field embedding J with
  V_h J = J(K_h tensor I_2) carries a conserved internal field spin together
  with the compact energy certificate E6 on the same state.
- **[Row 8, annihilation](PROOF.md#row-8-annihilation).** A signed pair meets
  and is removed by native contact. The contact state has a complete energy
  pointer with mean energy 11.051928582384567 in transfer units, and the
  contact current releases a record field that propagates by an exact
  discrete wave equation.
- **[Row 20, electric force](PROOF.md#row-20-electric-force).** Residence
  records of one unit, read from its native currents, give the far-field 1/r
  form: fitted exponent 1.035144 with 95% profile interval
  [0.983376, 1.086373]. The strength alpha is a separate question.
- **[Row 31, orbital levels](PROOF.md#row-31-orbital-levels).** Native
  first-hit feedback creates an attractive well, and its controlled limit has
  hydrogen-like standing shells e_n = -(lambda c_*)^2/(4n^2) with n^2 orbital
  states per shell.

```sh
python3 code/openwave_partial_records/verify.py
python3 tools/build_openwave_partial_records_claim.py --check
```

`verify.py` checks every archived receipt in `evidence/` against its original
sealed manifest ([`EVIDENCE.json`](EVIDENCE.json)), binds the sealed receipt
fields and runs [`recompute.py`](recompute.py), which rederives the results
from the carrier graph, the native rule and the archived records. The claim
ledger is
[`claims/openwave_partial_records.json`](../../claims/openwave_partial_records.json);
every row in it is partial and officially untested.
