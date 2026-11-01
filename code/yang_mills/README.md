# Yang--Mills collar-gap certificates

## Source-contract verifier

`verify_collar_gap_certificate.py` checks the declared finite row arithmetic used
by the issue-306 theorem witness. Both certificate interfaces use
`collar_gap_contract.py` for exact rational accounting. They require source
descriptions, positive declared rate bounds, refinement targets inside the
declared table, and normalized conditional-probability pairs with exact
total-variation distances below the supplied bounds. Every influence and its
integer multiplicity contributes to the row sum. The reported expression is

```text
gap_lower = c_floor * (1 - eta_upper).
```

Run the bundled contract witness with:

```bash
python3 code/yang_mills/verify_collar_gap_certificate.py \
  code/yang_mills/certificates/issue_306_theorem_contract_witness.json
python3 -m pytest code/yang_mills/test_collar_gap_certificate.py
```

The bundled JSON is marked `theorem_contract_witness`. It tests the arithmetic
and the explicit constant `3/8`. Source descriptions, rate bounds and the
completeness of conditional rows remain supplied premises: comparing one pair
does not establish that it maximizes influence over a source kernel. The
checker does not reconstruct a Gibbs measure or prove the rate, refinement or
continuum hypotheses of the paper's theorem.

Both interfaces reject `scope: physical_source_receipt`, including fixtures with
asserted provenance or continuum fields. Neither implements the source,
transfer, OS/noncollapse or continuum verifications needed for that status.
Accepted outputs have `physical_clay_receipt: false`.

## Finite calibration fixture

`finite_collar_gap_certificate.py` expands and checks a finite Ising
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

A compact family repeats the complete template with self-targets and independent
copies of its data. Supply either this family or an explicit table. A template
must omit generated IDs and targets; it cannot override them. Boolean or
fractional multiplicities, malformed source fields, and duplicate JSON keys
are rejected. Receipt replay checks every field and its JSON type, including
the full manifest and expanded-table hashes.

The retained audit under [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033)
reproduces a lost-influence defect at `de60560b`: compact bounds `3/4` and `1/4`
were reduced to `1/4`, so a noncontractive table returned a positive gap.
The witness interface also truncated fractional multiplicities and both
interfaces permitted unsupported physical promotion. Independent controls
maximize probability differences over all events, check heterogeneous tables,
and construct a two-spin heat-bath operator with a complete exact eigenbasis.
They verify attained gap bounds at both signs of the correlation, within
`2^-100` of the contractivity boundary, and at rate scales beyond binary64.
Three-spin controls derive every conditional row from the full supplied law
and check the entire weighted operator inequality by exact Schur elimination,
with unequal rates and two influences per site. A singular-support control
constructs a nonconstant centered zero-energy mode and requires rejection at
`eta = 1`; its zero gap is calculated from the conditional expectations.

```bash
python3 -m pytest -q code/yang_mills/tests/test_collar_certificate_accounting.py \
  code/yang_mills/tests/test_finite_collar_gap_certificate.py \
  code/yang_mills/test_collar_gap_certificate.py -W error
```

The 244-type receipt reproduces byte-for-byte; its `1/2` floor and the witness's
`3/8` floor are unchanged. The finite Z2 transfer diagnostic, frozen evidence,
paper and book claims, and registry payloads are unchanged. These corrections
establish no physical mass gap.

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
  fit leaves a 3.8% to 24% relative residual at `L = 3`;
* for the Kogut-Susskind Hamiltonian the single-flip form is exact with
  fiber-dependent rates `c_l(o) = lam (r + 1/r)`, `r = Omega(o)/Omega(X_l o)`,
  so the scalar cross-fiber equality fails (rate spread 1.04 to 2.05 at
  `L = 3`) while `c_l(o) >= 2 lam` remains an analytic positive floor;
* among the interacting tested points, the Dobrushin sum `eta_*` is below one
  at Wilson `beta_s = beta_t = 0.1` and Kogut-Susskind `lam = 2`.

### Numerical calculation and replay

The free Hamiltonian is evaluated analytically as
`H = (log(coth beta_t)/2) sum_l (I - X_l)`, with a uniform Perron state.
Only exact `beta_s == 0` takes this route. For interacting inputs, forming the
dense transfer matrix and then diagonalizing it loses small kinetic eigenvalues.
Instead, let `G` be the gauge-mask group and `C = G^perp` its electric cycle
space. With `n` orbits and `E` links,

```text
F[o,z] = (-1)^(rep_o dot z) / sqrt(n),  z in C,
K = (2 cosh beta_t)^E F diag(tanh(beta_t)^|z|) F^T,
A = diag(exp(beta_s P/2)) F diag(tanh(beta_t)^(|z|/2)),
T = (2 cosh beta_t)^E A A^T.
```

Character orthogonality gives `F^T F = I`: averaging over orbit representatives
is the full-configuration character sum divided by the common gauge-orbit size.
The singular vectors of `A` therefore give the transfer eigenvectors, and
`H = U diag(2 log(s_max/s_i)) U^T`. This computes the same Wilson operator;
it changes no couplings or tested points. A common row scale is removed and
restored in `lambda_max`.

The [LAPACK Jacobi SVD](https://www.netlib.org/lapack/explore-html/d8/d78/group__gejsv_gaca7ba7f1e8002c7a1d5bffa4ccbb541f.html)
is used with `JOBA=E`, whose relative singular-value accuracy is controlled by
the column-equilibrated factor. Here that factor is `diag(exp(beta_s P/2)) F`,
with condition number `exp(abs(beta_s) (max(P)-min(P))/2)`, independent of the
tiny kinetic eigenvalues. The source checks the spatial condition, retained
rank, solver status, scaled column reconstruction and Perron separation.
It also checks amplification by the Doob division: an entrywise Hamiltonian
uncertainty estimate `eps n ||H||_2` becomes a Frobenius estimate
`eps n ||H||_2 ||1/Omega||_2` for unit-norm `Omega`. The relative estimate
must fit the same resolution policy. For example, interacting `L=2, (3,3)`
has an accurate Hamiltonian but unresolved Doob entries and is refused;
every point in the original retained grid remains supported.
Row conservation is then checked relative to each row's own absolute mass.
An unresolved row causes refusal; the code does not repair its diagonal or
mask a small bad row with the scale of a larger one.
The Kogut-Susskind ground-state solve is scaled and checked against its ground
separation and smallest Perron component too. This support gate is conservative:
it can refuse an accurate small-coupling result, including `L=2, lam=.01`.
The retained grid and independently controlled `lam=.05` and `lam=1e4` pass.
Both interacting families also require the computed Dobrushin influence to
exceed a dimension-scaled conditional-probability contrast floor at the requested
relative resolution. Thus nearly uniform computed populations cannot turn a
lost interaction response into a claimed zero. Only the source-exact free law
bypasses that check. Masked operators and unused parameter names are rejected.
Narrowing an extended-precision coupling must preserve its value within ordinary
relative rounding; becoming zero or a badly rounded subnormal causes refusal.
These are numerical safeguards at a `1e-7` resolution policy, not certified
interval bounds or a relative-error guarantee for every derived scalar at every
accepted noncanonical coupling. Tiny cancellation-derived diagnostics remain
floating estimates; the retained grid is checked component by component across
platforms and against the independent controls below.
The generic rounded-matrix logarithm separately refuses an
unresolved bottom spectrum or Perron state. Near-zero positive couplings and
an unresolved leading eigengap are not replaced by an exact free model.

The defect was reproduced on `de60560b` before correction:

| Original input | Incorrect result | Independent control / correction |
| --- | --- | --- |
| Free `L=2, beta_t=.01` | Nonzero constant-rate residual and multi-flip mass | Exact free Hamiltonian, including every entry and rate |
| Free `L=2, beta_t=18` | `pi_max/pi_min` about 117, `eta_*` about 3.45 | Uniform law, `eta_* = 0`, heat-bath gap 2 |
| Interacting `L=2, beta_s=beta_t=.01` | About 0.4% full-H error | Original-input 60/90-digit transfer calculations |
| Interacting `L=3, beta_s=beta_t=.1` | Outside-single-flip mass about 6413.72 | Corrected value about 2152.32; fit residual .03777034 |

Floating defects depend on LAPACK and platform; the analytic controls do not.
The last row corrects live numerical evidence while preserving the finite-grid
constant-rate failure. The paper's rounded residuals `.038`, `.084`, `.24`
and the exact Kogut-Susskind rate theorem remain valid. No physical mass-gap
claim, claim-registry status, Lean theorem or frozen registration changes.

`verify_z2_finite_transfer_receipt.py` validates strict JSON, the complete grid,
source and run digests, scalar identities, and every replayed component.
Nonzero JSON lexemes that underflow to zero are rejected before binding checks.
Necessary normalization, least-squares, total-variation, projector-spectrum and
reversible fiber-mass bounds reject impossible summaries even after rehashing.
Substantive quantities use relative-only tolerance; absolute allowances apply
only to specified zero diagnostics and scale with the fresh calculation.
Replay is not an independent interacting solver. Independence comes from the
original-input high-precision, full-configuration, analytic near-free and
positive Perron-iteration controls in the tests. Ubuntu CI executes
these controls and the full grid, rather than merely collecting the tests.

This is a finite-grid diagnostic on a toy gauge system (`physical_clay_receipt:
false`), not a universal no-go or a compact-simple-gauge receipt. Anisotropic
couplings, the anisotropic Wilson-to-Hamiltonian limit, alternate transfer
objects, and fiber-dependent positive rates remain viable. The next proof
target is the quotient-space variable-rate approximate-tensorization lemma.

```bash
python3 code/yang_mills/z2_finite_transfer_receipt.py --L 2 3 \
  --output code/yang_mills/receipts/z2_finite_transfer_receipt.json
python3 code/yang_mills/kogut_susskind_fiber_rate_instances.py
python3 -m pytest code/yang_mills/tests/test_z2_finite_transfer_receipt.py
```

The second command refreshes the live exact-instance receipt's upstream byte
binding; its rational rates and gaps do not change. Both producers write LF
bytes on every platform. For a separate portable replay, write a fresh transfer
receipt to a temporary path and pass the retained and fresh paths to
`verify_z2_finite_transfer_receipt.py`.
