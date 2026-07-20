# Strict-one-loop Physical W/Z complex-pole map

Prepared: 2026-07-20

Integrated archive SHA-256:
`fb871c5ac810f8703fd49b9fdcf621096c1bd60bc0c426d34ba42dc2fd3bae8c`.

This source-of-truth port preserves the fail-closed checker hardening from the
audited handoff under `survival-proof-3/pro_strict_one_loop/`. The archive
digest and pristine upstream file manifest are retained here for provenance;
runtime verification depends only on files in this directory.

## Result

This package supplies the missing **pole-map kernel**. It turns a declared,
renormalized strict-one-loop input packet into the charged and neutral complex
poles without borrowing the experimental running-width convention and without
attaching a one-loop receipt to the archived order-2.5 row.

The map is

\[
w=\frac{g^2v_F^2}{4},\qquad
z=\frac{(g^2+g'^2)v_F^2}{4},
\]

with inverse-propagator convention

\[
\Gamma_W^T(s)=s-w-\Delta_{WW}^{(1)}(s)+O(h^2),
\]

\[
\Gamma_N^T(s)=
\begin{pmatrix}
s-\Delta_{AA}^{(1)}(s)&-\Delta_{AZ}^{(1)}(s)\\
-\Delta_{ZA}^{(1)}(s)&s-z-\Delta_{ZZ}^{(1)}(s)
\end{pmatrix}+O(h^2).
\]

Every `Delta` in the input already includes the one-loop factor,
counterterms, tadpoles and the declared retained-term mask. The strict poles
are therefore

\[
s_W^{[1]}=w+\Delta_{WW}^{(1)}(w),\qquad
s_Z^{[1]}=z+\Delta_{ZZ}^{(1)}(z).
\]

The neutral mixing product

\[
\Delta_{ZA}^{(1)}\Delta_{AZ}^{(1)}
\]

is retained in the neutral-matrix diagnostic but excluded from a strict
one-loop root because it has loop power two. The code records the two-loop
effective-inverse term as

\[
-\frac{\Delta_{ZA}^{(1)}(s)\Delta_{AZ}^{(1)}(s)}{s}.
\]

For each pole, the package emits two deliberately separate readouts:

1. the strict coefficients
   \[
   \delta M_V^{(1)}=\frac{\Re\Delta_V^{(1)}}{2m_{V,0}},\qquad
   \Gamma_V^{(1)}=-\frac{\Im\Delta_V^{(1)}}{m_{V,0}};
   \]
2. the exact coordinate transform of the one-loop-truncated complex pole,
   \[
   s_V=(M_V-i\Gamma_V/2)^2.
   \]

The latter contains harmless kinematic terms beyond strict one loop; it must
not be used in place of the strict coefficients for a finite-order Nielsen or
gauge-independence test.

## What is evaluated now

`data/conditional_smdr_order1_fixture.json` reconstructs the archived SMDR
v1.3 order-1 result at `Q=160 GeV`. It is a regression fixture, not an
independent self-energy producer. The evaluated output is:

| quantity | W | Z |
|---|---:|---:|
| tree mass (GeV) | 79.969486678981 | 90.863694261919 |
| strict one-loop mass readout (GeV) | 80.374161202712 | 90.680036075608 |
| strict one-loop width (GeV) | 2.007425074735 | 2.402420059845 |
| exact square-root readout of truncated pole, M (GeV) | 80.379345721647 | 90.687836666971 |
| exact square-root readout of truncated pole, Gamma (GeV) | 1.997189095430 | 2.407078720028 |

The reconstructed truncated poles are

\[
s_W^{[1]}=(6459.842027569383-160.532752773045i)\;\mathrm{GeV}^2,
\]

\[
s_Z^{[1]}=(8222.835212344102-218.292761806439i)\;\mathrm{GeV}^2.
\]

## Scientific status

The package closes the algebraic **pole-map kernel** layer, not the OPH-native
source layer. The supplied fixture fails closed because it has:

- no OPH-selected FJ renormalized VEV;
- no complete explicit-tadpole versus converted-tadpole-free equality receipt;
- no target-clean D10/D11 matching packet;
- no licensed source law/covariance;
- no independent general-gauge self-energy engine;
- no finite-order BRST, Slavnov--Taylor, Ward and Nielsen receipt;
- no physical source clock; and
- known target ancestry in D11.

Accordingly the generated status is

```text
CONDITIONAL_STRICT_1L_POLE_MAP_NOT_OPH_NATIVE_PHYSICAL
```

and `physical_promotion_allowed=false`.

The template `data/oph_fj_input_template.json` is the exact interface the
future evidence pipeline must fill. Its input booleans are only unverified
claims. This v1 producer and schema can never self-promote, even if every
caller-supplied boolean is true. A separately versioned aggregate verifier must
resolve every hash-bound evidence object, validate it independently, and bind
it to the exact numerical subject before any physical promotion is possible.

The upstream checker originally trusted those booleans and did not bind the
receipt's numerical subject to its fixture. A self-consistent unrelated pole
could therefore pass with physical promotion. The integrated checker now
binds the fixture, parameters, self-energies, masks, tolerances, covariance,
derived fields, neutral diagnostics, status, blockers, and promotion policy.
It also imposes checker-owned tolerance caps. The adversarial regression tests
freeze these repairs.

## Run

```bash
python3 run_all.py
```

The command regenerates the receipt, executes 16 unit/adversarial tests,
validates the JSON Schema and runs an independent checker that does not import
the map implementation.

## Files

- `src/wz_pole_map.py`: strict-one-loop map and order guard.
- `checker/check_receipt.py`: small independent recomputation checker.
- `schemas/physical_wz_strict_1l_pole_map_receipt.schema.json`: Draft 2020-12 schema.
- `data/conditional_smdr_order1_fixture.json`: conditional numerical fixture.
- `data/oph_fj_input_template.json`: required OPH+FJ producer interface.
- `outputs/conditional_strict_1l_pole_map_receipt.json`: generated receipt.
- `tests/test_wz_pole_map.py`: regression and adversarial tests.
- `THEOREMS_AND_PROOFS.md`: formal statements, proofs, limitations, and
  numerical derivation.
- `IMPLEMENTATION_CONTRACT.md`: ordered instructions for completing the physical producer.
