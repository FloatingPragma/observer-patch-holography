# Regulator Gluing Evidence Bundle

Machine gate for the quantum regulator gluing datum of the specialist derivation papers
(Definition `def:qregdatum`, Proposition `prop:regulatedrealization`,
Lemma `lem:qregdatum-realized`). It implements the acceptance test of
paper-audit finding 003 (issue #529):

- every overlap recharting map is displayed explicitly;
- every displayed map is verified as an invertible `*`-isomorphism on its
  declared chart-change domain (exact linear bijection, product preservation,
  adjoint preservation, unit preservation; for unitary-implemented maps,
  exact `U*U = UU* = I` and the `Ad_U` restriction checks);
- every triple-overlap composition law is verified, strictly or as a central
  unit-scalar defect, with the Cech 2-cocycle identity checked on quadruples;
- a bare finite interface projection is rejected with structured reasons.

All matrix entries are exact Gaussian rationals (`Fraction` real and
imaginary parts). No floating point enters any check, so there is no
tolerance policy to state.

## Model format

A bundle JSON contains `witnesses` (documents with
`model: "quantum_regulator_datum"`) and a `countermodel`
(`model: "classical_patch_net"`). A quantum regulator datum declares:

- `patches`: `{id, hilbert_dim}` per regulator cell; the cell algebra is
  `M_{hilbert_dim}(C)`;
- `overlaps`: per edge, a `chart_change_domain`
  (`full_matrix_algebra` or `generated_subalgebra` with explicit generator
  matrices) and a `recharting`, either
  `{type: "unitary", direction: [dst, src], matrix: [[..]]}` implementing
  `Ad_U`, or `{type: "algebra_isomorphism", generator_images: [..]}`;
- `triples`: `{patches: [i, j, k], law: "strict"}` requires
  `U_ij U_jk = U_ik` exactly; `law: "central"` requires
  `U_ij U_jk U_ik^{-1}` to equal the declared central unit scalar
  `expected_defect`;
- `quadruples`: `[i, j, k, l]` lists on which the alternating defect product
  `z_jkl z_ikl^{-1} z_ijl z_ijk^{-1} = 1` is verified.

Scalars serialize as `[re, im]` Fraction strings; matrices as row lists of
such pairs.

## Witnesses

- `strict_cocycle_witness`: patches A, B, C with `M_2(C)` and the strict
  unitary cocycle `Ad_X`, `Ad_Z`, `Ad_XZ` on the triple (A, B, C), plus
  patch D with `M_4(C)` glued through the declared proper subalgebra
  `M_2 (x) 1`, which exercises the chart-change-domain checks.
- `central_defect_witness`: patches P, Q, R, S with `M_2(C)` whose triples
  close up to the nontrivial central scalars `i, -i, -i, i`, with the
  2-cocycle identity verified on the quadruple.
- `bare_interface_projection_countermodel`: the issue #529 net
  `S_i = {0,1,2,3}`, `S_j = {0,1}`, one-bit interface, many-to-one
  projection. The gate rejects it with `NO_RECHARTING_MAP`,
  `STATE_SPACE_BIJECTION_IMPOSSIBLE` (4 versus 2),
  `ALGEBRA_DIMENSION_MISMATCH` (16 versus 4), and
  `MANY_TO_ONE_PROJECTION`.

## Run

From the repo root:

```bash
python3 code/regulator_gluing/verify_regulator_gluing_bundle.py            # exit 0 iff the gate holds
python3 code/regulator_gluing/verify_regulator_gluing_bundle.py --show-maps
python3 code/regulator_gluing/verify_regulator_gluing_bundle.py --bundle path.json
python3 -m pytest code/regulator_gluing/test_regulator_gluing_bundle.py
```

The verifier exits nonzero when any positive witness fails or when the
countermodel passes. The emitted artifact is
`code/regulator_gluing/runs/regulator_gluing_evidence_bundle_current.json`;
it stores the full model and the check report, including every displayed
overlap map and the structured countermodel rejection.

## Claim boundary

The gate certifies the declared finite-dimensional gluing datum and its
composition laws at fixed cutoff. It does not certify continuum limits,
channel-typed (non-isomorphism) transitions, or the genuinely noncentral
crossed-module branch, which is carried by the higher-gauge theorems of the
specialist derivation papers.
