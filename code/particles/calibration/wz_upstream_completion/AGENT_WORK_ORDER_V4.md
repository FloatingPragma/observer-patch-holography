# Agent work order v4 — build the missing W/Z source-to-pole bridge

## Global rule

Do not change the downstream strict pole map to accommodate an upstream defect. The frozen inverse-propagator sign, strict series, A–Z mixing order and square-root branch are part of the contract.

## Workstream A — Source/action producer

**Owner:** source-theory agent  
**Target:** `SM_EFT_ACTION_1`

1. Run with no target data mount.
2. Emit the complete field census and the exact action AST.
3. Freeze `gprime_convention = SM_hypercharge`, not GUT `g1`.
4. Emit every retained/excluded operator and its source reason.
5. Emit `v_chart` and `v_F` as different types unless an equality theorem is present.
6. Emit full `Yu,Yd,Ye`, basis rotations and CKM; a top-only packet fails.
7. Emit source DAG and canonical hashes.
8. Never emit a mass core or an inverse target adapter.

**Stop condition:** if the source cannot emit the action or full Yukawa packet, mark the lane `OPH_CHART_ONLY`; do not fabricate placeholders as zeros.

## Workstream B — EFT/RG matching

**Owner:** matching agent  
**Target:** `EFT_MATCHING_1`

1. Partition the scale axis into intervals.
2. For each interval emit active fields, representations, spin/statistics and operator basis.
3. Derive one-loop beta coefficients independently from the census.
4. Freeze schemes and every finite `DRbar ↔ MSbar` or decoupling map.
5. Carry interval Jacobians and deterministic remainders.
6. Reject MSSM beta coefficients in a pure-SM interval.
7. Output one canonical `SM_MSbar_FJ(Q)` vector.

## Workstream C — Feynman rules and diagram universe

**Owners:** rule agent A and independent rule agent B  
**Target:** `RULE_EQUIVALENCE_1`

- A differentiates the action AST.
- B uses a separately authored implementation and may not import A's rule tables.
- Canonicalize vertices by field multiset, momentum convention and tensor basis.
- Compare exact vertex hashes and coefficients.
- Enumerate the complete one-loop 1PI two-point universe independently of both engines.

## Workstream D — Counterterms and chiral restoration

**Owner:** algebraic-renormalization agent  
**Target:** `RENORMALIZATION_ST_1`

1. Use `d=4-2epsilon` and explicitly freeze BMHV gamma5 or a fully equivalent declared prescription.
2. Substitute bare maps into the complete invariant + GF + ghost + external-source action.
3. Expand to first order; do not hand-enter CT vertices.
4. Compute UV poles exactly.
5. Solve the linearized ST restoration equations for finite symmetry-restoring CTs.
6. Emit the counterterm-basis rank, null space and selected normalization conditions.

## Workstream E — Direct FJ engine

**Owner:** engine A  
**Target:** `FJ_DIRECT_1`

- Generate the nonlinear gauge and ghosts from the hashed action.
- Keep explicit tadpoles.
- Compute one-point, transverse, longitudinal, Goldstone and mixed blocks.
- Keep complex absorptive parts and exact UV poles.
- Emit per-diagram records and summed coefficient balls.

## Workstream F — Converted engine

**Owner:** engine B  
**Target:** `FJ_CONVERTED_1`

- Independently calculate in the declared tadpole-free chart.
- Implement the complete finite map `p_L = p_F + hbar Δp`, with the loop
  marker kept distinct from the Higgs fluctuation and Planck's constant.
- Transform fields, parameters, mass arguments, counterterms and self-energies.
- Re-expand before truncation.
- Use a separate scalar-integral implementation.
- Never call engine A.

## Workstream G — BRST/ST/Nielsen checker

**Owner:** small-checker agent  
**Target:** `WARD_ST_NIELSEN_1`

The checker may read coefficient records but may not import either generator. Recompute:

- BRST nilpotence;
- action invariance and `S_gf+gh=sPsi`;
- the complete QCD plus electroweak gauge-fixing/ghost sector (the 45-point
  electroweak grid may hold `xiS=1` only at strict electroweak one loop);
- UV cancellation;
- linearized ST identity;
- all generated two-point ST projections;
- photon Ward zeros;
- charged and neutral Nielsen identities;
- determinant Nielsen identity;
- direct/converted FJ equality;
- strict pole coefficients across the frozen gauge grid.

Use exact algebra first and complex balls second. A float-only threshold is not a proof.

## Workstream H — Analytic continuation and physical amplitude

**Owner:** pole-interpretation agent  
**Target:** `PHYSICAL_CURRENT_POLE_1`

1. Freeze cut conventions and second-sheet vector.
2. Isolate one zero by argument principle/Rouché balls.
3. Prove the derivative/nondegeneracy ball excludes zero.
4. Build at least one BRST-invariant physical scattering amplitude in a charged-current channel and one in a neutral-current channel.
5. Verify the same pole occurs with nonzero amplitude residue.
6. Do not demand positivity of the complex propagator residue.

## Workstream I — Source law/covariance

**Owner:** uncertainty agent  
**Target:** `SOURCE_COVARIANCE_1`

- Use `Cz=0` only after all source parents are mathematically unique and deterministic.
- Otherwise emit a justified joint law.
- Keep source covariance, imported nuisance covariance and nonstochastic theory errors separate.

## Workstream J — Clock and unit attachment

**Owner:** metrology agent  
**Target:** `SOURCE_CLOCK_1`

Until the source reference transition is complete, output only dimensionless pole coordinates. The existing GeV checksum is not a source clock.

## Workstream K — Target firewall

**Owner:** provenance agent

- Source jobs run in a container with no PDG/target files, no target environment variables and no writable comparison path.
- Comparison runs receive only immutable source receipts.
- Because W/Z targets are already known, classify this campaign `post_exposure_validation`.

## Machine receipt contracts

The current directory has nine checklist schemas for ten instances; W and Z
reuse the physical-pole schema. Do not edit a schema to accommodate a failed
producer. These v1 files cannot authorize promotion. Build proof-bearing v2
schemas and an independent resolver that consumes no producer code, resolves
every artifact, recomputes every digest, and derives rather than trusts result
flags.

## Required CI jobs

```text
schema
source_dag_blacklist
action_ast_hash
anomaly_arithmetic
beta_from_census
vertex_equivalence
diagram_universe_completeness
uv_pole_cancellation
st_master
ward_zero
nielsen_matrix
fj_direct_vs_converted
strict_pole_series
complex_ball_precision_nesting
physical_current_pole
mutation_suite
promotion_conjunction
```

## Delivery rule

A partial success is acceptable and should be reported precisely. Never set a parent gate true because a downstream number is numerically plausible.


## Frozen computational sidecars

- nonlinear gauge grid: `data/nonlinear_gauge_grid_v1.json`
- count: 45
- canonical SHA-256: `6e0265eda5d55e5def430548441777f8809671fbbddeeb1ea4c5f0ace588abfd`
- receipt dependency DAG: `data/receipt_dependency_dag_v4.json`
- DAG canonical SHA-256: `79a8ec98d77286e8063e57b2395ddc36f9507cdc4e0721eb7995b0b8f519dfca`
