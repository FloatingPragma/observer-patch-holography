# Agent instructions: turn the constructed map into an OPH-native physical map

The pole-map algebra is now fixed. Do not edit its sign convention, root
branch, strict-order rule or fail-closed v1 policy to make an input pass.
Replacing the conditional fixture is necessary but not sufficient: physical
promotion belongs to a separately versioned aggregate verifier that resolves
independent hash-bound evidence.

## 1. Freeze the claim class

Use one of these mutually exclusive classes:

```text
EXTERNAL_SM_EFT_VALIDATION
OPH_NATIVE_PHYSICAL
CONDITIONAL_POST_EXPOSURE_DIAGNOSTIC
TARGET_COMPARISON_ONLY
```

The current fixture is `CONDITIONAL_POST_EXPOSURE_DIAGNOSTIC`. A clean
imported Standard Model packet may validate the field-theory engine but cannot
be relabeled OPH-native.

The `evidence_gates` keys in the input template are declarations made by the
caller. V1 records them as `unverified_evidence_claims`; they never change its
fail-closed output gates.

## 2. Produce the OPH+FJ renormalized input packet

The source must emit a canonically normalized one-Higgs Standard Model action
and exactly one VEV identity. For the preferred branch:

\[
V(H)=m^2H^\dagger H+\lambda(H^\dagger H)^2,
\qquad m^2(Q)=-\lambda(Q)v_F(Q)^2,
\qquad v_F>0.
\]

The source declaration must distinguish this from the Landau effective-
potential minimum `v_L` and the Fermi coordinate `v_GF`. Add the following
receipt fields:

```text
potential_normalization
field_normalization
scheme = SM_MSbar_FJ
Q
m2, lambda, v_F
proof that m2 + lambda*v_F^2 = 0 in the declared interval
source DAG and no-target-ancestry result
```

Do not rename the existing OPH transmutation coordinate `v(P)` to `v_F`.
Construct the action-level bridge or retain the field as `v_chart`.

## 3. Implement both one-loop tadpole routes

Build two separately implemented producers that consume the same matching
packet:

1. **Direct FJ:** parameter-defined `v_F`, explicit tadpole diagrams and
   tadpole counterterms retained.
2. **Converted tadpole-free:** compute in the declared tadpole-free coordinate,
   apply the complete finite parameter transformation, then re-expand to
   strict one loop.

For every tree parameter `p_a`, apply

\[
O_F^{(1)}=O_L^{(1)}+\sum_a\delta p_a^{(1)}
\frac{\partial O^{(0)}}{\partial p_a}.
\]

The converter must transform masses appearing inside loop integrals,
counterterms, wave-function terms and tadpoles before truncation. Compare the
complex coefficients, separately:

```text
Re Delta_WW(w)
Im Delta_WW(w)
Re Delta_ZZ(z)
Im Delta_ZZ(z)
Delta_AA(z), Delta_AZ(z), Delta_ZA(z)
```

A match only after taking the nonlinear square root is insufficient.

## 4. Replace the D10/D11 split boundary with a source-closed EFT packet

Emit, for every scale interval:

```text
local action and active field census
spin/statistics, representations and multiplicities
kinetic and hypercharge normalization
renormalization scheme
beta functions with contribution-level monomial mask
threshold masses and threshold scales
finite decoupling and DRbar/MSbar conversion maps
matching remainder intervals
J_match
```

The checker must independently recompute one-loop beta coefficients from the
field census. A package using MSSM beta coefficients cannot be fed directly
into a pure-SM pole engine without the full threshold/decoupling map.

D11 must independently emit `lambda(Q)`, `yt(Q)` and `Q/E_star` from the same
source branch. Reject any seed whose core values or Jacobian descend from a
calibration/target surface, regardless of a local
`predictive_promotion_allowed=true` flag.

## 5. Build the independent self-energy engine

The production engine must not call SMDR or TSIL. It must generate the one-loop
renormalized complex coefficients in a general linear `R_xi` gauge augmented
by the five nonlinear gauge parameters. Retain:

```text
physical gauge bosons
Goldstones
Faddeev--Popov ghosts
longitudinal/mixed blocks
fermions and Higgs
counterterms
explicit FJ tadpoles
absorptive parts on a frozen analytic sheet
```

It must output the complete neutral matrix, not only `Delta_ZZ`:

\[
\Gamma_N^T=
\begin{pmatrix}
s-\Delta_{AA}&-\Delta_{AZ}\\
-\Delta_{ZA}&s-z-\Delta_{ZZ}
\end{pmatrix}.
\]

Use certified complex balls or intervals. Re-run at 128, 192 and 256 bits;
enclosures must nest and shrink. The final radius, not decimal agreement,
defines the numerical tolerance.

## 6. Verify finite-order gauge identities

For each frozen gauge point, check the strict coefficients before the square
root readout. Required checks are:

```text
extended BRST master identity
charged Slavnov--Taylor sector
neutral Slavnov--Taylor sector
Gamma_AA^T(0)=0
Gamma_AZ^T(0)=0
charged Nielsen identity
neutral matrix Nielsen identity
determinant Nielsen identity
UV-pole cancellation
```

A multi-gauge numerical grid is a mutation/stress test. The proof-bearing
object is the symbolic or certified-ball identity. Do not cite an all-order
Nielsen theorem as evidence that the finite implementation included the right
diagrams and counterterms.

## 7. Supply the source law and covariance

Use

\[
\mathrm{Law}(z)=\delta(z-z_0),\qquad C_z=0
\]

only after global root uniqueness, selector uniqueness and exact primitive
selection are proved. Otherwise provide a target-independent joint law or a
hash-pinned ensemble with source-justified weights. Never infer weights from
interval widths and never assign equal weights to unresolved branches.

Keep these separate from `C_z`:

```text
root enclosure
branch ambiguity
matching truncation
pole truncation
scale variation
rounding/solver error
```

Then compute

\[
C_{\rm pole,source}=J_{\rm pole}J_{\rm match}J_{\rm source}C_z
J_{\rm source}^TJ_{\rm match}^TJ_{\rm pole}^T.
\]

## 8. Keep the strict order guard

The map rejects `1.5`, `2`, `2.5` or any unspecified mixed label. Do not weaken
this guard. A future two-loop map must be a separately versioned schema and
must include at least

\[
-\Delta_{ZA}^{(1)}(s)\Delta_{AZ}^{(1)}(s)/s
\]

plus pole-expansion derivative terms and a contribution mask identical across
the VEV converter and pole engine.

## 9. Close the clock last

Until the target-free physical clock is proved, the native output is

\[
s_W/E_\star^2,\qquad s_Z/E_\star^2,
\]

not a mass in GeV. A GeV fixture may be used in
`EXTERNAL_SM_EFT_VALIDATION` or `CONDITIONAL_POST_EXPOSURE_DIAGNOSTIC`, with
that provenance visible in every receipt.

## 10. Algebraic acceptance command

After replacing `data/oph_fj_input_template.json` with a real packet, run:

```bash
python3 src/wz_pole_map.py data/oph_fj_input.json \
  --output outputs/oph_native_strict_1l_pole_map_receipt.json
python3 checker/check_receipt.py \
  --fixture data/oph_fj_input.json \
  --receipt outputs/oph_native_strict_1l_pole_map_receipt.json
```

In the future aggregate pipeline, promotion is licensed only when its checker reports `PASS` and
`physical_promotion_allowed=true`. A numerical pole with any false source or
field-theory gate remains conditional.

For the hardened v1 schema, `physical_promotion_allowed` is always false, so
the command above certifies only subject-bound pole algebra. The quoted
promotion condition is a requirement for the future aggregate verifier, not
an output this checker can produce. That verifier must load the referenced
source, matching, FJ, covariance, gauge/BRST, clock, and ancestry receipts;
recompute their hashes and schemas; require exact subject/order/mask identity;
and derive status rather than accept it from the pole input.
