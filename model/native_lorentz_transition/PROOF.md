# Twelve source directions determine a conformal Lorentz transition

Let `u_p` be the twelve exact spatial vectors returned by
[`common_source.slow_geometry.raw_frame`](../common_source/slow_geometry.py).
Their squared radius is `R² = φ + 2`, where `φ = (1 + √5)/2`. Introduce the
declared four-dimensional lifts and bilinear form

\[
k_p=(1,u_p),\qquad G=\operatorname{diag}(-R^2,1,1,1).
\]

The source determines the vectors and their Gram matrix. The time coordinate
and interpretation of these lifts as physical null displacements are
additional premises. In particular, a sphere of directions is not itself a
four-dimensional spacetime.

## The finite transition criterion

**Theorem.** Let `L` be any invertible real `4 × 4` matrix. If
`(L k_p)ᵀ G (L k_p) = 0` for every one of the twelve source ports, then there
is a unique positive real `λ` such that

\[
L^T G L=\lambda G.
\]

Consequently `Λ = L/√λ` is a Lorentz isometry in this raw coordinate chart.
The matrix `L` preserves the future cone exactly when `L₀₀ > 0`; it preserves
spacetime orientation exactly when `det L > 0`. Both tests are needed to
obtain the proper orthochronous component. Neither sign follows from the
twelve null equations alone.

**Proof.** The existing
[`null_tomography_certificate`](../common_source/null_tomography.py)
constructs the `12 × 10` evaluation matrix `A` on symmetric forms and an
exact `9 × 12` decoder `D`. Its identities over `Q(√5)` establish that
`rank A = 9` and `ker A = span(G)`. The identities hold over the reals as
well: a left inverse on the nine-dimensional `B₀₀ = 0` slice establishes
its real injectivity, and `G` is a nonzero kernel vector. This argument
does not extrapolate a finite numerical sample to a universal claim.

The symmetric form `B = Lᵀ G L` has zero evaluations by hypothesis, so
`B = λG`. Its `(0,0)` entry gives the unique value
`λ = B₀₀/G₀₀`. Since `L` is invertible, Sylvester's law of inertia gives
`B` one negative and three positive directions. The value zero is singular;
a negative `λ` would give three negative and one positive directions.
Therefore `λ > 0`, and division by `√λ` proves the isometry statement.

For the orientation test, `L e₀` is timelike because its squared norm is
`−λR²`; hence `L₀₀ ≠ 0`. A Lorentz isometry permutes the two connected
timelike-cone components, so the sign of the time coordinate of `L e₀`
selects its action on the whole future component and its nonzero causal
boundary. Equivalently, from `Λ G⁻¹ Λᵀ = G⁻¹`,
`Λ₀₀² − R² ∑ᵢ Λ₀ᵢ² = 1`. When `Λ₀₀ > 0` and
`t > 0, |x| ≤ Rt`, the transformed time satisfies
`t' ≥ t(Λ₀₀ − R sqrt(∑ᵢ Λ₀ᵢ²)) > 0`. The inverse has the same
time-orientation sign. Finally `det Λ = det L/λ²`; its sign is the
independent spacetime-orientation test. ∎

The Lean theorem
[`nine_null_directions_determine_mod_metric`](../../Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean)
proves the corresponding kernel result for a different, explicit
nine-direction coordinate frame. The twelve native-port kernel is the
executable exact algebraic certificate cited above. The transition theorem
here is a prose proof with exact executable controls, not an additional
Lean theorem.

## Composition, scale and finite-error readouts

If `L₁ᵀGL₁ = λ₁G` and `L₂ᵀGL₂ = λ₂G`, then

\[
(L_2L_1)^T G(L_2L_1)=\lambda_1\lambda_2 G.
\]

Thus the metric factor multiplies, and its positive square root multiplies
as a length scale. Normalized isometries compose in the same order. The
inverse factor is `λ⁻¹`. These identities apply to supplied composable
maps; they do not manufacture a chart-overlap cocycle from records.

For an arbitrary real symmetric form `B`, put
`λ₀ = B₀₀/G₀₀` and `H = B − λ₀G`. Then `H₀₀ = 0` and the same decoder
gives the other nine independent entries as

\[
H_{\rm slice}=D\,(k_p^T B k_p)_{p=1}^{12}.
\]

In particular, if all twelve readings have absolute value at most `ε`,
then

\[
\max_{i,j}|H_{ij}|\leq \lVert D\rVert_\infty\,\epsilon.
\]

The certificate computes `||D||∞ = 1 − √5/5` exactly from the source. This is an
entrywise metric-residual bound in the stated chart. It gives neither a
positive `λ₀` nor a distance from `L` to a Lorentz matrix without extra
invertibility and conditioning information. Arbitrary noisy twelve-tuples
need not lie in the evaluation image; the underlying decoder rejects
incompatible exact readings rather than silently fitting them. The bound
above concerns actual evaluations of a common form, or bounded errors
relative to such evaluations.

## Countercontrols and the source attachment

Spatial parity `diag(1,−1,1,1)` preserves `G` and the future cone, passes
all twelve tests, and has determinant `−1`. Time reversal also passes
every null test. A rank-one map `L = k₀ e₀ᵀ` sends every test vector to
the same null vector while being singular. These examples show why the
invertibility, future-orientation and properness gates are separate.
The executable controls include all three cases, an exact boost, a
rational rotation, a dilation and noncommuting compositions. A cyclic
spatial coordinate rotation also permutes the actual twelve source ports;
the certificate verifies that permutation directly.

The formal
[`OrientedLorentzEquiv`](../../Lean/Geometry/LorentzOverlapCocycle.lean)
contract assumes preservation of the Lorentz pairing and the future
null/timelike loci. It does not require determinant `+1`, so parity is
within that contract. [`CelestialSoldering`](../../Lean/Geometry/CelestialSoldering.lean)
then transports supplied future-null event displacements covariantly on
the quotient of null rays. This finite criterion can discharge the
pairing-preservation part after normalization if actual source chart
transitions are shown to satisfy its hypotheses.

No source observer/scheduler construction of a boost or general chart
transition is supplied here. The native coordinate permutation is a
finite source symmetry, not such a derivation. The null physical readout,
the linearity and invertibility of overlap maps, and a common clock or
volume rule fixing their conformal factors are separate requirements.
The continuous spinor correspondence involving `SL(2,C)`, `PSL(2,C)` and
the proper orthochronous Lorentz group is not formalized by this package;
the finite source `PSL(2,F5)` action is a different group. No propagation
law, interacting relativistic field theory or Einstein equation is
inferred from these transition tests.
