# C1 canonical Lorentz module handoff

Issue: `#689` (`[C1] Canonical Lorentz module`).

## Exact attained packet

Five issue-scoped modules provide the full intrinsic C1 target and its
exact coordinate handoff to the existing Einstein tensor stack.

| Module | Exact content |
|---|---|
| `Geometry/CanonicalLorentzModule.lean` | `Herm2 = ℝ × (Fin 3 → ℝ)`; the explicit Pauli realization `toMatrix`; Hermiticity; determinant identity `det_toMatrix`; additivity and real homogeneity; injectivity; unique representation of every Hermitian `2 × 2` complex matrix; `finrank_Herm2 = 4`; and the constructive diagonal inertia certificate `(+---)` with a positive scalar line and negative-definite three-space. |
| `Geometry/CelestialNullCone.lean` | The subtype of nonzero future null vectors, positive-scaling setoid, actual quotient `FutureNullRay`, normalization to the unit two-sphere, inverse representative `(1,n)`, and the exact equivalence `futureNullRayEquivCelestial`. |
| `Geometry/ObserverFrameHyperboloid.lean` | The future unit-timelike hyperboloid `FrameHyperboloid`, its standard point, the identity `t² = 1 + |x|²`, and `t ≥ 1`. |
| `Geometry/ObserverRestSpace.lean` | The Lorentz-orthogonal kernel at every frame; `finrank_restSpace = 3`; the sign-reversed Lorentz form; its additivity, real homogeneity, symmetry, strict positivity away from zero, and zero-norm characterization. |
| `Geometry/EinsteinTensorBridge.lean` | The explicit linear equivalence `pauliEinsteinChart : Herm2 ≃ₗ[ℝ] OPH.EinsteinBranch.V 3`; the exact sign identity `lorentzQ v = -quadOf (eta 3) (pauliEinsteinChart v)`; null-cone equivalence; and future-null transport with unchanged time orientation. |

The celestial result is a set-level equivalence. No topology or smooth
structure is bundled. The future-unit-timelike subtype is the algebraic
hyperboloid model underlying `H³`, and the orthogonal kernel is called a rest
space without asserting a formal manifold tangent-bundle identification.

The inertia result is deliberately a constructive diagonal certificate rather
than an invocation of a general-purpose signature API: in the canonical
scalar/traceless splitting, determinant is exactly one positive square minus
three squares, and the corresponding one- and three-dimensional summands are
proved positive and negative definite.

## Verification

The targeted Lean checks are:

```text
lake env lean Geometry/CanonicalLorentzModule.lean
lake env lean Geometry/CelestialNullCone.lean
lake env lean Geometry/ObserverFrameHyperboloid.lean
lake env lean Geometry/ObserverRestSpace.lean
lake env lean Geometry/EinsteinTensorBridge.lean
```

All five pass. The embedded `#print axioms` audit reports only Lean/mathlib's
standard `propext`, `Classical.choice`, and `Quot.sound`.  There is no
`sorryAx`, native-decide axiom, or new project axiom in this packet.

## Boundary and remaining gaps

This closes the canonical *internal* Lorentz-geometry target only.  It does
not prove any of the following:

- that observer-patch data select this Hermitian module or its determinant;
- that `FrameHyperboloid` is the frame bundle of physical spacetime;
- a soldering map, coincidence relation, affine translation sector, or
  Lorentz overlap cocycle;
- a physical clock, rods, causal dynamics, or operational calibration;
- a continuum/refinement limit or Einstein dynamics.

The finite algebraic portions of coincidence descent, affine/Lorentz overlap
transport, induced celestial/frame/rest covariance, and the candidate
rank-three source-frame bridge are now supplied by the separate bounded C2
contract in `C2_EVENT_FRAME_SOLDERING.md`. Source realization, population,
topology, physical cone attachment, causality, refinement, rods, clocks, and
physical spacetime remain downstream. None is imported into the epistemic
status of C1.
