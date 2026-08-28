# OphGap

Lean 4 / Mathlib corroboration of the OPH source-gap brick, folded into
`observer-patch-holography` so it shares this repository's dependency pin and
proof-debt checks.

This proves the gap for ONE frozen graph: the graph pinned by
`domain_freeze_sha256 = a0be6fc64aecf9ca375fd91c57315e8af5e5cf161c99611f4844ba8f452ae7ff`,
with 8,662 carriers and 11,816 seams. It does NOT prove a general theorem about
all graphs, all OPH data, or all source-gap claims. It also does NOT verify the
Python extraction that produced the frozen Lean literals; those literals are
the proof's input.

## Reproduce the check

The project uses the toolchain named in `Lean/lean-toolchain`
(`leanprover/lean4:v4.29.1`) and pins Mathlib to `v4.29.1` in
`Lean/lakefile.lean`. From the repository root, run:

```sh
cd Lean && lake build OphGap
```

`OphGap` is intentionally not a default Lake target. A dedicated Lean CI job
builds it when its sources, pinned Lean inputs, or CI wiring change; unrelated
Lean edits do not re-evaluate the frozen certificate.

`OphGap/SignedGap.lean` proves, with zero `sorry`:

1. `signed_kinetic_identity` — x ⬝ (Dᵀ D x) = Σ_e (σ_e x(tgt e) − x(src e))².
2. `gram_posDef_of_injective` — injective D ⇒ Dᵀ D positive definite.
3. `triangle_posDef` — the all-negative triangle matrix is positive definite (a CONTROL).
4. `kernel_trivial_of_negative_cycle` / `gram_posDef_of_negative_cycle` — connected +
   a negative signed cycle ⇒ ker D = 0 ⇒ Dᵀ D positive definite (signed-balance-to-injectivity).

Negative control: the all-negative square is shown NOT injective, NOT positive definite,
and to admit NO negative closed walk (`square_no_negative_cycle`).

## gapwitness (2026-08-20): the rule applied to OPH's actual graph

* `OphGap/Checker.lean` — a certificate checker `check n es cyc : Bool` (spanning tree rooted at
  vertex 0 given by edge order; negative closed walk given by edge indices) and its SOUNDNESS
  `posDef_of_check : check n es cyc = true → (Dᵀ D).PosDef`, zero `sorry`.
* `OphGap/OphData.lean` — the frozen OPH visible seam complex, 8,662 nodes / 11,816 edges, all
  signs −1, regenerated from `MAIN_CONFIG` (seed 20260751) and matched byte-for-byte to the
  receipt's `domain_freeze_sha256 = a0be6fc6…`; vertices BFS-relabelled from the witness start.
* `OphGap/OphWitness.lean` — `check_true` is decided by the Lean KERNEL (`decide +kernel`, no
  `native_decide`), giving `oph_gram_posDef`, `oph_eigenvalues_pos`, `oph_kernel_trivial` for the
  real graph. Negative control `tampered_refused`: flipping ONE sign on a witness edge makes the
  kernel reject the same certificate.

## gapinvariance (2026-08-20): the theorem now names OPH's ORIGINAL carrier ids

The gapwitness theorems were about a BFS-relabelled copy of the graph; the relabelling was
checked only in Python. That step is now proved in Lean (Ending A, a transport theorem), with
the transport's hypotheses decided by the kernel on the literal data.

* `OphGap/Transport.lean` — `transportCheck n newT oldT nodes origEs es : Bool` (vertex map
  injective on the carrier list via a kernel-checked left inverse; every edge endpoint a carrier;
  edge permutation + inverse carrying each canonicalised original edge onto the certified list,
  checked in one linear pass with `O(log n)` lookup tries) and its soundness
  `origGram_posDef_of_transportCheck : check … = true → transportCheck … = true →
  ((origD nodes origEs _)ᵀ * origD nodes origEs _).PosDef`, where `origD` is indexed by
  `{v // v ∈ nodes}` (OPH's carrier ids) and `Fin origEs.length` (OPH's edges, OPH's order and
  orientation). Vacuity controls on the gapless square: `sqRelabBij_refused`,
  `sqRelabCollapse_refused`, `sqToTri_refused`; positive control `triPerm_accepted`.
* `OphGap/OphOrigData.lean` — OPH's node list and edge list with ORIGINAL ids (pure input, emitted
  by `tools/gen_orig_data.py` from OPH's own `seam_complex` output), plus two lookup-tree hints.
* `OphGap/Slices/*.lean` — the certificate evaluated by the kernel in thirteen windows, one per
  module, import-chained so Lake builds them one at a time (a single `lean` process retains the
  memory of every finished `decide +kernel` and would cross the 16 GB guard); build with
  `LEANBUILD_JOBS=1 leanbuild build`.
* `OphGap/OphOrigWitness.lean` — `transport_true` (assembled from the slices by proved `_split`
  lemmas), `oph_original_gram_posDef`, `oph_original_eigenvalues_pos`, `oph_original_kernel_trivial`.
* `OphGap/OphOrigControls.lean` — Ending B measured (`endingB_direct_check_false`); negative controls
  `tampered_refused`, `collision_refused`, `edge_collision_refused`, `sign_mismatch_refused`,
  `wrong_vertex_set_refused`. The build covers all five refusal controls: malformed or tampered
  certificates are rejected by the checker rather than silently accepted.

Python now only PRODUCES data files; no step of the argument is justified outside Lean.
