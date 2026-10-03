# McKay / golden-field scope audit

Audit date: 2026-10-02
Target: `arithmon/oph-mckay-galois-bridge-0` at K7-Lean preregistration SHA `2553ed170d04db8777f84fb2f9ef2b66b386f7fb`
Comparison baselines captured at audit start: K7 `366f69eb0a5a4fdc53e0ff0337b5c1b671d9b70a`; OPH preregistration SHA `4ae2148a26ce15591adaac78ff408fb2cc32d3a2`. The shared K7 checkout later advanced to `de250c5178f1326aa428e45ad0f953320d4d2554`; that later head was not used.

## Verdict

**`CURRENT_ARITHMON_MCKAY_CLAIM_EXCEEDS_CERTIFIED_CONTENT`** for the claim that the current Arithmon theorem establishes a McKay route from E8 or the binary icosahedral group to the golden ratio.

This is a statement about what the checked-in Arithmon theorem proves. It is not a disproof of the classical McKay correspondence, nor a finding that an independent reconstruction cannot establish the golden character field.

## Recalculation and classification

In `GIFT/Relations/GoldenRatio.lean`, `phi_path_mckay` has exactly these conjuncts:

1. `rank_E8 + dim_G2 + rank_E8 = 30`, which reduces to `8 + 14 + 8 = 30`.
2. `30 = 6 * Weyl_factor`, which reduces to `30 = 6 * 5`.

Both are **ARITHMETIC_IDENTITY** propositions. The theorem body discharges them with `native_decide`. Its docstring's chain `E8 -> Binary Icosahedral -> Icosahedron -> Golden Ratio` is **INTERPRETIVE COMMENT**, not a theorem premise or conclusion. There is no group carrier, representation, character, fusion rule, graph, field, or real-embedding choice in the proposition. The accurate classification is **ARITHMETIC_IDENTITY_WITH_MCKAY_INTERPRETATION**.

Other audited files:

| File | Relevant certified content | Classification for this question |
|---|---|---|
| `GIFT/Foundations/GoldenRatio.lean` | Real definitions of `phi` and `psi`, algebraic identities, and finite Fibonacci identities | **ARITHMETIC_IDENTITY**; no McKay/group-theoretic statement |
| `GIFT/Foundations/E8Mathlib.lean` | Mathlib E8 Coxeter matrix reference, finite root-count arithmetic, and E8 dimension calculations | **GRAPH THEORY / ARITHMETIC IDENTITY**; no affine McKay graph or finite-group bridge |
| `GIFT/Foundations/E8Lattice.lean` | Euclidean-space and E8-lattice definitions with lattice lemmas | **ACTUAL LATTICE MATHEMATICS**, unrelated to a representation ring or this implication |
| `GIFT/Relations/KoideAssembly.lean` | Assembly and bounds for a mass-ratio expression using existing model inputs | Excluded by the runbook; it supplies no premise to this audit |

Repository-wide exact-term search found no K7-Lean construction or theorem for `SL(2,F5)`, `2I`, the binary icosahedral group, its irreducible characters, a representation ring, McKay fusion, or the affine E8 graph. Matches for `PSL(2,7)` are separate claims and do not supply this bridge. The import roots, `Verification/AllImports.lean`, declaration inventory, and generated documentation expose no hidden McKay declaration.

Mathlib in the pinned local dependency does contain the typed carrier `Matrix.SpecialLinearGroup (Fin 2) (ZMod 5)`. That makes a standard group carrier available for a future derivation; it does not establish this group's order, a faithful complex doublet, its character field, or its McKay graph in the current K7-Lean source.

## External comparison and trust boundary

No OPH result was used as a premise, no OPH data were copied into Arithmon, and no Koide or physical data were used. The existing OPH checkout remained at `ed657eb1f98dd0607212baa4b4ca9ed11d181f57`. I fetched the preregistered OPH commit `4ae2148a26ce15591adaac78ff408fb2cc32d3a2` from `upstream` and recorded the target producer path `code/a5_closure/sl2f5_mckay_e8_certificate.py` at blob `2938db085477737cab6eb09a3a87f7330c37f473`. The producer contents and generated receipt were not consumed: there is no independent Arithmon group/representation reconstruction to compare against, so the OPH cross-control was **NOT RUN**.

## Claim adjudication and next route

| Claim | Verdict | Reason / next evidence required |
|---|---|---|
| Current `phi_path_mckay` proves a McKay representation-theoretic route | **HOLD / unsupported by this theorem** | The proposition is only the two recalculated integer equalities above. |
| Current K7-Lean proves `Q(sqrt(5))` is a character field | **NOT ESTABLISHED** | No finite group character or trace-field construction is present. |
| Current K7-Lean proves the nontrivial Galois embedding is unselected by McKay | **NOT ESTABLISHED** | No conjugate representation or recomputed fusion graph is present. |
| A future independent `SL(2,F5)` derivation could establish the stronger result | **OPEN** | Requires the group, faithful irreducible doublet, independently derived irreducibles/fusion, explicit graph isomorphism, trace-field generation, and Galois-conjugate fusion computation. |

The next mathematical route is a separately reviewed workstream implementing those constructions from the typed `SL(2,F5)` carrier. Only after that evidence exists should K7-Lean add a certificate module or state the stronger `MCKAY_DERIVES_GOLDEN_CHARACTER_FIELD__GALOIS_EMBEDDING_UNSELECTED` outcome. This audit deliberately does not infer that stronger conclusion from the OPH comparator.

## Work status

- Lean theorem statements and executable certificates: unchanged.
- Source comment correction: labels the McKay bullet as motivation and scopes the existing golden-ratio module to its actual identities.
- No canonicals, observables, prediction tables, papers, or OPH files changed.
- `lake build`: passed (8,855 jobs); `lake build Verification`: passed (8,859 jobs).
- `scripts/test_proof_inventory.py`: 3 tests passed; `proof_inventory.py --check`: 150 Lean files, 14 axioms, 0 holes, 1,476 `native_decide` occurrences; `update_verification_imports.py --check`: passed.
- `scripts/check_blueprint_sync.sh`: passed (321 declarations, 0 mismatches); `git diff --check`: passed.
- `lake exe checkdecls blueprint/lean_decls`: invoked but not completed because `GIFT.olean` is absent; the repository workflow treats this as a warning (`|| echo`) rather than a blocking gate. Building the default target and `Verification` did not produce that root object.
