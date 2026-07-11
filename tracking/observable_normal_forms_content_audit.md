# Observation-Determined Normal Forms: current content and compatibility audit

Date: 2026-07-11

Status: audit of the current repo-native manuscript and standalone Lean artifact

## Scope boundary

This audit covers only the mathematics of finite and compact constraint/rewrite systems:

- quotient state spaces and consistent subsets;
- observable fibers and partial normalizers;
- computational terminals versus semantic consistency;
- residual/error-bound and inverse-observation moduli;
- receipt-conditioned endpoint comparison;
- approximate naturality, refinement towers, and inverse limits;
- local-repair and selection obstructions;
- ranked functional systems, linear observability, and formal proof receipts.

The physical \(\chi_\nu\) hypothesis is explicitly out of scope. This report does not assess any coherent-matter source claim, collar coefficient, force law, device, experimental protocol, G9/G10 calibration, or empirical prediction. References to the \(\chi_\nu\) corpus below are limited to its finite quotient, repair, boundary-fiber, normalizer, and proof-receipt interfaces.

## Audited current state

| Surface | State audited |
|---|---|
| Repo-native manuscript | `extra/observable_normal_forms.tex`, SHA-256 `ce9d9dda9febec0ee9cfedc6797724b01767f5c198be94720ff7a9e8a3371817` |
| Bibliography | `extra/observable_normal_forms.bib`, SHA-256 `a183a981f984a2c1674bca7fddbbda356cea05dca467ee0102ae853b3a1e2ad3` |
| Current PDF | `extra/observable_normal_forms.pdf`, SHA-256 `1bf54d9a189f1e876ae67737a4d943c4ecb386cdf1c33e08efbd10617de888ed`; rebuilt from the audited TeX/BibTeX sources with visible global release `r1527` |
| Standalone Lean artifact | `Lean/ObserverPatchHolography/Proofs/ObservableNormalForms/`, Lean and Mathlib `v4.29.1`, pinned by `lean-toolchain`, `lakefile.lean`, and `lake-manifest.json` |
| Independent Lean build in this audit | standalone build succeeded through target 8259; the merged parent build succeeded through target 8265; all new-paper and bridge theorem dependencies were subsets of `propext`, `Classical.choice`, and `Quot.sound` |
| Intended `chi_nu_test` comparison repo | Clean `main...origin/main` at `0f9e43b36386ad15e94947751500bf32ee9ccc58`, proof-chain v10 |
| Independent comparison build | `chi_nu_test/proof_chain/formal: lake build` succeeded through target 8287; linter warnings only |
| Mathematical \(\chi_\nu\) note | `extra/chi_nu_susceptibility_bounds.tex` |

The paper, dedicated artifact, merged Lean hierarchy, and consistency updates
are tracked on the current `main` line. This audit records the released `r1526`
state and its final bridge/hash verification.

## Executive verdict

The repo-native paper has implemented the material corrections identified in the original bundle audit. Its mathematical architecture is now coherent, and no counterexample was found to the revised main theorems. In particular:

- the manuscript is application-neutral and contains no OPH reference;
- computational terminal states are correctly separated from semantic consistency;
- refinement estimates operate on explicit realizable/reconstructing domains;
- same-level implementation agreement has no spurious refinement tail;
- the genuinely projective result now compares different levels through a common anchor;
- the cofinal result has the required compatibility, completeness, nonexpansiveness, summability, and vanishing-receipt hypotheses;
- the edge-case nonemptiness assumptions, stochastic endpoint scope,
  synchronous-depth scope, and zero-space singular-value convention are
  explicit;
- cross-source endpoint uniqueness modulo a silent equivalence is characterized
  exactly, rather than only given as a sufficient singleton-fiber condition;
- finite observation-fiber resampling is separated from representative
  selection and proved to be the weighted conditional-expectation projector.

The finite mathematical architecture is compatible with `chi_nu_test` T6/T27/T31/T32 once equality is read on the declared quotient and raw decode normal forms are not conflated with consistent states. The revised audited-terminal theorem now captures the formerly missing Route-A case directly.

The technical draft gate is now clear: the manuscript and bibliography compile,
the PDF carries the current title and authors, the formalization prose matches
the standalone artifact, and the Lean build and axiom receipt reproduce. The
remaining submission decisions are editorial rather than hidden proof claims:

1. affiliations, addresses, corresponding-author information, and ORCIDs need
   author-approved values;
2. important paper theorems remain only partially or not formalized, especially
   the cofinal/inverse-limit results, intrinsic moduli, sharpness, compact
   uniformity, linear certificates, and complexity results;
3. the candidate-new cross-level/refinement package still needs external
   specialist priority review.
4. venue-specific submission metadata and external specialist review remain
   editorial gates; the repository and public paper surfaces are synchronized
   at release `r1526`.

## Completed manuscript corrections

### 1. Authors, title, metadata, and neutrality

Completed:

- The title is now *Observation-Determined Normal Forms: Stability, Obstructions, and Refinement in Constraint and Rewrite Systems*.
- Bernhard Mueller, David Matscheko, and Jonathan Hill are named as authors in the TeX source.
- PDF metadata in the TeX source names the same authors and current title.
- The introduction explicitly distinguishes “observation-determined normal form” from the established control-theory phrase “observable normal form.”
- Searches of the current `.tex` and `.bib` find no `OPH`, `Observer-Patch`, or observer-overlap branding.
- The former OPH application section and OPH bibliography entries are gone.

Remaining:

- Add affiliations, postal addresses if required by the venue, email addresses, a corresponding-author designation, and ORCIDs.
- Obtain author approval for those metadata rather than inferring them from other submissions.

### 2. Constructive/classical status of the partial normalizer

Completed in the `Proof-carrying interpretation` remark:

- finite decidable data are explicitly separated from the ordinary classical construction;
- injectivity alone is no longer described as an algorithm for deciding fiber inhabitation;
- succinct realizability complexity is distinguished from constructive decidability.

The universal property remains correct in the paper's classical setting.

### 3. Audited-terminal theorem \(T\supseteq C\)

Completed as `thm:audited-terminal`.

The paper now defines:

- \(T=\operatorname{NF}(\to)\), the computational terminal set;
- \(C\subseteq T\), the semantically consistent terminals.

Under observation preservation, weak normalization of every state, and injectivity of \(B|_T\), it proves a unique reachable terminal \(t_x\), with

\[
t_x\in C
\quad\Longleftrightarrow\quad
C_{B(x)}\ne\varnothing.
\]

The proof is correct. It also states the proper limitation: weak normalization gives a unique reachable terminal but does not assert that every execution is strongly normalizing.

This resolves the main compatibility gap with `chi_nu_test` Route A, where every input has a unique quiescent endpoint but an endpoint is consistent exactly when its protected fiber is realizable.

Lean status: fully formalized as `existsUnique_auditedTerminal`, with no nonstandard axiom dependency.

### 4. Precise admissible domains for metric refinement

Completed at the start of `sec:refinement`.

The manuscript defines

\[
Q^{\mathrm{adm}}
=
\{q\in Q:\text{there exists a unique }c\in C\text{ with }B(c)=B(q)\}
\]

and states that all normalizers and coarse maps in the metric tower act on chosen admissible subspaces preserved by the restrictions. Obstructed and ambiguous fibers remain in the exact layer and are not assigned an undeclared finite metric normalizer.

The inverse-limit theorem also says explicitly that its \(Q_\infty\) is the inverse limit of these admissible domains, not automatically the full ambient inverse limit.

This fixes the earlier scope discontinuity between the partial exact normalizer and the total metric tower.

### 5. Same-level agreement versus genuine cross-level comparison

Completed correctly.

`prop:same-level-implementations` now proves, for two outputs at the same fine level and on the same input,

\[
d_n(\rho_{m,n}\widehat x_m,\rho_{m,n}\widehat y_m)
\le
K_{m,n}(e_m+e'_m),
\]

and explicitly notes that no refinement-defect term occurs. This is the sharp architecture identified in the audit: both outputs compare to the same exact fine normal form.

`thm:anchored-cross-level` is now genuinely projective. It compares outputs at arbitrary levels \(m\) and \(\ell\) after restriction to a common level \(n\), through an anchor \(h\), and includes:

- both solver receipt errors;
- both pathwise exact-normalizer defect segments;
- inverse-observation separation at the common anchor.

The proof's five triangle segments are correct, and minimizing over admissible anchors is legitimate.

`cor:nested-cross-level` correctly removes the first path segment and anchor mismatch when the inputs are nested/observation-compatible.

`cor:cofinal-projective` now has load-bearing hypotheses:

- a compatible admissible input in the inverse limit;
- complete coordinate spaces;
- nonexpansive restriction maps;
- summable pathwise normalizer defects;
- vanishing solver receipt errors.

Its Cauchy, implementation-independence, cofinal-subsequence, inverse-limit compatibility, and distance-to-\(N_nq_n\) conclusions follow correctly.

Lean status:

- same-level agreement: fully formalized;
- anchored and nested metric cores: formalized, while the full stage-indexed manuscript wrappers remain partial;
- cofinal Cauchy/completeness and inverse-limit conclusions: not yet formalized.

### 6. Nonempty hypotheses and empty-set countermodels

Completed:

- the finite product calculus requires a nonempty finite index family;
- the collar-section criterion requires a nonempty write space;
- the finite no-repair certificate repeats the nonempty-write qualification;
- the real-valued robust margin requires \(R\ne\varnothing\);
- the text explains the vacuous empty-domain repair and the convention issue for ordinary real-valued distance to the empty relation;
- ranked functional systems require nonempty boundary layer and nonempty finite site state spaces;
- tower spaces and uniform families carry the needed nonemptiness/boundedness qualifications.

Lean formalizes both load-bearing countermodels:

- `empty_write_space_counterexample`;
- `empty_relation_repairMargin_zero`.

These changes repair genuine edge-condition holes rather than merely editorial omissions.

### 7. Finite Markov drift receipt

Completed as `prop:markov-receipt`.

The proposition correctly iterates the affine drift inequality

\[
KV\le \kappa V+\xi,
\qquad 0\le\kappa<1,
\]

and applies Markov's inequality to obtain a one-time endpoint probability bound. It correctly requires positive-support transitions to preserve \(B\) before declaring a \((\delta,0)\) receipt, and it separately states that a schedule kernel must be supported on rewrite steps or stutters.

The caveat is correct and important: this is an endpoint bound, not a pathwise confinement theorem; excursions require another hypothesis.

Lean status: the finite Markov-operator drift iteration and one-time tail inequality are formalized. The complete kernel-law/event/positive-support wrapper remains only partially packaged.

### 8. Synchronous dependency-depth result

Completed as `cor:functional-synchronous`.

The induction is correct: after round \(k\), ranks through \(k\) agree with the generated extension, so depth \(m\) settles the system. The paper correctly limits this to the explicitly declared synchronous semantics and does not infer a compiler or hardware depth bound without a separate correspondence theorem.

Lean status: fully formalized as `RankedSynchronousSystem.synchronous_depth_settling`, with generated-extension uniqueness also formalized.

### 9. Singular-value and zero-space correction

Completed:

- \(\kappa_B=\|(B|_H)^\dagger\|\) is called the reciprocal of the smallest singular value only when \(H\ne\{0\}\);
- the zero-space pseudoinverse convention \(\kappa_B=0\) is explicit;
- the dynamic corollary explicitly assumes \(T(G)\subseteq G\) and \(S(G)=0\);
- the time-tube singular-value statement is restricted to \(G^\perp\ne\{0\}\);
- the static finite-field Rule-90 result is no longer incorrectly described as an instance of the real/complex time-tube singular-value corollary.

The linear stability proof remains correct.

Lean status: the linear quotient, rank, and singular-value results are not yet formalized in the standalone artifact.

### 10. Related-work and novelty discipline

Completed:

- the paper now cites and distinguishes approximate inverse-system/almost-commuting-diagram work;
- it states that the ambient inverse system is exact while normalizer/receipt compatibility is controlled metrically;
- it explicitly disclaims invention of approximate inverse systems;
- it cites formal rewriting libraries and distinguishes its artifact scope;
- the novelty statement is narrowed to an integrated certificate architecture, especially the anchored cross-level comparison and cofinal consequence.

Remaining:

- External review is still required from at least one rewriting/semantics specialist and one inverse-systems, numerical multilevel, or variational-analysis specialist.
- The anchor theorem is a transparent triangle/telescope construction. The paper should retain “candidate-new architecture/application” language unless a specialist search establishes a closer priority boundary.
- Approximate inverse systems are now covered, but adjacent multilevel numerical-analysis, approximate commuting-diagram, and robust observability literature should still be checked before submission.

### 11. Observer-confluence decomposition and finite projection

Completed:

- `thm:cross-source-modulo` proves that, under observation preservation
  and \(\operatorname{NF}=C\), identification of consistent states modulo an
  arbitrary silent relation is equivalent to agreement of normal endpoints
  reached from all equally observed sources;
- the silent-quotient corollary states the corresponding injectivity criterion
  on \(C/{\sim}\);
- the prose separates same-input confluence, cross-source endpoint
  identification, endpoint existence, and all-schedule liveness;
- the aligned two-bit system jointly witnesses preserved observation,
  completeness, normalization, same-source confluence, and cross-source
  uniqueness;
- `thm:fiber-conditional-expectation` proves finite full-support
  observation-fiber averaging is idempotent, weighted self-adjoint,
  \(L^2\)-contractive, and fixes exactly observation-measurable functions;
- Lean formalizes this finite algebraic characterization and the recognition
  theorem, but does not wrap the operator as Mathlib's measure-theoretic
  `condexp` onto \(\sigma(B)\); the proof index therefore marks the manuscript
  theorem as partially formalized;
- its R1--R3 matrix receipt must be checked against an independently extracted
  transition matrix, so it is not a projector-compared-to-itself certificate;
- the Markov section now includes stationary and expected-occupation bounds,
  and the ranked section includes an explicit acyclic Boolean-circuit
  compilation corollary.

## Current theorem and Lean coverage audit

The standalone project contains eleven library source files plus its Lake
configuration, zero `sorry`, zero `admit`, zero custom `axiom`
declarations, and 36 theorem-level `#print axioms` commands. Fresh
standalone and parent `lake build` runs succeeded during this audit. Every
printed theorem depended only on a subset of the standard Mathlib axioms
`propext`, `Classical.choice`, and `Quot.sound`; several
results were axiom-free.

### Fully formalized current-paper results

- proof-carrying partial-normalizer equivalence, main `(a) \leftrightarrow (b)` direction;
- reachable-normal-form fiber theorem;
- all three empty/singleton/confluence clauses;
- audited-terminal alternative;
- observation-leak counterexample;
- boundary-identification iff cross-source endpoint uniqueness modulo an
  arbitrary relation;
- weak-normalization endpoint existence and the fine/coarse two-bit
  confluence separation;
- heterogeneous and symmetric two-output estimates;
- approximate schedule endpoint comparison;
- one-step approximate naturality and exact naturality;
- same-level implementation agreement;
- collar-section criterion and no-repair certificate;
- robust margin and its empty-set audit counterexamples;
- the finite weighted resampling kernel, its measurable fixed space,
  idempotence, weighted self-adjointness, Pythagorean energy identity, and
  \(L^2\) contraction;
- the exact R1 fiber-support, R2 equal-row, and R3 weighted-detailed-balance
  matrix-recognition equivalence for conditional resampling;
- synchronous dependency-depth settling and generated-extension uniqueness;
- standalone width-three Rule-90 kernel, image, good/bad readouts, and no-reverse-repair results.

### Formalized mathematical cores with manuscript packaging remaining

- finite Markov drift and tail bound;
- sensor-enrichment certificate;
- arbitrary metric telescoping;
- receipt-to-exact specialization;
- fine-to-coarse solver estimate;
- anchored and nested cross-level metric cores.

### Current paper results not yet formalized

- the image-subtype bijection part of `thm:universal`;
- intrinsic finite moduli and minimality;
- rate transfer and sharpness;
- high-probability schedule wrapper;
- stationary and expected-occupation wrappers;
- compact-space and uniform-family results;
- exact product calculus;
- consistency-model perturbation;
- full stage-indexed tower/summability/modulus wrappers;
- cofinal common-limit and inverse-limit normalizer theorems;
- equivariant-section/stabilizer theorem;
- the silent-quotient packaging of the modulo-identification theorem;
- the acyclic Boolean-circuit compiler wrapper;
- asynchronous functional update count and path-gain estimate;
- linear quotient/rank/singular-value results;
- succinct NP/coNP/\(\Pi_2^P\) classifications.

This is substantial, honest partial formalization, not complete mechanization of the paper.

## Relationship to Jonathan Hill's existing Lean development

The dedicated artifact does **not** supersede
`Lean/ObserverPatchHolography/Source/` as a whole.  The two developments have
different roles:

- `ObserverPatchHolography` supplies the concrete carrier interface
  (`OPHCarrier`, dependent records, exposed overlap data, gauge equivalence),
  the H1--H4 asynchronous-repair architecture, a Newman's-lemma route to
  confluence, and the nontrivial Rule-90 carrier witnesses;
- `ObservableNormalForms` supplies a paper-neutral abstract normal-form layer,
  exact cross-source boundary identification modulo an arbitrary relation,
  quantitative stability/refinement results, finite weighted resampling, and
  an admission-free submission artifact.

The new theorem
`boundaryIdentifiesModulo_iff_observerEndpointUniqueModulo` broadens the
*generic logical packaging* around
`OPH.boundary_fiber_observer_unique`: it proves both directions of the
boundary-identification/cross-source-endpoint equivalence and does not bake in
one carrier or one gauge relation.  More precisely, Jonathan's H1--H3 laws
establish the neutral `CompleteFor` premise, boundary preservation gives
`ObservationPreserving`, and Hfib is the forward boundary-identification
condition; the old endpoint theorem is therefore recovered as the forward
direction of the new equivalence.  The converse is new generic packaging.  It
still does not prove that the concrete declared boundary of the older carrier
identifies its actual consistent quotient.  The old Rule-90 good/bad
information-set results remain the concrete finite witnesses, and the old
same-input confluence machinery remains useful for repair-law questions.

The merged hierarchy now records that correspondence explicitly in
`OPH.boundaryIdentifiesModulo_iff_observerEndpointUniqueModuloLR`, under
`Lean/ObserverPatchHolography/Source/ObserverPatchHolography/Bridges/`.
Its axiom audit contains only standard axioms.  The bridge deliberately leaves
the actual boundary map and its consistent-quotient injectivity as inputs, so
it clarifies rather than hides the remaining live #304 obligation.

The full parent build succeeds, but the older `Primitives.lean` intentionally
retains exactly three admissions: `localRepair`, `Repair`, and
`repair_respects_gauge`.  Its printed boundary-identification, termination,
completeness, confluence-demo, and Rule-90 theorems do not depend on `sorryAx`.
The new artifact neither imports nor conceals those declarations.  Therefore
there is no wholesale supersession: the new work is the paper's
submission-facing generic layer, while the old development remains the
concrete carrier/application layer and retains its separately useful
same-input confluence results and unfinished global repair operator.

## Manuscript/artifact synchronization completed

The current `sec:lean` and `sec:availability` now describe the standalone
artifact rather than the obsolete three-file snapshot. They record:

- the pinned Lean and Mathlib versions and pinned, checksummed Lake manifest;
- the independently successful parent and standalone builds;
- the theorem-level axiom audit and absence of `sorryAx`;
- the exact, partial, and unformalized coverage boundary from
  `Lean/ObserverPatchHolography/Proofs/ObservableNormalForms/PROOF_INDEX.md`;
- the build receipt, submission manifest, and relative hashes.

The manuscript does not claim complete mechanization. Its formalization table
matches the current proof index, including the remaining cofinal,
inverse-limit, linear, and complexity gaps.

## Mathematical compatibility with `chi_nu_test`

### Compatibility matrix

| Current paper surface | `chi_nu_test` mathematical anchor | Compatibility verdict |
|---|---|---|
| Quotient-first state space \(Q=\Sigma/\Gamma\) | `QuotientRepair.lean:80-125`, `QuotientRepairPresentation` | Direct match. Equality in the paper is equality on the quotient, not necessarily literal equality of presentations. |
| Boundary-preserving total normalizer | T6: `globalRepair_mem_CQ` (`:281`), `globalRepair_boundary` (`:286`), `globalRepair_idem` (`:304`) | Direct specialization on the globally realizable Route-B presentation. |
| Exact rewrite completeness \(\operatorname{NF}=C\) | T6: `normalForm_iff_CQ` (`:231`) | Direct specialization where `Hcomp` holds. |
| Schedule independence under confluence/descent | T6: `schedule_independence` (`:340`) | Compatible exact result. The paper's receipt bound is a different, conditional approximate endpoint theorem. |
| Confluence versus boundary determination | `Core/Primitives.lean`: `confluence_of_commute` (`:1253`), `boundary_fiber_observer_unique` (`:1191`), `demoCarrier_not_confluent` (`:1371`) | Strong conceptual match. Path confluence and injectivity of the protected boundary on consistent quotient states are correctly separated. |
| Gauge-singleton fibers | `CarrierBridge.rule90Cylinder_Hfib_tube_gauge` (`:216`), `hfib_strictly_weaker_than_informationSet` (`:382`) | Compatible after quotienting gauge. Raw-record information-set injectivity is stronger and must not be identified with quotient-fiber injectivity. |
| Empty/reconstructing/ambiguous trichotomy | T31: `Rule90Readout.readout_trichotomy` (`:301`), `tubeData_bijective_iff` (`:321`) | Direct conceptual match. T31 is an exact finite-field/cardinality instance of the paper's fiber trichotomy. |
| Unique terminal on unrealizable fibers | T27: `routeA_world_exists_unique` (`:1162`), `routeA_world_consistent_iff` (`:1175`), `routeA_observer_uniqueness` (`:1148`) | Now directly captured by `thm:audited-terminal` with \(T\supseteq C\). It is not an instance of the stronger \(\operatorname{NF}=C\) corollary. |
| Existence of unrealizable Route-A fibers | T31 landing: `exists_unrealizable_tube_iff` (`RouteA.lean:1614`) | Compatible. The paper's obstruction output and admissible-domain restriction preserve rather than hide these fibers. |
| Strong normalization of every decode schedule | T32: `decodeStep_wellFounded` (`:1781`), `routeA_universal_settlement` (`:1821`) | Compatible but logically independent. The paper's audited-terminal theorem assumes only weak normalization; T32 supplies a stronger carrier-specific fact. The approximate endpoint theorem does not reprove T32. |
| Strong local-repair obstruction | `Core/Rule90.rule90_no_frustrationFree_repair` (`:212`); Route-A cylinder no-go | Same obstruction theme, but definitions differ. The paper's strong repair is a one-map section into a local relation; Route A is a multi-step transactional decoder. No equivalence should be claimed without a bridge. |
| Ranked functional reconstruction | `LayeredCarrier.sweep_eq_extend` (`:231`), `hfib_singleton` (`:272`), `reconstruction_of_boundary_preserving_repair` (`:299`) | Direct specialization of the paper's functional-audit architecture. The new audited-terminal theorem handles failed global audits cleanly. |
| Linear quotient/information boundary | `Rule90Cylinder.tube_information_set_iff` (`:309`) and related exact finite-field theorems | Compatible at the quotient-rank level. The paper's real singular-value estimate is an additional conditioning result, not something already supplied by T9/T31. |
| Uniform inverse conditioning across refinement | No corresponding quantitative theorem in the audited `chi_nu_test` chain | A new obligation. Stagewise information-set theorems do not imply uniform inverse moduli. |
| Approximate naturality/cross-level tower | No full metric/Lipschitz counterpart in the audited `chi_nu_test` chain | Compatible proposed extension, not an already closed `chi_nu_test` result. |
| Finite Markov drift receipt | No matching proof-chain theorem | An optional mechanism for generating settling receipts; it does not alter T6/T27/T31/T32. |
| Width-three Rule-90 regression | `Core/Rule90.lean`: good/bad boundary and reverse-repair no-go | Statements match at width three. The paper artifact does not subsume the much larger v10 cylinder/worldline theorem family. |

### T6 versus T27 is now represented correctly

The paper contains two distinct exact interfaces:

1. `cor:empty-singleton` assumes \(\operatorname{NF}=C\). It matches the T6 `QuotientRepair` presentation, where quiescence and consistency coincide.
2. `thm:audited-terminal` assumes \(C\subseteq T=\operatorname{NF}\). It matches T27 Route A, where every fiber has one quiescent terminal but an unrealizable fiber's terminal is not consistent.

This distinction is essential. Applying the first theorem to raw Route-A decode normal forms would be false; applying the second is mathematically faithful.

### T31 and the admissible-domain restriction

T31 proves a sharp readout trichotomy:

- above the exact counting threshold, the readout is surjective but has ghosts;
- at equality, it is bijective;
- below it, it is injective but not surjective and unrealizable readings exist.

The paper's \(Q^{\mathrm{adm}}\) metric domain correctly retains only reconstructing, realizable fibers for total finite-valued normalizers while leaving T31-style unrealizable fibers visible in the exact layer. This is compatible and avoids silently inventing a metric value for \(\bot\).

### `H_fib` is not raw information-set injectivity

In `chi_nu_test`, `H_fib` may conclude `gaugeEquiv` on raw records, whereas an information set concludes literal raw-record equality. `CarrierBridge.hfib_strictly_weaker_than_informationSet` gives an explicit carrier where the former holds and the latter fails.

The current paper is quotient-first, so its \(B|_C\) injectivity corresponds naturally to `H_fib` only after declared gauge redundancy is removed. Do not describe it as equivalent to raw information-set injectivity unless the presentation-to-quotient map is injective on the chosen representatives.

### Mathematical interface in the \(\chi_\nu\) note

Within the permitted mathematical scope, `extra/chi_nu_susceptibility_bounds.tex` uses:

- a fixed-regulator quotient \(Q_r=\Sigma_r/\Gamma_r\);
- a quotient normal-form map \(n_r\) inside the activation \(q\oplus e\);
- a receipt requiring invariance under relabeling and accepted repair-schedule changes.

The current paper is compatible with those interfaces conditionally:

- a total \(n_r\) must act on a declared realizable/reconstructing domain or be independently supplied as a total activation structure;
- any presentation-level activation must descend through the quotient;
- schedule/relabeling invariance requires exact naturality or a quantitative defect receipt, not merely the name “normalizer.”

No statement here assesses the physical hypothesis in which those mathematical interfaces are used.

## Remaining submission and editorial decisions

### 1. Author affiliation block is incomplete

The names are correct, but `amsart` currently has no `\address`, `\email`, `\curraddr`, `\thanks`, corresponding-author marker, or ORCID metadata. Obtain the author-approved details rather than inventing them.

### 2. Formal coverage remains partial

The standalone artifact is strong and clean, but its `PROOF_INDEX.md` correctly records substantial gaps. The abstract and introduction may say “selected results are formalized”; they must not say “the paper is fully machine-checked.”

### 3. Novelty remains provisional

The manuscript now has appropriately cautious language, but a negative literature search is not a priority result. Preserve the theorem-by-theorem classical/supporting/candidate-new boundary and obtain external specialist review.

### 4. Optional mathematical strengthening

Useful but nonblocking additions would be:

- a general finite-field fiber-cardinality theorem, making the T31-style empty/unique/ghost trichotomy a reusable paper result rather than only a Rule-90 example;
- a short proposition relating raw information sets, quotient information sets, and gauge-singleton fibers;
- a partial-normalizer naturality theorem that tracks obstruction values across refinement instead of restricting all metric results to \(Q^{\mathrm{adm}}\).

## Build and PDF verification

The corrected bibliography now contains ordinary LaTeX accent sequences and no
control characters. A fresh isolated Tectonic build succeeded with no undefined
citations, undefined references, TeX errors, or overfull boxes. The resulting
PDF is 28 pages, carries the current title and all three authors in its metadata,
and embeds all fonts. The repository build helper reproduced the same current
paper surface. Release `r1526` rebuilt the full paper bundle, regenerated the
release manifest and Markdown exports, and synchronized the public paper
surface.

## Final recommendation

The mathematical content has moved from “major revision required” to a
coherent, application-neutral technical draft with synchronized source, PDF,
and a substantial standalone Lean artifact. Before external submission:

1. add author-approved affiliations, correspondence details, and identifiers;
2. obtain specialist review of the priority boundary and the cofinal
   refinement package;
3. freeze theorem names and numbering, rerun the Lean receipt and hashes, and
   zip the standalone artifact;
4. keep claims of complete mechanization and established novelty out of the
   paper;
5. retain this tracking report as the external `chi_nu_test`/OPH
   mathematical crosswalk rather than reintroducing OPH into the manuscript.

With those editorial steps, the paper can remain fully OPH-neutral while
serving as a precise mathematical foundation for the quotient,
observable-fiber, repair, refinement, and receipt obligations tracked elsewhere
in the workspace.
