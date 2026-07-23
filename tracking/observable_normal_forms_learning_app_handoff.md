# Learning-app handoff: observation-determined normal forms

Date: 2026-07-11

Target: the coding agent maintaining the OPH Textbooks application at
[learn.floatingpragma.io](https://learn.floatingpragma.io/)

Canonical mathematics:
[Observation-Determined Normal Forms: Stability, Obstructions, and Refinement in Constraint and Rewrite Systems](../extra/observable_normal_forms.pdf)

Canonical source:
[observable_normal_forms.tex](../extra/observable_normal_forms.tex)

Formal artifact:
[ObservableNormalForms](../Lean/ObservableNormalForms/)

## Purpose

Update every learning-app explanation of consensus, confluence, normal forms,
observer agreement, and repair so that four logically separate proof
obligations are no longer collapsed:

1. same-source confluence;
2. cross-source identification from protected observations;
3. normalization and all-schedule liveness;
4. constructibility by the declared local repair support.

The standalone paper provides the application-neutral mathematical interface.
It contains no OPH-specific claim. The learning app may explain how the generic
theorems constrain an OPH application, but it must not rewrite the standalone
paper as an OPH derivation.

## Repository boundary

The source of `learn.floatingpragma.io` is not present in the
`oph-meta/WebProjects` checkout and was not found among the visible GitHub
repositories at the time of this handoff. The current deployment is a Lovable
single-page application.

The instructions below therefore use exact visible-copy search anchors rather
than private source paths. Apply them to every source record or component that
emits the matching text. Do not patch only the compiled JavaScript bundle.

The closest local implementation is `WebProjects/oph-lab`. It has a separate
consensus page and is not the source of `learn.floatingpragma.io`.

## Required mathematical model

Introduce the following objects before making claims about observer agreement:

- `Q`: states after quotienting representation changes that are declared
  silent;
- `C \subseteq Q`: semantically consistent states;
- `B : Q \to \mathcal B`: protected observation, boundary, input, sector, or
  record data;
- `\to`: accepted rewrite or repair relation;
- `\approx`: the declared silent or gauge equivalence on consistent states;
- `\operatorname{NF}(\to)`: normal forms of the rewrite relation.

For an observation `b`, the consistent observation fiber is

~~~math
C_b = C \cap B^{-1}(b).
~~~

Teach the fiber trichotomy explicitly:

| Fiber | Meaning |
|---|---|
| `C_b = \varnothing` | The protected observation is unrealizable. |
| `C_b` has one element modulo `\approx` | The observation reconstructs a unique consistent state modulo the declared silence. |
| `C_b` has multiple inequivalent elements | The observation is ambiguous. |

The trichotomy is semantic. A failed search is not automatically a proof that a
fiber is empty.

## The four obligations

### 1. Same-source confluence

Confluence compares repair paths that start from one source:

~~~math
x \to^* y,\qquad x \to^* z
\quad\Longrightarrow\quad
y \downarrow z.
~~~

When termination and the required local-diamond condition hold, this gives a
schedule-independent normal form `\operatorname{nf}(x)` for that fixed source
`x`.

Approved language:

> Different accepted repair orders from one fixed initial quotient state reach
> the same normal form.

Do not write:

> Confluence means that all initial states, observers, or boundary-compatible
> descriptions end at one global state.

The latter is a cross-source claim and needs another theorem.

### 2. Cross-source identification

Under the two explicit premises

1. accepted steps preserve `B`; and
2. `\operatorname{NF}(\to)=C`,

the standalone theorem proves the equivalence

~~~math
\bigl[
  c,d\in C \land B(c)=B(d)
  \Longrightarrow c\approx d
\bigr]
\quad\Longleftrightarrow\quad
\bigl[
  B(x)=B(y)
  \Longrightarrow n_x\approx n_y
\bigr],
~~~

where `n_x` and `n_y` are any reachable normal endpoints of `x` and `y`.

Approved language:

> Endpoints reached from different sources with the same protected observation
> agree modulo the declared silent equivalence exactly when the observation
> identifies the consistent quotient.

The theorem is an equivalence, not a proof that an arbitrary physical boundary
map has the required injectivity. The application must still prove

~~~math
u,v\in C_{\mathrm{application}}
\;\land\;
B_{\mathrm{application}}(u)=B_{\mathrm{application}}(v)
\quad\Longrightarrow\quad
u\approx_{\mathrm{gauge}}v.
~~~

Do not report the generic theorem as an application-specific injectivity
certificate.

### 3. Normalization and liveness

Keep the following levels distinct:

- weak normalization: at least one normal endpoint exists;
- strong normalization: no infinite rewrite sequence exists;
- fairness or another scheduling condition: enabled work is not postponed
  forever;
- all-schedule settlement: every allowed execution reaches a terminal result.

Confluence alone does not supply endpoint existence. Weak normalization plus
endpoint uniqueness does not prove that every possible execution settles.

Approved language:

> Weak normalization supplies at least one endpoint. Strong normalization,
> fairness, or another liveness certificate is needed before every allowed
> schedule may be said to settle.

### 4. Local repairability

Existence and uniqueness of a consistent extension do not show that a
prescribed local write support can construct it.

In the finite collar model, let `X_W` be the writable interior, `X_D` the
protected collar, and `R\subseteq X_W\times X_D` the locally consistent
relation. A total exact collar-preserving repair exists exactly when

~~~math
\pi_D : R \to X_D
~~~

is surjective, assuming the stated nonempty-domain conditions.

Approved language:

> A unique consistent state may still be unreachable by the declared local
> repair support. Repairability needs its own support or collar-projection
> certificate.

## Exact examples for an interactive lesson

### Negative example: confluence without reconstruction

Use two states:

~~~text
Q = C = {a, b}
rewrite relation = empty
B(a) = B(b) = 0
~~~

Every state is already a normal form, so the relation is confluent. The single
observation fiber contains two distinct consistent states. Confluence therefore
does not imply observation-based reconstruction.

Suggested interaction:

- let the learner select source `a` or `b`;
- show that each source has no branching repair path;
- then reveal that both expose the same `B` value but remain different;
- ask why the boundary does not determine the state.

### Positive example: a protected two-bit repair

Use

~~~text
Q = F_2 x F_2
B(u, v) = u
C = {(0, 0), (1, 0)}
repair: (u, 1) -> (u, 0)
~~~

The repair preserves `B`, normal forms are exactly `C`, and every source
`(u,v)` terminates at `(u,0)`. Sources with the same protected bit therefore
have the same endpoint.

Suggested interaction:

- allow the learner to toggle `u` and `v`;
- animate only `v` being repaired;
- visibly lock `u` as the protected observation;
- compare two sources with the same `u` and different `v`;
- then compare sources with different `u` and explain why different endpoints
  are expected.

### Observation-leak example

Use

~~~text
Q = {x, c, d}
C = {c, d}
B(x) = B(c) = 0
B(d) = 1
repair: x -> d
~~~

The consistent fiber over `0` is the singleton `{c}`, but `x` repairs to `d`
because the step changes the protected observation. This makes the
observation-preservation premise visible.

## Current live copy that must be replaced

Search every textbook, chapter variant, summary, glossary, quiz, diagram,
metadata record, and generated page for these anchors:

- `Chapter 4: Reality as a Consensus Protocol`
- `The Fixed-Point Theorem: Objectivity from Confluence`
- `Confluence. A system is confluent if`
- `no matter what choices you make along the way`
- `they all reach the same lake`
- `The theorem says` followed by `every state has a unique normal form`
- `Apply repairs in any order`
- `Unique normal form | Objectivity: physics is the same for all observers`
- `the recording provides a unique resolution`
- `negotiation converges to a unique result`
- `same endpoint regardless of repair order`
- diagram labels `initial states` and `unique normal form`
- chapter summaries that map `Consensus` directly to `Objectivity` without an
  observation-fiber premise.

The current app contains more than one Chapter 4 presentation. Fix every
duplicate, not only the first search hit.

## Required chapter structure

Refactor the Chapter 4 theorem presentation into this order:

1. local repair and the inconsistency potential;
2. termination or the exact liveness premise actually available;
3. local diamond and same-source confluence;
4. source-indexed normal form `\operatorname{nf}(x)`;
5. observation fibers and the empty/singleton/ambiguous trichotomy;
6. cross-source endpoint identification;
7. silent or gauge quotient;
8. local repairability;
9. holonomy or global-extension obstructions;
10. application status and open certificates.

Do not call `\operatorname{nf}(x)` one observer-independent world until both the
same-source theorem and the application-specific observation-identification
premise are displayed.

## Required diagram changes

Any Newman/confluence funnel must show one source branching into several repair
schedules:

~~~text
                 schedule A
              /--------------> nf(x)
source x ----<
              \--------------> nf(x)
                 schedule B
~~~

Do not label several unrelated dots as `initial states` and funnel them to
`nf(s)` merely because confluence holds.

If the app wants to visualize cross-source identification, use a second,
explicitly conditional diagram:

~~~text
x --repair--> n_x
|             |
B(x)=B(y)     n_x ≈ n_y   only if B identifies C modulo ≈
|             |
y --repair--> n_y
~~~

Diagram alt text must preserve the same distinction.

## Suggested UI cards

Render four compact cards with these titles and questions:

| Card | Question |
|---|---|
| Same-source confluence | Do different repair orders from this source agree? |
| Cross-source identification | Do equally observed sources agree modulo silence or gauge? |
| Liveness | Does every allowed schedule actually settle? |
| Local repairability | Can the declared write support construct the extension? |

Each card should expose:

- required hypotheses;
- theorem output;
- what it does not imply;
- evidence status;
- a link to the relevant paper section or formal artifact.

Do not encode these four statuses in one Boolean such as `isConfluent`.

Suggested data shape:

~~~ts
type ProofObligation = {
  id: 'same-source' | 'cross-source' | 'liveness' | 'repairability';
  question: string;
  requires: string[];
  establishes: string;
  doesNotEstablish: string[];
  evidence: 'paper' | 'lean-backed' | 'application-open';
};
~~~

## Robust and refinement layer

Add an advanced callout after the exact finite lesson.

Exact injectivity is not quantitative stability. The paper separates:

- an error-bound modulus `\eta`, carrying residual error to distance from the
  consistent set;
- an inverse-observation modulus `\omega`, carrying observation discrepancy
  between consistent states to state discrepancy.

The learning point is:

~~~text
small residual
  -> near the consistent set       controlled by eta
  -> close consistent states       controlled by omega
~~~

Stagewise injectivity at every finite resolution does not guarantee uniform
stability under refinement. The app must not summarize the refinement theorem
as “finite uniqueness automatically survives the continuum limit.”

The projective result also has explicit premises: admissible compatible input,
complete coordinate spaces, nonexpansive restrictions, summable pathwise
defects, and vanishing solver-receipt errors.

## Conditional-expectation callout

For a finite full-support weighted state space, averaging within observation
fibers is the conditional expectation onto observation-measurable functions.
It is:

- idempotent;
- weighted self-adjoint;
- weighted `L^2`-contractive;
- fixed exactly on functions constant within each observation fiber.

The matrix recognition receipt checks:

- R1: no transition between different observation fibers;
- R2: equal rows for states in the same fiber;
- R3: weighted detailed balance.

State two exclusions prominently:

1. fiber averaging does not select one representative from an ambiguous fiber;
2. being a projector does not by itself prove a spectral gap or convergence
   rate for another repair dynamics.

## Evidence labels

Use these labels consistently:

### Lean-backed finite core

- partial normalizer and exact observable-fiber results;
- reachable normal forms remain in the protected fiber;
- cross-source identification equivalence modulo an arbitrary relation;
- same-source versus cross-source counterexamples;
- observational output estimates;
- selected refinement metric cores;
- finite collar repair criterion and no-repair margin;
- conditional-resampling projector and R1--R3 recognition;
- synchronous ranked-system settling;
- finite Markov drift iteration core;
- width-three Rule 90 regression results.

### Paper proof, not fully machine checked

- intrinsic-modulus minimality and sharpness;
- compact and uniform-family extensions;
- full cofinal and inverse-limit constructions;
- equivariant selector obstruction;
- linear-algebraic certificates;
- succinct complexity classifications;
- full manuscript wrappers around some partially formalized metric cores.

### Application-specific certificate

- the actual physical observation or boundary map identifies the actual
  consistent quotient;
- the actual repair transitions instantiate the abstract matrix receipt;
- the actual scheduler satisfies the liveness premise;
- the actual local write support satisfies the repairability criterion.

Never promote an application-specific item to “Lean-backed” merely because the
generic theorem it would consume is formalized.

## Glossary updates

Add or revise these terms:

- Confluence: joining of paths from one source.
- Normal form: a source-indexed terminal state.
- Observation fiber: consistent states with one protected observation.
- Observation-determined normal form: endpoint determined by observation data
  modulo declared silence.
- Cross-source identification: agreement of endpoints from equally observed
  sources.
- Weak normalization: existence of at least one normal endpoint.
- Strong normalization: absence of infinite rewrite paths.
- Repairability: constructibility by the declared write support.
- Silent equivalence: differences intentionally hidden from the observation.

## Required comprehension checks

At minimum, add these questions:

1. Does confluence imply that two different states with the same boundary have
   the same endpoint?
   - Correct answer: No. That needs identification of consistent states by the
     boundary, plus observation preservation and normal-form completeness.
2. Does a decreasing Lyapunov score prove that every schedule settles?
   - Correct answer: Only with the stated discreteness, well-foundedness, or
     scheduling hypotheses; liveness is a separate obligation.
3. Does a unique consistent extension prove that a one-site repair can build
   it?
   - Correct answer: No. The declared write support needs a repairability
     certificate.
4. What does an empty observation fiber mean?
   - Correct answer: The observation is inconsistent with the declared
     consistency model. It is not merely an algorithm timing out.
5. Does conditional expectation choose one hidden state inside an ambiguous
   fiber?
   - Correct answer: No. It averages fiber-measurable information.

## Link targets

Use stable repository links until the first-party public paper page is
available:

- PDF:
  `https://github.com/FloatingPragma/observer-patch-holography/blob/main/extra/observable_normal_forms.pdf`
- TeX:
  `https://github.com/FloatingPragma/observer-patch-holography/blob/main/extra/observable_normal_forms.tex`
- Lean artifact:
  `https://github.com/FloatingPragma/observer-patch-holography/tree/main/Lean/ObservableNormalForms`
- Application discussion:
  `https://github.com/FloatingPragma/observer-patch-holography/issues/304`

When the first-party paper page is published, prefer that HTML page for
reader-facing links while retaining the repository artifact link for audit.

## Non-goals

Do not:

- introduce OPH terminology or citations into the standalone paper;
- claim that the standalone theorem selects the correct physical boundary map;
- turn same-source confluence into a global uniqueness theorem;
- turn weak normalization into all-schedule liveness;
- turn conditional expectation into representative selection;
- infer a spectral gap from idempotence;
- infer a continuum theory from an inverse-limit normalizer;
- infer a hardware implementation from an abstract repair map;
- use Rule 90 as a selected physical law rather than a regression example;
- import the physical chi-nu hypothesis into this lesson.

## Acceptance checklist

- [ ] Every duplicated Chapter 4 source has been audited.
- [ ] No learner-facing definition says confluence compares arbitrary initial
      states.
- [ ] Every normal form in the confluence lesson is indexed by its source.
- [ ] Cross-source agreement displays the observation-identification premise.
- [ ] Observation preservation and `\operatorname{NF}=C` are visible where the
      equivalence is invoked.
- [ ] Weak normalization, strong normalization, fairness, and all-schedule
      settlement are not conflated.
- [ ] Local repairability is presented as a fourth obligation.
- [ ] Empty, singleton, and ambiguous fibers are all shown.
- [ ] The negative two-state and positive two-bit examples behave correctly.
- [ ] The observation-leak example shows why preservation is load-bearing.
- [ ] Newman/confluence diagrams branch from one source.
- [ ] Cross-source diagrams are explicitly conditional.
- [ ] Glossary, summary tables, quizzes, chapter recaps, diagram alt text, and
      SEO snippets use the corrected language.
- [ ] The standalone paper is listed as a mathematical companion, not as a
      seventh core OPH paper.
- [ ] Lean-backed, paper-only, and application-open evidence labels are
      distinct.
- [ ] Links resolve to the paper, source, and formal artifact.
- [ ] Type checking, linting, unit tests, and the production build pass.
- [ ] The deployed `learn.floatingpragma.io` bundle is inspected after publish,
      not only the source preview.

## One-paragraph replacement summary

Use this when a compact chapter recap is needed:

> Local repair can be schedule-independent without making a protected boundary
> globally reconstructing. Confluence compares repair paths from one fixed
> source. Agreement between different sources with the same observation is a
> separate property: under observation preservation and exact normal-form
> completeness, it holds modulo the declared silent equivalence exactly when
> the observation identifies the consistent quotient. Endpoint existence and
> all-schedule liveness require their own normalization or fairness premises,
> and a unique extension may still be unreachable by the declared local write
> support. These four obligations must remain separate in every theorem,
> diagram, example, and status label.
