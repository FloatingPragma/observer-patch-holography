# Established mathematics for observer agreement

This guide governs theorem reuse in the
[executable model](../model/candidate_laws/contrast_memory_spherical_v1/README.md) and its
[specification](../model/candidate_laws/contrast_memory_spherical_v1/SPEC.md). The target is one physical
process with compatible carrier, response, record and clock readouts.
Established results about distributed recovery and sheaf diffusion are tools
for that construction. They do not select its physical law.

The [working terms](INTERPRETATION_AND_TERMS.md) fix the distinctions between
the proposed ontology, a theorem's model, its states and its uses of time.

For each proposed result, identify the nearest theorem below, map its objects
and hypotheses to the actual candidate, and isolate the remaining obligation.
Reuse the conclusion once those hypotheses hold. A renamed consensus result
or a longer simulation of it is not an independent theoretical contribution.

The executable model requires source-specific mappings for these theorems.
Its contrast-memory law is a candidate native repair law; its seam memory is
private working memory of the two port observers of a seam, and its reader
admits a record only when both observers' copies agree. A theorem about
distributed recovery applies to it once its objects and hypotheses are mapped
to these seams, memories and records.

## One structure: complete history and internal computation

OPH proposes two descriptions of the same structure: a timeless solution to
self-reference and a world that computes itself through its internal local
dependencies. States, records, read/write relations and responses to admissible
interventions constitute its operation. The complete history can contain
change without an external host or a further time in which the whole history
is played. This is the closure and intrinsic-realization hypothesis of the
[basis](BASIS.md#the-closure-hypothesis).

An application must construct the correspondence between those descriptions,
preserving the observer's accessible records and interventions. A fixed-point
equation alone does not supply that implementation. The
[consensus paper's Introduction](../paper/reality_as_consensus_protocol.tex)
distinguishes a solver's repair schedule from the observer-facing record
history it selects. Internal computation can have local event order without
a global clock; a solver index or an authenticated event order supplies no
physical duration by itself.

Observer-local time is a further reading within that same structure. On the
stated geometric modular-flow branch, the observer's algebra-state pair
supplies dimensionless modular ordering; physical time additionally needs
a readable transition, event correspondence and calibrated clock instrument.
See the [synthesis paper](../paper/observers_are_all_you_need.tex) and the
[observer-clock certificate](../paper/tex_fragments/OBSERVERS_TECHNICAL_SUPPLEMENT_PORT.tex).
Identifying this ordering with experienced time is an interpretation.
Clock compatibility requires local readings tied to the declared observable
evolution and consistent comparisons at shared events. It does not require
identical clock rates, preferred global time or universal simultaneity.
Any shared time parameter used by a restricted benchmark remains a declared
model choice; an auxiliary simulation index gains no physical status from it.
The dual description alone establishes no final instantaneous cosmic endpoint.

## 1. Recovery is a property of a system and its execution assumptions

[Dijkstra (1974), *Self-stabilizing systems in spite of distributed control*](https://www.cs.utexas.edu/~EWD/transcriptions/EWD04xx/EWD426.html)
provides local finite-state constructions that enter a legitimate set from
arbitrary initial states and remain legitimate. His central daemon selects
one enabled move at a time; the guarantee covers every such selection.
A distributed daemon is explicitly outside that paper's scope. Legitimate
states have enabled moves, so legitimate operation continues after recovery.
This is not a theorem that arbitrary local rules stabilize, or that every
legitimate state is a fixed point.

**OPH application contract.** Before invoking recovery language, supply:

- The complete state and an independently specified legitimate set `L`.
  Defining `L` as whatever the chosen solver eventually reaches is insufficient.
- The initial domain, permitted perturbations and conserved sectors. Identify
  a fault span closed under the modeled fault and protocol moves. If faults
  can alter supposedly protected records, a proof that assumes those records
  are intact covers a narrower fault model.
- The enabled moves and execution class: central or concurrent scheduling,
  atomicity, stale reads, delivery, crashes, random choices and any fairness.
- **Closure:** each allowed fault-free move from `L` stays in `L`.
- **Recovery:** every maximal fault-free execution from the admitted fault span
  reaches `L`, with deterministic, almost-sure or expected guarantees stated
  separately. State the unit of a recovery bound.

This is the programme's contract for recovery after a finite disturbance.
Continued forcing requires a tracking or noise-tube result. Recovery does not
promise safety throughout the recovery interval or restoration of information
that the allowed disturbance erased. Fair scheduling alone gives no wall-clock
bound; a time interpretation needs additional timing premises.

**Immediate use:** apply this contract to the executable model's loads, its
private seam memories and its agreed seam records separately. Keep ensemble
convergence distinct from convergence of individual recorded trajectories; do
not assert Dijkstra's ring theorem for this different transition system.

## 2. Termination, confluence and stabilization answer different questions

[Newman (1942)](https://doi.org/10.2307/1968867) supplies the established
termination-plus-local-confluence route to confluence. The
[finite consensus paper](../paper/reality_as_consensus_protocol.tex),
`thm:quotient-repair-normal-form`, uses that route on its accepted quotient
relation, with separate completeness and protected-boundary premises.

Continue to verify those concrete premises rather than reproving the abstract
lemma. A unique terminal state from one initial quotient does not imply that
different histories share one terminal state. That stronger statement needs
the boundary-fibre uniqueness condition. Conversely, a legitimate operating
set may contain cycles and need no terminal normal form. Preserve the
distinction when specifying an observer's continued record-producing activity.

The existing noisy fair-block branch also has its own contraction, stopping-time
and occupation assumptions. It is not a corollary of the word “asynchronous.”

## 3. Use cellular sheaves for a declared linear overlap model

For a finite graph, let each vertex carry a finite-dimensional inner-product
space `F(v)`, each edge a comparison space `F(e)`, and each incidence a fixed
linear restriction `R_ve`. With an orientation on each edge `e=(u,v)`, set

```text
(delta x)_e = R_ve x_v - R_ue x_u,
L_F = delta* delta,
E(x) = ||delta x||^2 / 2.
```

[Hansen and Ghrist, *Opinion Dynamics on Discourse Sheaves*](https://arxiv.org/html/2005.12798),
§§2 and 4, gives `ker L_F = H^0(G;F) = ker delta` and convergence of
`dx/dt = -alpha L_F x`, `alpha>0`, to the orthogonal projection of the initial
state onto that space. Agreement compares restricted expressions, so it need
not make distinct vertex states identical. The theorem assumes the fixed
linear sheaf and this evolution equation; it does not establish arbitrary
nonlinear, asynchronous, adaptive or quantum-channel convergence.

**Immediate use:** factor the existing weighted scalar disagreement into its
actual restriction maps and adjoints. Record the mass/inner-product convention,
nullspace, protected coordinates and boundary data before comparing spectra.
If a native quantum candidate admits a suitable operator-space description,
prove that identification and its positivity, invariant-sector and dynamical
conditions first. Calling density matrices “stalks” does not discharge them.
Keep conservative Hamiltonian evolution distinct from dissipative heat flow.

Changing the restrictions can reduce the reported discrepancy by changing what
is compared. Therefore an adaptive proposal must preserve a declared observation
map or expose and test its changed meaning. A zero residual under a modified
readout cannot validate the original physical interpretation.

## 4. Boundary uniqueness does not imply exact agreement

[Hansen and Ghrist, *Toward a Spectral Theory of Cellular Sheaves* (2019)](https://arxiv.org/html/1808.01513),
§4.1, establishes harmonic extension results. For 0-cochains with values pinned
on vertices `B`, an extension minimizing `||delta x||^2` always exists; it is
unique exactly when `H^0(G,B;F)=0`. With interior vertices `I`, it solves

```text
L_II x_I = -L_IB x_B.
```

This interior harmonicity is weaker than `delta x=0`. Boundary data must lie
in the image of `H^0(G;F) -> C^0(B;F)` to admit a global section. When the
boundary data are a section on a boundary subcomplex, the connecting map into
relative first cohomology expresses the extension obstruction. For the
programme's graph computations, explicitly testing the linear image is enough.

For example, an ordinary three-vertex path pinned to 0 and 1 at its endpoints
has the unique harmonic interior value `1/2` and positive disagreement.
Do not label it complete agreement. Compute both boundary uniqueness and the
remaining residual. On homogeneous unconstrained linear sheaves, the zero
section always exists; nonzero `H^1` is not a blanket assertion that no global
section exists.

**Immediate use:** before adding a boundary selector, test whether the current
protected data extend and whether a nonzero global section vanishes on the
boundary. This separates an inconsistent target from a free interior mode or
a failing solver. It also tests the hypothesis needed by OPH's existing
boundary-conditioned uniqueness result.

The same source, §4.2, also limits boundary reduction: eliminating an internal
vector-valued node can induce `z=x+y`, which no pairwise graph on the unchanged
boundary stalks can enforce. Every pair of these boundary coordinates is freely
variable. Retain the mediator or justify a richer representation before claiming
an exact reduction. The favorable scalar case requires every original vertex
stalk to have dimension at most one, including eliminated vertices.

## 5. Keep relational, probabilistic and quantum extension separate

[Abramsky and Carù (2019)](https://arxiv.org/html/1911.03521), Eq. (4),
Theorem 4.1 and Propositions 5.1/6.2, formalize local agreement and global
extension. For relations, a nonempty natural join supplies a satisfying
assignment, while an exact global relation must recover every local relation
under projection. Probabilistic extension additionally requires the specified
weights. Full local support can hide a Bell obstruction. The natural-join
criterion for relations does not justify multiplying probability marginals
and calling the product their extension.

**Immediate use:** name the valuation type of each candidate's agreement
target before choosing a solver. Check probability compatibility separately
from support compatibility. A positive joint law for all counterfactual
measurement outcomes is not the same object as a global quantum density
operator with context-dependent measurements. A classical extension failure
rejects that classical interpretation; it does not reject every quantum
realization. A2 itself
[does not imply global extension](AXIOM_REFERENCE.md#a2-observer-agreement).

### Finite-lattice methods for record constraints

[Ghrist and Riess, *Cellular Sheaves of Lattices and the Tarski Laplacian*](https://arxiv.org/html/2007.04099),
§§3.2–3.4, identifies sections with `Fix(id ∧ L)` and gives descending-chain
termination conditions. A concrete reusable case has a fixed finite graph,
finite lattice stalks, and join-preserving incidence maps `f_ie` with Galois
right adjoints `f_ie^+`. The adjoints here are order-theoretic residuals.

[Riess and Ghrist, *Diffusion of Information on Networked Lattices by Gossip*](https://arxiv.org/html/2204.00167),
§§III–IV, supplies the corresponding update for broadcasting vertices `tau_t`:

```text
x_i[t+1] = x_i[t] ∧ meet_{j in N(i) ∩ tau_t} f_ie^+(f_je(x_j[t])),
e = {i,j}; the empty meet is the top lattice element.
```

The stated model requires each vertex to broadcast infinitely often, receipt
and computation within the step, and fault-free simultaneous broadcasts.
Finite descending state space and liveness give finite settling on each such
execution. Arbitrarily long firing gaps preclude a uniform clock bound.
This theorem does not cover arbitrary stale messages or crashes.

Each update is monotone and fixes every section. Consequently the settled
state is the greatest section below the initial state: any such section stays
below every iterate. This proves independence from the live firing order for
that initial state, without making different initial states equivalent.

**Immediate use:** test a finite-record encoding as a constraint-propagation
baseline, declaring its order, maps and admissible information loss. With
powersets and direct-image overlap maps, a section is a family of supports with
matching projections. For binary variables, the supports of `x=y`, `y=z` and
`x≠z` all project to `{0,1}` on their overlaps, yet admit no joint assignment.
Empty supports can also satisfy the section equations. Require a separate
nonvacuity and global-extension test; neither probability weights nor quantum
realizability follow from this support computation.

## 6. Apply distributed impossibility results only to matching models

[Fischer, Lynch and Paterson (1985), Theorem 1](https://groups.csail.mit.edu/tds/papers/Lynch/jacm85.pdf)
rules out guaranteed termination of deterministic nontrivial binary consensus
in its fully asynchronous message-passing model with a possible crash, even
with eventual reliable delivery to nonfaulty processes. It is a relevant
obstruction if an OPH implementation makes those commitments. It is not a
no-go theorem for fault-free finite rewriting, asymptotic diffusion, or a
protocol with genuinely supplied additional synchrony assumptions.

Before invoking this result, map processes, decisions, validity, communication,
failures and admissible schedules. Before claiming to evade it, identify the
changed hypothesis and its implementation cost. The
[BFT appendix](../paper/appendix_B_bft_qecc_extensions.tex) separates
protocol safety from the partial synchrony needed by its liveness branch.
No network theorem supplies a physical clock merely by naming an event a tick.

## 7. Work products and scientific contribution

The following bounded outputs use the existing issue ownership. They add no
new issue queue or mechanism and retain the existing candidate evidence.
Reuse a completed certificate when the mapping is unchanged; calculate only
what the actual candidate's missing premise or changed interface requires.

| Work | Required output | Existing owner |
| --- | --- | --- |
| Specify native recovery | One closure/recovery contract with fault span, scheduler and protected records; distinguish deterministic and probabilistic claims. | [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076), [#1026](https://github.com/FloatingPragma/observer-patch-holography/issues/1026) |
| Reuse the linear theory | Actual stalks/restrictions, inner products, `delta`, nullspace, boundary extension residual and uniqueness test for the scalar comparison. | [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076) |
| Preserve one implementation | Accepted operation grammar and same-source preparation, evolution and record history; unchanged reference checks plus the relevant theorem premises. | [#1077](https://github.com/FloatingPragma/observer-patch-holography/issues/1077) |
| Discriminate physical laws | Two viable laws with the same preparation, effects, boundary and clock; a predicted transient or response difference with declared calibration freedom. | [#1025](https://github.com/FloatingPragma/observer-patch-holography/issues/1025), [#1076](https://github.com/FloatingPragma/observer-patch-holography/issues/1076) |
| Review the inference | Nearest primary theorem, premise mapping, unchanged conclusion, candidate-specific remaining work, and scope of any obstruction. | [#1033](https://github.com/FloatingPragma/observer-patch-holography/issues/1033), [#1078](https://github.com/FloatingPragma/observer-patch-holography/issues/1078) |

Record each reusable result as **applicable with checked premises**,
**conditional on named missing premises**, or **not applicable to this model**.
These are research-workflow labels, not changes to the scientific claim registry.
Keep links to the supplied premise evidence; a citation alone is not that evidence.

Recovery, linear diffusion, harmonic extension, confluence and the distinction
between compatibility and extension are established mathematics. An OPH-specific
contribution must identify the common physical realization or a necessary
additional result absent from the cited theorem, and ultimately supply a
discriminating observable relation. Shared equilibria do not choose transition
rates, instruments, reference states or physical clocks. The existing geometry,
gauge, continuum and empirical obligations are listed in the
[premise register](registers/PREMISE_REGISTER_V3.md).

The sources above have different conventions and scopes. In particular, the
arXiv HTML renderings can renumber theorem counters relative to the PDFs;
the linked sections and displayed equations identify the intended results.
