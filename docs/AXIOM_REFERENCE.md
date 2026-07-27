# The OPH Axiom Reference

Canonical reader-facing reference for the three-axiom basis of Observer Patch
Holography. The machine-readable normative source is
`claims/axiom_registry.yaml`; papers that state the full basis include
`paper/tex_fragments/THREE_AXIOM_BASIS.tex`. A surface that only consumes an
axiom names the exact dependency and cites this reference.

OPH rests on exactly three core axioms. Everything else in the corpus is a
theorem, an exact result inside a named finite realization, a declared open
interface, an independence result with countermodels, a physical
identification, or a withdrawn claim. The seven-value status vocabulary is
`axiom_forced`, `exact_named_realization`, `discovery_only`,
`conditional_open_interface`, `independence_limited`,
`physical_identification`, and `withdrawn`.

---

## A1. Oriented twelve-port observer screen

**Plain language.** There exists an observer patch net on an oriented
spherical screen. At every finite resolution, each local carrier has twelve
primitive boundary ports. The ports form the vertices of an oriented
triangular boundary with 30 edges and 20 faces: five edges meet each port,
the neighbors of each port form a five-cycle, and every edge belongs to two
oppositely oriented faces, so the boundary is combinatorially the boundary of
an icosahedron. Carriers join through typed seams and coherent triple
overlaps, refine to an oriented spherical support, and expose local state,
readback, records, repair moves, and checkpoints.

**Formal statement.** Let \((R,\preceq)\) be a directed regulator system. For
every \(r\in R\), an observer-screen object exists:

\[
\mathfrak N_r =
(\mathcal P_r,\mathcal A_r,\mathcal R_r,\mathcal I_r,
 \mathcal U_r,\mathcal C_r,N_r,S_r,b_r),
\]

with the following typed components.

1. \(\mathcal P_r\) is a finite patch and overlap category, closed under the
   intersections needed by the theory.
2. \(\mathcal A_r\) is an isotone local algebra net. Disjoint visible regions
   obey the declared locality law.
3. \(\mathcal R_r\) contains accessible records, ports, boundary observables,
   and checkpoints. \(\mathcal I_r\) contains the visible restriction,
   translation, and readback interfaces. \(\mathcal U_r\) contains the
   allowed update and repair instruments.
4. \(\mathcal C_r\) is a federation of finite local carriers with finite
   carrier index set \(\mathrm{Car}_r\) and a typed map
   \(c_r:\mathrm{Car}_r\to\operatorname{Ob}(\mathcal P_r)\). The accessible
   algebra of carrier \(i\) is that of the patch \(c_r(i)\). A local carrier
   \(i\) has the quotient-visible tuple
   \[
   \mathcal C_{r,i}=
   \left(
   \mathcal A_{r,i},
   \mathcal R_{r,i},
   \{(\mathcal I_{r,i,p},\pi_{r,i,p},e_{r,i,p})\}_{p\in P_{r,i}},
   \mathcal M_{r,i},
   \mathcal U_{r,i},
   \operatorname{Chk}_{r,i}
   \right),
   \]
   where the \(\mathcal A_{r,i}\) are finite accessible algebras, the
   \(\mathcal R_{r,i}\subseteq Z(\mathcal A_{r,i})\) are central record
   algebras, and the \(e_{r,i,p}\) are twelve primitive pairwise-orthogonal
   central port projections with \(\sum_p e_{r,i,p}=1\). The
   \(\mathcal I_{r,i,p}\) are interface algebras with visible readout maps
   \(\pi_{r,i,p}:\mathcal A_{r,i}\to\mathcal I_{r,i,p}\); the
   \(\mathcal M_{r,i}\) contain mismatch scores, the \(\mathcal U_{r,i}\)
   bounded local updates and repairs, and \(\operatorname{Chk}_{r,i}\)
   continuation data.
5. Each carrier has a boundary packet
   \[
   K_{r,i}=(P_{r,i},E_{r,i},F_{r,i},o_{r,i}),
   \qquad
   |P_{r,i}|=12,\quad |E_{r,i}|=30,\quad |F_{r,i}|=20,
   \]
   whose edges and faces are incidence relations, not additional ports.
   Every port has degree five and a five-cycle link; every edge occurs in
   exactly two faces with opposite induced orientations under \(o_{r,i}\).
   The resulting oriented simplicial complex is combinatorially isomorphic to
   the boundary complex of a regular icosahedron; "regular" identifies the
   combinatorial reference complex and adds no Euclidean placement. No
   preferred vertex labels, embedding coordinates, edge lengths, gauge
   labels, or particle labels are part of the packet. A single unlabeled
   oriented template \(K^\circ\) represents this isomorphism class, and every
   carrier supplies an orientation-preserving incidence isomorphism
   \(\kappa_{r,i}:K_{r,i}\overset{\sim}{\to}K^\circ\), understood modulo
   orientation-preserving automorphisms.
6. The federation has a seam algebra \(\mathcal B_{r,ij}\) for every declared
   overlap with unital restriction maps
   \(\rho_{r,i\to ij}:\mathcal A_{r,i}\to\mathcal B_{r,ij}\); declared
   nonempty triple overlaps have algebras \(\mathcal B_{r,ijk}\), typed
   restrictions from the pairwise seams, and coherent restriction cocycles.
   The resulting overlap nerve is \(N_r\). A connected subfederation realizes
   complete availability of the readback, record, feedback or repair,
   prediction, and checkpoint interfaces required of an operational observer.
   Availability does not imply repair termination or confluence.
7. \(S_r\) is the declared finite normalized edge-midpoint refinement
   \(\operatorname{Ref}_{n(r)}(K^\circ)\) with its seed map from
   \(K^\circ\), carrying an orientation-preserving geometric homeomorphism
   \(\iota_r:|S_r|\to S^2\) with mesh tending to zero. A source-bound
   simplicial map \(b_r:N_r\to S_r\) connects the federation nerve to this
   support, and the nerve contains a designated oriented two-cycle \(z_r\)
   with \((b_r)_*[z_r]=[S_r]\in H_2(S_r;\mathbb Z)\). This degree-one
   condition prevents a constant or collapsing bridge and forces support
   coverage. For \(r\preceq s\) the tower supplies nerve and support maps
   \(n_{s\to r}\) and \(s_{s\to r}\).
8. Routing, quotienting, incidence, records, checkpoint data, support
   projection, and refinement commute: \(b_r\circ n_{s\to r}
   = s_{s\to r}\circ b_s\), the support realizations obey the controlled
   compatibility estimate
   \(\sup_x d_{S^2}(\iota_s(x),\iota_r(|s_{s\to r}|x))\le\varepsilon_r\)
   with \(\varepsilon_r\to0\), and refinement carries \(z_s\) to \(z_r\).
   The specification states the compatible finite-stage realizations and
   their mesh convergence; the analytic \(S^2\) refinement limit is carried
   by executable receipts on the declared tower rather than asserted as a
   separate axiom clause.
9. An admissible presentation equivalence preserves accessible algebras,
   interface maps, oriented incidence, records, the response interface, the
   repair interface, checkpoint continuation, and refinement. Hidden labels,
   embedding coordinates, materials, and memory layouts do not enter
   protected outputs.

The local carrier boundary, the federation overlap nerve, and the global
\(S^2\) support are distinct objects connected by explicit maps. A1 asserts
the existence and architecture of the observer screen.

**A1 constrains.** The carrier boundary combinatorics, the typed
carrier-federation-support architecture with its explicit maps, the
operational interface surface, and the refinement tower with controlled
spherical support realization.

**A1 does not imply.** Semantic agreement, repair termination, confluence,
global state extension, any state selection, equal state weights
\(\rho(e_p)=1/12\), the integer load fiber and quadratic readback, the
perturb/readback response law, the inverse-port involution and six axes, the
\(A_5\) rotation group as axiom content, any selection between the
Galois-conjugate rank-three Gram frames, the compact current algebra, collar
recovery, generalized-entropy structure, or any laboratory identification.
The listed incidence consequences (inverse port, axes, \(A_5\), frame pair,
antipode polynomial) are theorems about A1 carriers, proved separately; the
response law and current algebra are exact results inside a named
realization; equal state weights and the load fiber are derivation targets.
No gauge, particle, coupling, measured-target, or fitted coordinate may
appear in the source packet used for any of these derivations.

**Architecture firewall.** Every load-bearing A1 field is quotient-visible;
has an operational role in state access, readback, repair, records,
checkpoints, routing, or refinement; is invariant under complete relabeling
of the carrier presentation; is realizable and independently testable at
finite cutoff; was fixed without gauge, particle, mass, coupling, or
measured-target data; is deletion-tested against a clearly identified
observer operation; and keeps same-A1 alternative completions and off-branch
carrier alternatives as separate controls. A field failing any test leaves
A1 and becomes a derivation target, an open physical attachment, or a
withdrawn claim.

**Finite realization.** The A1 class is inhabited: the released carrier
packet realizes the boundary complex exactly, the incidence-nerve federation
realizes twelve charts, thirty seams, and twenty coherent triple overlaps
with an operational observer, and the geodesic tower realizes the oriented
spherical support with mesh convergence and refinement naturality
(`code/a5_closure/` receipts and the simulator federation and bridge
producers). These receipts verify architecture and implementation; they do
not derive A1 from A2 or A3.

---

## A2. Observer agreement

**Plain language.** Observers operating on the screen agree on the meaning
of the data they jointly interpret.

**Formal statement.** For each regulator \(r\), let \(\mathsf{Data}_r\) be
the category of observer-accessible data and \(\mathsf{Meaning}_r\) the
category of operational interpretations, including probabilities, update
effects, and public record values. A2 states that the interpretation map

\[
\mathcal J_r:\mathsf{Data}_r\longrightarrow\mathsf{Meaning}_r
\]

is natural with respect to every visible overlap restriction, recharting,
seam translation, higher-overlap map, federation map, and refinement map:
every declared data-access diagram for accepted public data continues to
commute after interpretation. For two patches \(P,Q\) sharing \(O\),
accepted data satisfy

\[
\mathcal J_O\!\left(\operatorname{res}_{P\to O}d_P\right)
=
\mathcal J_O\!\left(
  \tau_{Q\to P}\operatorname{res}_{Q\to O}d_Q
\right),
\]

where \(\tau_{Q\to P}\) is the declared translation between the two observer
interfaces, and across resolutions
\(\mathcal J_r\circ c_{s\to r}=C_{s\to r}\circ\mathcal J_s\), where
\(c_{s\to r}\) transports accessible data and \(C_{s\to r}\) transports its
operational meaning.

**A2 constrains.** Naturality of operational meaning on accepted shared
data, across pairwise overlaps, higher overlaps, rechartings, seam
translations, and refinement.

**A2 does not imply.** Global extendability, termination, confluence, unique
public normal forms, Byzantine safety, record durability, or exact closure
of an exponential family. A2 applies to accepted shared data; raw records,
private coordinates, and repairable mismatches may differ before acceptance.
Local states are specified first: a global state cannot be assumed and then
used to make overlap agreement tautological. A1 (or an explicit
implementation interface) types the admissible or accepted record domain; A2
tests naturality on that domain and cannot define the domain over which it
quantifies. Temporal record persistence belongs to an A1 continuation
theorem, not to semantic agreement. A2 naturality does not imply that
independently optimized A3 states push forward to one another; optimizer
compatibility is a separate theorem.

**Countermodel boundary.** Federations satisfying A1 and A2 exist in which
accepted repair fails to terminate, in which distinct schedules reach
distinct raw configurations before acceptance, and in which no single global
state extends the compatible local family. These are retained as permanent
controls; the corresponding strengthenings are consensus interfaces, not
axioms.

---

## A3. Conditional maximum randomness

**Plain language.** Everything that observer agreement leaves unconstrained
is maximally random.

**Formal statement.** At finite regulator \(r\), a state is a compatible
family of local normalized states
\(\rho_r=(\rho_{r,P})_{P\in\mathcal P_r}\) on the accessible algebra net;
this typing does not presume a single global state extension. Let
\(\mathcal K_r\) be the nonempty convex set of such families satisfying the
finite observer-visible constraints supplied by A1 and A2. A valid
specification includes an A1-generated observable and constraint grammar, an
enumeration or factorization theorem showing that every A2-visible
constraint passes through that grammar, and omitted-constraint mutation
tests. The constraint family is target-independent, quotient-visible, and
finite in the declared regulator scheme.

A3 includes an exact reference and aggregation rule: a compatible local
reference family \(\tau_r\), a finite A1-generated observer cover
\(\mathcal G_r\subseteq\mathcal P_r\), and strictly positive exact weights
\(w_{r,P}\) from quotient-visible A1 data. The cover restriction map
\[
\rho\longmapsto(\rho_{r,P})_{P\in\mathcal G_r}
\]
is injective on \(\mathcal K_r\). Equivalently at smooth points, no nonzero
feasible tangent is invisible to every member of the cover. The aggregation
rule states its representation or base-measure dependence, cover and weight
normalization, and behavior under refinement, and it is natural under
admissible presentation equivalence. For compatible local families,

\[
\mathcal D_r(\rho\Vert\tau_r)
=
\sum_{P\in\mathcal G_r}
w_{r,P}\,D(\rho_{r,P}\Vert\tau_{r,P}),
\]

with finite-algebra Umegaki relative entropy in the declared representation
and trace convention. The realized state is the information projection

\[
\rho_r
=
\operatorname*{arg\,min}_{\rho\in\mathcal K_r}
\mathcal D_r(\rho\Vert\tau_r).
\]

When every local reference density is proportional to the identity relative
to its declared faithful matrix trace (or the weighted expectation of
\(\log\tau_{r,P}\) is constant on \(\mathcal K_r\)), this is equivalent to
weighted local entropy maximization. A general tracial state on a reducible
algebra may carry unequal central-block weights and does not establish this
equivalence. Invariance alone does not prove uniqueness on a reducible
algebra: if distinct invariant reference, cover, or weight rules survive,
the A3 specification for that use remains open and no dependent selection is
promoted.

Existence, uniqueness, support, and faithfulness conditions are stated for
the selected regulator model. A continuum use of A3 is defined through a
controlled finite-regulator or relative-modular limit and never relies on an
undefined von Neumann entropy for a type-III algebra.

**The meaning of "random".** Least informative relative to the declared
reference and the complete observer-agreement constraint set. It does not
mean independent noise, equal probabilities in every coordinate system, or
global thermal equilibrium.

**A3 object-type gate.** Every use of A3 names one of three optimizer
types: an ontic state on a fixed accessible algebra; an observer-inference
state on a fixed evidence algebra; or a transition or schedule distribution
on a fixed finite move simplex. A3 does not directly select a policy, repair
coefficient, field list, multiplicity, source character, response law, or
capacity ontology. Such a claim requires an A1-fixed state space, complete
constraints, and a proved map from the unique optimizer to that output;
otherwise it is `conditional_open_interface` or `discovery_only`.

**A3 model-space boundary.** A3 selects a state inside one A1-fixed feasible
space. It can compare port weights, record distributions, repair schedules,
sector occupancy, or scalar channels only when A1 supplies one common
ambient algebra, the A1 grammar generates the observables, and every
A2-visible constraint is proved to factor through that grammar. A3 cannot
compare unrelated Hilbert spaces, arbitrary field lists, or different
ontology classes without an independently derived common source space. This
prevents maximum randomness from becoming a renamed economy rule.

**A3 does not imply.** Optimizer pushforward across refinement, state
alignment, collar recovery, Markovity, mixing, maximally mixed edge
degeneracies, a \(\log d\) entropy term without a proved factorized
identity-proportional edge reference, exponential-family closure under
coarse-graining, or any selection among unrelated model spaces. The
identity-channel and collar countermodels delimit these nonimplications and
remain permanent controls.

---

## How the axioms combine

A1 supplies the observer-screen architecture. A2 turns compatible
interpretation into public meaning, which together with A1 fixes the
feasible family of public states. A3 selects the least informative state in
that family relative to the declared reference. Every gravity, gauge,
matter, flavor, and quantitative conclusion then displays the additional
structure it actually uses: a theorem covering every model of the applicable
axioms is `axiom_forced`; an exhaustively certified result inside one named
finite realization is `exact_named_realization`; explicit mathematical or
implementation premises are `conditional_open_interface`; countermodel-
limited claims are `independence_limited`; empirical attachments are
`physical_identification`; and simulation patterns without coverage theorems
are `discovery_only`.

Two scopes never share one label. An exact enumeration proves uniqueness
only inside its declared menu until a grammar-completeness theorem covers
all relevant three-axiom models; a named-realization witness proves
existence and exactness inside that realization, not that the axioms select
it.

## Glossary

- **Patch**: a finite system with local state, boundary, records, readback,
  and repair moves.
- **Observer**: a connected subfederation with complete readback, record,
  feedback or repair, prediction, and checkpoint interfaces.
- **Local carrier**: one A1 carrier with the twelve-port oriented boundary
  packet.
- **Federation screen (nerve)**: the overlap nerve \(N_r\) of the carrier
  federation with its seam and triple-overlap algebras.
- **Support screen**: the oriented spherical support \(S_r\) with its
  refinement tower and realization \(\iota_r\).
- **Public record**: accepted shared data in the typed record domain, the
  objects A2 quantifies over.
- **Meaning**: the operational interpretation of data: probabilities, update
  effects, public record values.
- **Constraint grammar**: the A1-generated family of observables and
  constraints through which every A2-visible constraint factors.
- **Reference state**: the exact A3 reference family \(\tau_r\) with its
  cover and weights.
- **\(A_5\)**: the alternating group on five letters, isomorphic to the
  orientation-preserving incidence automorphism group of the carrier
  boundary. The symbol names a group, never an axiom.

## Dependency rules

1. A surface that states or counts the basis gives exactly three axioms and
   pairs each informal sentence with its formal statement or a citation to
   this reference.
2. A surface that consumes an axiom names the exact dependency
   (`A1 observer_screen_architecture`, `A2 observer_meaning_agreement`,
   `A3 conditional_maximum_randomness`) and any additional interfaces from
   the registry catalogue.
3. Untyped phrases of the form "the axioms give recovery" or "the framework
   selects the matter content" are rejected by the registry checks.
4. Retired principles (the former recovery bundle and the former economy
   axiom) never appear as premises of active claims; their surviving content
   is carried by classified interfaces with countermodels.
