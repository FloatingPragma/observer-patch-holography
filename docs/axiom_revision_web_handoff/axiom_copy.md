# Axiom Copy: Canonical Renderings

Exact short, medium, long, and formal renderings of the three OPH axioms. `axiom_copy.json` mirrors this file field for field; sites consume the JSON and never paraphrase the formal text. Authority order: `claims/axiom_registry.yaml`, then `docs/AXIOM_REFERENCE.md`, then this bundle. Formal renderings are quoted verbatim from the axiom reference in KaTeX-safe notation: `\( \)` inline, `\[ \]` display, no custom macros.

---

## A1. Oriented twelve-port observer screen

Registry key: `observer_screen_architecture`

### Short

Each local observer carrier meets the screen through twelve primitive boundary ports arranged as the oriented boundary of an icosahedron, and carriers federate along typed seams into an oriented spherical support.

### Medium

There exists an observer patch net on an oriented spherical screen. At every finite resolution, each local carrier has twelve primitive boundary ports. The ports form the vertices of an oriented triangular boundary with 30 edges and 20 faces, combinatorially the boundary of an icosahedron: five edges meet each port, the neighbors of each port form a five-cycle, and every edge belongs to two oppositely oriented faces. Carriers join through typed seams and coherent triple overlaps, refine to an oriented spherical support, and expose local state, readback, records, repair moves, and checkpoints.

### Long

There exists an observer patch net on an oriented spherical screen. At every finite resolution, each local carrier has twelve primitive boundary ports. The ports form the vertices of an oriented triangular boundary with 30 edges and 20 faces, combinatorially the boundary of an icosahedron: five edges meet each port, the neighbors of each port form a five-cycle, and every edge belongs to two oppositely oriented faces. Carriers join through typed seams and coherent triple overlaps, refine to an oriented spherical support, and expose local state, readback, records, repair moves, and checkpoints.

A1 is the architecture axiom. It fixes the carrier boundary combinatorics (twelve ports, thirty edges, twenty faces, degree-five ports with five-cycle links, coherently oriented faces), the typed carrier-federation-support structure with its explicit maps, the operational interface surface (state access, readback, records, repair, checkpoints), and the refinement tower with a controlled spherical support realization. It names no gauge group, particle, coupling, or measured number. The orientation-preserving symmetry group of the boundary is the alternating group on five letters; that is a theorem about A1 carriers, never an axiom. Semantic agreement and state selection belong to A2 and A3.

### Formal

Let \((R,\preceq)\) be a directed regulator system. For
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

### Constrains

The local carrier boundary combinatorics, the typed carrier-federation-support architecture, the interface and repair type surface, and the refinement tower with its controlled spherical support realization.

### Does not imply

- semantic agreement (A2)
- any state selection (A3)
- repair termination, confluence, or global state extension
- equal state weights rho(e_p) = 1/12
- the integer load fiber or quadratic readback law
- the perturb/readback response law or the inverse-port response
- the inverse-port involution, six axes, A5 rotation group, or Gram frame pair as axiom content (these are derived incidence consequences)
- any physical selection between the Galois-conjugate frame pair
- the compact current algebra or any gauge, particle, coupling, or metric content
- collar recovery, generalized-entropy identities, or gravity structure
- a laboratory identification of any carrier object

---

## A2. Observer agreement

Registry key: `observer_meaning_agreement`

### Short

Observers operating on the screen agree on the meaning of the data they jointly interpret.

### Medium

For each finite resolution, interpretation maps observer-accessible data to operational meanings: probabilities, update effects, public record values. A2 states that this map commutes with every declared way of moving data around the screen: overlap restriction, recharting, seam translation, higher-overlap maps, federation maps, and refinement. Two observers who accept the same shared data assign it the same meaning.

### Long

For each finite resolution, interpretation maps observer-accessible data to operational meanings: probabilities, update effects, public record values. A2 states that this map commutes with every declared way of moving data around the screen: overlap restriction, recharting, seam translation, higher-overlap maps, federation maps, and refinement. Two observers who accept the same shared data assign it the same meaning.

A2 quantifies over accepted public data only. Raw records, private coordinates, and repairable mismatches may differ before acceptance. Local states are specified first; a global state is never assumed, which keeps overlap agreement from becoming tautological. A1 types the record domain that A2 tests; A2 cannot define the domain it quantifies over. Agreement implies neither repair termination, nor confluence, nor a single global state: countermodels for each are permanent controls, and the corresponding strengthenings are consensus interfaces.

### Formal

For each regulator \(r\), let \(\mathsf{Data}_r\) be
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

### Constrains

Naturality of operational meaning on accepted shared data, across pairwise overlaps, higher overlaps, rechartings, seam translations, and refinement.

### Does not imply

- global state extension
- repair termination or confluence
- unique public normal forms
- Byzantine safety, dissemination, or record durability
- closure of an exponential family under coarse-graining
- commutation of the A3 optimization with refinement (optimizer pushforward is a separate theorem)
- the domain of accepted data (A1 or an explicit implementation interface types it; A2 cannot define the domain it quantifies over)

---

## A3. Conditional maximum randomness

Registry key: `conditional_maximum_randomness`

### Short

Everything that observer agreement leaves unconstrained is maximally random.

### Medium

At each finite resolution, a state is a compatible family of local states on the accessible algebra net. A1 and A2 carve out the convex set of families satisfying every observer-visible constraint. A3 selects from that set the least informative family relative to a declared reference: the information projection under a weighted sum of local relative entropies over a declared, state-determining observer cover with positive exact weights. When every local reference density is identity-proportional in its declared trace, this equals weighted local entropy maximization.

### Long

At each finite resolution, a state is a compatible family of local states on the accessible algebra net. A1 and A2 carve out the convex set of families satisfying every observer-visible constraint. A3 selects from that set the least informative family relative to a declared reference: the information projection under a weighted sum of local relative entropies over a declared, state-determining observer cover with positive exact weights. When every local reference density is identity-proportional in its declared trace, this equals weighted local entropy maximization.

Random means least informative relative to the declared reference and the complete observer-agreement constraint set; it does not mean independent noise, equal probabilities in every coordinate system, or global thermal equilibrium. Every use of A3 names its optimizer type: an ontic state on a fixed accessible algebra, an observer-inference state on a fixed evidence algebra, or a transition distribution on a fixed finite move simplex. A3 selects a state inside one A1-fixed feasible space and cannot compare unrelated Hilbert spaces, field lists, or ontology classes; this model-space boundary keeps maximum randomness from becoming a renamed economy rule. Commutation of the optimizer with refinement is a separate theorem.

### Formal

At finite regulator \(r\), a state is a compatible
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

### Constrains

One state-selection principle: the least informative compatible local state family relative to the declared reference and the complete observer-agreement constraint set.

### Does not imply

- independent noise or equal probabilities in every coordinate system
- global thermal equilibrium
- exact closure of an exponential family under arbitrary coarse-graining
- optimizer pushforward across refinement without a separate compatibility theorem
- collar recovery, state alignment, Markovity, or mixing
- maximally mixed edge degeneracies or a log d entropy term without a proved factorized identity-proportional edge reference
- selection among unrelated Hilbert spaces, field lists, or ontology classes (the model-space boundary)
- any policy, repair coefficient, field content, multiplicity, source character, response law, or capacity ontology without an A1-fixed state space, complete constraints, and a proved optimizer-to-output map
