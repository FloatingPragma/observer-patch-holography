# Three-Axiom Learning-Material Handoff

Self-contained source packet for external authors of OPH learning material.
The normative sources are `claims/axiom_registry.yaml` (machine-readable) and
`docs/AXIOM_REFERENCE.md` (reader-facing reference); where this packet and
those sources differ, the sources win. Style rules for produced material are
binding and live in `docs/STYLE_GUIDE.md`. Formal statements in this packet
are quoted verbatim from the axiom reference in KaTeX-safe notation:
`\( \)` delimits inline math, `\[ \]` delimits display math, and no custom
macros are used.

Teaching copy produced from this packet never carries repository paths,
receipt names, issue numbers, or project-internal identifiers. The axiom
names A1, A2, A3, the group name \(A_5\), and the seven status-class words
are public vocabulary and may appear.

---

## 1. One-sentence description

Observer Patch Holography (OPH) is a research program that treats finite
observers and their agreements as the primitive layer of physics and derives
physical structure from three axioms: an observer-screen architecture, an
observer-agreement law, and a conditional maximum-randomness principle.

---

## 2. The three axioms in one page

Physics texts usually begin with fields on a spacetime background and
introduce observers at the end, as consumers of predictions. OPH inverts the
order of presentation and takes the observers as primitive. The basic
objects are patches: finite systems with local state, boundary ports,
records, readback, and repair moves. Three axioms govern how patches form a
world.

A1, the oriented twelve-port observer screen, is the architecture axiom.
Every local carrier, at every finite resolution, has twelve primitive
boundary ports. The ports sit as the vertices of a triangulated boundary
with thirty edges and twenty faces: five edges meet each port, the
neighbors of each port form a five-cycle, and every edge lies in two faces
whose orientations disagree along it. Those counting rules force the
boundary to be, combinatorially, the surface of an icosahedron. Carriers
join along typed seams, seams cohere on triple overlaps, and the resulting
federation refines toward an oriented spherical support. A1 states what an
observer screen is. It names no particle, no coupling, and no measured
number.

A2, observer agreement, is the semantic axiom. When two observers hold data
about a shared region and each interprets it, translating the data first
and interpreting second gives the same result as interpreting first and
translating the meaning. The axiom quantifies over accepted public data
only. Two observers may hold conflicting raw records, and repair may run
for a long time; A2 is silent about all of that. It requires only that once
data is accepted as public, its meaning cannot depend on who reads it or
through which chart.

A3, conditional maximum randomness, is the selection axiom. A1 fixes which
observables exist and A2 fixes which constraints bind the public state.
Among all states compatible with those constraints, A3 selects the least
informative one relative to a declared reference state: the information
projection. Random here has a precise meaning: least informative given the
reference and the complete constraint set. It does not mean uniform
weights, independent noise, or thermal equilibrium.

Everything beyond these three statements is derived or explicitly
classified. Some consequences are theorems that hold in every model of the
axioms. Some are exact results certified inside one named finite
realization. Some are open interfaces with stated premises, some are
limited by countermodels, some are empirical identifications, and some are
withdrawn. The icosahedral boundary has the alternating group on five
letters as its orientation-preserving symmetry group, which is why the
symbol \(A_5\) appears throughout the corpus. The symbol names a group,
never an axiom.

---

## 3. Axiom cards

Each axiom has one informal card and one formal card. Informal cards quote
the plain-language statements; formal cards quote the formal statements
from the axiom reference verbatim.

### A1. Oriented twelve-port observer screen

**Informal card.** There exists an observer patch net on an oriented
spherical screen. At every finite resolution, each local carrier has twelve
primitive boundary ports. The ports form the vertices of an oriented
triangular boundary with 30 edges and 20 faces, combinatorially the
boundary of an icosahedron: five edges meet each port, the neighbors of
each port form a five-cycle, and every edge belongs to two oppositely
oriented faces. Carriers join through typed seams and coherent triple
overlaps, refine to an oriented spherical support, and expose local state,
readback, records, repair moves, and checkpoints.

**Formal card.**

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

### A2. Observer agreement

**Informal card.** Observers operating on the screen agree on the meaning
of the data they jointly interpret.

**Formal card.**

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

### A3. Conditional maximum randomness

**Informal card.** Everything that observer agreement leaves unconstrained
is maximally random.

**Formal card.**

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

---

## 4. Constrains / does not imply

One table per axiom. The "does not imply" rows are the canonical
nonimplication lists; each is backed by a countermodel or an explicit scope
rule in the axiom reference.

### A1

| Class | Content |
|---|---|
| Constrains | The local carrier boundary combinatorics, the typed carrier-federation-support architecture, the interface and repair type surface, and the refinement tower with its controlled spherical support realization. |
| Does not imply | Semantic agreement (that is A2). |
| Does not imply | Any state selection (that is A3). |
| Does not imply | Repair termination, confluence, or global state extension. |
| Does not imply | Equal state weights \(\rho(e_p)=1/12\). |
| Does not imply | The integer load fiber or quadratic readback law. |
| Does not imply | The perturb/readback response law or the inverse-port response. |
| Does not imply | The inverse-port involution, six axes, the \(A_5\) rotation group, or the Gram frame pair as axiom content; these are derived incidence consequences, proved separately. |
| Does not imply | Any physical selection between the Galois-conjugate frame pair. |
| Does not imply | The compact current algebra or any gauge, particle, coupling, or metric content. |
| Does not imply | Collar recovery, generalized-entropy identities, or gravity structure. |
| Does not imply | A laboratory identification of any carrier object. |

### A2

| Class | Content |
|---|---|
| Constrains | Naturality of operational meaning on accepted shared data, across pairwise overlaps, higher overlaps, rechartings, seam translations, and refinement. |
| Does not imply | Global state extension. |
| Does not imply | Repair termination or confluence. |
| Does not imply | Unique public normal forms. |
| Does not imply | Byzantine safety, dissemination, or record durability. |
| Does not imply | Closure of an exponential family under coarse-graining. |
| Does not imply | Commutation of the A3 optimization with refinement; optimizer pushforward is a separate theorem. |
| Does not imply | The domain of accepted data; A1 or an explicit implementation interface types that domain, and A2 cannot define the domain it quantifies over. |

### A3

| Class | Content |
|---|---|
| Constrains | One state-selection principle: the least informative compatible local state family relative to the declared reference and the complete observer-agreement constraint set. |
| Does not imply | Independent noise or equal probabilities in every coordinate system. |
| Does not imply | Global thermal equilibrium. |
| Does not imply | Exact closure of an exponential family under arbitrary coarse-graining. |
| Does not imply | Optimizer pushforward across refinement without a separate compatibility theorem. |
| Does not imply | Collar recovery, state alignment, Markovity, or mixing. |
| Does not imply | Maximally mixed edge degeneracies or a \(\log d\) entropy term without a proved factorized identity-proportional edge reference. |
| Does not imply | Selection among unrelated Hilbert spaces, field lists, or ontology classes; this is the model-space boundary. |
| Does not imply | Any policy, repair coefficient, field content, multiplicity, source character, response law, or capacity ontology without an A1-fixed state space, complete constraints, and a proved optimizer-to-output map. |

---

## 5. Worked finite toy models

### 5.1 A1: the icosahedral boundary and a two-carrier seam

Take twelve ports labeled \(1,\dots,12\). Place port 1 as an upper pole,
ports \(2,3,4,5,6\) as an upper pentagon, ports \(7,8,9,10,11\) as a lower
pentagon offset by half a step, and port 12 as a lower pole. The edge set:

- five edges from port 1 to each of \(2,\dots,6\);
- five pentagon edges \(2\text-3, 3\text-4, 4\text-5, 5\text-6, 6\text-2\);
- ten cross edges, each upper-pentagon port joined to the two nearest
  lower-pentagon ports (for example \(2\text-7\), \(2\text-8\), \(3\text-8\),
  \(3\text-9\), and so on around);
- five pentagon edges \(7\text-8, 8\text-9, 9\text-10, 10\text-11, 11\text-7\);
- five edges from port 12 to each of \(7,\dots,11\).

That is \(5+5+10+5+5=30\) edges. The faces are the five upper-cap triangles
around port 1, ten band triangles alternating between the pentagons, and
five lower-cap triangles around port 12, for 20 faces. Checks a reader can
run by hand: the Euler characteristic is \(12-30+20=2\), the sphere value;
the handshake identity gives average degree \(2\cdot30/12=5\), and every
port has degree exactly five; the link of port 1 is the five-cycle
\((2,3,4,5,6)\). Orient every face so that walking its boundary keeps the
face on the left. The edge \(2\text-3\) then appears in the cap face
\((1,2,3)\) traversed \(2\to3\) and in the band face containing \(2,3\)
traversed \(3\to2\): two faces, opposite induced orientations. These
incidence rules pin the complex down as the boundary of an icosahedron. No
coordinates, lengths, gauge labels, or particle labels enter; the packet is
pure incidence plus orientation.

For the seam example, take two carriers \(a\) and \(b\), each with its own
twelve-port boundary packet and its own accessible algebra
\(\mathcal A_a\), \(\mathcal A_b\). A declared overlap between them carries
a seam algebra \(\mathcal B_{ab}\) with unital restriction maps
\(\rho_{a\to ab}:\mathcal A_a\to\mathcal B_{ab}\) and
\(\rho_{b\to ab}:\mathcal A_b\to\mathcal B_{ab}\). Concretely, let the seam
expose two shared readout channels: the observables that carrier \(a\)
reads through its ports 2 and 3 restrict to the same seam elements as the
observables carrier \(b\) reads through its ports 7 and 8. The seam is an
incidence-typed shared boundary region with its own algebra; it adds no
thirteenth port to either carrier, and it does not identify the two
carriers' hidden presentations. When a third carrier joins and the three
overlaps have a common region, A1 requires a triple-overlap algebra with
typed restrictions from all three pairwise seams and a coherence condition
on the restriction cocycles. The two-carrier picture is where meaning
transport (A2) gets its stage.

### 5.2 A2: two patches translating a shared record

Patches \(P\) and \(Q\) share an overlap region \(O\). The overlap carries
a two-outcome readout. Patch \(P\) encodes the outcomes as \(\{0,1\}\);
patch \(Q\) encodes the same physical outcomes as \(\{-,+\}\). The declared
translation is \(\tau_{Q\to P}(+) = 1\), \(\tau_{Q\to P}(-) = 0\).
Interpretation \(\mathcal J\) sends an accepted data record to its
operational meaning; in this toy, the meaning of a record is the empirical
frequency it assigns to each outcome.

Accepted data: \(d_P\) is a record whose restriction to \(O\) reads
\((1,0,1)\); \(d_Q\) is a record whose restriction to \(O\) reads
\((+,-,+)\). The A2 naturality square:

```
   d_Q  --- restrict to O --->  (+, -, +)
    |                               |
    |  translate tau_{Q->P}         |  translate tau_{Q->P}
    v                               v
  (translated record)  --->     (1, 0, 1)
    |                               |
    |  interpret J_O                |  interpret J_O
    v                               v
  meaning              =        frequency 2/3 for outcome 1
```

Both paths land on the same meaning: outcome 1 (equivalently \(+\)) has
frequency \(2/3\). In symbols, the accepted data satisfy

\[
\mathcal J_O\!\left(\operatorname{res}_{P\to O} d_P\right)
=
\mathcal J_O\!\left(\tau_{Q\to P}\operatorname{res}_{Q\to O} d_Q\right).
\]

The scope rule matters as much as the square. Suppose carrier \(Q\)'s raw
log reads \((+,+,+)\) before repair, while \(P\)'s reads \((1,0,1)\). A2
says nothing about this pair, because neither is accepted shared data. The
mismatch is a repair problem, and A1 supplies the repair interfaces.
Only after acceptance does the square bind. This is also why A2 cannot be
made tautological by positing one global state and restricting it: local
states come first, and agreement is a condition on their translations.

### 5.3 A3: a two-outcome information projection, computed

Outcome space \(\{a,b\}\); a state is a probability pair
\((p_a, p_b)\) with \(p_a+p_b=1\). Declared reference: the uniform state
\(\tau=(1/2,1/2)\). The relative entropy is

\[
D(p\Vert\tau) = p_a\ln\frac{p_a}{1/2} + p_b\ln\frac{p_b}{1/2}
             = \ln 2 - H(p),
\]

with \(H\) the Shannon entropy in nats. Suppose observer agreement
certifies one constraint: for the agreed observable \(X\) with \(X(a)=1\),
\(X(b)=0\), accepted public records force
\(\mathbb E_p[X] = p_a \ge 0.7\). The feasible set is
\(\mathcal K = \{p : p_a\ge 0.7\}\), a closed convex set. The information
projection is

\[
\rho = \operatorname*{arg\,min}_{p\in\mathcal K} D(p\Vert\tau).
\]

The unconstrained minimizer is \(\tau\) itself, which violates the
constraint, so the minimizer sits on the boundary \(p_a=0.7\):
\(\rho=(0.7,\,0.3)\), with

\[
D(\rho\Vert\tau) = \ln 2 - H(0.7) \approx 0.6931 - 0.6109 = 0.0823
\ \text{nats}.
\]

With the equality constraint \(\mathbb E_p[X]=0.7\) instead, the same
projection arises in exponential-family form
\(p_a = e^{\lambda}/(e^{\lambda}+1)\) with
\(\lambda=\ln(0.7/0.3)\approx 0.847\): the projection tilts the reference
by exactly the constraint and no more.

Two teaching points follow directly from the computation. First, if the
constraint set does not bind (say the record only certifies
\(p_a \ge 0.4\)), the projection is the reference itself: maximum
randomness leaves everything at the reference. Second, "random" is
reference-relative. Repeat the computation with the non-uniform reference
\(\tau'=(3/4,1/4)\) and no binding constraint: the projection is
\((3/4,1/4)\), and equal weights appear nowhere. A3's "maximally random"
means least informative relative to the declared reference and the
complete constraint set, and the full axiom states exactly how the
reference family, the observer cover, and the weights are supplied.

---

## 6. Counterexamples: what the axioms do not force

Each paragraph names a property that the three axioms do not imply and
sketches the shape of the countermodel that proves it. These countermodels
are permanent controls; learning material places them next to the claims
they limit.

**Global extension.** Take three patches with pairwise overlaps, each
overlap carrying a \(\pm1\) readout, and local states asserting perfect
anticorrelation across each of the three overlaps. Every pairwise marginal
is consistent, so A1 and A2 are satisfied on accepted overlap data. A
single global assignment would need three \(\pm1\) values whose three
pairwise products are all \(-1\); the product of the three pairwise
products is a perfect square, hence \(+1\), while three factors of \(-1\)
multiply to \(-1\). No global state extends the compatible local family.
Global extension is a consensus interface, never a consequence of the
axioms.

**Confluence.** Build a federation whose repair system has one mismatch
configuration with two applicable repair moves, leading to two distinct
terminal raw configurations, with no further move joining them. Both runs
satisfy every A1 interface and every A2 naturality square on accepted
data, because acceptance happens after a schedule is chosen. Distinct
schedules reach distinct raw configurations; the diamond property is an
added interface, and so is termination.

**Exponential-family closure.** At a fine resolution the A3 projection
lies in an exponential family generated by the constraint observables.
Coarse-grain by a many-to-one outcome map that merges outcomes across the
constraint level sets. The image of the exponential family is a mixture
family that fails closure with respect to the pushed-forward statistics,
and the coarse-level projection differs from the pushforward of the
fine-level projection. Neither A2 naturality nor A3 forces these to agree;
optimizer compatibility across coarse-graining and refinement is a
separate theorem with its own premises.

**Collar recovery.** Take a carrier whose collar channel acts as the
identity on one subalgebra and erases its complement. Every A1 interface
is present, accepted data satisfy A2, and the A3 projection exists and is
unique on the visible algebra. Yet no channel recovers the erased
complement: the identity-channel countermodel shows that recoverability of
collar data is an added gravity-side interface. A related alignment
countermodel (two carriers sharing an entangled pair whose local states
match while no aligning isomorphism exists) delimits state alignment the
same way. Recovery statements always carry their interface premises.

**Matter selection.** Take any model of the three axioms and adjoin a
sterile sector: extra degrees of freedom invisible to every declared
interface, record, and constraint. The constraint set is unchanged, the
information projection on the visible algebra is unchanged, and every
axiom remains satisfied. The axioms therefore do not select a unique
matter package, do not fix a generation count, and do not exclude extra
light sectors. Exact matter results live inside named finite realizations
with declared premises, and completeness of any menu is a separate open
question.

---

## 7. Dependency map: from three axioms to the corpus

The map below states, for each downstream area, what the axioms supply,
what additional named structure the area consumes, and the status class of
the result. Status classes are the seven-value vocabulary:
`axiom_forced`, `exact_named_realization`, `discovery_only`,
`conditional_open_interface`, `independence_limited`,
`physical_identification`, `withdrawn`.

| Area | What the three axioms supply | Additional named structure | Status class of the headline result |
|---|---|---|---|
| Public records and the event algebra | A1 types records and checkpoints; A2 makes accepted meaning chart-independent. Total atom readouts glue to finite compatible global sections, and their functions form a commutative event algebra. | Stated gluing premises of the record theorem. | Theorem under its stated premises; the physical-universe capacity attachment is `conditional_open_interface`. |
| Lorentz kinematics | A1 supplies the oriented spherical support with its refinement tower. On the typed global-support branch, the orientation-preserving conformal group of the two-sphere is isomorphic to the proper orthochronous Lorentz group, events organize as 3+1 records, and observer space is hyperbolic three-space. | The typed global-support branch: the finite-cap certificate clauses and the complete modular algebra-state package on the same refinement tower. | `conditional_open_interface` (conditional theorem on the named branch). |
| Einstein branch | A2-accepted repair bookkeeping and the A3 projection feed an entropy-stationarity argument yielding the Einstein equation with one metric-proportional term per connected component; coupling universality is theorem-grade for every icosahedrally equivariant law on the declared family. | The named gravity interfaces listed below, plus two pending physical identifications. | `conditional_open_interface` (conditional composition theorem; every clause instrumented). |
| Finite gauge chain | A1's boundary incidence alone proves the \(1+3+3'+5\) coefficient decomposition, the unique central graph involution \(J\), and \(10J=A^3-4A^2-5A+10I\); a target-blind impulse/readback producer derives the response \(R=-J\); the current certificate constructs \(\mathfrak{su}(3)+\mathfrak{su}(2)+\mathfrak u(1)\); determinant balance and tensor-kernel descent compute the common \(\mathbb Z_6\) kernel. | The named carrier realization and its declared premises: a full-rank physical current with inner \(A_5\) action, the conjugate matter-projector pair. No matter-content selection rule enters. | `exact_named_realization` for the chain inside the named realization; laboratory current identification is `physical_identification`, attachment pending. |
| Matter witnesses | Inside the named realization, the selected exterior package \(\Lambda^2V\oplus\Lambda^4V\) over a trace-balanced carrier of dimensions three and two branches to exactly one 15-state chiral generation (the rank-15 witness), cancels the listed anomalies with hypercharge balance, and carries the \(\mathbb Z_6\) kernel. | The declared exterior-module menu (exhaustively scanned; menu completeness open); one common multiplicity object must precede any generation number. | `exact_named_realization` for the witnesses. The generation count is open inside the conditional window \(3\le N_g\le 5\); \(N_g=3\) is a declared completion, never a derivation. |
| Quantitative tests | A3's projection plus the record capacity yield declared closure maps whose fixed points are interval-certified: each declared map has exactly one fixed point on its stated domain. | Physical transport and endpoint maps, scheme maps, and identification premises, each named per test. | `physical_identification` with attachment pending for the measured comparisons; the arithmetic fixed-point certificates are exact on their stated domains. |

Named interfaces of the Einstein branch (each is a catalogued premise, and
each carries its own countermodels or controls):

| Interface | Plain description | Class |
|---|---|---|
| `edge_center_decomposition` | The cut algebra decomposes into edge and center parts. | `conditional_open_interface` |
| `central_interface_state_alignment` | Alignment of central interface states across a cut; delimited by Bell-pair and identity-channel countermodels. | `independence_limited` |
| `exact_conditional_markov_split` | The state satisfies an exact conditional Markov split across the collar. | `conditional_open_interface` |
| `approximate_collar_recovery_rate` | A quantitative recovery rate for collar data. | `conditional_open_interface` |
| `edge_reference_factorization_and_traciality` | The edge reference factorizes and is tracial. | `conditional_open_interface` |
| `edge_entropy_split_and_cut_additivity` | Edge entropy splits and adds across cuts. | `conditional_open_interface` |
| `finite_generalized_entropy_identification` | The finite generalized-entropy functional is identified. | `conditional_open_interface` |
| `fixed_cap_entropy_stationarity` | Entropy stationarity at fixed cap. | `conditional_open_interface` |
| `support_visible_bw_modular_normalization` | The support-visible modular normalization clause. | `conditional_open_interface` |
| `repair_conservation_and_stress_attachment` | Repair bookkeeping conserves and attaches to a stress tensor. | `conditional_open_interface` |
| `fixed_volume_small_ball_variation` | The small-ball variation at fixed volume. | `conditional_open_interface` |
| `ward_bianchi_coupling_equality` | Equality of the Ward and Bianchi couplings. | `conditional_open_interface` |
| `quantum_focusing_extension` | Optional focusing extension; no registered active theorem consumes it. | `conditional_open_interface` |
| `vacuum_reference_and_cosmological_constant` | The vacuum reference and cosmological-constant attachment. | `physical_identification`, attachment pending |
| `absolute_metric_and_newton_scale` | The absolute metric normalization and Newton-scale attachment. | `physical_identification`, attachment pending |

---

## 8. Glossary

- **Patch**: a finite system with local state, boundary, records, readback,
  and repair moves.
- **Observer**: a connected subfederation with complete readback, record,
  feedback or repair, prediction, and checkpoint interfaces. One carrier is
  never automatically one observer.
- **Local carrier**: one A1 carrier with the twelve-port oriented boundary
  packet.
- **Federation screen (nerve)**: the overlap nerve of the carrier
  federation with its seam and triple-overlap algebras.
- **Support screen**: the oriented spherical support with its refinement
  tower and geometric realization.
- **Public record**: accepted shared data in the typed record domain, the
  objects A2 quantifies over.
- **Meaning**: the operational interpretation of data: probabilities,
  update effects, public record values.
- **Constraint grammar**: the A1-generated family of observables and
  constraints through which every A2-visible constraint factors.
- **Reference state**: the exact A3 reference family with its observer
  cover and weights.
- **\(A_5\)**: the alternating group on five letters, isomorphic to the
  orientation-preserving incidence automorphism group of the carrier
  boundary. The symbol names a group, never an axiom.

---

## 9. Lesson sequence

Five lessons, each with objectives, exercises, and an answer key. The
sequence assumes comfort with basic probability and no physics background.

### Lesson 1: Reading the screen

Objectives: state the A1 boundary counting rules; verify the combinatorial
facts by hand; distinguish the three screen objects (carrier boundary,
federation nerve, support screen).

Exercises.

1.1 Compute the Euler characteristic of the carrier boundary from its
counts and say what surface it indicates.

1.2 Use the handshake identity to find the average port degree, and state
the extra A1 clause that makes every port reach that degree.

1.3 The orientation-preserving incidence automorphism group of the
boundary has order 60. Which named group is it, and what is its status in
the theory: axiom or theorem?

Answer key.

1.1 \(12-30+20=2\), the Euler characteristic of the sphere.

1.2 \(2E/V = 60/12 = 5\); A1 requires every port to have degree exactly
five with a five-cycle link, so the average is attained at every port.

1.3 The alternating group \(A_5\). Its role is a theorem about A1
carriers, proved from the incidence rules; it is never an axiom, and the
axioms do not mention the group.

### Lesson 2: Carriers, seams, observers

Objectives: define patch, local carrier, federation, and observer; explain
seam algebras and triple-overlap coherence; explain why interface
availability differs from repair termination.

Exercises.

2.1 In the two-carrier seam example of the toy models, list the typed data
the overlap must carry.

2.2 Explain why a single carrier is never automatically an observer.

2.3 True or false: A1 guarantees that accepted repair terminates.

Answer key.

2.1 A seam algebra for the overlap, unital restriction maps from each
carrier's accessible algebra into it, and, when a third carrier shares the
region, a triple-overlap algebra with typed restrictions from the pairwise
seams and coherent restriction cocycles.

2.2 An observer is a connected subfederation with complete availability of
readback, record, feedback or repair, prediction, and checkpoint
interfaces. A single carrier may lack part of that interface set.

2.3 False. A1 asserts interface availability; termination and confluence
are separate interfaces with countermodels.

### Lesson 3: Agreement as naturality

Objectives: draw and read the A2 naturality square; apply it to a
translation example; state the acceptance scope rule.

Exercises.

3.1 Two patches encode a shared three-symbol record as \((1,0,1)\) and
\((+,-,+)\) with translation \(+\mapsto1\), \(-\mapsto0\). Verify the
naturality square for the frequency interpretation.

3.2 Give two kinds of data that A2 does not constrain.

3.3 Does A2 imply that a single global state underlies all patches?
Sketch the countermodel.

Answer key.

3.1 Translating then interpreting gives frequency \(2/3\) for outcome 1;
interpreting the first record directly gives the same \(2/3\); the square
commutes.

3.2 Raw records before acceptance, private coordinates, and repairable
mismatches.

3.3 No. Three patches with pairwise perfect anticorrelation on their
overlaps have compatible pairwise data, while the parity of the three
pairwise products obstructs any single global assignment.

### Lesson 4: Conditional maximum randomness

Objectives: define the information projection; compute it in a two-outcome
model; explain what "random" does and does not mean.

Exercises.

4.1 With uniform reference and the constraint \(p_a\ge0.7\), find the
projection and its relative entropy in nats.

4.2 With reference \((3/4,1/4)\) and no binding constraint, find the
projection, and state what this shows about equal weights.

4.3 With uniform reference and the equality constraint
\(\mathbb E[X]=0.7\), write the projection in exponential-family form and
give \(\lambda\).

Answer key.

4.1 \(\rho=(0.7,0.3)\); \(D=\ln2-H(0.7)\approx0.0823\) nats.

4.2 The projection is the reference \((3/4,1/4)\) itself. Maximum
randomness is reference-relative; it does not mean equal probabilities.

4.3 \(p_a=e^\lambda/(e^\lambda+1)\) with
\(\lambda=\ln(0.7/0.3)\approx0.847\).

### Lesson 5: From axioms to physics: the status vocabulary

Objectives: name the seven status classes; classify sample claims; state
the forbidden promotions.

Exercises.

5.1 Classify: "the orientation-preserving symmetry group of the carrier
boundary is \(A_5\)."

5.2 Classify: "the screen-current chain yields
\(\mathfrak{su}(3)+\mathfrak{su}(2)+\mathfrak u(1)\)."

5.3 A draft page says "OPH derives three generations of matter." Repair
the sentence.

Answer key.

5.1 A theorem about every A1 carrier, proved from incidence; on a public
surface it reads as an `axiom_forced` consequence of A1.

5.2 `exact_named_realization`: an exact result inside the named carrier
realization under its declared premises (a full-rank physical current with
inner \(A_5\) action and the conjugate matter-projector pair). It is never
stated as a consequence of the axioms alone.

5.3 "Inside the named realization, one exact 15-state generation witness
exists, and the generation count is open inside the conditional window
\(3\le N_g\le5\); the value three is a declared completion, and its family
attachment is an open interface."

---

## 10. FAQ and objection responses

**How many axioms does OPH have?** Exactly three: A1 the oriented
twelve-port observer screen, A2 observer agreement, A3 conditional maximum
randomness. Everything else in the corpus is a theorem, an exact result
inside a named finite realization, a declared open interface, an
independence result with countermodels, a physical identification, or a
withdrawn claim.

**Is recovery an axiom?** No. Collar recoverability, generalized edge
entropy, fixed-cap stationarity, quantum focusing, and durable records are
each a classified gravity derivation target or a physical identification
with its own countermodels; durable records are typed as observer-screen
continuation contracts. A retired recovery bundle never appears as a
premise of an active claim.

**Does maximum randomness imply recovery?** No. The identity-channel
countermodel satisfies all three axioms while erased collar data stays
unrecoverable, and the Bell-pair countermodel delimits state alignment the
same way. Every recovery statement carries its interface premises
explicitly.

**Does OPH derive three generations?** No. One exact 15-state generation
witness exists inside the named realization. The generation count is open
inside the conditional window \(3\le N_g\le5\), and the value three is a
declared completion pending the family-attachment derivation. Economy-based
generation counting is withdrawn.

**Does the gauge theorem use a selection rule?** No. The finite chain from
boundary incidence to \(\mathfrak{su}(3)+\mathfrak{su}(2)+\mathfrak u(1)\)
and the \(\mathbb Z_6\) kernel is exact inside the named realization under
its declared premises, and no matter-content selection rule enters it.
Economy selection is retired as a global axiom; finite order definitions
survive only as diagnostics and select no physics.

**What does \(A_5\) mean?** The alternating group on five letters. It is
isomorphic to the orientation-preserving incidence automorphism group of
the twelve-port carrier boundary, which is why it appears throughout the
corpus. The symbol names a group, never an axiom, and the basis contains
no axiom numbered five.

**Does A2 secretly assume a global state?** No. Local states are specified
first, and agreement is a naturality condition on their translations. A
global state cannot be assumed and then used to make overlap agreement
tautological, and a three-patch parity countermodel shows compatible local
families without any global extension.

**Is "maximally random" a physical claim about noise?** No. It means least
informative relative to the declared reference and the complete
observer-agreement constraint set. It does not mean independent noise,
equal probabilities in every coordinate system, or global thermal
equilibrium.

---

## 11. Website copy blocks

### Short (one paragraph)

Observer Patch Holography rests on three axioms. The first fixes the
architecture: finite observers meet the world through carriers with twelve
boundary ports arranged, combinatorially, as the surface of an icosahedron,
joined by typed seams into a spherical screen. The second fixes semantics:
observers agree on the meaning of the data they jointly accept. The third
fixes the state: whatever agreement leaves unconstrained is maximally
random, in the precise sense of an information projection. Everything else
is a theorem, an exact result inside a named finite model, a declared open
interface, or a withdrawn claim, and each result carries its label.

### Medium (three paragraphs)

Observer Patch Holography starts where most physics ends: with the
observers. Its primitive objects are patches, finite systems with local
state, records, readouts, and repair moves. Three axioms govern them. The
architecture axiom states that every local carrier has twelve primitive
boundary ports forming an oriented triangulated boundary with thirty edges
and twenty faces, the combinatorics of an icosahedron, and that carriers
federate along typed seams into an oriented spherical support.

The agreement axiom states that observers agree on the meaning of the data
they jointly interpret: translating accepted data between observer
interfaces and then interpreting gives the same result as interpreting and
then translating the meaning. The randomness axiom states that everything
agreement leaves unconstrained is maximally random, realized as the least
informative state compatible with all observer-visible constraints,
relative to a declared reference.

The discipline of the program lies in what the axioms do not claim. The
symmetry group of the carrier boundary, the alternating group on five
letters, is a theorem, never an axiom. The chain from boundary incidence
to the gauge algebra of the Standard Model is exact inside one named
finite realization under declared premises. The generation count is open
inside a stated window, with three as a declared completion. Recovery,
confluence, and global states are interfaces with countermodels. Every
public claim carries one of seven status labels, and countermodels are
published next to the claims they limit.

### Long (one page)

Most presentations of fundamental physics begin with fields on a spacetime
background and introduce observers at the end. Observer Patch Holography
inverts that order. The primitive layer is a network of finite observers:
patches with local state, boundary ports, records, readback, and repair
moves. Three axioms, and only three, govern this layer.

The first axiom fixes the architecture of the observer screen. At every
finite resolution, each local carrier meets the world through twelve
primitive boundary ports. The ports form the vertices of an oriented
triangulated boundary with thirty edges and twenty faces: five edges meet
each port, the neighbors of each port form a five-cycle, and every edge
lies in two faces with opposite orientations. These counting rules force
the boundary to be, combinatorially, the surface of an icosahedron.
Carriers join along typed seams, coherent on triple overlaps, and the
federation refines toward an oriented spherical support. The axiom names
no particle, no coupling, and no measured number.

The second axiom fixes semantics. Observers operating on the screen agree
on the meaning of the data they jointly interpret: for accepted public
data, translating between observer interfaces and then interpreting gives
the same result as interpreting and then translating the meaning. The
axiom is silent about raw records, private coordinates, and mismatches
under repair; it binds only after acceptance.

The third axiom fixes the state. The first two axioms determine which
observables exist and which constraints bind the public state. Among all
compatible states, the realized one is the least informative relative to a
declared reference: the information projection. Random here is a technical
term, relative to reference and constraints; it does not mean uniform,
independent, or thermal.

Everything else is derived or classified under a seven-value status
vocabulary. The alternating group on five letters appears as the symmetry
group of the carrier boundary: a theorem, never an axiom. From boundary
incidence, a finite chain constructs the gauge algebra of the Standard
Model and its \(\mathbb Z_6\) global structure, exact inside one named
finite realization under declared premises, with the laboratory
identification tracked separately. One exact 15-state matter generation
witness exists inside that realization; the generation count is open
inside the window three to five, with three as a declared completion. On
the gravity side, an entropy-stationarity argument yields the Einstein
equation as a conditional composition theorem whose named interface
premises are instrumented end to end. Quantitative comparisons are run
under a falsification program with certified arithmetic and frozen
verdicts. Where the axioms do not force a property (global states,
confluence, recovery, matter selection), a countermodel is published next
to the claim it limits. The program's claim discipline is the point: three
axioms, explicit interfaces, and labels that move only by dated artifact.

---

## 12. Term migration table

Internal mapping from retired terminology to the terminology of this
packet. Learning material uses only the replacement terms; the retired
terms appear solely in scan lists and redirects.

| Retired term | Replacement | Disposition |
|---|---|---|
| Five axioms (also "OPH5") | Three axioms: A1 observer screen architecture, A2 observer agreement, A3 conditional maximum randomness | The basis contains exactly three core axioms; any surface counting five is stale. |
| Recoverable generalized entropy axiom (the recovery bundle) | The gravity interface catalogue | Each component (collar recoverability, generalized edge entropy, fixed-cap stationarity, quantum focusing, durable records) is a classified derivation target or physical identification with countermodels; durable records are observer-screen continuation contracts. |
| Minimal admissible realization / economy rule | Declared completions with open status | Economy selection is retired as a global axiom; finite order definitions survive only as diagnostics and select no physics. Economy-based generation counting, one-Higgs selection, and no-extra-light-sector claims are withdrawn. |
| MaxEnt with refinement stability | Conditional maximum randomness (A3), with optimizer compatibility as a separate theorem | A3 selects the information projection at each finite regulator; commutation of the optimizer with refinement (optimizer pushforward) is a catalogued interface, never part of the axiom. |
