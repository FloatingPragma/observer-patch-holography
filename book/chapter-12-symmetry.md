# Chapter 12: Symmetry on the Sphere

## 12.1 The Intuitive Picture: Symmetries Are Aesthetic Choices

Start with the aesthetic picture of symmetry.

Symmetries are aesthetic preferences. The universe could have been asymmetric,
lopsided, or irregular, but it happens to be symmetric in certain ways.
Physicists chose to study symmetric systems because they're easier to analyze
and more beautiful. Symmetry is a convenience, not a necessity.

This view treats symmetry as a happy accident or an unexplained gift. The laws of physics happen to look the same in all directions (rotational symmetry). They happen to be stable from one day to the next (time translation symmetry). But there's no deeper reason for this. The universe could have been otherwise.

Noether broke that picture.

## 12.2 The Surprising Hint: Symmetries Imply Conservation Laws

In 1918, Emmy Noether proved one of the most important theorems in physics.

### Noether's Revolution

Noether had come to Göttingen in 1915 at the invitation of Hilbert and Klein, to help them understand energy conservation in general relativity. What she discovered was far more general.

Her position in that story matters. Noether was one of the strongest
algebraists in Europe, yet for years her Göttingen courses had to be
advertised under Hilbert's name, because the university would not appoint a
woman. When the faculty resisted her habilitation, Hilbert lost patience: "We
are a university," he told them, "not a bathhouse." The theorem that
came out of that period became one of the load-bearing beams of modern
physics. It is a reminder that the edifice was built by many hands, including
people whose institutions did not always know how to recognize them.

**Noether's Theorem**: Every continuous symmetry of the action corresponds to a conserved quantity.

Time-translation symmetry gives
conservation of energy. Space-translation symmetry gives conservation of
momentum. Rotation symmetry gives conservation of angular momentum. Gauge
symmetry gives conservation of charge.

Conservation laws are geometric consequences of symmetry, with nothing arbitrary about them.

This is where physics stops looking like a cabinet full of separate rules. Energy conservation, momentum conservation, and charge conservation are the quantities that remain fixed when the same action can be read from shifted, rotated, or phase-twisted points of view.

Once that connection lands, symmetry stops being decorative. It becomes the reason repeated measurements made by different observers can be stitched into one account without inventing conservation laws by decree.

Symmetries are connected to the deepest physical laws. The "stuff" of physics (energy, momentum, charge) is really just "geometry" (symmetry). If symmetry were optional, conservation would be optional. But conservation laws are among the most precisely tested facts in all of science.

## 12.3 The First-Principles Reframing: Symmetries Are Consistency Requirements

The deeper question is why symmetry keeps showing up as law, not decoration.

### Symmetry Enables Agreement

Recall our thesis: reality is the process of making observations consistent between observers.

Consider two astronomers observing the same galaxy. One measures energy in her reference frame. The other measures energy in his frame, moving at a different velocity. Their numbers are different.

They are compatible because they are related by a Lorentz transformation. On
the certified screen branch, this symmetry is the conformal motion of the
refined sphere. Independently supplied modular flow must match that motion and
its normalization. The Lorentz transformation tells the astronomers how to
translate between their observations.

**Symmetry is the grammar of consistency.** Without symmetry, different observers could not compare notes. Their measurements would be incommensurable.

### The Overlap Algebra

Observers have patches with algebras of observables. When patches overlap, they
must agree on the overlap region.

Conservation laws are the simplest form of this agreement. If I measure total energy in my region and you measure total energy in your region, and our regions overlap, then we must agree on the energy in the overlap, because energy is conserved.

**Symmetry provides the translation manual that makes different viewpoints compatible.**

## 12.4 Why Symmetry Lives on the Screen

In the symmetric screen chart, an observer-accessible cut is represented by
the holographic screen $S^2$. The natural angular symmetry group is **SO(3)**.

SO(3) is the group of ordinary rotations in three-dimensional space. Calling it
a group only means that rotations can be composed, undone, and compared in a
consistent way.

This has immediate consequences. Whatever physics lives on the screen must organize itself into **representations** of SO(3), which specify how fields transform under rotations.

The representations are labeled by angular momentum $l=0,1,2,\ldots$. The
scalar mode $l=0$ stays unchanged under rotation. The vector mode $l=1$
transforms like an arrow and carries three components. The tensor mode $l=2$
transforms like a stress matrix and carries five.

This explains part of the angular-momentum structure: fields on the sphere
decompose into discrete angular modes because spherical harmonics are labeled
by integers. Intrinsic spin is a separate representation-theoretic input, which
for fermions enters through the spinor and double-cover structure discussed
next.

## 12.5 The Spinor Mystery

But electrons have spin 1/2. There's no l = 1/2 representation of SO(3).

If you rotate an electron by 360 degrees, it doesn't return to its original state. It picks up a minus sign. You must rotate by 720 degrees to get back.

### The Double Cover

Electrons transform under **SU(2)**, the double cover of SO(3). Every rotation in SO(3) corresponds to two elements in SU(2), differing by a sign.

Objects transforming under SU(2) are called **spinors**. They have half-integer spin.

### The Dirac Belt Trick

Human shoulders are the wrong prop for this lesson.

Use a belt, a ribbon, or a plate attached to ribbons. Hold one end fixed and
rotate the other end by 360 degrees. The belt carries a twist that cannot be
removed while both ends keep the same orientation. Rotate the end by another
360 degrees, and the doubled twist can be combed away without rotating the end
again.

The belt is not a spinor. It shows the topology spinors obey. A 360-degree
rotation lands on the other lift in the double cover. A 720-degree rotation
returns to the starting lift.

### Why Half-Integers Exist

Quantum mechanics allows **projective representations**. Physical states are rays in Hilbert space, vectors defined only up to an overall phase. This phase freedom permits the double cover SU(2).

A ray is a direction, not one particular arrow. Multiplying a quantum
state by an overall phase changes the vector but not the physical state. That
small freedom is what lets spinors carry the minus sign after a 360-degree
rotation without changing observable probabilities.

Half-integer-spin matter sectors become possible because quantum mechanics allows projective representations of the screen's symmetry group.

## 12.6 Wigner's Classification

In 1939, Eugene Wigner classified all possible elementary particles.

A particle is a representation of the Poincare group, the symmetry group of special relativity.

The Poincare group collects the basic moves that leave special relativity
unchanged: translations in space and time, rotations, and Lorentz boosts between
observers moving at constant relative velocity.

Irreducible representations are labeled by two numbers only: mass $m$, which
is continuous and non-negative, and spin $s$, which comes in the familiar
discrete ladder $0, 1/2, 1, 3/2, 2, \ldots$.

Mass and spin are the only quantum numbers that follow from spacetime symmetry.
The complete catalogue entry for an elementary particle is two numbers. As
catalogue entries go, this is admirably brief.

**Particles are representations of symmetries.** Spacetime symmetry fixes the mass-and-spin labels, while the realized internal charges and matter content require additional structure.

This changes what a particle is. An allowed transformation pattern replaces
the image of a tiny marble with a fixed identity tag. Mass tells you how the
excitation sits with time translations. Spin tells you how it sits with
rotations.

The spare label set matters. Once the symmetry group is fixed, only a limited menu of stable transformation patterns is left. The particle table starts to look less like a box of arbitrary ingredients and more like a list of admissible roles.

## 12.7 The Standard Model Gauge Groups

One usually writes the Standard Model gauge group as:

$$G_{SM} = SU(3) \times SU(2) \times U(1)$$

Start with the word **group**. A group is a set of moves that can be followed
by other moves, can be undone, and includes a do-nothing move. Rotations form a
group: rotate a cup, rotate it again, and the result remains a rotation. Every
rotation has an inverse rotation that takes you back.

In physics, the moves are often not visible rotations of ordinary objects. They
are transformations of fields. If the allowed transformations form a group,
then physicists can ask how particles, forces, and charges respond to those
moves.

The letters name the kind of transformation:

**U** means **unitary**. A unitary transformation preserves quantum
probabilities, the way an ordinary rotation preserves length.

**S** means **special**. For these matrix groups, it removes a shared overall
phase. Mathematically this is the determinant-one condition. In practical
terms, it leaves the nontrivial internal rotation while discarding a redundant
common twist.

The number says how many complex components the transformation acts on.
$U(1)$ acts on one complex phase. It is a circle of possible phase rotations,
written $e^{i\theta}$. $SU(2)$ acts on two-component objects, which is why it
naturally organizes weak doublets. $SU(3)$ acts on three-component objects,
which is why it naturally organizes the three color labels of quarks.

$U(1)$ is abelian, which means the order of two transformations does not matter.
$SU(2)$ and $SU(3)$ are non-abelian, which means the order can matter. That is
why the weak and strong interactions have richer self-interactions than plain
electromagnetism.

The multiplication sign names three independent internal transformation systems
running side by side, a tensor-product bookkeeping rather than ordinary
numerical multiplication. A particle can carry color, weak isospin, and
hypercharge labels at the same time.

A **gauge group** is a group of transformations that change the mathematical
description without changing the physical situation. The word gauge adds one
extra feature: the transformation can be chosen locally. Different observers,
or different points in spacetime, may use different internal bookkeeping
choices, and the predictions must agree. Gauge fields are what make that
local agreement possible.

$SU(3)$ carries the strong-force color bookkeeping. $SU(2)$ carries the weak
interaction before symmetry breaking. $U(1)$ carries hypercharge and later
feeds electromagnetism through its mixing with $SU(2)$.

Where do these internal symmetries come from?

The useful picture is practical. These groups are the accounting systems that specify which transformations count as physically equivalent in the strong, weak, and electromagnetic sectors. The real question of the chapter, and of the later Standard Model chapter, is why nature settles on exactly this trio.

### Extra Dimensions

Maybe the screen is $S^2 \times K$, where K is a tiny internal manifold.

If K is a circle, you get U(1). If K is more complex (like a Calabi-Yau space), you can get non-Abelian groups like SU(3).

### Boundary Currents

AdS/CFT provides another route. If the boundary theory has a global symmetry, the bulk has a corresponding gauge field.

*Global symmetry on boundary corresponds to gauge symmetry in bulk.*

A conserved current on the screen creates a gauge boson in the bulk.

### Our Route: Gauge Group from Gluing

In this book we take a different route. The gauge group is not assumed in
advance. We look at what happens when observer patches are glued together. The
edge charges seen at a fixed resolution combine according to a fusion table, a
finite chart of which charges can merge into which, and as the resolution is
refined that fusion data carries through consistently. A
reconstruction theorem then works backward from the surviving charge data to the
symmetry group behind it.

This is classification before selection. Transport and refinement reconstruct
a compact group under their stated hypotheses. The economy rule
then selects the smallest admissible one-Higgs matter package, including
three colors and a conditional minimum of three
families. Attaching that family count to physical chiral matter is a separate
construction. On the selected matter package, six shared center
transformations act trivially. Counting those duplicates once gives
$S(U(3)\times U(2))\cong(SU(3)\times SU(2)\times U(1))/\mathbb Z_6$.
Here $\mathbb Z_6$ is a six-step cyclic overlap between the three factors;
the quotient says those six shared moves are one move.

The notation looks forbidding, but the roles are practical. $SU(3)$ is the
color accounting system for quarks. $SU(2)$ is the weak doublet accounting
system. $U(1)$ is the single continuous charge direction that feeds ordinary
electromagnetism after symmetry breaking.

### The Icosahedral Finite Symmetry

The twelve-port screen supplies an independent route to the same Lie-algebra
shape. The proper rotations of an icosahedron form $A_5$, a group with sixty
elements. It is the same abstract group as the even permutations of five
objects. In the carrier picture it acts on twelve vertices, twenty faces, and
thirty edges. Their stabilizers have orders five, three, and two. One small
group therefore organizes the ports, the triangular overlap neighborhoods, and
the edge-sector and boundary bookkeeping without choosing a preferred
direction.

In the declared echosahedral lineage, the twelve ports arrange themselves
under the construction rules; the antipodal pairing and the icosahedral
symmetry come out rather than being put in. The result is scoped to that
architecture.

The twelve real port readings split under $A_5$ into four bundles:

$$
12=1+3+3'+5.
$$

The $1$ is the uniform reading. The $5$ is the traceless quadrupole-like part.
The two threes are inequivalent icosahedral vector bundles. Because none of
these pieces repeats, an exactly $A_5$-invariant linear response cannot mix
them arbitrarily. It acts as one scalar on each irreducible block. The
symmetry therefore fixes degeneracy patterns while leaving the actual
frequencies and amplitudes to the dynamics.

The same rigidity appears in angular harmonics. The icosahedral vertices form
a spherical five-design, so angular averages through degree five reproduce the
round-sphere averages. The five quadrupole modes stay together, and the first
nonconstant $A_5$-invariant spherical harmonic appears at degree six. This is
why the icosahedron is useful as a finite screen: it hides the preferred axes
from low-order angular probes unusually well. A realized state can be
aligned or anisotropic. The symmetry constrains the allowed block structure,
not the state chosen inside it.

The golden ratio appears here as geometry. Icosahedral coordinates contain
$\sqrt5$ and $\varphi=(1+\sqrt5)/2$, and the adjacency spectrum contains
$\pm\sqrt5$. This is the shape arithmetic of the port frame. The later
golden-ratio balance coordinate plays a separate self-similar role in the
quantitative construction.

With outward orientation and a compact noncommuting current law, the same
$1+3+3'+5$ coefficient space carries
$\mathfrak u(1)\oplus\mathfrak{su}(2)\oplus\mathfrak{su}(3)$. Chapter 14 shows
how the bracket works. The finite register decomposition does not supply that
law. The result concerns coefficient and current directions; the public record
projectors themselves commute. A physically constructed, refinement-stable current
realization is the separate bridge into the gauge sector.

$A_5$ is a finite regulator symmetry, neither the Standard Model gauge group
nor the Lorentz group; continuous Lorentz symmetry and the physical matter
package arrive by the separate routes already described. Tying the antipodal
$\mathbb Z_6$ residue to its action on matter is work in progress.

### Yang-Mills and the Gap

Once a compact gauge group is reconstructed, the theory has to select the field
theory carried by that group. In the OPH route, compact gauge data curve around
small loops, the four-dimensional scaling limit makes that curvature smooth,
and the local maximum-entropy rule leaves the curvature-squared term as the
relevant pure-gauge action:

$$
S_E[A]=\frac{1}{4g^2}\int_{\mathbb R^4}\langle F_{\mu\nu},F_{\mu\nu}\rangle\,d^4x.
$$

The field strength $F$ measures curvature of the gauge connection. The formula
says that the action is built from curvature squared, integrated over
four-dimensional Euclidean space. In OPH this is the continuum form of
compact-gauge patch bookkeeping.

The mass-gap proposal has a simple image. At finite resolution, local repair
removes inconsistencies from the gauge data and returns the system toward its
repaired vacuum. The Yang-Mills gap is then the repair gap: the least
nontrivial energy needed to leave the repaired vacuum survives as the least
nontrivial gauge excitation in the continuum,

$$
\Delta_{\mathrm{YM}}=\Delta_{\mathrm{rep}}.
$$

Behind the image sits an exact finite theorem with a separate physical gate.
For repair to enforce a common energy floor, every allowed kind of local
repair must retain some minimum strength, and neighboring repairs must not
frustrate one another more and more as the screen grows or its boundary
changes. Once the finite source lists every boundary-layer type and prints
one uniform influence bound below the instability threshold, the theorem
gives one positive gap that works at every location, size, boundary, and
refinement stage. Locality or mixing alone cannot do this; small finite
counterexamples show why both controls matter. The actual compact-gauge
source still has to produce that list and those bounds, and the finite
dynamics still has to survive as a nontrivial four-dimensional quantum
theory, before the accounting becomes literal.

A single finite Ising test case calibrates the arithmetic of that argument;
it does not stand in for the compact-gauge construction itself.

## 12.8 Symmetry Breaking

The universe has beautiful symmetries. But the symmetries are also hidden.

In the minimal electroweak action and its ordinary vacuum, the electromagnetic
quadratic field has no hard mass term while the W and Z fields do. Why?

### The Mexican Hat

The Higgs potential:

$$V(\phi) = -\mu^2 |\phi|^2 + \lambda |\phi|^4$$

has rotational symmetry. But the minimum is in a circular valley, not at the center.

$V$ is potential energy for the Higgs field $\phi$. The parameter $\mu^2$ sets
the scale of the unstable central point, and $\lambda$ controls how steeply the
potential rises at large field values. The squared magnitude $|\phi|^2$ says
that the energy depends on distance from the origin in field space, not on a
particular direction. That is why the equation is symmetric even though the
chosen ground state is not.

The system picks a point in the valley. The symmetry is **spontaneously broken**. The equations are symmetric; the state is not.

### The Higgs Mechanism

When the Higgs field settles to a non-zero value, the would-be Goldstone modes
are absorbed by the gauge bosons, the $W$ and $Z$ become massive, the Higgs
boson remains as the physical excitation, and fermion masses are fed through
their Higgs couplings. The underlying symmetry $SU(2)\times U(1)$ narrows to
$U(1)_{\mathrm{em}}$.

"Absorbed" is physicists' shorthand. The would-be massless Goldstone degrees of
freedom do not appear as separate particles. They become the extra polarization
states needed by massive $W$ and $Z$ bosons.

In the quantized electroweak theory, the surviving
$U(1)_{\mathrm{em}}$ direction carries the Maxwell kinetic term in the
deconfined vacuum. The photon shows up as an exactly massless messenger,
while the broken directions carry the massive $W$ and $Z$ excitations.

Symmetry breaking corresponds to the screen "freezing" into a specific
configuration. We live in a frozen shard of a more symmetric world.

This section describes the ordinary electroweak mechanism with its gauge,
scalar, source, and quantum premises taken as given; building those
ingredients from observer-like sources inside OPH is work in progress.

## 12.9 CPT: The Unbreakable Symmetry

Most symmetries can be broken. But one cannot: **CPT**.

$C$ swaps particles with antiparticles. $P$ reflects the world in a mirror.
$T$ runs the movie backward.

The **CPT theorem**: Any Lorentz-invariant local quantum field theory is invariant under CPT.

You can break C, P, T, CP, CT, PT individually. But if you apply all three together, physics must look the same.

The consequences are famously sharp. Every particle has an antiparticle with
exactly the same mass, and particle and antiparticle lifetimes are identical.

On the screen, CPT combines charge conjugation with orientation reversal and
reversal of the record flow. In the Lorentzian field-theory limit this becomes
the usual combined charge, parity, and time-reversal symmetry.

CPT is the immune system of reality, a consistency check that can never be bypassed.

## 12.10 Noether's Theorem: The Calculation

Consider a field theory with action:

$$S = \int d^4x \, \mathcal{L}(\phi, \partial_\mu\phi)$$

Here $S$ is the action, the quantity whose stationary points give the field
equations. The integral $\int d^4x$ means "add contributions over spacetime."
$\mathcal{L}$ is the Lagrangian density, a local rule built from the field
$\phi$ and its spacetime derivatives $\partial_\mu\phi$. The Greek index $\mu$
labels spacetime directions.

Suppose an infinitesimal change $\phi \to \phi + \epsilon\,\delta\phi$ leaves
the action unchanged. Then a current is conserved:

$$\partial_\mu J^\mu = 0$$

where the conserved current is:

$$J^\mu = \frac{\partial\mathcal{L}}{\partial(\partial_\mu\phi)}\delta\phi$$

$J^\mu$ is the current associated with the symmetry. The equation
$\partial_\mu J^\mu=0$ is a continuity equation: what flows out of one region
must enter another. The variation $\delta\phi$ is the infinitesimal change of
the field under the symmetry, and the factor
$\partial\mathcal{L}/\partial(\partial_\mu\phi)$ measures how sensitively the
local rule responds to changes in the field's rate of change. Noether's
theorem says that if changing the field this way leaves the action fixed, a
current must be conserved.

For time translation, $\delta\phi = \partial_t\phi$. The conserved current is energy density.

For space translation, $\delta\phi = \partial_i\phi$. The conserved current is momentum density.

Together, these form the **stress-energy tensor**:

$$T^{\mu\nu} = \frac{\partial\mathcal{L}}{\partial(\partial_\mu\phi)}\partial^\nu\phi - \eta^{\mu\nu}\mathcal{L}$$

$T^{\mu\nu}$ is the stress-energy tensor. It records energy density, momentum
density, pressure, and stress in one object. The symbol $\eta^{\mu\nu}$ is the
flat-spacetime metric used in special relativity. This is the compact form of
the sentence that symmetry under spacetime translations gives conservation of
energy and momentum.

This is the precise sense in which conserved "stuff" (energy, momentum) is tied to symmetry.

The calculation earns its keep here. It shows that a conservation law is not an
extra commandment stapled onto the theory after the fact. The conserved current
is the shadow cast by an allowed infinitesimal transformation. If the action
does not change when you slide in time, rotate, or shift phase, a current must
exist whose flow is preserved. The chapter therefore treats symmetry as
operational structure, not decoration.

![Noether's theorem turns an allowed transformation that leaves the action fixed into a conserved current.](../assets/book_diagrams/symmetry-to-conservation.svg){width=80%}

## 12.11 What Symmetry Predicts

Symmetry earns its place in physics because it leaves hard fingerprints.
Noether ties symmetry to conservation. The sphere carries rotational
structure, angular harmonics, and the conformal bridge to Lorentz symmetry.
Spinors fit naturally on that sphere. Wigner shows that once relativity is in
place, particles are classified by how they transform.

In the detailed construction, this bridge is explicit. Sky directions, cap
boundaries, and the two sides of a cap all transform together under a Lorentz
motion. The book-level point is that the same symmetry moves the screen view and
the recovered spacetime view.

The observed world obeys the script. Conservation laws hold, CPT remains
intact, and physical quantum fields obey spin-statistics. OPH has to
reproduce that last connection; its spin, exchange, and quantum-field
construction is work in progress. Symmetry is one of the mechanisms by which
physics keeps itself coherent.

### Noether's Human Lesson

Noether arrived at her theorem through a problem that looked technical, the
status of conservation laws in general relativity, and left physics holding
one of its central pillars.

The theorem's lesson is simple enough to say without the machinery: if a
physical description can be changed in a certain way without changing the
action, then something must be conserved. The pairings listed at the start of
this chapter (time and energy, space and momentum, rotation and angular
momentum, gauge and charge) all instantiate the one rule. Each conserved
quantity is a public invariant, something different observers can carry
through their calculations without losing agreement.

This turns symmetry from beauty into bookkeeping. A symmetry is a rule for
translating descriptions while preserving what can be checked. In OPH language,
it is a compatibility rule
for patches. Two observers may use different coordinates, phases, frames, or
local bases. If their translation rule is a true symmetry, they agree
on the shared content.

Section 12.7 unpacked the Standard Model gauge group
$SU(3)\times SU(2)\times U(1)$ factor by factor. Read through Noether's
theorem, each factor is a translation rule with a conserved charge attached,
and the quotient by a shared center, which returns in Chapter 14, records
which transformations that look separate are actually identified on all
physical states.

Conservation laws are the durable threads that let observers compare the same
physical story from different cuts.

## 12.12 Reverse Engineering Summary

The old intuition treated symmetry as a kind of cosmic taste. The deeper
picture is harsher. Symmetry is the translation manual that lets different
observers describe one world without contradiction. Rotational symmetry keeps
descriptions compatible across direction. Time-translation symmetry keeps them
compatible across repeated comparison. Gauge symmetry keeps them compatible
across local descriptions of charge. Conservation laws are the public record
of that agreement.

Consistency has dictated grammar. What it says about size is the next
question.

---

The translation rules need an arena that carries them. Our universe expands,
accelerates, and hides information behind a cosmological horizon. The arena
itself has thermodynamics.

**Chapter 13: The de Sitter Patch** examines that question.
