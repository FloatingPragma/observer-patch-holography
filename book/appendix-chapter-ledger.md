# Appendix: Chapter-by-Chapter Field Notes and Symbol Guide

This appendix is a companion trail through the book. It has one purpose: make
the technical vocabulary and the human lineage harder to miss. The main
chapters tell the argument in narrative order. This appendix walks the same
path again, chapter by chapter, and asks three questions each time. What is
the chapter trying to reverse engineer? Which symbols or equations should the
reader keep steady? Which communities of builders made that chapter possible?

The guide is deliberately plain. It is a reader's workbench for the chapters.
If a symbol appears
earlier than its full technical life in the book, the note here gives the
plain-language meaning. If a discovery is associated with one famous name, the
note widens the frame enough to show the relay team behind it. OPH itself is
presented in the book as a synthesis, and a synthesis has to honor the tools
it inherits.

## Prologue: The Reverse-Engineering Posture

The prologue sets the method before it sets the physics. A reverse engineer
does not begin with the source code. A reverse engineer begins with behavior:
outputs, crashes, invariants, strange edge cases, and repeated symptoms. The
prologue applies that habit to reality. It says that modern physics has
produced a pile of symptoms that do not fit the naive story of objects sitting
in an observer-independent container. Relativity denies a preferred frame.
Quantum mechanics denies a complete spreadsheet of pre-existing measurement
answers. Black-hole thermodynamics denies naive volume counting. Holography
suggests that boundary information can encode bulk physics. OPH then asks
what architecture would make those symptoms natural.

The symbols in the prologue are light, but the conceptual terms are important.
An **observer patch** is the finite portion of the world to which an observer
has operational access. A **record** is physical information that can be
checked later or by someone else. An **overlap** is a shared region, shared
interface, or shared algebra on which two descriptions can be compared. A
**fixed point** is a stable result of an update or repair process: apply the
consistency rule again and nothing relevant changes. A **consensus process**
is a repair dynamic rather than a vote about taste: incompatible descriptions
are rejected or adjusted until the overlap-visible records agree.

The prologue also introduces the local pixel ratio $P$. In the detailed OPH
program, $P$ is the small dimensionless ratio that sets the effective local
grain of the screen description. It is the local fixed point at which the outside screen geometry and the inside
electromagnetic readout are required to match. The theory defines $P$ by that
matching equation and feeds it no measured number. The working value used by
the chains is located from the measured fine-structure constant under the
constant-identification axiom, and the requirement that the two readings match
exactly is the open closure test. Later chapters ask how far
that same fixed point can travel through the weak sector, the electromagnetic
coupling, the Higgs-top surface, quark structure, neutrinos, and the
gravity-facing side of the framework.

The human lineage behind the prologue is broader than any one field. Reverse
engineering belongs to computing and security culture, but physics has always
had a similar discipline. Galileo inferred laws from falling bodies and
projectiles. Newton inferred a universal gravitational structure from motions
on Earth and in the sky. Maxwell inferred fields from electric and magnetic
regularities. Einstein inferred new spacetime structure from clocks, rods,
light, and the failure of ether reasoning. Quantum theory inferred its rules
from spectra, heat radiation, scattering, and detector clicks. OPH inherits
that posture. The theory is the architecture that makes those symptoms belong
to one observer-consistency system.

The removed introductory paragraph used market-like labels, and the book is
stronger without it. The prologue should not sound as if it is optimizing for
search phrases. It should sound like a working physicist and programmer
looking at a system whose behavior cannot be explained by the official mental
model.

## Chapter 1: Consistency First

Chapter 1 names the central inversion: objectivity is not assumed first and
explained later. Objectivity is reconstructed from consistency across finite
perspectives. The chapter's most important question is therefore simple:
what must be true if many bounded observers can nevertheless inhabit one
shared public world?

The chapter uses the word **observer** in an operational sense. An observer is
not necessarily a human mind. It is any physical system capable of receiving
signals, forming records, updating internal state, and participating in a
network of comparisons. A detector, a laboratory, a biological nervous
system, and a future engineered patch-federation can all instantiate parts of
that role at different scales. The aim is to identify the formal surface they
share: bounded access,
record-making, and overlap comparison.

The **normal form** language matters. In computation, a normal form is a
canonical representative reached after applying rewriting rules. If two
different paths of repair land in the same normal form, the system has a
strong kind of coherence. Chapter 1 uses that image for physical reality.
Local descriptions may begin mismatched. The public world is what remains
after admissible repairs remove overlap-visible disagreement. A **Lyapunov
function**, introduced more formally elsewhere in OPH, is a quantity that
decreases along accepted repair steps. It proves that the repair process is
not wandering forever when the state space is finite. It does not by itself
prove that every repair order lands in the same result. That stronger claim is
confluence, and OPH treats it as a separate local-diamond plus repair-completeness
condition on the observer-facing quotient.
The theorem is relative to the fixed initial quotient state. A stronger
"same boundary data gives the same answer" statement also needs a unique
consistent extension in that boundary or sector fiber. The generic equivalence
between this boundary-identification property and cross-source endpoint
agreement, together with its Lean bridge to the finite repair architecture, is
proved in the companion paper *Observation-Determined Normal Forms*. Injectivity
of the declared physical boundary map on the full consistent quotient remains an
application-specific premise.

The five axiom groups are not meant as arbitrary declarations. They bundle
lessons from several mature fields. Finite screen capacity comes from
black-hole thermodynamics and holography. Local algebras come from quantum
theory and algebraic quantum field theory. Overlap gluing comes from sheaf
logic, quantum marginal problems, and operational comparison. Generalized
entropy and recoverability come from semiclassical gravity and quantum
information. The economy principle belongs to the selection part of the
program, where the admissible low-energy sectors are narrowed.

The chapter's diagram, the consensus funnel, is a repair picture. Many local
descriptions enter. Some mismatch. Some can be repaired.
Some are rejected. The final public structure is not the arithmetic average
of opinions. It is the stable normal form that survives all allowed overlap
checks.

The human chain is long. Bell, Kochen, Specker, Gleason, von Neumann,
Everett, Bohr, Wigner, Haag, Kastler, Haag again through algebraic QFT,
Wheeler, Zurek, Bekenstein, Hawking, 't Hooft, Susskind, Maldacena, Ryu,
Takayanagi, Preskill, Hayden, and many others contribute pieces. OPH's first
chapter asks readers to see those pieces as symptoms of one architecture.

## Chapter 2: The Original Hackers

Chapter 2 is the philosophical prehistory of the same architecture. Careful
thinkers found structural weaknesses in the naive picture long before the
laboratory made those weaknesses unavoidable; the chapter does not turn ancient
philosophers into disguised modern physicists.

Plato's cave supplies the projection motif. A lower-dimensional record can
carry information about a source the prisoners never see directly. The
relevant modern equation is

$$S_{max}=\frac{A}{4\ell_P^2}.$$

Here $S_{max}$ is maximum entropy or maximum information capacity, $A$ is the
boundary area, and $\ell_P$ is the Planck length. The formula does not turn
Plato into a physicist. It shows that the projection intuition later acquired
a quantitative gravitational form: capacity scales with boundary area.

Zeno supplies pressure against naive continuity. The paradoxes do not prove a
Planck lattice. They show that infinite divisibility is not innocent. Modern
physics adds several reasons to distrust the naive continuum as fundamental:
finite horizon entropy, ultraviolet divergences, Planck-scale dimensional
analysis, and quantum discreteness in many observables.

The Skeptics supply context. A property detached from the conditions of
observation is not the kind of thing quantum mechanics later lets us keep.
Bohr's complementarity, the Kochen-Specker theorem, and contextuality results
make that lesson technical. Descartes supplies the unavoidable observer.
Kant supplies the idea that space might be a form of representation rather
than a pre-built box. Godel and Hofstadter supply self-reference.

The chapter's new cave diagram is a reminder that projection is not enough.
The shadows become a world only because observers compare and reconstruct. A
single shadow is ambiguous. Several partial records, cross-checked through
overlap, can begin to constrain a shared model.

The human chain predates physics departments. It includes oral argument,
geometry, astronomy, logic, theology, optics, mechanics, and experimental
science. A modern reader should not treat philosophy as a failed
version of physics. Philosophy was the first debugging environment for
assumptions beyond the instruments of the time.

## Chapter 3: The Holographic Screen

Chapter 3 turns the philosophical projection hint into a physical storage
question. If gravitational systems have finite information capacity and that
capacity scales like area, what observer-facing chart exposes that capacity?
The chapter's answer is the holographic screen, especially the spherical screen
naturally associated with symmetric light-cone and horizon constructions.

The most important symbol is $S\leq A/(4\ell_P^2)$. The inequality form says
that a region's entropy is bounded by its boundary area in Planck units. $S$
is entropy. $A$ is area. $\ell_P^2$ is Planck area. If the book later writes
the same idea with $G$, $\hbar$, $c$, and $k_B$ restored, it is the same
physics with units made explicit. Natural units suppress those constants so
the structure can be seen.

The symbol $S^2$ means the two-sphere, the ordinary surface of a ball
abstracted as a mathematical object. It is two-dimensional because two
coordinates, like latitude and longitude, label points on it. In OPH this
sphere is a chart for observer-facing support data, not a literal ball that
stores the universe on its surface. The standard round metric on $S^2$ is not
the whole OPH theory. It is the symmetric starting geometry used to organize
patches, caps, and overlaps.

The chapter's screen diagram shows two observer patches and the
lens-shaped overlap where their descriptions can be compared. That overlap is
the operational hinge of the book, not decorative shading. There is one
shared screen net idealization, many observer-facing access cuts, and many
finite patches. A private patch can carry more than the shared lens. Public
reality is built from what many such lenses force to agree.

The historical builders include Bekenstein, who argued that black holes have
entropy, Hawking, who found black-hole radiation and fixed the temperature,
't Hooft, who framed the holographic principle, Susskind, who sharpened it,
and Maldacena, who gave the AdS/CFT duality as a controlled realization. But
the screen idea also rests on earlier work: Riemannian geometry, causal
structure in relativity, quantum field theory, statistical mechanics, and
information theory. The screen is a meeting point, not an isolated invention.

The practical lesson is this: once capacity is finite, continuum language
becomes effective language. Smooth geometry remains useful, just as fluid
dynamics is useful although water is molecular. The screen chapter asks what
kind of microscopic or finite-resolution bookkeeping could make smooth space
appear.

## Chapter 4: Entropy on the Edge

Chapter 4 explains why records have a direction. The basic puzzle is that many
microscopic laws can be run backward, while ordinary life cannot. The answer
is not a new force called time's arrow. The answer is a low-entropy record
boundary plus overwhelming counting.

Carnot's formula

$$\eta_{max}=1-\frac{T_{cold}}{T_{hot}}$$

uses $\eta$ for efficiency. $T_{hot}$ and $T_{cold}$ are absolute
temperatures. The formula says that useful work can be extracted from a
temperature difference, and that no engine can beat the ratio. Boltzmann's
formula

$$S=k_B\ln W$$

uses $S$ for entropy, $k_B$ for Boltzmann's constant, and $W$ for the number
of microstates compatible with a macrostate. The logarithm turns
multiplication of independent possibilities into addition of entropies.
Shannon's entropy

$$H=-\sum_i p_i\log_2 p_i$$

uses $p_i$ for probabilities of possible messages or outcomes. It measures
uncertainty in bits when the logarithm is base two.

The chapter's entropy-arrow diagram shows the physical price of memory. A
record does not float above thermodynamics. It consumes low-entropy resources
and exports heat or waste entropy. Landauer's principle, $k_BT\ln 2$ per
erased bit at temperature $T$, is the cleanest expression of that price.

The human chain is industrial as much as theoretical. Steam engines forced
Carnot and Clausius to think about heat. Boltzmann and Gibbs built the
statistical picture. Maxwell's demon exposed the information problem.
Szilard, Shannon, Brillouin, Landauer, Bennett, and many later workers tied
information to physical cost. Bekenstein and Hawking moved the same logic to
horizons. OPH needs all of that because a consensus world requires records,
and records require entropy gradients.

## Chapter 5: The Algebra of Questions

Chapter 5 asks what kind of mathematical object a measurement question is.
The classical answer is too simple: all properties sit together and can be
read in any order. Quantum mechanics says no. Some questions do not commute.

The core notation is

$$[X,P]=XP-PX=i\hbar.$$

$X$ is the position operator. $P$ is the momentum operator. $XP$ and $PX$ are
two different products, corresponding to two different orders of applying the
operators. The bracket $[X,P]$ is the commutator, the algebraic measure of
order-dependence. $i$ is the imaginary unit. $\hbar$ is Planck's constant
divided by $2\pi$. The uncertainty relation

$$\Delta X\,\Delta P\geq \frac{\hbar}{2}$$

uses $\Delta X$ and $\Delta P$ for spreads of possible outcomes in a state.
The statement concerns the state and the algebra, rather than poor
instruments.

The chapter's question-order diagram is meant to remove the mystery from the
notation. Asking A then B can be physically different from asking B then A.
That is what non-commutativity means.

The modular-flow formula

$$\sigma_t(A)=\Delta^{it}A\Delta^{-it}$$

uses $A$ for an observable, $t$ for a flow parameter, and $\Delta$ for the
modular operator associated with an algebra-state pair. The KMS condition is
the thermal-equilibrium signature of that flow. This prepares the reader for
the later claim that time-like structure can be read internally from a
restricted state.

The builders are Planck, Einstein, Bohr, Heisenberg, Born, Jordan,
Schrodinger, Dirac, von Neumann, Stone, Gelfand, Segal, Haag, Kastler,
Tomita, Takesaki, Connes, Rovelli, and many more. The algebra of questions is
what survived when the old orbit picture failed.

## Chapter 6: Overlap Consistency

Chapter 6 makes objectivity operational. Two observers do not need identical
private descriptions. They need matching predictions on the observables both
can access. In notation, if observers $i$ and $j$ assign states $\omega_i$
and $\omega_j$, then agreement on the shared algebra is written

$$\omega_i|_{\mathcal A(P_i\cap P_j)}
=
\omega_j|_{\mathcal A(P_i\cap P_j)}.$$

$P_i$ and $P_j$ are patches. $P_i\cap P_j$ is their intersection. The symbol
$\mathcal A$ denotes an algebra of observables. The vertical bar means
"restrict to." In plain language: ignore the private parts and compare only
the shared questions.

The chapter also uses Bell's CHSH expression. The classical bound is
$|S|\leq 2$, while quantum mechanics permits $|S|\leq 2\sqrt2$. Here $S$ is
a correlation combination built from four measurement settings, not entropy.
The Bell lesson is that overlap agreement cannot be explained by local hidden
instruction sheets.

The reduced-state notation $\rho_A=\mathrm{Tr}_B\rho_{AB}$ says that the
state on subsystem $A$ is obtained by tracing out subsystem $B$. $\rho$ is a
density matrix, a quantum state that may include classical uncertainty and
entanglement. $\mathrm{Tr}_B$ is the partial trace. This is the formal tool
for asking what a local patch sees when the total state contains more than
the patch can access.

The overlap-consistency diagram shows the simplest case: two patches with a
shared region and matching state assignments there. The deeper problem is
that many pairwise overlaps do not automatically glue to one global state.
That is why the quantum marginal problem matters. Compatibility is a theorem
to earn, not an assumption.

The human lineage includes Bell, Clauser, Horne, Shimony, Holt, Aspect,
Zeilinger, Kochen, Specker, Fine, Foulis, Randall, Fawzi, Renner, Petz,
Lieb, Ruskai, and many others. OPH reads their lesson as architecture:
public facts live where local quantum descriptions can be glued without
contradiction.

## Chapter 7: Recovery

Chapter 7 asks why the overlap web does not fall apart. Quantum information
cannot be copied freely, noise is everywhere, and horizons hide data. Yet the
world has durable records. The answer is recovery: information can survive in
encoded correlations even when local access is damaged.

The central information quantity is CMI:

$$I(A:C|B)=S(AB)+S(BC)-S(B)-S(ABC).$$

$A$, $B$, and $C$ are subsystems. $S(AB)$ is the entropy of the joint system
$AB$, and similarly for the other terms. The quantity measures how much
correlation remains between $A$ and $C$ once $B$ is known. If it is zero, the
state has a quantum Markov property: $B$ screens off $A$ from $C$ in the
right sense. If it is small, recovery is approximate, but exact Markov normal
forms require exact zero CMI or a controlled collar replacement limit.

The Fawzi-Renner theorem says, roughly, that small CMI implies the existence of a recovery map. Petz gave an earlier
canonical recovery map in the exact Markov setting. This recovery map supplies
a recovered comparison state, not a dimension-free one-shot projection onto
the exact Markov set. The book does not require the reader to compute the map,
but it does require the conceptual lesson: lost-looking local data can be
reconstructible from surrounding correlations.

The collar-tripartition diagram shows the split used repeatedly in
holographic recovery: a cap, a collar, and an exterior. The collar is the
buffer region that can make given-data independence possible. It is the
boundary information that lets inside and outside fit back together.

The historical chain begins with Shannon's noisy-channel problem, passes
through no-cloning, quantum error correction, strong subadditivity, Petz
recovery, Fawzi-Renner recovery, black-hole information, the Page curve,
Hayden-Preskill decoding, and entanglement wedge reconstruction. It is one of
the clearest examples of the book's larger point: practical communication
engineering, abstract inequality theory, and quantum gravity turn out to be
talking about the same survival problem.

## Chapter 8: Holography

Chapter 8 broadens the screen idea into a holographic reconstruction program.
The old intuition says a volume contains the independent degrees of freedom
of that volume. Black holes say otherwise. Holography says boundary data can
carry the physics of a bulk.

The Bekenstein-Hawking formula is the recurring anchor:

$$S_{BH}=\frac{k_B A}{4\ell_P^2}.$$

$S_{BH}$ is black-hole entropy, $k_B$ is Boltzmann's constant, $A$ is horizon
area, and $\ell_P$ is the Planck length. In units where $k_B=1$, the same
formula looks like $S=A/(4\ell_P^2)$. The Bekenstein bound, often written
$S\leq 2\pi ER/\hbar c$, says that entropy in a region is limited by energy
$E$ and radius $R$. The constants $\hbar$ and $c$ restore the quantum and
relativistic units.

AdS/CFT introduces a sharper dictionary. Boundary operators $\mathcal O$ have
scaling dimensions $\Delta$. Bulk fields $\phi$ approach boundary values that
act as sources for those operators. The GKPW relation schematically says that
a bulk partition function with boundary source equals a generating functional
for boundary correlation functions. The details are technical, but the
meaning is direct: boundary data can compute bulk physics.

OPH does not claim dS/CFT. Chapter 8 is careful about this. Anti-de
Sitter space has a global boundary and negative cosmological constant. Our
universe is closer to de Sitter, with positive cosmological constant and
observer-dependent horizons. OPH takes the boundary-encoding lesson but
rebuilds it around static patches, finite screen capacity, and overlap
agreement, not a single global CFT at infinity.

The human lineage runs through Bekenstein, Hawking, Gibbons, Perry, 't Hooft,
Susskind, Maldacena, Witten, Gubser, Klebanov, Polyakov, Ryu, Takayanagi,
Hubeny, Rangamani, Hamilton, Kabat, Lifschytz, Lowe, and many others. The
chapter's job is to show why boundary-first physics is an essential clue. It
is one of the most successful clues quantum gravity has produced.

## Chapter 9: Entanglement

Chapter 9 explains how a boundary can feel like a bulk. The answer is
entanglement. Correlations are not decoration on top of space. In holographic
settings, the pattern of correlations helps define the geometry.

The RT formula is the chapter's core bridge:

$$S(A)=\frac{\mathrm{Area}(\gamma_A)}{4G_N}.$$

$S(A)$ is entanglement entropy of boundary region $A$. $\gamma_A$ is the
bulk extremal surface anchored to the boundary of $A$. $G_N$ is Newton's
constant. The equation says that entropy of a boundary region is measured by
an area in the emergent bulk.

Bell states give the simplest entanglement example:

$$|\Phi^+\rangle=\frac{|00\rangle+|11\rangle}{\sqrt2}.$$

The kets $|00\rangle$ and $|11\rangle$ are two-qubit basis states. The
factor $1/\sqrt2$ normalizes the state so total probability is one. If one
qubit is considered alone, its reduced density matrix looks maximally mixed,
even though the pair as a whole is pure. Entanglement is joint-state
correlation, not ignorance about one hidden classical value.

The monogamy inequality, written with a tangle $\tau$, says that entanglement
cannot be shared freely with everyone at once. If $A$ is maximally entangled
with $B$, it has no entanglement budget left for $C$. This is one reason
entanglement networks can support locality. Not everything can be equally
near everything else.

The entanglement-wedge diagram in the chapter shows two boundary regions
whose reconstructed bulk wedges overlap. That shared wedge is the geometric
version of observer overlap. If two boundary descriptions reconstruct the
same bulk operator, consistency requires agreement there.

The builders include Einstein, Podolsky, Rosen, Schrodinger, Bell, Clauser,
Aspect, Zeilinger, Bekenstein, Hawking, Maldacena, Ryu, Takayanagi, Van
Raamsdonk, Swingle, Pastawski, Hayden, Preskill, Yoshida, Penington, and many
others. Entanglement began as a paradox and became a construction material.

## Chapter 10: Error Correction

Chapter 10 turns stability into a coding problem. The naive view says
information is protected by isolation or copying. Quantum theory blocks both
as general strategies. Unknown quantum states cannot be copied, and isolated
systems are rare. The deeper strategy is encoded redundancy.

The logical qubit notation

$$|\psi_L\rangle=\alpha|000\rangle+\beta|111\rangle$$

uses $|\psi_L\rangle$ for the encoded logical state. $\alpha$ and $\beta$ are
complex amplitudes. The basis states $|000\rangle$ and $|111\rangle$ are
physical three-qubit states. The code has not made three copies of an unknown
qubit. It has stored one logical state in correlations.

The Knill-Laflamme condition

$$P E_a^\dagger E_bP=\alpha_{ab}P$$

states when errors are correctable. $P$ projects onto the code space. $E_a$
and $E_b$ are possible errors. The dagger is the adjoint. The numbers
$\alpha_{ab}$ form syndrome data. The environment may learn which error
happened, but it must not learn the protected logical state.

The OPH boundary is explicit: a bare overlap graph is only a finite constraint
code. Code distance, min-cut resilience, and corrupted-carrier correction
require the additional code space, logical-operator, error-family, and recovery
certificate.

Landauer's cost $k_BT\ln2$ returns at the end of the chapter. $k_B$ is
Boltzmann's constant, $T$ is temperature, and $\ln2$ appears for one erased
bit. Error correction is therefore not free. Syndrome extraction, decoding,
resetting ancillas, and maintaining records all require physical work.

The error-correction diagram shows a logical state spread across physical
carriers. One carrier can be damaged, yet the logical pattern can survive if
the code structure remains intact. Holographic error correction translates
the same idea to bulk reconstruction from boundary regions.

The human chain includes Wootters and Zurek on no-cloning, Shor and Steane on
early quantum codes, Calderbank, Shor, and Steane on CSS codes, Gottesman on
stabilizer formalism, Knill and Laflamme on the correction condition, Kitaev
on topological codes, Preskill and many experimental communities on fault
tolerance, plus the holographic-code work of Pastawski, Yoshida, Harlow,
Almheiri, Dong, and others. The chapter is where quantum computing and
spacetime stop looking like separate subjects.

## Chapter 11: MaxEnt and the Arrow

Chapter 11 asks where time comes from when no external clock is allowed to
stand outside the universe. The answer is a stack: entropy,
records, restricted access, maximum entropy inference, and modular flow.

The MaxEnt distribution is usually written

$$p_i=\frac{e^{-\beta E_i}}{Z}.$$

$p_i$ is the probability of state $i$. $E_i$ is its energy. $\beta$ is inverse
temperature, often $1/(k_BT)$. $Z$ is the partition function that normalizes
the probabilities so they sum to one. Jaynes taught physicists to read this
as principled inference, not as a heat-bath trick: choose the least
biased distribution compatible with known constraints.

The Wheeler-DeWitt equation is often summarized as

$$\hat H\Psi=0.$$

$\hat H$ is the Hamiltonian constraint operator and $\Psi$ is the state of
the universe. The zero on the right is the source of the "problem of time":
the universal state does not evolve with respect to an external time
parameter in the ordinary Schrodinger way.

Modular theory supplies a different route. The flow
$\sigma_t(A)=\Delta^{it}A\Delta^{-it}$ is internal to an algebra-state pair.
The KMS condition gives it thermal character. The Unruh formula
$T=\hbar a/(2\pi ck_B)$ then shows a deep relation among acceleration
$a$, temperature $T$, light speed $c$, Boltzmann's constant $k_B$, and
Planck's constant $\hbar$.

The builders include Boltzmann, Gibbs, Einstein, Wheeler, DeWitt, Jaynes,
Tomita, Takesaki, Bisognano, Wichmann, Unruh, Connes, Rovelli, and Jacobson.
The chapter uses their work to make one point: time can be read as the
inside ordering of records and restricted states, not as a river outside the
system.

## Chapter 12: Symmetry

Chapter 12 treats symmetry as the translation manual for observers. If two
descriptions differ by an allowed transformation but describe the same
physics, something is preserved. Noether's theorem is the formal bridge.

The action

$$S=\int d^4x\,\mathcal L$$

is a spacetime integral of the Lagrangian density $\mathcal L$. The symbol
$d^4x$ means integrate over four spacetime coordinates. A field is written
$\phi$, and $\partial_\mu\phi$ is its derivative along spacetime direction
$\mu$. If an infinitesimal transformation $\delta\phi$ leaves the action
unchanged, there is a conserved current

$$J^\mu=\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}\delta\phi$$

with conservation equation

$$\partial_\mu J^\mu=0.$$

$J^\mu$ is the current. The equation says the current has no source or sink.
The stress-energy tensor

$$T^{\mu\nu}
=
\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}\partial^\nu\phi
-\eta^{\mu\nu}\mathcal L
$$

packages energy density, momentum density, pressure, and stress. The metric
$\eta^{\mu\nu}$ is the flat spacetime metric.

The gauge group $SU(3)\times SU(2)\times U(1)$ names the Standard Model's
strong, weak, and hypercharge symmetries before electroweak symmetry breaking.
Later quotient structure identifies transformations that act the same on
physical states.

The compact-gauge branch also supplies the Euclidean Yang-Mills action in four
dimensions on its controlled branch, under its declared compact-gauge
branch assumptions. The corresponding mass-gap statement identifies the first
nonzero Yang-Mills energy with the first nonzero repair eigenvalue:

$$\Delta_{\mathrm{YM}}=\Delta_{\mathrm{rep}}.$$

The Noether diagram in the chapter shows the sequence: symmetry, fixed
action, conserved current. The human story centers on Emmy Noether, but it
also includes Hilbert, Klein, Weyl, Wigner, Yang, Mills, Glashow, Salam,
Weinberg, and the experimental communities that turned symmetries into
measured particle physics. The OPH reading is that symmetry preserves the
possibility of shared description across patches.

## Chapter 13: The de Sitter Patch

Chapter 13 gives the cosmological arena. The standard late-time fit is close to
a positive cosmological constant. That gives observer-dependent horizons and a
finite entropy capacity on the de Sitter branch. The relevant
radius is

$$r_{dS}=\sqrt{\frac{3}{\Lambda}}.$$

$r_{dS}$ is the de Sitter horizon radius and $\Lambda$ is the cosmological
constant. The associated temperature is

$$T_{dS}=\frac{\hbar c}{2\pi k_B r_{dS}},$$

where $\hbar$ is Planck's constant divided by $2\pi$, $c$ is the speed of
light, and $k_B$ is Boltzmann's constant. The entropy is

$$S_{dS}=\frac{\pi r_{dS}^2}{\ell_P^2}.$$

The factor of $\pi$ appears because the horizon area is $4\pi r_{dS}^2$ and
the entropy formula divides by $4\ell_P^2$.

The chapter's capacity numbers are enormous but finite, around
$10^{122}$ to $10^{123}$ depending on convention. At the cosmic record fixed
point OPH reads dimensionless cosmological-constant products as global capacity
data: they fix the screen capacity on which finite observer-patch physics is
organized, while the SI curvature scale also uses the selected OPH scale certificate. The
input-free global closure is
$N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}})$, with
$\log|\Omega^{\mathrm{sc}}_N|-N$ as its finite-count representation.

The modular-anomaly branch introduces an effective dark component. The
benchmark acceleration

$$a_0^{(\mathrm{OPH})}
=
\frac{15}{8\pi^2}\frac{c^2}{r_{dS}}
$$

uses the de Sitter radius to set a deep-infrared scale. The proximity to the
empirical MOND scale is the reason the chapter treats galaxy-scale anomalies
as an extension branch outside the recovered-core theorem package.

The builders include de Sitter, Friedmann, Lemaitre, Hubble, Slipher,
Gamow, Penzias, Wilson, Guth, Starobinsky, Linde, Riess, Perlmutter, Schmidt,
Gibbons, Hawking, Bousso, Banks, Fischler, Verlinde, Milgrom, and many
observational teams. Cosmology is the most public of sciences: no one can
move the universe into a laboratory, so agreement depends on many telescopes,
surveys, calibrations, and cross-checks.

## Chapter 14: The Standard Model

Chapter 14 is the longest technical reconstruction chapter because particle
physics is full of structure. The Standard Model bundles a gauge group,
representation assignments, anomaly cancellations, generations, mixing, masses,
symmetry breaking, and measured couplings. OPH asks how much
of that structure can be organized by consistency.

The gauge group is

$$SU(3)\times SU(2)\times U(1).$$

$SU(3)$ is the color symmetry of the strong interaction. $SU(2)$ is weak
isospin. $U(1)$ is hypercharge. Fermions come in representations of this
group. A representation tells how a field transforms under the symmetry.

In OPH, the charge sectors that survive a consistency test around every loop, and
that keep surviving as the description is refined, reconstruct the compact gauge
group they carry. Choosing the smallest single-Higgs matter content that works is
the later step.

Hypercharge is written $Y$. Electric charge is related by

$$Q=T_3+\frac{Y}{2},$$

where $T_3$ is the third component of weak isospin. Anomaly cancellation
means that certain quantum inconsistencies vanish when all fields in a
generation are counted together. The chapter's Tannaka-Krein diagram shows a
deep reconstruction idea: a group can be read from its representation
category, meaning from the way its charged sectors transform and fuse.

The Higgs potential is often written

$$V(H)=-\mu^2|H|^2+\lambda|H|^4.$$

$H$ is the Higgs field, $\mu$ and $\lambda$ are parameters, and the sign
structure makes the symmetric point unstable so electroweak symmetry breaks.
Yukawa couplings connect fermions to the Higgs and generate masses after the
Higgs gets a vacuum expectation value.

The generation-count diagram marks a proved OPH statement: the window begins
at three for intrinsic CP capability and closes above five from weak-sector
ultraviolet consistency. The chapter is careful about which entries are proved,
which are benchmark checks, which use an empirical anchor, and which require
external data.

The builders are too many for a short list, but the relay includes Dirac,
Pauli, Fermi, Yang, Mills, Gell-Mann, Zweig, Glashow, Salam, Weinberg,
Higgs, Englert, Brout, Kibble, 't Hooft, Veltman, Kobayashi, Maskawa,
Cabibbo, Gross, Wilczek, Politzer, Fritzsch, Nambu, Goldstone, Lederman,
Perl, and the enormous detector collaborations at CERN, Fermilab, SLAC,
Brookhaven, KEK, DESY, and elsewhere. The Standard Model is a civilization
scale measurement artifact.

## Chapter 15: Relativity

Chapter 15 reconstructs spacetime behavior from modular and screen
constraints. It moves from light cones and boosts to the branch-qualified
Einstein equation and the Newtonian limit.

The Lorentz factor

$$\gamma=\frac{1}{\sqrt{1-v^2/c^2}}$$

uses $v$ for relative speed and $c$ for light speed. It tells how time and
length convert between inertial observers. The Einstein equation

$$G_{\mu\nu}+\Lambda g_{\mu\nu}=8\pi G T_{\mu\nu}$$

uses $G_{\mu\nu}$ for the Einstein tensor, $g_{\mu\nu}$ for the metric,
$\Lambda$ for the cosmological constant, $G$ for Newton's constant, and
$T_{\mu\nu}$ for stress-energy. The equation says geometry and energy-momentum
must fit together.

The generalized entropy

$$S_{gen}=\frac{A}{4G\hbar}+S_{out}$$

combines an area term with outside quantum entropy. In Jacobson-style and
entanglement-equilibrium reasoning, variations of this quantity help connect
thermodynamics to gravitational dynamics. The modular Hamiltonian $K$ is
defined by $\rho\propto e^{-K}$ for a reduced state $\rho$. It can differ from
ordinary energy; in special wedge or cap limits it behaves like a boost
generator.

In OPH, the finite consensus theorem supplies quotient normal forms for visible
records. The Einstein row additionally uses the recovered-core bridge:
geometric readout, BW modular covariance, the null-stress bridge, the
bounded-interval kernel, fixed-cap stationarity, the small-ball area identity,
remainder control, and the timelike tensor upgrade.

The modular-flow, null-blowup, and Newton-limit diagrams in the chapter are
three steps of the same story. Smooth modular flow gives a local clock.
Near a boundary point, curved screen geometry straightens into a null ray.
In the weak-field slow-motion limit, Einstein's equation reduces to Poisson's
equation and then Newton's force law.

The builders include Galileo, Newton, Maxwell, Michelson, Morley, Lorentz,
Poincare, Einstein, Minkowski, Riemann, Grossmann, Hilbert, Noether,
Eddington, Schwarzschild, Friedmann, Lemaitre, Penrose, Hawking, Unruh,
Bisognano, Wichmann, Jacobson, Wald, and many gravitational-wave and
astronomical teams. Relativity is often narrated as Einstein's miracle, but
the working structure is a network of mathematics, experiment, and
observation.

## Chapter 16: Matter, Motion, and Classical Physics

Chapter 16 explains why the quantum screen can look like ordinary matter in
ordinary space. The answer is stabilization. Some excitations survive
transport, symmetry constraints, decay channels, decoherence, and redundant
recording. Those are the patterns we call matter.

The relativistic energy relation

$$E^2=p^2c^2+m^2c^4$$

uses $E$ for energy, $p$ for momentum, $m$ for mass, and $c$ for light speed.
In natural units, where $c=1$, it becomes $E^2=p^2+m^2$. The action principle
uses

$$S=\int L\,dt,$$

with $L$ the Lagrangian and $t$ time. Stationary action means small changes to
the path do not change $S$ to first order. The stationary path can be a saddle
or maximum rather than a literal minimum.

The chapter lists particle masses and couplings, and the notes distinguish
their support level. A benchmark check tests the framework against known
values. An empirically anchored witness uses a measured input. A
first-principles prediction depends only on declared source inputs. A hadron
calculation also needs QCD binding, which dominates hadron mass.

The hierarchy/naturality relation has its own support level. The weak hierarchy is carried
by the selected local/global bridge between the pixel scale and the screen
capacity, with twelve screen ports, a 24-slot oriented repair register, zero
bridge residual, and zero Higgs naturality defect. That strengthens the particle
prediction record without upgrading the weak-boson benchmark pair, charged
lepton absolute masses, first-principles hadrons, Strong CP, or the SI gravity
clock stack.

The matter-stability diagram shows the ladder: screen excitations,
particles, atoms, chemistry, and classical objects. A macroscopic object is
not fundamental because it is large. It is public because environmental
records make it durably sampleable.

The builders include Dalton, Mendeleev, Thomson, Rutherford, Bohr, Moseley,
Chadwick, de Broglie, Schrodinger, Heisenberg, Dirac, Pauli, Fermi, Yukawa,
Anderson, Gell-Mann, Zweig, Feynman, Schwinger, Tomonaga, Dyson, Wilson,
Gross, Wilczek, Politzer, Higgs-sector experimental teams, lattice QCD
groups, and metrology collaborations. Matter is the public face of a deeply
quantum, deeply collective construction.

## Chapter 17: Darwin's Laws

Chapter 17 introduces selection as a way to think about why these laws rather
than arbitrary alternatives. It does not claim that equations reproduce like
organisms. It claims that candidate structures must pass filters before they
can belong to a stable public world.

The selection filters are finite capacity, overlap consistency, record
stability, compression, recoverability, and observer support. A law that
cannot fit the screen budget, cannot stabilize records, or cannot let
observers compare notes is not a viable public law. This is why the chapter
compares laws to protocols. A protocol survives because it enables reliable
coordination.

The entropy expressions $A/(4G)$ or $A/(4\ell_P^2)$ appear again as capacity
language. The area $A$ is the selection environment. The constants $G$ and
$\ell_P$ set the gravitational unit system. Compression is the information
theoretic demand that a law summarize observations more efficiently than a
lookup table. Quantum Darwinism supplies the microphysical cousin: pointer
states become classical because the environment copies them redundantly.

The selection-filters diagram shows a candidate pattern passing through
consistency, record, compression, and observer-support gates. The output is
not "truth by popularity." It is public physics: structure that can be
checked and used by many observers.

The human lineage includes Darwin, Wallace, Mendel, Fisher, Haldane, Wright,
Mayr, Dobzhansky, Dawkins, Zurek, Wheeler, Smolin, Susskind, Weinberg,
Tegmark, and many critics of anthropic reasoning as well. The chapter is not
asking readers to accept every selection story. It is asking them to notice
that modern physics repeatedly turns existence questions into constraint and
filter questions.

## Chapter 18: The Strange Loop

Chapter 18 turns the self-reference theme into physics. It starts where
mathematics first made self-reference precise, with Gödel's self-describing
sentence, Tarski on truth, Turing on halting, and Hofstadter's strange loops,
then asks what happens when the self-referential object is the whole universe.
The reverse-engineering target is the forcing: how the demand that a world read
itself removes freedom after freedom until only two numbers remain, and how the
loop then writes one equation for each of them.

The symbols to keep steady are the two closure equations. The local one is
$P = \varphi + \sqrt\pi/A_T(P)$: the pixel ratio $P$ has to return itself after
a full pass through the forward chain $A_T$, with $\varphi$ the golden-ratio
balance and $\sqrt\pi$ the boundary width. The global one is $N = F(N)$: the
record capacity $N$ has to equal the capacity the inside observers read back
through the map $F$. The chapter's claim is that these two equations come
straight from what self-reading requires, once geometry, signature, dimension,
the gravity law, and the gauge menu are forced. Chapter 18 tells that
subtraction as a story, and Chapter 19 folds the same chain into the synthesis.

The chapter is careful about status. The uniqueness of $P$ on the declared
domain is a certified result. The capacity fixed point is conditional on stated
premises. Both equations are open at the last step by a named quantity, the
low-energy hadronic transport term on the local side and a counting theorem on
the global side. The no-landscape corollary is the payoff: one loop, one
$(P, N)$, no vacuum to relocate into, and therefore a theory that a single blind
number can end.

The human lineage runs through Gödel, Tarski, Turing, Quine, and Hofstadter on
the mathematics of self-reference, Escher for the picture, Wheeler for the
participatory universe and it-from-bit, Chew and the later conformal-bootstrap
community for consistency-fixes-the-theory as a working method, and Tomita and
Takesaki, Jacobson, Tannaka and Krein, and the holographic-code builders for the
standard theorems that each turn out to be a small loop.

## Chapter 19: Synthesis

Chapter 19 gathers the machinery into one sentence: reality is the
consistency of observations across overlapping perspectives. The sentence is
short because the earlier chapters did the work. It includes finite
screens, algebras, overlaps, entropy, recovery, holography, entanglement,
codes, modular time, symmetry, de Sitter capacity, particle structure,
classical matter, and selection. The selection chain the chapter tells as a
story lists each requirement, the freedom it removes, its mathematical carrier,
and its current status. The full proof strategy runs in two parts, step by step,
with every open gap named.

The symbol guide is mostly referential. The Standard Model
quotient by $\mathbb Z_6$ says that a shared discrete center is identified
across the gauge factors. The pixel ratio $P$ is the local fixed point. The
fine-structure constant $\alpha$ measures electromagnetic coupling, and
$\alpha^{-1}$ is its inverse. The process $e^+e^-\to\mathrm{hadrons}$ is an
electron-positron annihilation channel whose data constrain hadronic
spectral payloads. Strong CP is encoded in a QCD angle often written
$\theta$ or $\bar\theta$, and the book is careful that OPH does not derive its
vanishing.

The synthesis chapter is also where support language matters most. A
reconstruction can be impressive without being uniform in theorem support.
Some parts are recovered-core structural consequences. Some are implemented
Phase-II closures with conditional input-free theorem targets. Some are
assumption-dependent continuations. Some are
empirical validations. Some depend on external payloads. A serious book should
not flatten those categories.

The human lineage here is the full relay. No one builds a theory of this
scope alone. The synthesis is made possible by thermodynamicists, quantum
founders, relativists, particle physicists, cosmologists, information
theorists, category theorists, condensed-matter physicists, experimental
collaborations, instrument builders, data analysts, and skeptical critics.
Every consistency check in the book has the same moral form as science
itself: partial perspectives are forced to meet.

## Chapter 20: Metaphysics

Chapter 20 asks what the physics means for experience, objectivity, and
existence. Its discipline is to keep metaphysics downstream of the technical
structure. Experience is read as the inside of observer patches, rather than a
ghostly substance added afterward. Objectivity is the overlap-stable public
record, with no God's-eye inventory underneath it.

The mathematical metaphor of a sheaf is central. Local data are assigned to
regions. If local data agree on overlaps, they may glue to a global section.
If not, the obstruction matters. In OPH, the situation is richer because the
local data are quantum states on algebras, but the sheaf idea captures the
logic of objectivity: public world is successful gluing.

The chapter also uses $\varphi$, $P$, $\sqrt\pi$, and $\alpha^{-1}$ in the
self-reference and pixel-fixed-point discussion. $\varphi$ is the golden
ratio when it appears. $P$ is the local pixel ratio. $\sqrt\pi$ appears as a
geometric normalization in the screen-side story. $\alpha^{-1}$ is the
inverse fine-structure constant. Each symbol belongs to the fixed-point
bookkeeping, not numerology. The claim stands or falls by the
declared equations and numerical checks.

The observer-loop diagram shows the metaphysical closure: world, observers,
records, and models feed back into one self-description. The loop is
structural, not ordinary backward-in-time causation.

The human builders include Nagel, Kant, Husserl, James, Peirce, Godel,
Turing, Hofstadter, Wheeler, Everett, Zurek, Rovelli, and many philosophers
of mind and science. The chapter's wager is that philosophy improves when it
does not float away from physics, and physics improves when it notices that
observers were never outside the system.

## Epilogue: Continuation and Restoration

The epilogue gives the continuation architecture. Observer patterns are
structural, identifiable through interfaces, and recoverable under the
restoration surface. The important distinction is between backup and
continuation. A backup is an external record. Continuation asks whether the
restored pattern carries the same internal flow of experience.

The key terms are **boundary-sector label**, **interface-relative interior state**,
and **interface**. A boundary-sector label tells how the observer pattern
glues to its environment. An interface-relative interior state is the inside pattern
specified relative to that interface. Given-data independence means that,
once the relevant boundary data are fixed, inside and outside do not need
extra direct information about each other to make compatible predictions.

The recovery language places immortality on the engineering surface. OPH's
operational surface contains exact or approximate restoration statements for
accessible checkpoint data under controlled interface conditions. That changes
the category of the question. Continuation becomes an engineering and identity
problem with mathematical boundaries.

The human chain includes memory research, neuroscience, quantum information,
cybernetics, computer science, philosophy of personal identity, cryonics
debates, and hardware design. It is the point where the book's
reverse-engineering posture turns back toward action.

## Shared Symbol Glossary

This final glossary collects symbols that recur across chapters. It is meant
for quick orientation, not formal completeness.

**$A$** usually means area. In entropy formulas it is often a horizon or
boundary area. When the book says capacity scales with $A$, it means that
gravity counts independent information by boundary area, not bulk
volume.

**$\ell_P$** is the Planck length. It is built from $G$, $\hbar$, and $c$.
The square $\ell_P^2$ is the Planck area. In gravitational entropy formulas,
area is measured in units of $\ell_P^2$.

**$S$** can mean entropy or a Bell correlation combination. In expressions
like $S(A)$, $S_{BH}$, or $S=k_B\ln W$, it means entropy. In CHSH contexts,
where bounds like $|S|\leq2$ appear, it means a correlation statistic.

**$k_B$** is Boltzmann's constant. It converts microscopic statistical
counting into thermodynamic units. It appears in entropy, temperature, and
Landauer-cost formulas.

**$\hbar$** is Planck's constant divided by $2\pi$. It sets the scale of
quantum action. Commutators, uncertainty relations, Unruh temperature, and
many quantum formulas carry it.

**$c$** is the speed of light. In relativity it is the invariant speed that
converts time and space units. In natural units physicists often set $c=1$.

**$G$** is Newton's gravitational constant. It controls the strength of
gravity and appears in Einstein's equation and gravitational entropy.

**$\Lambda$** is the cosmological constant. In the de Sitter chapters it sets
the horizon radius and therefore the global screen capacity at the cosmic
record-capacity fixed point \(N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}})\).

**$P$** is the OPH local pixel ratio. It is a dimensionless fixed-point value
linking the local screen grain to the electromagnetic readout side of the
program.

**$\alpha$** is the fine-structure constant, the dimensionless strength of
electromagnetism at a specified scale. Its inverse $\alpha^{-1}$ is commonly
quoted because the low-energy value is about 137.

**$\rho$** is a density matrix, the quantum state assigned to a system or
subsystem. It can describe pure states, mixtures, and reduced states obtained
by tracing out inaccessible degrees of freedom.

**$\omega$** is often a state as a functional on an algebra. Written this way,
it assigns expectation values to observables without needing a matrix display.

**$\mathcal A(P)$** is the algebra of observables associated with patch $P$.
The algebra records which questions can be asked and how they combine.

**$[A,B]$** is a commutator, $AB-BA$. If it is zero, the two operators commute.
If it is not zero, order matters.

**$\Delta$** can mean a spread, as in $\Delta X$, or the modular operator in
modular theory. Context decides. In uncertainty formulas it is a statistical
spread. In $\sigma_t(A)=\Delta^{it}A\Delta^{-it}$ it is the modular operator.

**$K$** is often a modular Hamiltonian, defined by a relation like
$\rho\propto e^{-K}$. It generates modular flow for a restricted state.

**$T_{\mu\nu}$** is the stress-energy tensor. It packages energy density,
momentum density, pressure, and stress, and it sources curvature in Einstein's
equation.

**$G_{\mu\nu}$** is the Einstein tensor, a curvature object built from the
metric. In Einstein's equation it stands on the geometry side.

**$g_{\mu\nu}$** is the spacetime metric. It determines distances, times,
angles, and causal structure.

**$SU(3)$, $SU(2)$, and $U(1)$** are symmetry groups. In the Standard Model
they organize color, weak isospin, and hypercharge.

**$Y$** is hypercharge. Together with weak isospin component $T_3$, it gives
electric charge through $Q=T_3+Y/2$.

**$H$** can mean a Hamiltonian, a Hilbert space, Shannon entropy depending on
context, or the Higgs field. The chapter context should always tell the
reader which one is meant.

**$\mu,\nu$** are spacetime indices. They label time and space directions in
relativistic formulas.

**$\gamma$** can be the Lorentz factor or a surface such as $\gamma_A$ in the
RT formula. Again, context decides. $\gamma_A$ is a geometric surface;
$\gamma=1/\sqrt{1-v^2/c^2}$ is a boost factor.

**$\theta_{\mathrm{QCD}}$ and $\bar\theta$** are strong-CP angle parameters.
The book is careful that OPH does not derive the vanishing of the physical
strong-CP phase.

## Three Worked Reading Examples

**Example 1: Reading an entropy formula.** Suppose the book writes
$S=A/(4\ell_P^2)$. First identify $S$ as entropy, not Bell $S$. Then identify
$A$ as area and $\ell_P^2$ as Planck area. The formula says the number of
independent distinguishable states is controlled by boundary area in Planck
units. It does not say the boundary is literally painted with classical bits.
It says the gravitational capacity has area scaling.

**Example 2: Reading an overlap equation.** Suppose the book writes
$\omega_i|_{\mathcal A(P_i\cap P_j)}
=\omega_j|_{\mathcal A(P_i\cap P_j)}$. Read from the inside out. $P_i$ and
$P_j$ are patches. Their intersection is the overlap. $\mathcal A$ turns that
overlap into a menu of observables. $\omega_i$ and $\omega_j$ are state
assignments. The vertical bar restricts each state to the shared menu. The
equation says both observers make the same predictions for shared questions.

**Example 3: Reading support level.** Suppose a chapter says a value is a
benchmark check, an empirically anchored witness, a first-principles prediction,
or an assumption-dependent continuation. These are not stylistic labels. A
benchmark check tests known physics. An empirically anchored witness uses a
measured target as part of the setup. A first-principles prediction is stronger
because it uses only declared source inputs. An assumption-dependent
continuation is a plausible extension whose assumptions are not all theorem
statements. Keeping those labels visible keeps the reader oriented.

The book's ambition is large, but its reading discipline is ordinary. Track
the symbols. Track the support level of each value. Track which part of the human
chain supplied each tool. Then ask whether the architecture makes the clues
cohere better than the naive picture it replaces.
