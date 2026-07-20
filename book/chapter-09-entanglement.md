# Chapter 9: Entanglement Builds Space

## 9.1 The Intuitive Picture: Space Is a Stage

Start with the old geometric intuition.

Space is a container. It's the stage on which physics happens. Objects exist at
locations in space. The distance between two objects is a property of that
stage, a fixed backdrop that exists independently of what occupies it.

This is Newton's absolute space. It's the intuition behind graph paper, GPS coordinates, and every map ever drawn. Space is geometry waiting to be filled. It exists whether or not anything is in it. Two points are close or far based on how the stage is built, not on any relationship between the things at those points.

The vacuum is what its name promises: empty. Nothing there. A container with nothing inside.

Quantum field theory broke that picture.

## 9.2 The Surprising Hint: The Vacuum Is Not Empty

### The Scissors of the Vacuum

Imagine you have a pair of quantum scissors and decide to cut the vacuum itself. You draw a boundary around a spherical region (nothing inside, just empty space) and snip.

In classical physics, this is boring. Space is just coordinates. You label one side A and the other side B. Nothing changes.

In quantum physics, the vacuum is full of fluctuating fields. Virtual particles
appear and disappear. When a pair appears near your cut, one half can end up
inside your sphere and the other outside. That pair is entangled. Your cut
separates the regions by severing a web of correlations between them.

### Experimental Evidence

You can see hints of this in the **Casimir effect**. Place two metal plates a fraction of a micron apart and they feel a tiny force pushing them together. This force comes from the vacuum modes restricted between the plates. The plates change which vacuum fluctuations can exist, and that changes the energy. The vacuum has structure, and that structure depends on boundaries.

Another hint is the **Unruh effect**. An accelerated observer sees the vacuum as a warm bath of particles. An inertial observer sees nothing. How can they disagree about whether particles exist? Because acceleration limits the accelerated observer's access to spacetime. There are regions they can't see: events behind their acceleration horizon. The loss of that information makes the vacuum look thermal.

### The Area Law

The deepest hint came from studying entanglement entropy. Take a region of space in its ground state. Draw a boundary. Compute the entanglement between inside and outside.

You might expect the entropy to scale with volume. Bigger regions have more stuff.

Instead, for ground states of local systems, the entropy scales with the **boundary area**:

$$S(A) \propto |\partial A|$$

This is the **area law** for entanglement entropy. Only degrees of freedom near
the boundary, within a correlation length of the cut, contribute to the
entanglement.

Read the notation literally. $A$ is the region being studied. $\partial A$ is
its boundary, the cut between inside and outside. $S(A)$ is the entropy seen
by an observer who has access to $A$ but not to its complement. The symbol
$\propto$ means "is proportional to." The equation says that the dominant
entropy comes from the cut surface, not from the volume of stuff enclosed by
the cut.

Space is woven from quantum correlations rather than supplied as a passive container. The vacuum is entangled across every boundary you can draw. Cut the entanglement, and you cut the connectivity of space itself.

## 9.3 The First-Principles Reframing: Space Emerges from Entanglement

The deeper question is why space looks geometric if the vacuum is really woven
from correlations.

### The Consistency Imperative

Recall our core thesis: reality is the process of making observations between observers consistent.

If there were no correlations across your cut, the vacuum wouldn't glue itself together. You couldn't walk from A to B without noticing a seam where observations fail to match.

**Space is the pattern of correlations that enables observer agreement.** The stage picture had the dependence backwards: the correlations come first, and the arena is what they add up to.

Two regions are "close" when they share many quantum correlations, so
observations in one constrain observations in the other. They are "far" when
they share few correlations and are nearly independent.

Distance emerges from the entanglement structure of the vacuum state; nothing about it is primitive.

### The Ryu-Takayanagi Formula

We introduced the RT formula in Chapter 8: entanglement entropy of a boundary region equals the area of the minimal bulk surface anchored on that region's boundary, divided by 4G. This looks exactly like the Bekenstein-Hawking formula for black hole entropy, except the same structure applies to any region.

The deep implication: **geometry encodes entanglement**.

That sentence is easy to repeat and easy to misunderstand. The precise claim is stronger: the amount of quantum correlation across a cut determines the size of the bulk surface associated with that cut. Entropy is doing geometric work. If the boundary state ties two regions together strongly, the bulk description between them is correspondingly thick and connected.

This idea did not arrive as one person's finished vision. Bekenstein and
Hawking made black holes count states. Srednicki and many others showed how
ordinary quantum fields produce area-law entanglement. Maldacena gave the
AdS/CFT setting in which a boundary theory and a gravitational bulk could be
compared. Ryu and Takayanagi found the sharp area formula. Hubeny, Rangamani,
Takayanagi, Faulkner, Lewkowycz, Maldacena, Engelhardt, Wall, Harlow, and many
others then turned the slogan into a working toolkit. OPH uses that inherited
toolkit and asks what common observer-first architecture those clues are
pointing toward.

### A Simple Example

Consider a 2D CFT on an interval of length L. The entanglement entropy is:

$$S = \frac{c}{3}\ln\frac{L}{\epsilon}$$

where c is the central charge and epsilon is a UV cutoff.

The central charge counts, roughly, how many independent quantum degrees of
freedom the CFT has. The UV cutoff is the shortest distance the calculation
allows. Without it, the field theory would keep counting arbitrarily tiny
correlations across the boundary of the interval.

$S$ is the entanglement entropy of the interval. $L$ is the interval length.
$\epsilon$ is the short-distance cutoff, the smallest distance the effective
field theory is allowed to resolve. The logarithm appears because a
two-dimensional conformal field theory organizes correlations scale by scale.
The formula is compact, but every symbol carries a physical role.

In AdS_3, the minimal "surface" is a geodesic, the shortest path through the bulk. Compute its length using the AdS metric. Divide by 4G.

**With the standard cutoff identification, they match.** A quantum field theory
calculation and a geometric calculation give the same answer.

## 9.4 Bell's Theorem: The Reality of Entanglement

Bell experiments force entanglement on us.

For suitable two-wing experiments, any local hidden-variable account obeys

$$|S| \leq 2.$$

Quantum mechanics allows a stronger pattern and reaches the Tsirelson limit

$$|S| \leq 2\sqrt{2}.$$

Here $S$ means the Bell-CHSH correlation combination, a number built from four
correlation measurements chosen by the two experimenters. The vertical bars
mean absolute value. Classical local hidden-variable models keep
that number at or below 2. Quantum mechanics permits a larger value, but not an
arbitrary one. The ceiling $2\sqrt{2}$ is the quantum limit.

That stronger pattern has been observed. In 2015, Ronald Hanson's group in
Delft ran the first loophole-free test, with entangled electron spins held in
two diamonds 1.3 kilometers apart on opposite sides of the campus; experiments
in Vienna and at NIST followed within weeks. Together they closed the major
loopholes at the same time and ruled out the simple local hidden-variable
models Einstein hoped would survive.

The Bell pair makes the structure vivid:

$$|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle).$$

The total state is pure, while each qubit by itself is maximally mixed:

$$\rho_A = \frac{1}{2}|0\rangle\langle0| + \frac{1}{2}|1\rangle\langle1|.$$

The whole is ordered in a way the parts are not. That is the operational signature of entanglement.

The ket $|\Phi^+\rangle$ names one Bell state. The symbols $|00\rangle$ and
$|11\rangle$ mean that the two qubits are both 0 or both 1. The factor
$1/\sqrt{2}$ normalizes the state so total probability is 1. The matrix
$\rho_A$ is the state seen by an observer who can access only qubit $A$.
Because that observer has ignored the other qubit, the local state looks like
a fair coin even though the two-qubit state is perfectly ordered.

Chapter 6 gives the overlap-consistency reading of Bell correlations. The lesson needed here is narrower. Space cannot be built out of classical bookkeeping alone. The vacuum correlations belong to the same experimentally compulsory quantum structure.

## 9.5 ER = EPR: Wormholes Are Entanglement

Chapter 8 flagged the coincidence: wormholes and entanglement entered physics in the same year, 1935, in two papers that both carry Einstein's name, and for eighty years nobody connected them.

In 2013, Juan Maldacena and Leonard Susskind made a bold proposal: **ER = EPR**.

The proposal is that Einstein-Rosen bridges (wormholes) and Einstein-Podolsky-Rosen correlations (entanglement) are deeply linked: two languages for the same underlying connectivity in the right regimes.

### The Thermofield Double

The strongest evidence comes from the **thermofield double state**:

$$|\text{TFD}\rangle = \sum_n e^{-\beta E_n/2} |n\rangle_L |n\rangle_R$$

This state lives on two copies of a system. It is an entangled purification of a thermal state at temperature T = 1/beta.

The sum runs over energy states labeled by $n$. $E_n$ is the energy of state
$n$. The subscripts $L$ and $R$ label the left and right copies of the system.
$\beta$ is inverse temperature, equal to $1/(k_B T)$ when Boltzmann's constant
is kept explicit. The exponential weights higher-energy states less strongly,
just as in ordinary thermal physics.

"Purification" means that a mixed, thermal-looking state can be treated as part
of a larger pure state if we include a second auxiliary copy. The thermal
uncertainty is then reinterpreted as entanglement with that second copy.

In AdS/CFT, the thermofield double is dual to an **eternal two-sided black hole**. The two boundaries correspond to two copies of the CFT. They're connected by a smooth wormhole through the interior.

Break the entanglement and the wormhole collapses. Maintain the entanglement and
the connection holds.

### Traversable Wormholes

In 2017, Gao, Jafferis, and Wall showed that with a small coupling between the two boundaries, the wormhole becomes **traversable**. You can send a message from one side to the other.

In the dual setting, the same protocol can be read in quantum-information language as **quantum teleportation** and in bulk language as sending a signal through the wormhole.

## 9.6 Bit Threads: A Flow Picture

The RT formula uses minimal surfaces. In 2016, Freedman and Headrick introduced an equivalent picture: **bit threads**.

Draw threads: imaginary lines carrying entanglement. The density of threads can't exceed 1/4G at any point. Subject to this constraint, maximize the number of threads connecting region A to its complement.

The maximum number equals the RT entropy.

This is a **max-flow, min-cut theorem** in a gravitational setting. The minimal surface is the bottleneck where thread density is maximized.

In the language of this book, threads are the links that let observers compare notes. The more threads between two regions, the more they can agree about shared observations.

## 9.7 Tensor Networks: Circuits for Spacetime

The RT formula tells you the answer. Tensor networks give you the mechanism.

A **tensor network** builds a large quantum state from small pieces. Each tensor is a multi-index array. The connections between tensors represent entanglement.

### MERA: Building in Scale

The **Multi-scale Entanglement Renormalization Ansatz** (MERA) handles critical systems by building in scale. Layer by layer, you move to larger scales. The network grows upward into a new dimension.

In 2012, Brian Swingle noticed that the geometry of a MERA network is
**hyperbolic**, just like AdS space. The depth in the network plays the role of
the radial direction in AdS.

MERA is more than a numerical trick. It is one influential discrete model
showing how entanglement can organize geometry in an AdS-like way.

### The HaPPY Code

In 2015, Hayden, Pastawski, Preskill, Nezami, and Yoshida built a toy model called the **HaPPY code**.

They tiled a hyperbolic disk with perfect tensors, and two things happened at
once. The RT formula became exact, with boundary entropy equal to the number of
legs cut by a minimal path. Bulk operators also became recoverable from
different boundary regions.

In that toy-model setting, this redundancy is quantum error correction. The bulk exists because it
is encoded in boundary degrees of freedom with a supplied code structure.

## 9.8 Monogamy: Why Space Is Local

If entanglement builds space, why does space look local? Why can't you step from New York to Tokyo in one move?

**Monogamy of entanglement** supplies the answer.

Quantum entanglement is jealous. If system A is maximally entangled with system B, it can't be entangled with system C at all:

$$\tau_{A:BC} \geq \tau_{A:B} + \tau_{A:C}$$

This forces the entanglement network to be sparse. You can't make a complete graph where everything is equally close to everything else. You're pushed toward a lattice-like structure with modest connectivity.

The symbol $\tau$ is the tangle, a measure of how much entanglement one system
shares with another. The label $A:BC$ means "system $A$ compared with the joint
system made from $B$ and $C$." The inequality says that entanglement cannot be
freely duplicated. If $A$ uses up its entanglement budget with $B$, less is
available for $C$.

That is one reason locality can emerge. Things can only be near a limited
number of other things. Geometry can then inherit sparse structure from the
constraints of entanglement monogamy.

## 9.9 Entanglement Wedges and Reconstruction

The RT surface divides the bulk into pieces. The region between a boundary region A and its RT surface is called the **entanglement wedge** of A.

**Subregion duality**: The physics inside the entanglement wedge can be reconstructed from boundary region A alone.

![Two boundary regions can reconstruct overlapping entanglement wedges, forcing agreement on the shared bulk region.](../assets/book_diagrams/entanglement-wedge-overlap.svg){width=78%}

### Overlapping Wedges

Consider two observers with access to different boundary regions. If their entanglement wedges overlap, they can both reconstruct the same bulk physics. That overlap is where their observations must agree.

This is consistency made geometric. The structure of entanglement forces their reconstructions to match in the overlap.

### Black Holes and Islands

In AdS/CFT and related semiclassical holographic models, Hawking radiation can
produce a late-time effect: the radiation's entanglement wedge can include a
region **inside** the black hole, an "island."

In those models, the island formula reproduces the Page curve and shows how the radiation can encode information that semiclassical bulk reasoning seemed to lose.

The island formula shows what holographic encoding can do. Information that
looks lost in one description can remain stored in a nonlocal way and reappear
when the right radiation data are assembled. That is the lesson the book needs.
The black-hole interior is encoded inside the same quantum description. It is
no separate hidden storage room.

## 9.10 From Entanglement to the Classical World

If everything is entangled, why does the world look classical?

The answer involves **decoherence** and **quantum Darwinism**.

When a quantum system interacts with its environment, certain "pointer states" become stable: states that can be copied into the environment without being destroyed. The environment measures them repeatedly, storing redundant records.

Classical facts are quantum information that got copied everywhere. You look at a chair. I look at the same chair. We agree because we're both sampling redundant records in the environment.

This is error correction as a law of physics. Reality stabilizes itself through redundancy.

## 9.11 What Entanglement Predicts

If geometry is built from entanglement, several things have to happen. In the
holographic regime, boundary entropy has to track extremal surfaces. Low-energy
states have to care about boundary area more than bulk volume. Entropy
inequalities have to act like geometric constraints. Bulk regions have to be
reconstructible from the right boundary data. Black-hole evaporation has to
respect unitarity through encoded interior information.

That is the direction the evidence points. Ryu-Takayanagi works in the
settings where it should. Area-law scaling is widespread in the regimes that
matter here. Entanglement wedge reconstruction works in explicit examples.
Black-hole information is read through encoded radiation, not simple deletion,
in the semiclassical holographic models where the island/Page-curve picture is
under control.

### Entanglement as a Discovery Chain

Entanglement has one of the strangest histories in physics because it began
as a complaint. Einstein, Podolsky, and Rosen used it in 1935 to argue that
quantum mechanics could not be complete. Schrödinger answered by naming the
phenomenon and rating it, in his 1935 reply, as "the characteristic trait of
quantum mechanics, the one that enforces its entire departure from classical
lines of thought." For decades the issue looked
philosophical. Then John Bell found the inequality that moved the argument
from taste to experiment. Clauser, Freedman, Aspect, Zeilinger, and many
others turned the test into a laboratory program. Modern loophole-free Bell
experiments make the point hard to evade: the world does not carry local
classical instruction sheets that determine all possible measurement results.

The later holographic story is just as communal. Bekenstein and Hawking made
horizons thermodynamic. 't Hooft and Susskind drew the boundary-first lesson.
Maldacena supplied the AdS/CFT laboratory. Ryu and Takayanagi wrote the
minimal-area formula that made entanglement look geometric. Hubeny,
Rangamani, and Takayanagi extended it to time-dependent cases. Van Raamsdonk
emphasized that changing entanglement can change connectivity. Swingle,
Pastawski, Hayden, Preskill, Yoshida, Almheiri, Engelhardt, Penington, and
many others connected tensor networks, error correction, wedges, and islands.

That history matters because the formula

$$S(A)=\frac{\mathrm{Area}(\gamma_A)}{4G_N}$$

is easy to admire without digesting. $S(A)$ is the entanglement entropy of
boundary region $A$. The surface $\gamma_A$ is the bulk extremal surface
anchored to the boundary of $A$. $G_N$ is Newton's gravitational constant.
The equation says that a quantum-information quantity on a boundary is
measured by a geometric area in the bulk. If area changes, entropy changes.
If entropy changes, the emergent geometry changes.

The Bell inequality in the same chapter uses the symbol $S$ in a different
way, for the CHSH combination of correlations,
usually written as a sum and difference of expectation values for four
measurement settings. The classical limit is $|S|\leq 2$. Quantum theory can
reach $2\sqrt 2$. The repeated letter is unfortunate but common. The book
uses context to distinguish them: entropy $S(A)$ belongs to regions and
areas, while Bell $S$ belongs to correlation tests.

OPH needs both meanings. Bell-type entanglement says local classical
bookkeeping is too small. RT-type entanglement says quantum bookkeeping can
be geometric. Entanglement wedges then say that two observers may reconstruct
the same interior region from different boundary data. If they do, agreement
is forced. The shared wedge is the geometric version of an overlap.

For finite evidence bundles, that overlap language has a strict neutral-bulk
version. Observer features from different shards are comparable only after
gauge and port presentation labels have been quotiented, terminal normal forms
have been reached, and interface transport has glued the charts without hidden
cycle holonomy. A distance built from the descended features is first a
pseudometric. It becomes a metric only after feature collisions are quotiented
or the features are proved jointly separating.

## 9.12 Reverse Engineering Summary

Space acts less like a passive backdrop than a web of quantum correlations,
with the structure of that web becoming geometry. Area-law scaling,
Ryu-Takayanagi, entanglement wedges, tensor networks, and the ER=EPR intuition
all point in the same direction. Distance is a measure of shared
correlations. Cut enough entanglement and you cut space itself.

---

Space emerges from entanglement. But why is this structure stable? Why doesn't the entanglement web unravel?

In the next chapter, this picture connects to quantum error correction. In controlled holographic
settings, spacetime is encoded in a structure that can protect information. The OPH overlap
network begins more modestly as a finite constraint code, and it becomes a quantum
error-correcting claim only when the code subspace, logical operators, error model, and recovery
map are supplied.
