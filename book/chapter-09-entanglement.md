# Chapter 9: Entanglement Builds Space

## 9.1 The Intuitive Picture: Space Is a Stage

Start with the old geometric intuition.

Space is a container. It's the stage on which physics happens. Objects exist at
locations in space. The distance between two objects is a property of that
stage, a fixed backdrop that exists independently of what occupies it.

This is Newton's absolute space. It's the intuition behind graph paper, GPS coordinates, and every map ever drawn. Space is geometry waiting to be filled. It exists whether or not anything is in it. Two points are close or far based on how the stage is built, not on any relationship between the things at those points.

The vacuum-empty space-is simply... empty. Nothing there. A container with nothing inside.

Quantum field theory broke that picture.

## 9.2 The Surprising Hint: The Vacuum Is Not Empty

### The Scissors of the Vacuum

Imagine you have a pair of quantum scissors and decide to cut the vacuum itself. You draw a boundary around a spherical region-nothing inside, just empty space-and snip.

In classical physics, this is boring. Space is just coordinates. You label one side A and the other side B. Nothing changes.

In quantum physics, the vacuum is anything but empty. Fields fluctuate. Virtual particles pop in and out of existence. When a pair appears near your cut, one half can end up inside your sphere and the other outside. That pair is entangled. Your cut doesn't just separate two regions-it severs a web of correlations that tied them together.

### Experimental Evidence

You can see hints of this in the **Casimir effect**. Place two metal plates close together-just a fraction of a micron apart-and they feel a tiny force pushing them together. This force comes from the vacuum modes restricted between the plates. The plates change which vacuum fluctuations can exist, and that changes the energy. The vacuum has structure, and that structure depends on boundaries.

Another hint is the **Unruh effect**. An accelerated observer sees the vacuum as a warm bath of particles. An inertial observer sees nothing. How can they disagree about whether particles exist? Because acceleration limits the accelerated observer's access to spacetime. There are regions they can't see-events behind their acceleration horizon. The loss of that information makes the vacuum look thermal.

### The Area Law

The deepest hint came from studying entanglement entropy. Take a region of space in its ground state. Draw a boundary. Compute the entanglement between inside and outside.

You might expect the entropy to scale with volume. Bigger regions have more stuff.

Instead, for ground states of local systems, the entropy scales with the **boundary area**:

$$S(A) \propto |\partial A|$$

This is the **area law** for entanglement entropy. Only degrees of freedom near the boundary-within a correlation length of the cut-contribute to the entanglement.

Space is not a passive container. It's woven from quantum correlations. The vacuum is entangled across every boundary you can draw. Cut the entanglement, and you cut the connectivity of space itself.

## 9.3 The First-Principles Reframing: Space Emerges from Entanglement

The deeper question is why space looks geometric if the vacuum is really woven
from correlations.

### The Consistency Imperative

Recall our core thesis: reality is the process of making observations between observers consistent.

If there were no correlations across your cut, the vacuum wouldn't glue itself together. You couldn't walk from A to B without noticing a seam-a glitch where observations would fail to match.

**Space is not a stage that matter lives on. Space is the pattern of correlations that enables observer agreement.**

Two regions are "close" when they share many quantum correlations-when observations in one region constrain observations in the other. Two regions are "far" when they share few correlations-when they are nearly independent.

Distance is not a primitive. It emerges from the entanglement structure of the vacuum state.

### The Ryu-Takayanagi Formula

We introduced the RT formula in Chapter 8: entanglement entropy of a boundary region equals the area of the minimal bulk surface anchored on that region's boundary, divided by 4G. This looks exactly like the Bekenstein-Hawking formula for black hole entropy, except the same structure applies to any region.

The deep implication: **geometry encodes entanglement**.

That sentence is easy to repeat and easy to misunderstand. The claim is not that geometry and entanglement vaguely resemble each other. The claim is that the amount of quantum correlation across a cut determines the size of the bulk surface associated with that cut. Entropy is doing geometric work. If the boundary state ties two regions together strongly, the bulk description between them is correspondingly thick and connected.

### A Simple Example

Consider a 2D CFT on an interval of length L. The entanglement entropy is:

$$S = \frac{c}{3}\ln\frac{L}{\epsilon}$$

where c is the central charge and epsilon is a UV cutoff.

In AdS_3, the minimal "surface" is a geodesic-a shortest path through the bulk. Compute its length using the AdS metric. Divide by 4G.

**With the standard cutoff identification, they match.** Two completely different calculations-one from quantum field theory, one from geometry-give the same answer.

## 9.4 Bell's Theorem: The Reality of Entanglement

Entanglement is not a decorative idea. Bell experiments force it on us.

For suitable two-wing experiments, any local hidden-variable account obeys

$$|S| \leq 2.$$

Quantum mechanics allows a stronger pattern and reaches the Tsirelson limit

$$|S| \leq 2\sqrt{2}.$$

That stronger pattern has been observed. The 2015 loophole-free Bell tests
closed the major loopholes at the same time and ruled out the simple local
hidden-variable models Einstein hoped would survive.

The Bell pair makes the structure vivid:

$$|\Phi^+\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle).$$

The total state is pure, while each qubit by itself is maximally mixed:

$$\rho_A = \frac{1}{2}|0\rangle\langle0| + \frac{1}{2}|1\rangle\langle1|.$$

The whole is ordered in a way the parts are not. That is the operational signature of entanglement.

Chapter 6 gives the overlap-consistency reading of Bell correlations. The lesson needed here is narrower. Space cannot be built out of classical bookkeeping alone. The vacuum correlations used in this chapter belong to the same experimentally compulsory quantum structure.

## 9.5 ER = EPR: Wormholes Are Entanglement

Einstein and Rosen wrote about wormholes in 1935. Einstein, Podolsky, and Rosen wrote about entanglement the same year. For eighty years, no one connected them.

In 2013, Juan Maldacena and Leonard Susskind made a bold proposal: **ER = EPR**.

The proposal is that Einstein-Rosen bridges (wormholes) and Einstein-Podolsky-Rosen correlations (entanglement) are deeply linked-two languages for the same underlying connectivity in the right regimes.

### The Thermofield Double

The strongest evidence comes from the **thermofield double state**:

$$|\text{TFD}\rangle = \sum_n e^{-\beta E_n/2} |n\rangle_L |n\rangle_R$$

This state lives on two copies of a system. It is an entangled purification of a thermal state at temperature T = 1/beta.

In AdS/CFT, the thermofield double is dual to an **eternal two-sided black hole**. The two boundaries correspond to two copies of the CFT. They're connected by a smooth wormhole through the interior.

Break the entanglement and the wormhole collapses. Maintain the entanglement and
the connection holds.

### Traversable Wormholes

In 2017, Gao, Jafferis, and Wall showed that with a small coupling between the two boundaries, the wormhole becomes **traversable**. You can send a message from one side to the other.

In the dual setting, the same protocol can be read in quantum-information language as **quantum teleportation** and in bulk language as sending a signal through the wormhole.

## 9.6 Bit Threads: A Flow Picture

The RT formula uses minimal surfaces. In 2016, Freedman and Headrick introduced an equivalent picture: **bit threads**.

Instead of drawing a surface, draw threads-imaginary lines carrying entanglement. The density of threads can't exceed 1/4G at any point. Subject to this constraint, maximize the number of threads connecting region A to its complement.

The maximum number equals the RT entropy.

This is a **max-flow, min-cut theorem** in a gravitational setting. The minimal surface is where thread density is maximized-the bottleneck.

In the language of this book, threads are the links that let observers compare notes. The more threads between two regions, the more they can agree about shared observations.

## 9.7 Tensor Networks: Circuits for Spacetime

The RT formula tells you the answer. Tensor networks give you the mechanism.

A **tensor network** builds a large quantum state from small pieces. Each tensor is a multi-index array. The connections between tensors represent entanglement.

### MERA: Building in Scale

The **Multi-scale Entanglement Renormalization Ansatz** (MERA) handles critical systems by building in scale. Layer by layer, you move to larger scales. The network grows upward into a new dimension.

In 2012, Brian Swingle noticed something striking: the geometry of a MERA network is **hyperbolic**-just like AdS space. The depth in the network plays the role of the radial direction in AdS.

MERA isn't just a numerical trick. It's a discrete version of AdS/CFT-the first concrete circuit that turns entanglement into geometry.

### The HaPPY Code

In 2015, Hayden, Pastawski, Preskill, Nezami, and Yoshida built a toy model called the **HaPPY code**.

They tiled a hyperbolic disk with perfect tensors, and two things happened at
once. The RT formula became exact, with boundary entropy equal to the number of
legs cut by a minimal path. Bulk operators also became recoverable from
different boundary regions.

This redundancy is quantum error correction. The bulk exists because it's the error-corrected version of the boundary.

## 9.8 Monogamy: Why Space Is Local

If entanglement builds space, why does space look local? Why can't you step from New York to Tokyo in one move?

The answer is **monogamy of entanglement**.

Quantum entanglement is jealous. If system A is maximally entangled with system B, it can't be entangled with system C at all:

$$\tau_{A:BC} \geq \tau_{A:B} + \tau_{A:C}$$

This forces the entanglement network to be sparse. You can't make a complete graph where everything is equally close to everything else. You're pushed toward a lattice-like structure with modest connectivity.

That's what locality means. Things can only be near a limited number of other things. Geometry emerges from the constraints of entanglement monogamy.

## 9.9 Entanglement Wedges and Reconstruction

The RT surface divides the bulk into pieces. The region between a boundary region A and its RT surface is called the **entanglement wedge** of A.

**Subregion duality**: The physics inside the entanglement wedge can be reconstructed from boundary region A alone.

### Overlapping Wedges

Consider two observers with access to different boundary regions. If their entanglement wedges overlap, they can both reconstruct the same bulk physics. That overlap is where their observations must agree.

This is consistency made geometric. The structure of entanglement forces their reconstructions to match in the overlap.

### Black Holes and Islands

In AdS/CFT and related semiclassical holographic models, there is a striking late-time effect: as Hawking radiation accumulates, the radiation's entanglement wedge can include a region **inside** the black hole-an "island."

In those models, the island formula reproduces the Page curve and shows how the radiation can encode information that semiclassical bulk reasoning seemed to lose.

This matters because it shows what holographic encoding can do. Information
that looks lost in one description can remain stored in a nonlocal way and
reappear when the right radiation data are assembled. That is the lesson the
book needs. The black-hole interior is not a second hidden storage room. It is
an encoded part of the same quantum description.

## 9.10 From Entanglement to the Classical World

If everything is entangled, why does the world look classical?

The answer involves **decoherence** and **quantum Darwinism**.

When a quantum system interacts with its environment, certain "pointer states" become stable-states that can be copied into the environment without being destroyed. The environment measures them repeatedly, storing redundant records.

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
Black-hole information is read through encoded radiation rather than deletion.

## 9.12 Reverse Engineering Summary

Space is not a passive backdrop. The vacuum is a web of quantum correlations,
and the structure of that web is what becomes geometry. Area-law scaling,
Ryu-Takayanagi, entanglement wedges, tensor networks, and the ER=EPR intuition
all point in the same direction. Distance is a measure of shared
correlations. Cut enough entanglement and you cut space itself.

---

We've seen that space emerges from entanglement. But why is this structure stable? Why doesn't the entanglement web unravel?

In the next chapter, we'll see how this picture connects to quantum error correction. Spacetime isn't just entanglement-it's a code that protects information. The bulk exists because it's the error-corrected version of the boundary. And this connection explains why spacetime is stable: the same mechanisms that protect quantum computers protect reality itself.
