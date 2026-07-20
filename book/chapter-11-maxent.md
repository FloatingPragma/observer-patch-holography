# Chapter 11: MaxEnt and the Arrow

## 11.1 The Intuitive Picture: Time Is Fundamental

Start with the Newtonian picture of time.

Time is a fundamental external parameter. It flows from past to future,
independent of anything in the universe. Events happen in time, just as objects
exist in space. The clock ticks whether or not anything is happening. Time is
the stage; physics is the play.

This is Newton's absolute time: "Absolute, true, and mathematical time, of itself, and from its own nature, flows equably without relation to anything external."

The arrow of time feels fundamental in this picture. We remember yesterday. We do not remember tomorrow. Eggs break. They do not unbreak. Time has a direction, built into its very nature.

General relativity and quantum theory broke that picture.

![Known constraints select a least-biased state, and the restricted algebra-state pair carries an internal modular clock.](../assets/book_diagrams/maxent-clock.svg){width=78%}

## 11.2 The Surprising Hint: Time Is Not Fundamental

### The Scandal of the Second Law

Physics has a scandal.

Almost all our fundamental laws are time-reversible. Newton's F = ma works the same forward and backward. Maxwell's equations are reversible. Schrödinger's equation is reversible. Einstein's General Relativity is reversible.

Film a planet orbiting a star and play it backward-it looks perfectly physical. But film an egg breaking and play it backward? Absurd.

This is the **Arrow of Time**. Where does it come from? The microscopic laws
are the wrong place to look.

### No Preferred Time in GR

In general relativity, there's no preferred time coordinate. Different observers slice spacetime differently; none is privileged.

The Wheeler-DeWitt equation-the analog of Schrödinger's equation for the universe-is:

$$H\Psi = 0$$

The Hamiltonian acting on the wavefunction of the universe gives zero. There is no explicit external time derivative in this formalism, so the universe can look *frozen* at the fundamental level.

$H$ is the Hamiltonian constraint, the operator that would normally generate
time evolution. $\Psi$ is the wavefunction of the universe in this formal
setting. The equation does not say that nothing happens in experience. It says
that the fundamental constraint equation has no outside time parameter built
into it.

This is the **problem of time** in quantum gravity. If the fundamental description has no time, where does time come from?

Time enters this story without a fundamental external clock. The microscopic laws are time-symmetric. Something else must generate the arrow of time we experience.

## 11.3 The First-Principles Reframing: Time Emerges from Modular Flow

The deeper question is why we experience time at all if the fundamental
description has no preferred clock.

### The Thermal Time principle

In the 1990s, Alain Connes and Carlo Rovelli proposed a stark idea. Time can be
read from incomplete knowledge. Start with the observer's limited state
$\rho$. At finite cutoff, or in a special type-I representation, one writes a
modular Hamiltonian $K=-\ln\rho$. In the general operator-algebraic setting the
more fundamental object is the modular automorphism group itself; an inner
density-matrix generator need not exist. The thermal-time proposal reads that
flow as the time the observer actually experiences.

Here $\rho$ is the observer's density matrix, the quantum bookkeeping object for
what the observer can know. The modular Hamiltonian is not necessarily the
ordinary energy of the whole universe but the operator that generates the
natural evolution of the restricted state the observer actually has access to
when such an operator representative exists.

This is a strange move the first time one sees it. In ordinary mechanics, the
Hamiltonian is given first and time evolution follows. Here the restricted
state itself furnishes the clock. Time is tied to access, ignorance, and coarse
graining.

### Tomita-Takesaki Theory

The deeper theorem behind that proposal comes from operator algebra. Once an
observer has a rich enough algebra of questions and a state that probes it
fully, the pair carries a preferred internal flow whether or not anyone inserts
an external master clock. The formal machinery is called
**Tomita-Takesaki theory**.

An automorphism is a reshuffling of the allowed questions that preserves their
algebraic rules. Modular flow is a continuous family of such reshufflings,
indexed by a parameter that behaves like time. The concrete content is that the
observer's restricted state tells the question menu how to move.

An observer with partial access does not sit in a timeless fog. The
restriction itself orders experience into before and after. The flow depends on
the algebra-state pair, which is why different observers can inherit different
clocks from different access conditions. It also carries the thermal
equilibrium structure that makes temperature and time appear together in the
same move.

Modular flow matters here because it turns time from a background stage into
something earned by an observer's horizon and state. Partial knowledge has its
own dynamics.

### The Rindler Wedge

This abstract mathematics connects to reality through the Unruh effect.

First recall what a Lorentz boost is. In special relativity, two observers
moving at constant velocity do not split space and time in the same way. A
Lorentz boost is the transformation that converts one observer's space-time
coordinates into the other's. It is like a rotation, but a rotation in
spacetime: it tilts the time axis and one space axis while preserving the light
cone and the spacetime interval.

The word "generator" means the infinitesimal version of that transformation.
Just as angular momentum generates ordinary rotations, the boost generator
generates changes of inertial frame. A steadily accelerating observer can be
thought of as passing through a sequence of nearby inertial frames. Step by
step, their time evolution is built from tiny Lorentz boosts.

An observer accelerating uniformly sees only the **Rindler wedge**, the part of
spacetime from which signals can reach that observer. A horizon forms
behind them. For the vacuum state restricted to this wedge, the
Bisognano-Wichmann theorem shows that the modular Hamiltonian is exactly the
generator of Lorentz boosts.

For an accelerating observer, a Lorentz boost *is* time translation. The modular flow equals ordinary time evolution.

The modular temperature works out to:

$$T_{Unruh} = \frac{\hbar a}{2\pi k_B c}$$

The Unruh effect is Tomita-Takesaki theory applied to spacetime. The "time"
experienced by an observer is determined by their restricted access to the
quantum state.

$T_{Unruh}$ is the temperature seen by the accelerating observer. $a$ is the
observer's proper acceleration. The constants are the same ones used in the
black-hole temperature formulas. The larger the acceleration, the hotter the
restricted vacuum appears.

This is the first major payoff of the chapter. The mathematics does not stay abstract. Restrict the vacuum to what one observer can access, and the restriction behaves like a thermal state with its own clock.

This is also the point where modular theory stops sounding like rarefied operator algebra and starts sounding like lived physics. Restriction generates both a temperature and a time flow. Losing access to part of the state has thermodynamic and temporal consequences.

## 11.4 The Arrow of Time

In Chapter 4, we saw Boltzmann's insight: entropy $S = k \ln W$ measures the number of microstates compatible with a macrostate, and entropy increases because high-entropy states vastly outnumber low-entropy ones.

But why does the accessible cosmological record have a low-entropy side?

### The Past Hypothesis

The deeper answer to the arrow of time is the **Past Hypothesis**: the record
we inhabit is anchored on a state of extraordinarily low entropy.

Standard cosmology describes the hot dense side of our branch as smooth, with
matter spread almost uniformly. That uniformity is low gravitational entropy.

Why is that side low entropy? Standard physics treats this as an unexplained
boundary condition. OPH gives the condition a role in record formation.

**The Past Hypothesis as a consistency requirement**: For observers to exist at all, they must be able to form and compare records. Records require entropy gradients-you can only write information by pushing entropy elsewhere. A universe in thermal equilibrium contains no observers, no records, no consistency-checking.

The MaxEnt principle says: assign the maximum-entropy state consistent with your constraints. But one constraint is that someone must exist to apply MaxEnt. This rules out equilibrium. The very existence of observers selecting MaxEnt states presupposes a universe far from equilibrium.

The specific numerical entropy of the hot dense record belongs to physical
cosmology. The OPH point is the consistency role: observers who compare records
require a low-entropy side, and the arrow of time points in the direction that
allows records to be made.

## 11.5 Jaynes: Entropy as Ignorance

Edwin Jaynes rewrote statistical mechanics in information-theoretic terms.

**Entropy is not a property of the gas. Entropy is a property of our knowledge about the gas.**

### The Maximum Entropy Principle

Suppose you know only the average energy. What probability distribution should you assign?

Choose the distribution that maximizes Shannon entropy subject to your constraints:

$$S = -\sum_i p_i \ln p_i$$

MaxEnt gives the Boltzmann distribution:

$$P(x) = \frac{1}{Z} e^{-\beta E(x)}$$

Thermal states are ubiquitous because they're the unique states of maximum ignorance given energy constraints.

In the entropy formula, $p_i$ is the probability of outcome $i$. In the
Boltzmann distribution, $P(x)$ is the probability of state $x$, $E(x)$ is its
energy, $\beta$ is inverse temperature, and $Z$ is the partition function that
normalizes all probabilities so they add to 1.

## 11.6 Time on the Holographic Screen

In the simplest finite type-I illustration, a support region $P$ on the screen
chart is cut from a global pure state. The restriction gives a density matrix:

$$\rho_P = \text{Tr}_{\bar{P}} |\Psi\rangle \langle \Psi|$$

This density matrix defines a modular Hamiltonian:

$$K_P = -\ln \rho_P$$

which generates the finite modular flow labeled by $t_P$. That flow supplies
the observer's internal time.

$\bar P$ means the complement of support region $P$, everything outside the
region in this finite illustration. The trace over $\bar P$ discards
inaccessible degrees of freedom and leaves the state available to the observer.
The logarithm then turns that restricted state into the modular generator.

This density-matrix formula is the finite illustration. The general observer
patch is described by its accessible algebra-state pair. That pair carries a
modular automorphism group, which supplies the internal clock even when no
density matrix lives inside the continuum algebra.

### Consistency of Clocks

For a shared clock, two overlapping observers have compatible modular actions
on their shared operational content. That compatibility supports a common
causal structure.

### Cosmic Time

Why do we all agree on a "cosmic time"?

As local modular flows synchronize across accessible records, a shared
coarse-grained cosmic time emerges as their collective clock rather than a
second fundamental time parameter. It is the time read by the reconstructed
world from within.

### From Modular Time to Gravity

The chain is clean once the pieces are visible. Recovery structure from Chapter
7 makes the time generator local near patch boundaries. The order and
orientation of the finite caps, their modular convergence, and their thermal
normalization turn that local flow into a standard geometric motion on the
sphere. Geometric time flow gives Lorentz kinematics on the screen. Fixed-cap
variation, local energy, the small-ball area identity, and the
scalar-to-tensor upgrade then give the tensor first-variation form of the
Einstein relation. A vacuum reference on the same common-domain
tower evaluates the integration tensor and gives the absolute equation.

The internal clock feeds into gravity through this chain of compatible
readouts.

## 11.7 Jacobson's Derivation

In 1995, Ted Jacobson performed one of the most beautiful derivations in theoretical physics.

He started with thermodynamics-the first law:

$$\delta Q = T \, dS$$

He then made three linked identifications. Entropy scaled with boundary area.
Heat became energy flux across a local horizon. Temperature became Unruh
temperature, proportional to surface gravity.

He demanded the relation hold for all local horizons.

Einstein's field equations are the geometric form of that requirement:

$$R_{\mu\nu} - \frac{1}{2}R g_{\mu\nu} = 8\pi G T_{\mu\nu}$$

This displayed version omits the cosmological-constant term. In OPH, the local
first-variation relation fixes the curvature response to stress. A vacuum
reference evaluates the metric residue, and global screen capacity separately
proposes a numerical value for the resulting cosmological term.

Finite repair first gives quotient normal forms. One common family of repaired
records then supplies the geometric readout of caps and diamonds, modular flow,
null stress, bounded-interval response, fixed-cap stationarity, small-ball area
variation, continuum control, and the tensor upgrade across local observer
directions. These linked steps turn local thermodynamic agreement into
spacetime dynamics when that common family also carries universal coupling, a
vacuum reference, and independent scale readouts. Construction and
certification of one inhabited family with all these readouts are work in
progress.

Jacobson inverted the logic of physics. Usually we think of gravity as fundamental, implying thermodynamic properties for horizons. Jacobson showed the reverse: **if you assume thermodynamics is fundamental, gravity is derived.**

**On Jacobson's thermodynamic reading, gravity is local thermodynamic equilibrium written geometrically; the usual force-law picture is secondary.**

The force of the argument lies in its austerity. Jacobson does not start with planets tracing curves through a manifold. He starts with heat flow, horizon entropy, and the insistence that the same thermodynamic accounting must work in every infinitesimal causal patch. Einstein's equation is what that insistence looks like when written geometrically.

In plain language, gravity becomes horizon bookkeeping done consistently
everywhere. If every tiny causal patch has to balance heat, entropy, and
temperature in the same way, the spacetime metric has to bend so that the
bookkeeping closes. Curvature is the public face of that accounting rule.

## 11.8 Complexity and the Growth of Interiors

For an eternal black hole in AdS/CFT, the boundary state is thermal and time-independent. But the bulk geometry is not static-the wormhole interior keeps growing.

What dual quantity is growing?

Leonard Susskind proposed: **computational complexity**.

Entropy measures how many states are consistent with observations. Complexity measures how hard it is to prepare a state-how many quantum gates you need.

Complexity keeps growing long after entropy saturates. This gives the interior-growth story a computational reading: the hidden work of preparing the state can keep increasing even when coarse entropy has stopped changing.

## 11.9 Special Relativity from Modular Structure

The Bisognano-Wichmann theorem contains a stunning implication: Lorentz symmetry-the foundation of special relativity-can be tied to the modular structure of the vacuum.

### The Unruh Effect: Where It Begins

In 1976, William Unruh discovered that an accelerating observer sees the vacuum differently. An observer accelerating through empty space sees thermal radiation-a bath of particles at temperature:

$$T_U = \frac{\hbar a}{2\pi c k_B}$$

where a is the acceleration. An inertial observer sees vacuum. An accelerating observer sees heat.

This is an exact result of quantum field theory. The vacuum looks different depending on your state of motion.

Why? Acceleration creates a **Rindler horizon**-a boundary beyond which signals can never reach the accelerating observer. This horizon has thermodynamic properties identical to a black hole horizon. The temperature comes from quantum fluctuations near this horizon.

### The Bisognano-Wichmann Theorem

In 1975-1976, Bisognano and Wichmann proved something deeper. Consider the vacuum state of a quantum field theory. Restrict attention to a Rindler wedge-the region accessible to a forever-accelerating observer.

The wedge modular automorphism is geometric. In a cutoff or pedagogical
type-I representation this is often written as a thermal density matrix:

$$\rho_R = \frac{e^{-2\pi K}}{Z}$$

where K is the Lorentz boost generator. In that representation the modular
Hamiltonian-which generates "time evolution" within the wedge-is proportional
to the boost:

$$H_{mod} = 2\pi K$$

In this wedge case, **the modular automorphism group is the Lorentz boost
flow**.

The wedge example shows how an abstract modular flow becomes a familiar
geometric transformation from relativity. Ordinary local clocks arise by
synchronizing this algebraic ordering across overlapping observer regions; the
theorem supplies the bridge between operator algebra and spacetime motion.

$$\Delta^{it} = e^{-2\pi i K t}$$

The natural modular evolution of the wedge state is exactly a Lorentz
transformation.

That means the same structure that tells the observer "this restricted state is thermal" also tells the observer how boosts and clocks fit together. Thermal language and relativistic geometry are two descriptions of one modular fact.

One structure is doing two jobs at once. Read algebraically, it is the modular evolution of a restricted state. Read geometrically, it is the boost symmetry of the wedge. The two readings agree because the observer's horizon cuts the vacuum in exactly the right way.

### Boosts from Thermal Structure

Start with thermal structure. Ask: what is the natural notion of time evolution? In the wedge setting, the answer is Lorentz boosts.

This reverses the usual logic in QFT. We do not postulate Lorentz symmetry and then discover thermal horizons; the Bisognano-Wichmann theorem shows the boost structure is encoded in modular flow.

That modular-boost link recovers Lorentz kinematics and a universal causal
speed on the screen.

### Connection to OPH

A support cap carries the relevant thermal modular data. Finite screen cells
regulate that chart, while the finite patch federation remains the carrier.
The construction extracts the geometric cap pair, transports its modular data,
and preserves the support behavior and thermal normalization. The cap modular
automorphism becomes geometric on the sphere, and that geometric action gives
the Lorentz symmetry.

### The Speed of Light

Why is there a maximum speed, and why is it the same for everyone?

The Unruh formula T = ℏa/(2πck_B) contains c. For the thermal-to-boost correspondence to work, there must be a universal velocity relating acceleration to temperature.

The conformal geometry of $S^2$ supplies the Lorentz group, whose invariant
speed is $c$ once physical units and the clock calibration are chosen. The
numerical value of $c$ is the conversion between those units; its universality
comes from the shared causal structure. Quantum no-signaling fits that
structure because entanglement alone cannot transmit a controllable message
outside the light cone.

### The Causal Structure

The light-cone structure of spacetime, the question of which events can
influence which, emerges from entanglement. Spacelike-separated regions can be
correlated without signaling. Timelike-separated events can have causal
influence. Null separation marks the dividing line between those two regimes.

The modular flow provides the time direction. Entanglement provides
correlations. No-signaling prevents faster-than-light communication. Together
they reproduce Minkowski causal structure.

### Why This Matters

Einstein discovered special relativity in 1905 by thinking about light and
motion. QFT gives the same structure another reading: Lorentz boosts are tied
to horizon thermodynamics via the Bisognano-Wichmann theorem. In OPH the
Lorentz group appears as the geometry of modular automorphisms on extracted
caps.

The laws of physics look the same to all inertial observers because thermal
states on wedge-shaped regions naturally evolve via boosts. Carrying that
modular-boost structure over to the screen and reading it back into bulk
kinematics gives the universal speed.

## 11.10 What Time Predicts

The thermal-time picture is anchored in physics. Tomita-Takesaki says an
algebra-state pair carries its own flow. The KMS condition gives that flow the
structure of thermal equilibrium. Bisognano-Wichmann turns modular time into a
Lorentz boost, and the corresponding cap flow gives the screen its geometric
Lorentz action. Boltzmann explains why irreversible records emerge out of
reversible microscopic laws.

The physical world fits this picture with surprising loyalty. Accelerating
observers inherit Unruh temperature from the same horizon logic that produces
Hawking radiation. Jacobson's thermodynamic route points toward Einstein's
equation. The microscopic laws are largely time-symmetric, while the arrow of
time rides on a low-entropy record boundary and the thermodynamics of record
making.

---

## 11.11 Memory and Records

Why do we remember the past but not the future?

A **memory** is a physical record-a low-entropy structure correlated with a past event. Creating a record requires work-you must push entropy somewhere else.

When you remember something, you're consulting a present record created at the cost of increasing entropy elsewhere. The record only makes sense if entropy was lower when the recorded event happened.

The arrow of time is the arrow of record-keeping. Time flows in the direction we can make and preserve consistent records.

## 11.12 Reverse Engineering Summary

Time need not be laid down as a primitive external river. General relativity
removes any preferred slicing. Quantum gravity sharpens that loss.
OPH reads time from the inside, through the modular flow attached to a
restricted state. The arrow points in the direction records can be made and
kept. Boltzmann explains why entropy rises. Jaynes explains why ignorance has
structure. Tomita-Takesaki supplies the clock. Bisognano-Wichmann ties that
clock to relativity. Jacobson turns the same thermodynamic language into
gravity.

---

We have located a source of time without putting time in by hand. Incomplete knowledge, restricted access, and record-building are enough to generate clocks and an arrow.

The harder question concerns translation. Different observers inherit different local clocks, different horizons, and different cuts through the state. Why do the conversion rules between their descriptions lock into the rigid form of symmetry and conservation law, with no case-by-case negotiation?

That is where **Chapter 12: Symmetry on the Sphere** begins.
