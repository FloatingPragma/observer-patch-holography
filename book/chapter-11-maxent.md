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

![Known constraints select a least-biased state, and the restricted algebra-state pair carries an intrinsic modular ordering.](../assets/book_diagrams/maxent-clock.svg){width=78%}

## 11.2 The Surprising Hint: Time Is Not Fundamental

### The Scandal of the Second Law

Physics has a scandal.

Almost all our fundamental laws are time-reversible. Newton's F = ma works the same forward and backward, Maxwell's equations run happily in either direction, and so do Schrödinger's equation and Einstein's general relativity.

Film a planet orbiting a star and play it backward; it looks perfectly physical. Film an egg breaking and play it backward, and the result is absurd.

This is the **Arrow of Time**. Where does it come from? The microscopic laws
are the wrong place to look.

### No Preferred Time in GR

In general relativity, there's no preferred time coordinate. Different observers slice spacetime differently; none is privileged.

The Wheeler-DeWitt equation, the analog of Schrödinger's equation for the universe, is:

$$H\Psi = 0$$

The Hamiltonian acting on the wavefunction of the universe gives zero. There is no explicit external time derivative in this formalism, so the universe can look *frozen* at the fundamental level.

$H$ is the Hamiltonian constraint, the operator that would normally generate
time evolution. $\Psi$ is the wavefunction of the universe in this formal
setting. Nothing in the equation forbids things from happening in experience.
It says only that the fundamental constraint carries no outside time
parameter.

This is the **problem of time** in quantum gravity. If the fundamental description has no time, where does time come from?

Time enters this story without a fundamental external clock. The microscopic laws are time-symmetric. Something else must generate the arrow of time we experience.

## 11.3 The First-Principles Reframing: From Modular Flow to Physical Time

The deeper question is why we experience time at all if the fundamental
description has no preferred clock.

### The Thermal Time principle

In 1994, Alain Connes, a French mathematician who had won the Fields Medal
for his work on operator algebras, and Carlo Rovelli, an Italian
quantum-gravity theorist, proposed a stark idea: time can be read from
incomplete knowledge.

The intuition rests on one fact of ordinary thermodynamics. Boltzmann's
distribution says that in a thermal state, probability falls off
exponentially with energy: cheap states are common, expensive states are
rare. Run that backward. Hand someone only the probabilities, and by
taking a logarithm they can recover the energy ladder that produced them. And
in quantum mechanics, energy is precisely the thing that generates time
evolution. So a state of partial knowledge secretly contains an energy
ranking, and an energy ranking contains a flow.

The Connes-Rovelli recipe is that inversion written in one line. Start with
the observer's limited state $\rho$, the quantum bookkeeping object for what
the observer can know, and write the modular Hamiltonian $K=-\ln\rho$. Take
the logarithm of what you know, and what comes out generates a flow. The
thermal-time proposal reads that flow as time; OPH treats it first as
intrinsic ordering. The modular Hamiltonian belongs to the observer's
restricted state, and it need not equal the ordinary energy of the whole
universe.

This is a strange move the first time one sees it. In ordinary mechanics, the
Hamiltonian is given first and time evolution follows. Here the restricted
state furnishes a preferred dimensionless ordering. A physical clock requires
an observer-readable transition, an event correspondence, and a calibration to
turn that ordering into a duration. Time is therefore tied to access and
coarse graining without being reduced to mere ignorance.

### Tomita-Takesaki Theory

The deeper theorem behind that proposal comes from operator algebra. Once an
observer has a rich enough algebra of questions and a state that probes it
fully, the pair carries a preferred internal flow whether or not anyone inserts
an external master clock. The formal machinery is called
**Tomita-Takesaki theory**. Minoru Tomita announced the core theorem in 1967,
in notes so compressed that much of the field could not follow them;
Masamichi Takesaki's 1970 lecture notes rebuilt the proofs and made the
theory usable.

The name is heavier than the picture. Give an observer a complete menu of
accessible questions and a sufficiently informative state. The pair contains
its own natural ordering of those questions. No clock has appeared, though the
raw material for one has.

An automorphism is a reshuffling of the allowed questions that preserves their
algebraic rules. Modular flow is a continuous family of such reshufflings,
indexed by a parameter that behaves like time. The concrete content is that the
observer's restricted state tells the question menu how to move.

An observer with partial access does not sit in an algebraic fog. The
restriction supplies an order of possible changes. The flow depends on the
algebra-state pair, which is why different observers can inherit different
modular orderings from different access conditions. It also carries the thermal
equilibrium structure that links temperature to the flow parameter.

Modular flow matters here because it supplies an internal source of temporal
ordering. The observer's horizon and state determine that flow.

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

In 1976, William Unruh discovered that an accelerating observer sees the
vacuum differently. An observer accelerating through empty space sees a
thermal bath at temperature:

$$T_U = \frac{\hbar a}{2\pi c k_B}$$

$T_U$ is the temperature seen by the accelerating observer. $a$ is the
observer's proper acceleration. The constants are the same ones used in the
black-hole temperature formulas. The larger the acceleration, the hotter the
restricted vacuum appears. An inertial observer sees vacuum. An accelerating
observer sees heat. This is an exact result of quantum field theory: the
vacuum looks different depending on your state of motion.

Acceleration creates a **Rindler horizon**, a boundary beyond which signals
can never reach the accelerating observer. The region visible to that
observer is the **Rindler wedge**. The horizon has thermodynamic properties
identical to a black hole horizon, and the temperature comes from quantum
fluctuations near it.

In 1975 and 1976, Joseph Bisognano and Eyvind Wichmann, working at Berkeley,
proved something deeper about this wedge. For the vacuum state restricted to
the wedge, the modular automorphism is geometric. In the simplified version
where the wedge has a density matrix (a regularized description, good enough
for intuition), the state is thermal:

$$\rho_R = \frac{e^{-2\pi K}}{Z}$$

where $K$ is the Lorentz boost generator. In that simplified form the modular
Hamiltonian, which generates "time evolution" within the wedge, is
proportional to the boost:

$$H_{mod} = 2\pi K$$

For a uniformly accelerating observer, boost motion supplies the relevant
time translation, and under the theorem's relativistic hypotheses the modular
flow follows that motion. In this wedge case, **the modular automorphism
group is the Lorentz boost flow**.

The Unruh effect is Tomita-Takesaki theory applied inside relativistic
spacetime. The restricted vacuum state carries the boost flow seen by the
accelerating observer. Limited access therefore has thermodynamic
consequences, and it supplies the ordering used by the clock construction.

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

**The Past Hypothesis as a consistency requirement**: For observers to exist at all, they must be able to form and compare records. Records require entropy gradients; writing information pushes entropy elsewhere. A universe in thermal equilibrium contains no observers and no records for them to compare.

The specific numerical entropy of the hot dense record belongs to physical
cosmology. The OPH point is the consistency role: observers who compare records
require a low-entropy side, and the arrow of time points in the direction that
allows records to be made.

## 11.5 Jaynes: Entropy as Ignorance

In 1957, Edwin Jaynes published a pair of papers titled "Information Theory
and Statistical Mechanics" and rewrote the subject.

**Entropy measures our knowledge about the gas. The gas itself carries no such number.**

### The Maximum Entropy Principle

Suppose you know only the average energy. What probability distribution should you assign?

Choose the distribution that maximizes Shannon entropy subject to your constraints:

$$S = -\sum_i p_i \ln p_i$$

Here $p_i$ is the probability of outcome $i$.

MaxEnt gives the Boltzmann distribution:

$$P(x) = \frac{1}{Z} e^{-\beta E(x)}$$

Thermal states are ubiquitous because they're the unique states of maximum ignorance given energy constraints.

In the Boltzmann distribution, $P(x)$ is the probability of state $x$, $E(x)$
is its energy, $\beta$ is inverse temperature, and $Z$ is the partition
function that normalizes all probabilities so they add to 1.

Jaynes liked to demonstrate the principle on a die. Told only that the
average roll is 4.5 instead of the fair 3.5, MaxEnt assigns each face the
least dramatic bias consistent with that single fact, and nothing more. No
story about how the die was loaded, no extra assumptions smuggled in: just
the flattest distribution the evidence permits.

### The Exact OPH Rule

In ordinary language, the first two OPH axioms draw the boundary of what the
observer can consistently know. The third chooses the least opinionated state
inside that boundary. Agreement supplies the shape; maximum randomness fills
the room without hiding extra furniture in it.

OPH applies Jaynes locally. At one finite resolution, a state is a
compatible family of patch states, and the first two axioms fix the set of
families the observer's constraints allow. The third axiom picks, from that
set, the family closest to a fixed reference: the realized state minimizes

$$
\mathcal D(\rho\Vert\tau)
=
\sum_{P}w_P\,D(\rho_P\Vert\tau_P).
$$

Here $D$ is relative entropy, a standard measure of how far one probability
assignment sits from another, computed patch by patch with weights $w_P$ over
the observer's patches. When the reference is flat, this is plain local
entropy maximization: the least opinionated state the constraints allow. It
is Jaynes's die, played across a federation of patches at once.

The finite matrix formalization has an exact first step and a sharp warning.
For any projective record partition, erasing cross-record coherences is the
uniform average over all independently signed block reflections. But a
matrix logarithm that simply declares $\log 0=0$ cannot define the required
support-aware relative entropy: it assigns zero to two orthogonal pure states
even though one state's support is absent from the other. The full spectral
majorization, extended-entropy, and information-projection package therefore
needs an explicit support layer; it cannot be obtained by silently
totalizing the logarithm.

### Where the Four Laws Come From

Chapter 4 treated thermodynamics as a hint, a set of laws discovered from
engines and gases and then found to be about information. At one finite
regulator stage, the third axiom yields a conditional theorem package only
after the weighted local objective has a global representation, one faithful
source reference is shared by the state and transition problems, and the
active source collar realizes the complete repaired-visible fibre. Physical
energy and clock calibration and refinement-uniform low-temperature control
are separate receipts, not consequences of that finite package.

Read the axiom on states with its supplied faithful reference and it selects
the Gibbs family, the exponential
weighting by energy that a thermodynamics course normally writes down as an
assumption. Read the same axiom on transitions, on the repair step that
carries one round of consensus into the next, and it selects a single repair
kernel: a map that leaves the reference alone, changes nothing when applied
twice, and preserves the quantities the patches agree on.

The laws follow from those two readings. Distance from the reference, measured
as relative entropy, shrinks under repair and never grows. That contraction is
the Second Law, and the Clausius inequality relating entropy change to heat
falls out of it, as does Landauer's price for erasing a bit. Two systems
placed in contact settle at a common inverse temperature, which is the Zeroth
Law. The energy bookkeeping splits exactly into a heat term and a work term,
which is the First Law. The weight the construction can place on excited
states is bounded, so entropy has a floor set by the number of ground states,
and no finite sequence of repair steps reaches it. That is the Third Law.

The contraction needs stationarity rather than microscopic reversibility. A
lazy walk around a three-stop one-way loop gives an exact counterexample: its
uniform distribution stays fixed and relative entropy cannot increase, even
though the forward and reverse equilibrium fluxes differ. Entropy descent
therefore carries less information than detailed fluctuation symmetry or
Onsager reciprocity. Those stronger symmetries need their own physical
evidence.

The source-counted repair table was checked across all fifteen nonempty
field-subset projections of its four recorded fields. One eight-state
repair-load description forms an irreducible
stationary chain and moves toward its stationary distribution, while retaining
a measurable directional imbalance. The finer description has only one closed
class, a single freezeout state. Its record-family label takes one value, so
the apparent conservation test is empty. This table supports the finite
stationary-kernel calculation. It supplies neither a common state-and-transition
reference nor a physical energy or clock calibration.

The separately preregistered state table and recurrent transition chain cannot
be joined by the hoped-for common-object map. The state-side resampling step is
idempotent, while the transition chain has a nonconstant mode multiplied by
$665437/726948$, strictly between zero and one. Any map intertwining these two
steps must erase that mode. Independently, the transition stationary mass
$7155/61511$ falls between two adjacent masses available from the $16{,}384$
equally weighted state observations, so no deterministic regrouping of that
table produces it. An invented random coupling would not be evidence from the
run. This is a final negative result for the current artifact, not a no-go for
all possible source dynamics.

The finite theorem package therefore does not import thermodynamic identities
as an extra axiom, but its physical interpretation is conditional. Five
requirements remain explicit for any replacement source object: a global
representation of the weighted objective;
one source-derived reference shared by both optimizers; realization of the
conditional-resampling kernel on the complete repaired-visible fibre by the
source-collar transition matrix; physical energy and clock calibration; and
refinement-uniform low-temperature control for a continuum third law.
Stochasticity, stationarity, and charge preservation alone establish only the
weaker stationary H-theorem branch and do not close the third receipt. Until
all five are supplied, this is not a thermometer reading.

## 11.6 Time on the Holographic Screen

In the simplest finite illustration, a support region $P$ on the screen
chart is cut from a global pure state. The restriction gives a density matrix:

$$\rho_P = \text{Tr}_{\bar{P}} |\Psi\rangle \langle \Psi|$$

This density matrix defines a modular Hamiltonian:

$$K_P = -\ln \rho_P$$

which generates the finite modular flow labeled by $t_P$. That flow supplies
an intrinsic ordering parameter for the observer's accessible algebra.

$\bar P$ means the complement of support region $P$, everything outside the
region in this finite illustration. The trace over $\bar P$ discards
inaccessible degrees of freedom and leaves the state available to the observer.
The logarithm then turns that restricted state into the modular generator.

This density-matrix formula is the finite illustration. The general observer
patch is described by its accessible algebra-state pair, and that pair
carries a modular flow even when no density matrix exists (in the full
continuum theory, none does; the algebra-state pair is the object that
survives, which is why Tomita-Takesaki earns its keep). The flow orders
algebraic change; the clock construction turns that order into time.

### Consistency of Clocks

For a shared clock, two overlapping observers need compatible modular flows
on what they share, calibrated instruments, and an agreed correspondence
between their events. When those conditions hold, a common causal structure
follows.

### Cosmic Time

Why do we all agree on a "cosmic time"?

On a cosmological branch where those clock conditions hold across observers,
local modular flows can support a shared coarse-grained cosmic time. It is a
collective clock read by the reconstructed world from within, rather than a
second fundamental time parameter.

### From Modular Time to Gravity

The full argument is a chain of cross-checks. The finite record structure
must independently hand over its geometry, its flow, and its energy
bookkeeping, and a single theorem then consumes all three and identifies the
internal flow with a motion on the sphere. From there the path to Einstein's
equation runs exactly as Jacobson mapped it.

## 11.7 Jacobson's Derivation

In 1995, Ted Jacobson published a short paper whose title gives the game away: "Thermodynamics of Spacetime: The Einstein Equation of State." It contains one of the most beautiful derivations in theoretical physics.

He started with the first law of thermodynamics:

$$\delta Q = T \, dS$$

In words, the heat flowing in equals temperature times the entropy change it
buys.

He then made three linked identifications. Entropy scaled with boundary area.
Heat became energy flux across a local horizon. Temperature became Unruh
temperature, proportional to surface gravity.

He demanded the relation hold for all local horizons.

Einstein's field equations are the geometric form of that requirement:

$$R_{\mu\nu} - \frac{1}{2}R g_{\mu\nu} = 8\pi G T_{\mu\nu}$$

The left side is curvature, the right side is matter and energy; the equation
says one is the ledger of the other.

This displayed version omits the cosmological-constant term. In OPH, the same
local argument fixes how curvature responds to stress, and the screen's
global capacity supplies the scale of the cosmological term.

Jacobson inverted the logic of physics. Usually we think of gravity as fundamental, implying thermodynamic properties for horizons. Jacobson showed the reverse: **if you assume thermodynamics is fundamental, gravity is derived.**

The force of the argument lies in its austerity. Jacobson does not start with planets tracing curves through a manifold. He starts with heat flow, horizon entropy, and the insistence that the same thermodynamic accounting must work in every infinitesimal causal patch. Einstein's equation is what that insistence looks like when written geometrically.

In plain language, gravity becomes horizon bookkeeping done consistently
everywhere. If every tiny causal patch has to balance heat, entropy, and
temperature in the same way, the spacetime metric has to bend so that the
bookkeeping closes. Curvature is the public face of that accounting rule.

In the OPH version, one common family of repaired records has to supply every
readout the argument needs: the geometry, the modular flow, and the energy
bookkeeping. When that family also carries a universal coupling and a vacuum
reference, local thermodynamic agreement turns into spacetime dynamics.
Building one family with all these readouts in place is work in progress.

## 11.8 Complexity and the Growth of Interiors

For an eternal black hole in AdS/CFT, the boundary state is thermal and
time-independent, yet the bulk geometry keeps changing as the wormhole
interior grows. What dual quantity is growing? Leonard Susskind proposed
computational complexity, the number of quantum gates needed to prepare the
state. Complexity keeps growing long after entropy saturates. In record
terms, the hidden work of preparing a state can keep accumulating even when
its coarse ledger has stopped changing, which is what the cost of
record-making looks like from the inside.

## 11.9 Special Relativity from Modular Structure

The Bisognano-Wichmann theorem contains a stunning implication: Lorentz
symmetry, the foundation of special relativity, can be tied to the modular
structure of the vacuum.

Section 11.3 laid out the wedge result: the natural modular evolution of the
vacuum restricted to a Rindler wedge is exactly a Lorentz boost. In ordinary
quantum field theory the theorem is proved inside a relativistic theory, so
it cannot double as a derivation of relativity; the circularity would be
visible from orbit. OPH instead builds a finite cap state without assuming
the geometry, then checks that it carries the same modular structure.
Ordinary local clocks then require calibrated transports across overlapping
observer regions.

One structure is doing two jobs at once. Read algebraically, it is the modular evolution of a restricted state; read geometrically, it is the boost symmetry of the wedge. The same fact that tells the observer "this restricted state is thermal" also tells the observer how boosts and clocks fit together, because the horizon cuts the vacuum in exactly the right way.

### Boosts from Thermal Structure

Start from thermal structure and ask what time evolution naturally means; in
the wedge setting, the answer is a Lorentz boost. The ordinary theorem reads
that as boost structure encoded in modular flow. The OPH reconstruction
reverses the reading: build the cap states first, without assuming
relativity, and the modular-boost link then supplies Lorentz kinematics and a
universal causal structure on the screen.

### Connection to OPH

A cap on the screen carries thermal modular data. Finite screen cells chart
it, while the finite patch federation remains the carrier. Transport that
modular data across the sphere and it becomes geometric: the cap's modular
flow acts as a motion on the sphere, and that geometric action is the
Lorentz symmetry.

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

The light-cone structure of spacetime answers which events can influence which.
In established relativistic theory, spacelike-separated regions can be
correlated without signaling, timelike-separated events can have causal
influence, and null separation marks the boundary. OPH has to reconstruct that
cone rather than borrow it.

The matched modular flow supplies an oriented cap motion. Entanglement
supplies correlations, and no-signaling forbids faster-than-light
communication. With the event and cone construction, these ingredients yield
Minkowski causal structure.

### Why This Matters

Einstein discovered special relativity in 1905 by thinking about light and
motion. Quantum field theory gives the same structure a second reading:
boosts tied to horizon thermodynamics through the Bisognano-Wichmann theorem.
OPH gives it a third, as the geometry of matched modular flows on the sphere,
with one universal causal speed forced by the shared structure.

## 11.10 Memory and Records

Why do we remember the past but not the future?

A **memory** is a physical record, a low-entropy structure correlated with a past event. Creating a record requires work and pushes entropy somewhere else.

When you remember something, you're consulting a present record created at the cost of increasing entropy elsewhere. The record only makes sense if entropy was lower when the recorded event happened.

The arrow of time is the arrow of record-keeping. Time flows in the direction we can make and preserve consistent records.

A record, once made, is a thing in its own right. It can outlast everything
about the moment that wrote it.

## 11.11 Reverse Engineering Summary

Time need not be laid down as a primitive external river. General relativity
removes any preferred slicing, and quantum gravity sharpens the loss into the
frozen equation $H\Psi = 0$. OPH builds time from the inside instead: an
observer's restricted state carries its own flow (Tomita-Takesaki), the flow
carries a temperature (the KMS condition, the algebraic signature of thermal
equilibrium), and the clock construction turns flow into physical time, with
the arrow pointing in the direction records can be made and kept. Boltzmann
explains why entropy rises, Jaynes explains why ignorance has structure,
Bisognano-Wichmann matches modular flow to Lorentz motion, and Jacobson turns
the same thermodynamic language into gravity. The physical world fits this
picture with surprising loyalty: accelerating observers inherit Unruh
temperature from the same horizon logic that produces Hawking radiation.

---

We have located an internal ordering without putting an external time
parameter in by hand. Restricted access and record-building orient that
ordering and give it an arrow, and the clock construction turns it into
physical time.

The harder question concerns translation. Different observers inherit different local clocks, different horizons, and different cuts through the state. Why do the conversion rules between their descriptions lock into the rigid form of symmetry and conservation law, with no case-by-case negotiation?

That is where **Chapter 12: Symmetry on the Sphere** begins.
