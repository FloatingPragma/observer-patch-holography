# Chapter 11: MaxEnt and the Arrow

## 11.1 The Intuitive Picture: Time Is Fundamental

Start with the Newtonian picture of time.

Time is a fundamental external parameter. It flows from past to future,
independent of anything in the universe. Events happen in time, just as objects
exist in space. The clock ticks whether or not anything is happening. Time is
the stage; physics is the play.

This is Newton's absolute time: "Absolute, true, and mathematical time, of itself, and from its own nature, flows equably without relation to anything external."

The arrow of time, the fact that we remember yesterday but not tomorrow, that eggs break but don't unbreak, seems fundamental to this picture. Time has a direction, built into its very nature.

General relativity and quantum theory broke that picture.

## 11.2 The Surprising Hint: Time Is Not Fundamental

### The Scandal of the Second Law

Physics has a scandal.

Almost all our fundamental laws are time-reversible. Newton's F = ma works the same forward and backward. Maxwell's equations are reversible. Schroedinger's equation is reversible. Einstein's General Relativity is reversible.

Film a planet orbiting a star and play it backward-it looks perfectly physical. But film an egg breaking and play it backward? Absurd.

This is the **Arrow of Time**. Where does it come from? It's not in the microscopic laws.

### No Preferred Time in GR

In general relativity, there's no preferred time coordinate. Different observers slice spacetime differently; none is privileged.

The Wheeler-DeWitt equation-the analog of Schroedinger's equation for the universe-is:

$$H\Psi = 0$$

The Hamiltonian acting on the wavefunction of the universe gives zero. There is no explicit external time derivative in this formalism, so the universe can look *frozen* at the fundamental level.

This is the **problem of time** in quantum gravity. If the fundamental description has no time, where does time come from?

Time is not a fundamental external parameter. The microscopic laws are time-symmetric. Something else must generate the arrow of time we experience.

## 11.3 The First-Principles Reframing: Observer-Relative Time from Modular Flow

The deeper question is why we experience time at all if the fundamental
description has no preferred clock.

### The Thermal Time principle

In the 1990s, Alain Connes and Carlo Rovelli proposed a stark idea. Time can be
read from incomplete knowledge. Start with the observer's limited state
$\rho$. From it one forms the modular Hamiltonian, $K=-\ln\rho$. That operator
generates a flow, and the thermal-time proposal reads that flow as a candidate
observer-relative time parameter.

This is a strange move the first time one sees it. In ordinary mechanics, the
Hamiltonian is given first and time evolution follows. Here the restricted
state itself furnishes a candidate internal clock. Time is tied to access, ignorance, and coarse
graining.

### Tomita-Takesaki Theory

The deeper theorem behind that proposal comes from operator algebra. Once an
observer has a rich enough algebra of questions and a state that probes it
fully, the pair carries a preferred modular flow whether or not anyone inserts
an external master clock. The formal machinery is called
**Tomita-Takesaki theory**.

An observer with partial access does not sit in a timeless fog. In the thermal-time
reading, the restriction itself orders experience into before and after. The flow depends on
the algebra-state pair, which is why different observers can inherit different
candidate local parameters from different access conditions. It also carries the thermal
equilibrium structure that makes temperature and time appear together in the
same move.

Modular flow matters here because it turns observer-relative time from a background stage into
something read from an observer's horizon and state. Partial knowledge has its own dynamics.

### The Rindler Wedge

This abstract mathematics connects to reality through the Unruh effect.

An observer accelerating uniformly sees only the **Rindler wedge**-part of spacetime. For the vacuum state restricted to this region, the Bisognano-Wichmann theorem shows that the modular Hamiltonian is exactly the generator of Lorentz boosts.

For an accelerating observer, a Lorentz boost *is* time translation. The modular flow equals ordinary time evolution.

The modular temperature works out to:

$$T_{Unruh} = \frac{\hbar a}{2\pi k_B c}$$

The Unruh effect isn't a separate phenomenon-it's Tomita-Takesaki theory applied to spacetime. In the thermal-time reading, the observer's restricted access selects a candidate observer-relative time parameter.

This is the first major payoff of the chapter. The mathematics does not stay abstract. Restrict the vacuum to what one observer can access, and the restriction behaves like a thermal state with its own clock.

This is also the point where modular theory stops sounding like rarefied operator algebra and starts sounding like lived physics. Restriction generates a temperature and a candidate modular flow; reading that flow as time is the conditional thermal-time/BW-branch interpretation. Losing access to part of the state has thermodynamic and temporal consequences.

## 11.4 The Arrow of Time

In Chapter 4, we saw Boltzmann's insight: entropy $S = k \ln W$ measures the number of microstates compatible with a macrostate, and entropy increases because high-entropy states vastly outnumber low-entropy ones.

But why did the universe start with low entropy in the first place?

### The Past principle

The deeper answer to the arrow of time is the **Past principle**: the universe began in a state of extraordinarily low entropy.

We're not riding a random fluctuation. We're riding the expansion from a very special initial condition-the Big Bang. The early universe was hot but smooth, with matter spread almost uniformly. That uniformity is low gravitational entropy.

Why was the Big Bang low entropy? Standard physics treats this as an unexplained initial condition. But our model offers a different perspective.

**The Past principle as a consistency requirement**: For observers to exist at all, they must be able to form and compare records. Records require entropy gradients-you can only write information by pushing entropy elsewhere. A universe in thermal equilibrium contains no observers, no records, no consistency-checking.

The MaxEnt principle says: assign the maximum-entropy state consistent with your constraints. But one constraint is that someone must exist to apply MaxEnt. This rules out equilibrium. The very existence of observers selecting MaxEnt states presupposes a universe far from equilibrium.

This doesn't derive the specific numerical entropy of the Big Bang. But it reframes the question: the Past principle isn't an arbitrary input to be explained by some deeper theory. It can be read as a consistency requirement. A universe containing observers who check for consistency appears to require a low-entropy beginning. The arrow of time points in the direction that allows records to be made.

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

## 11.6 Time on the Holographic Screen

Each observer has a patch $P$ on the holographic screen. The global state
restricts to a density matrix:

$$\rho_P = \text{Tr}_{\bar{P}} |\Psi\rangle \langle \Psi|$$

This density matrix defines a modular Hamiltonian:

$$K_P = -\ln \rho_P$$

which generates a modular parameter $t_P$ for that observer. In OPH this parameter is read as
observer-relative time only on the conditional BW geometric branch.

**Every observer has a candidate local modular parameter, not a standalone operational-clock theorem.**

### Consistency of Clocks

If two observers' patches overlap, their BW-branch observer-relative modular parameters must be compatible on the overlap. This is a strong constraint. Reality hangs together because the modular flows mesh.

### Cosmic Time

Why do we all agree on a "cosmic time"?

If the global state is highly entangled in a particular pattern, the modular flows of local patches are synchronized. In this conditional reading, cosmic time is an effective synchronization of local BW-branch modular times rather than a new fundamental clock.

### Roadmap: From Modular Time to Gravity

The chain is clean once the pieces are visible. Recovery structure from Chapter
7 makes the time generator local near patch boundaries. A key theorem then
identifies that local flow with a standard geometric motion on the sphere and
fixes its normalization. BW-branch geometric modular flow gives Lorentz kinematics on the
screen, and entanglement equilibrium together with the local energy bridge
yields Einstein's equation.

This chapter builds the time ingredient. The next sections show how it feeds into gravity.

## 11.7 Jacobson's Derivation

In 1995, Ted Jacobson performed one of the most beautiful derivations in theoretical physics.

He started with thermodynamics-the first law:

$$\delta Q = T \, dS$$

He then made three linked identifications. Entropy scaled with boundary area.
Heat became energy flux across a local horizon. Temperature became Unruh
temperature, proportional to surface gravity.

He demanded the relation hold for all local horizons.

Out popped **Einstein's field equations**:

$$R_{\mu\nu} - \frac{1}{2}R g_{\mu\nu} = 8\pi G T_{\mu\nu}$$

Jacobson inverted the logic of physics. Usually we think of gravity as fundamental, implying thermodynamic properties for horizons. Jacobson showed the reverse: **if you assume thermodynamics is fundamental, gravity is derived.**

**On Jacobson's thermodynamic reading, gravity is not fundamental in the usual force-law sense; it is what local thermodynamic equilibrium looks like geometrically.**

The force of the argument lies in its austerity. Jacobson does not start with planets tracing curves through a manifold. He starts with heat flow, horizon entropy, and the insistence that the same thermodynamic accounting must work in every infinitesimal causal patch. Einstein's equation is what that insistence looks like when written geometrically.

Put less formally, gravity becomes horizon bookkeeping done consistently everywhere. If every tiny causal patch has to balance heat, entropy, and temperature in the same way, the spacetime metric has to bend so that the bookkeeping closes. Curvature is the public face of that accounting rule.

## 11.8 Complexity and the Growth of Interiors

For an eternal black hole in AdS/CFT, the boundary state is thermal and time-independent. But the bulk geometry is not static-the wormhole interior keeps growing.

What dual quantity is growing?

Leonard Susskind proposed: **computational complexity**.

Entropy measures how many states are consistent with observations. Complexity measures how hard it is to prepare a state-how many quantum gates you need.

Complexity keeps growing long after entropy saturates. One continuation program relates interior growth to quantum complexity. This is a suggestive conjectural bridge, not a settled derivation of cosmic expansion.

## 11.9 Special Relativity from Modular Structure

The Bisognano-Wichmann theorem contains a stunning implication: Lorentz symmetry-the foundation of special relativity-can be tied to the modular structure of the vacuum.

### The Unruh Effect: Where It Begins

In 1976, William Unruh discovered that an accelerating observer sees the vacuum differently. An observer accelerating through empty space sees thermal radiation-a bath of particles at temperature:

$$T_U = \frac{\hbar a}{2\pi c k_B}$$

where a is the acceleration. An inertial observer sees vacuum. An accelerating observer sees heat.

This isn't a quirk or approximation. It's an exact result of quantum field theory. The vacuum looks different depending on your state of motion.

Why? Acceleration creates a **Rindler horizon**-a boundary beyond which signals can never reach the accelerating observer. This horizon has thermodynamic properties identical to a black hole horizon. The temperature comes from quantum fluctuations near this horizon.

### The Bisognano-Wichmann Theorem

In 1975-1976, Bisognano and Wichmann proved something deeper. Consider the vacuum state of a quantum field theory. Restrict attention to a Rindler wedge-the region accessible to a forever-accelerating observer.

The reduced density matrix on this wedge turns out to be thermal:

$$\rho_R = \frac{e^{-2\pi K}}{Z}$$

where K is the Lorentz boost generator. The modular Hamiltonian-which generates "time evolution" within the wedge-is proportional to the boost:

$$H_{mod} = 2\pi K$$

Here's the punchline: **modular flow IS Lorentz boost** (in QFT wedges).

$$\Delta^{it} = e^{-2\pi i K t}$$

The natural time evolution of a thermal state in a wedge-shaped region is exactly a Lorentz transformation.

That means the same structure that tells the observer "this restricted state is thermal" also tells the observer how boosts and clocks fit together. Thermal language and relativistic geometry are two descriptions of one modular fact.

One structure is doing two jobs at once. Read algebraically, it is the modular evolution of a restricted state. Read geometrically, it is the boost symmetry of the wedge. The two readings agree because the observer's horizon cuts the vacuum in exactly the right way.

### Boosts from Thermal Structure

Start with thermal structure. Ask: what is the natural notion of time evolution? The answer is Lorentz boosts.

This reverses the usual logic in QFT. We don't postulate Lorentz symmetry and then discover thermal horizons; the BW theorem shows the boost structure is encoded in modular flow.

That modular-boost link is the route by which Lorentz kinematics and a
universal light speed are recovered on the screen.

### Connection to Our Framework

Each observer's patch has a boundary, that boundary carries a horizon
temperature, and the modular flow of the horizon state supplies a candidate
observer-relative evolution parameter. Carried over from wedges in ordinary spacetime to
caps on the holographic screen, that flow becomes an actual geometric motion on
the sphere. Once that happens, the conformal symmetry of the sphere reproduces
Lorentz symmetry.

### The Speed of Light

Why is there a maximum speed, and why is it the same for everyone?

The Unruh formula T = ℏa/(2πck_B) contains c. For the thermal-to-boost correspondence to work, there must be a universal velocity relating acceleration to temperature.

From the boundary perspective: information propagates on the S² screen at a maximum rate determined by the entanglement structure. This rate, translated to the bulk, becomes c. The no-signaling theorem of quantum mechanics (entanglement can't transmit information) becomes, in the bulk, the statement that nothing travels faster than light.

### The Causal Structure

The light-cone structure of spacetime, the question of which events can
influence which, emerges from entanglement. Spacelike-separated regions can be
correlated without signaling. Timelike-separated events can have causal
influence. Null separation marks the dividing line between those two regimes.

On the BW geometric branch, modular flow provides the observer-relative time direction. Entanglement provides correlations. No-signaling prevents faster-than-light communication. Taken together, these ingredients reproduce the Minkowski-style causal structure targeted by the program.

### Why This Matters

Einstein discovered special relativity in 1905 by thinking about light and motion. Over a century later, we see it differently: in QFT, Lorentz boosts are tied to horizon thermodynamics via the Bisognano-Wichmann theorem. In our model the same pattern appears when the screen reaches its smooth geometric limit, so the Lorentz group shows up as the geometry of modular flow on caps.

The laws of physics look the same to all inertial observers because thermal states on wedge-shaped regions naturally evolve via boosts. In the OPH program, the universal speed emerges when that modular-boost structure is carried over to the screen and then read back into bulk kinematics.

## 11.10 What Time Predicts

The thermal-time picture does not float free of physics. Tomita-Takesaki says
an algebra-state pair carries its own flow. The KMS condition gives that flow
the structure of thermal equilibrium. Bisognano-Wichmann shows that modular
time becomes an actual Lorentz boost in the wedge setting. Boltzmann explains
why irreversible records emerge out of reversible microscopic laws.

The physical world fits this picture with surprising loyalty. Accelerating
observers inherit Unruh temperature from the same horizon logic that produces
Hawking radiation. Jacobson's thermodynamic route points toward Einstein's
equation. The microscopic laws are largely time-symmetric, while the arrow of
time rides on low-entropy initial conditions and the thermodynamics of record
making.

---

## 11.11 Memory and Records

Why do we remember the past but not the future?

A **memory** is a physical record-a low-entropy structure correlated with a past event. Creating a record requires work-you must push entropy somewhere else.

When you remember something, you're consulting a present record created at the cost of increasing entropy elsewhere. The record only makes sense if entropy was lower when the recorded event happened.

The arrow of time is the arrow of record-keeping. Time flows in the direction we can make and preserve consistent records.

## 11.12 Reverse Engineering Summary

Time does not need to be laid down as a primitive external river. General
relativity removes any preferred slicing. Quantum gravity sharpens that loss.
OPH reads observer-relative time from the inside, through the BW-branch modular flow attached to a
restricted state. The arrow points in the direction records can be made and
kept. Boltzmann explains why entropy rises. Jaynes explains why ignorance has
structure. Tomita-Takesaki supplies the modular flow. Bisognano-Wichmann ties that
flow to relativity on the geometric branch. Jacobson shows how the same thermodynamic language leans
toward gravity.

---

We’ve located a conditional source of observer-relative time without putting time in by hand. Incomplete knowledge, restricted access, and record-building supply candidate modular parameters and an arrow; they do not by themselves prove a universal operational clock.

The harder question concerns translation. Different observers inherit different local modular parameters, different horizons, and different cuts through the state. Why do the conversion rules between their descriptions lock into the rigid form of symmetry and conservation law instead of dissolving into case-by-case negotiation?

That is where **Chapter 12: Symmetry on the Sphere** begins.
