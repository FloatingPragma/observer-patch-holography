# Chapter 11: MaxEnt and the Arrow

## 11.1 The Intuitive Picture: Time Is Fundamental

Before we examine what physics discovered, let's articulate what seemed obvious for millennia.

**The intuitive picture**: Time is a fundamental external parameter. It flows from past to future, independent of anything in the universe. Events happen in time, just as objects exist in space. The clock ticks whether or not anything is happening. Time is the stage; physics is the play.

This is Newton's absolute time: "Absolute, true, and mathematical time, of itself, and from its own nature, flows equably without relation to anything external."

The arrow of time, the fact that we remember yesterday but not tomorrow, that eggs break but don't unbreak, seems fundamental to this picture. Time has a direction, built into its very nature.

And yet, nature gave us hints that shattered this picture.

## 11.2 The Surprising Hint: Time Is Not Fundamental

### The Scandal of the Second Law

Physics has a scandal.

Almost all our fundamental laws are time-reversible. Newton's F = ma works the same forward and backward. Maxwell's equations are reversible. Schrodinger's equation is reversible. Einstein's General Relativity is reversible.

Film a planet orbiting a star and play it backward-it looks perfectly physical. But film an egg breaking and play it backward? Absurd.

This is the **Arrow of Time**. Where does it come from? It's not in the microscopic laws.

### No Preferred Time in GR

In general relativity, there's no preferred time coordinate. Different observers slice spacetime differently; none is privileged.

The Wheeler-DeWitt equation-the analog of Schrodinger's equation for the universe-is:

$$H\Psi = 0$$

The Hamiltonian acting on the wavefunction of the universe gives zero. There's no time derivative. The universe, at the fundamental level, is *frozen*.

This is the **problem of time** in quantum gravity. If the fundamental description has no time, where does time come from?

**This is the hint**: Time is not a fundamental external parameter. The microscopic laws are time-symmetric. Something else must generate the arrow of time we experience.

## 11.3 The First-Principles Reframing: Time Emerges from Modular Flow

Now we reverse engineer. Why do we experience time if it's not fundamental?

### The Thermal Time principle

In the 1990s, Alain Connes and Carlo Rovelli proposed that time emerges from statistical mechanics-from our incomplete knowledge of the quantum state.

The logic:
1. We have a quantum system described by an algebra of observables
2. We have a state rho (a density matrix representing our knowledge)
3. Any density matrix defines a **modular Hamiltonian**: K = -ln rho

What is a modular Hamiltonian? In ordinary quantum mechanics, the Hamiltonian H generates time evolution via $e^{-iHt}$. The modular Hamiltonian instead generates a state-dependent modular flow. If you know the density matrix $\rho$, you can take its logarithm and get an operator $K = -\ln\rho$ whose modular evolution acts like an internal clock for that state.

4. This Hamiltonian generates a flow: sigma_s(A) = e^{iKs} A e^{-iKs}
5. **The Thermal Time principle**: This flow is proposed as what we experience as time.

Time isn't a coordinate on a manifold. Time is the modular flow of the statistical state.

Here is the reframing: **Time flows because we are in a state of incomplete knowledge.** The thermal-time program ties experienced time to state restriction and coarse-graining, rather than treating time as a primitive external parameter.

### Tomita-Takesaki Theory

The mathematical foundation is **Tomita-Takesaki theory**.

Let M be a von Neumann algebra and |Omega> a cyclic and separating vector. Tomita-Takesaki theory constructs, from M and |Omega> alone, a one-parameter group of automorphisms:

$$\sigma_t(A) = \Delta^{it} A \Delta^{-it}$$

Even without specifying a Hamiltonian, even without putting time in by hand, the algebra-state pair *generates its own time evolution*.

Key properties:
1. **KMS Condition**: The state satisfies thermal equilibrium at "temperature" beta = 1 with respect to modular time
2. **State dependence**: The modular flow is fixed by the chosen algebra-state pair; different faithful states generally give different flows

This theorem says: given any quantum system and any state of incomplete knowledge, there's a natural notion of time evolution.

### The Rindler Wedge

This abstract mathematics connects to reality through the Unruh effect.

An observer accelerating uniformly sees only the **Rindler wedge**-part of spacetime. For the vacuum state restricted to this region, the Bisognano-Wichmann theorem shows that the modular Hamiltonian is exactly the generator of Lorentz boosts.

For an accelerating observer, a Lorentz boost *is* time translation. The modular flow equals ordinary time evolution.

The modular temperature works out to:

$$T_{Unruh} = \frac{\hbar a}{2\pi k_B c}$$

The Unruh effect isn't a separate phenomenon-it's Tomita-Takesaki theory applied to spacetime. The "time" experienced by an observer is determined by their restricted access to the quantum state.

## 11.4 The Arrow of Time

In Chapter 4, we saw Boltzmann's insight: entropy $S = k \ln W$ measures the number of microstates compatible with a macrostate, and entropy increases because high-entropy states vastly outnumber low-entropy ones.

But why did the universe start with low entropy in the first place?

### The Past principle

The deeper answer to the arrow of time is the **Past principle**: the universe began in a state of extraordinarily low entropy.

We're not riding a random fluctuation. We're riding the expansion from a very special initial condition-the Big Bang. The early universe was hot but smooth, with matter spread almost uniformly. That uniformity is low gravitational entropy.

Why was the Big Bang low entropy? Standard physics treats this as an unexplained initial condition. But our model offers a different perspective.

**The Past principle as a consistency requirement**: For observers to exist at all, they must be able to form and compare records. Records require entropy gradients-you can only write information by pushing entropy elsewhere. A universe in thermal equilibrium contains no observers, no records, no consistency-checking.

The MaxEnt principle says: assign the maximum-entropy state consistent with your constraints. But one constraint is that someone must exist to apply MaxEnt. This rules out equilibrium. The very existence of observers selecting MaxEnt states presupposes a universe far from equilibrium.

This doesn't derive the specific numerical entropy of the Big Bang. But it reframes the question: the Past principle isn't an arbitrary input to be explained by some deeper theory. It's a consistency requirement. A universe containing observers who check for consistency must have started with low entropy. The arrow of time points in the direction that allows records to be made.

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

In OPH, each observer has a patch P on the holographic screen. The global state restricts to a density matrix:

$$\rho_P = \text{Tr}_{\bar{P}} |\Psi\rangle \langle \Psi|$$

This density matrix defines a modular Hamiltonian:

$$K_P = -\ln \rho_P$$

which generates modular time \(t_P\) for that observer.

**Every observer has their own emergent clock.**

### Consistency of Clocks

If two observers' patches overlap, their modular times must be compatible on the overlap. This is a strong constraint. Reality hangs together because the modular flows mesh.

### Cosmic Time

Why do we all agree on a "cosmic time"?

If the global state is highly entangled in a particular pattern, the modular flows of local patches are synchronized. Cosmic time emerges as the "center of mass" of all local modular times.

### Roadmap: From Modular Time to Gravity

Here is the high-level chain we use later:

1. **Recovery structure** from earlier chapters makes the time-generator local near patch boundaries.
2. **A key theorem** identifies that local time-flow with a standard geometric transformation on the sphere, fixing its normalization.
3. Geometric time-flow gives **Lorentz kinematics** on the screen.
4. **Entanglement equilibrium** plus a way to identify local energy yields Einstein's equation as an output.

This chapter builds the time ingredient. The next sections show how it feeds into gravity.

## 11.7 Jacobson's Derivation

In 1995, Ted Jacobson performed one of the most beautiful derivations in theoretical physics.

He started with thermodynamics-the first law:

$$\delta Q = T \, dS$$

Then made three assumptions:
1. **Entropy is area**: S proportional to boundary area
2. **Heat is energy flux**: delta Q is stress-energy integrated over a local horizon
3. **Temperature is Unruh temperature**: T proportional to surface gravity

He demanded the relation hold for all local horizons.

Out popped **Einstein's field equations**:

$$R_{\mu\nu} - \frac{1}{2}R g_{\mu\nu} = 8\pi G T_{\mu\nu}$$

Jacobson inverted the logic of physics. Usually we think of gravity as fundamental, implying thermodynamic properties for horizons. Jacobson showed the reverse: **if you assume thermodynamics is fundamental, gravity is derived.**

**Gravity is not a fundamental force. Gravity is what happens when the vacuum tries to maintain local thermodynamic equilibrium.**

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

### What This Means

Start with thermal structure. Ask: what is the natural notion of time evolution? The answer is Lorentz boosts.

This reverses the usual logic in QFT. We don't postulate Lorentz symmetry and then discover thermal horizons; the BW theorem shows the boost structure is already encoded in modular flow.

In the OPH program, the modular/boost link is part of the route by which Lorentz kinematics and a universal light speed are recovered. That stronger claim is program-level and should not be presented as a standalone theorem of Unruh/Bisognano-Wichmann alone.

### Connection to Our Framework

In our model:
1. Each observer's patch has a boundary
2. This boundary is a horizon with Gibbons-Hawking temperature
3. The modular flow of the horizon state generates time evolution

In [*Observers Are All You Need*](../paper/observers_are_all_you_need.pdf), this idea is carried over from wedges in ordinary spacetime to caps on the holographic screen. The natural time flow of a patch becomes an actual geometric motion on the sphere, and once that happens the conformal symmetry of the sphere reproduces Lorentz symmetry.

### The Speed of Light

Why is there a maximum speed, and why is it the same for everyone?

The Unruh formula T = ℏa/(2πck_B) contains c. For the thermal-to-boost correspondence to work, there must be a universal velocity relating acceleration to temperature.

From the boundary perspective: information propagates on the S² screen at a maximum rate determined by the entanglement structure. This rate, translated to the bulk, becomes c. The no-signaling theorem of quantum mechanics (entanglement can't transmit information) becomes, in the bulk, the statement that nothing travels faster than light.

### The Causal Structure

The light cone structure of spacetime-which events can influence which-emerges from entanglement:

- **Spacelike separation**: Regions can be correlated (entangled) but cannot signal
- **Timelike separation**: Events can have causal influence
- **Null separation**: The boundary between these regimes

The modular flow provides the time direction. Entanglement provides correlations. No-signaling prevents faster-than-light communication. Taken together, these ingredients reproduce the Minkowski-style causal structure targeted by the program.

### Why This Matters

Einstein discovered special relativity in 1905 by thinking about light and motion. Over a century later, we see it differently: in QFT, Lorentz boosts are tied to horizon thermodynamics via the Bisognano-Wichmann theorem. In our model the same pattern appears when the screen reaches its smooth geometric limit, so the Lorentz group shows up as the geometry of modular flow on caps.

The laws of physics look the same to all inertial observers because thermal states on wedge-shaped regions naturally evolve via boosts. In the OPH program, the universal speed emerges when that modular-boost structure is carried over to the screen and then read back into bulk kinematics.

## 11.10 Testable Predictions and Verified Results

The emergent time model includes both rigorous mathematical results and testable predictions:

**Rigorous results (mathematical theorems)**:

**1. Tomita-Takesaki theorem**: Once you specify both the observables available to an observer and the state they are in, there is a built-in notion of time flow. The algebra-state pair carries its own natural dynamics.

**2. KMS condition**: That natural time flow behaves exactly like thermal equilibrium. In other words, modular time is not an arbitrary parameter; it has the precise structure of heat physics.

**3. Bisognano-Wichmann theorem**: For a wedge-shaped region of spacetime, the natural modular time flow is literally a Lorentz boost. That is the deep reason thermal structure and relativity fit together so tightly. In OPH, the claim is that spherical patches on the screen behave in the same spirit.

**4. Boltzmann's H-theorem**: Under standard coarse-graining assumptions, entropy almost always rises. This is the classic bridge from reversible microscopic laws to an irreversible macroscopic arrow of time.

**Testable predictions**:

**1. Unruh effect**: Accelerating observers see thermal radiation at T = ℏa/(2πk_B c). While direct detection is beyond current technology (requires acceleration ~10²⁰ m/s²), the Unruh effect is closely related to Hawking radiation through the same horizon-thermodynamic structure, and the mathematics is well established.

**2. Jacobson's derivation**: If entropy ∝ area and temperature ∝ surface gravity, then Einstein's equations can be recovered under specific thermodynamic assumptions. This is strong evidence for the thermodynamic route, but not a proof that every horizon-thermodynamics construction must reduce uniquely to GR.

**3. Microscopic laws are largely time-symmetric**: Electromagnetic, strong, and gravitational dynamics are time-reversal symmetric in their standard formulations, while the weak interaction contains the known small CP/T-violating effects.

**4. Arrow of time from Past principle**: Given low-entropy initial conditions, the Second Law follows statistically. Confirmed by the entire edifice of thermodynamics and cosmology.

**Empirical validation signatures**:
- Microscopic laws with fundamental time asymmetry (beyond tiny CP violation)
- Modular flow failing to generate consistent time evolution
- Unruh temperature having wrong dependence on acceleration
- Jacobson's derivation failing for some horizon type

None of these contradicting observations has ever been made.

---

## 11.11 Memory and Records

Why do we remember the past but not the future?

A **memory** is a physical record-a low-entropy structure correlated with a past event. Creating a record requires work-you must push entropy somewhere else.

When you remember something, you're consulting a present record created at the cost of increasing entropy elsewhere. The record only makes sense if entropy was lower when the recorded event happened.

The arrow of time is the arrow of record-keeping. Time flows in the direction we can make and preserve consistent records.

## 11.12 Reverse Engineering Summary

Recap:

| Intuitive Picture | Surprising Hint | First-Principles Reframing |
|---|---|---|
| Time is a fundamental external parameter flowing from past to future | No preferred time in GR; the Wheeler-DeWitt equation \(H\Psi = 0\) points toward a fundamentally timeless description; microscopic laws are largely time-symmetric | Time emerges from modular flow of thermal states; the arrow of time is the direction of increasing entropy, driven by incomplete knowledge |

**The key reverse engineering insight**: We started with the intuition that time flows fundamentally, carrying the universe from past to future. General relativity showed that there is no preferred time; different observers slice spacetime differently. Quantum gravity shocked us further with the Wheeler-DeWitt equation, showing the universe is fundamentally frozen. Our model explains how time emerges: through the modular flow of density matrices. Time is what incomplete knowledge looks like. The arrow of time points in the direction of consistency-building, the direction where records can be made and compared.

**Additional lessons**:

1. **Boltzmann**: Entropy measures the number of microstates compatible with a macrostate. Entropy increases because high-entropy states vastly outnumber low-entropy states.

2. **Past principle**: The arrow of time exists because the Big Bang was a low-entropy state. Our model shows this isn't an arbitrary input; it's a consistency requirement. Observers need entropy gradients to form records, so a universe with observers must start far from equilibrium.

3. **Jaynes**: Entropy measures ignorance. In the Jaynes program, MaxEnt gives the least-biased probability distribution consistent with the stated constraints.

4. **Thermal Time principle**: Time is the modular flow of our statistical state.

5. **Tomita-Takesaki**: Any algebra-state pair generates its own time evolution.

6. **Jacobson**: Einstein's equations follow from thermodynamics. Gravity is an equation of state.

7. **Complexity**: The growth of wormhole interiors tracks computational complexity. Time flows because the universe is still computing.

8. **Records**: We remember the past because records require entropy flow from a low-entropy origin.

9. **Bisognano-Wichmann**: In QFT wedges, Lorentz boosts are modular flow. Our screen analog follows once the recovery and MaxEnt structure is in place (Chapter 15).

---

We've found the "engine" of reality: time emerges from incomplete knowledge, flowing in the direction of consistency-building.

Now we ask: why does the machine have these particular parts? Why these particles, these forces, these symmetries?

The answer lies in the geometry of the screen. That's the story of **Chapter 12: Symmetry on the Sphere**.
