# Chapter 3: The Screen and the Sphere

## 3.1 The Volume Hint

Here's what seems obvious about information: more space should hold more data.

A bigger hard drive stores more files. A bigger warehouse holds more boxes. A bigger brain should hold more memories. The amount of stuff you can fit into a container should scale with its volume.

This is the **intuitive picture**: information content scales with volume.

$$\text{Information} \propto V$$

If you have a box and you divide it in half, each half should hold half the information. If you double the size of a room, you should be able to fit twice as many things in it.

This seems so obvious that nobody questioned it for most of physics history.

And it's wrong.

The universe gave us a hint-a spectacular, unexpected hint-that information doesn't work this way at all. The hint came from the strangest objects in the cosmos: black holes.

## 3.2 The Teacup Problem: The Hint

In 1972, a graduate student named Jacob Bekenstein walked into John Wheeler's office at Princeton with a simple thought experiment.

Imagine a cup of hot tea. The tea has entropy-it is hot and messy, with many microscopic arrangements of molecules that produce the same macroscopic state.

Now lower the cup into a black hole.

The tea crosses the event horizon and vanishes. No one outside can ever see it again. If the tea is gone, so is its entropy. The total entropy of the observable universe has decreased.

But wait. The Second Law of Thermodynamics says toital entropy never decreases. The Second Law is the rule that makes time flow in a direction. It tells you why broken glasses don't unbreak, why scrambled eggs don't unscramble, why we remember the past but not the future.

If a black hole can erase entropy, the Second Law is wrong.

### Bekenstein's Bold Response

Bekenstein proposed that black holes must have entropy. When the tea falls in, the entropy doesn't disappear-it shows up as an increase in the black hole's own entropy.

But where could a black hole's entropy hide?

Black holes are supposed to be simple. In general relativity, a black hole is fully described by just three numbers: its mass, its electric charge, and its spin. Wheeler called this the "no-hair theorem"-black holes have no distinguishing features.

So where are the microstates? Where is the internal structure that entropy requires?

Bekenstein looked at the only thing that changes when you throw stuff in: the size of the event horizon. He made a guess-an educated guess, constrained by dimensional analysis and theoretical consistency-that the entropy is proportional to the **area** of the horizon:

$$S \propto A$$

Not the volume. The area.

### Hawking Confirms It

Stephen Hawking was skeptical. He set out to prove Bekenstein wrong by showing black holes have no temperature.

He studied quantum fields near a black hole horizon. What he found shocked him.

The vacuum of quantum field theory seethes with virtual particle pairs that pop into existence and annihilate. Near a horizon, one particle can fall in while the other escapes. To a distant observer, the black hole emits radiation-**Hawking radiation**.

Hawking calculated the temperature:

$$T_H = \frac{\hbar c^3}{8\pi G M k_B}$$

Once a black hole has temperature, it must have entropy. From thermodynamics, Hawking derived:

$$S_{BH} = \frac{A}{4 \ell_P^2}$$

where $\ell_P = \sqrt{\hbar G/c^3} \approx 1.6 \times 10^{-35}$ m is the Planck length.

The entropy of a black hole is proportional to its surface area, measured in Planck units.

### The Surprising Conclusion

**The hint**: Information scales with area, not volume.

**The lesson**: The intuitive picture-that information content scales with the size of a container-is fundamentally wrong. The boundary of a region, not its interior, is where the information lives.

**The first-principles reframing**: The 3D world we experience is not the fundamental level. The fundamental degrees of freedom live on 2D surfaces. The bulk is emergent, reconstructed from boundary data.

## 3.3 Entropy for Normal People

Before we go further, let's make sure we understand entropy. It sounds abstract but has a very concrete meaning.

### Counting Possibilities

Entropy is a count of possibilities. More precisely, it's the logarithm of the number of microscopic states (microstates) consistent with the same macroscopic description (macrostate).

Picture a row of ten coins. Each coin can be heads (H) or tails (T). There are $2^{10} = 1024$ possible arrangements.

Now suppose I tell you: "Five coins are heads and five are tails." How many arrangements fit this description? The answer is "10 choose 5" = 252.

But if I tell you "all ten are heads," there's only one arrangement: HHHHHHHHHH.

The five-and-five case has more entropy because there are more microstates matching the same macrostate. You have less information about the exact arrangement.

Boltzmann turned this into a formula:

$$S = k_B \ln W$$

where $W$ is the number of microstates and $k_B$ is Boltzmann's constant.

**Entropy measures ignorance.** A system with more possible microstates has higher entropy because there is more you don't know about its exact state.

### The Second Law

The Second Law-entropy increases-becomes almost obvious once you think this way.

If you have a gas in the left half of a box and you open a partition, the gas will spread to fill the whole box. Not because there's a force pushing it, but because there are vastly more arrangements with gas spread out than with gas confined.

The system wanders randomly through its possible states and inevitably ends up in high-entropy configurations because there are more of them.

### Information and Entropy: Two Names for One Thing

In 1948, Claude Shannon created information theory. He defined the **bit**-the fundamental unit of information-and derived a measure of uncertainty:

$$H = -\sum_i p_i \log_2 p_i$$

This is the **Shannon entropy**, measured in bits.

Shannon's entropy has the exact same mathematical form as Boltzmann's entropy. This is not a coincidence. Both measure missing information.

And in 1961, Rolf Landauer showed that erasing information costs energy-at least $k_B T \ln 2$ per bit. Information is physical. Bits are thermodynamic objects.

This is why the Bekenstein-Hawking formula is so important. It connects information (entropy) to geometry (area). It tells us where the data lives: on the boundary.

## 3.4 The Holographic Principle

Now let's follow the hint to its logical conclusion.

### The Bekenstein Bound

Bekenstein asked: what is the maximum entropy you can pack into a region of size $R$ containing energy $E$?

If you try to pack more and more entropy into a fixed region, you need more energy. But if you pack enough energy, gravity becomes strong. Eventually, the region collapses into a black hole.

A black hole of radius $R$ has maximum entropy $S = \pi R^2 / \ell_P^2$. You cannot pack more entropy into that region without forming a larger black hole.

This leads to the **Bekenstein bound**:

$$S \leq \frac{2\pi R E}{\hbar c}$$

The maximum entropy of a region is bounded by its surface area, not its volume.

### What This Means

This is weird. It means the interior of a region is somehow redundant. All the independent degrees of freedom can be counted by looking at the boundary.

In the early 1990s, Gerard 't Hooft and Leonard Susskind took this reasoning to its logical conclusion.

If the maximum information in a region scales with surface area, then **the fundamental degrees of freedom must be two-dimensional**. The three-dimensional interior is not fundamental-it is emergent, reconstructed from boundary data.

't Hooft called this the **holographic principle**, by analogy with holograms. A hologram is a two-dimensional film that encodes a three-dimensional image. When you illuminate it, you see depth that isn't really there-the depth is computed from the 2D pattern.

The holographic principle says the universe works the same way. The fundamental data lives on a 2D surface. The 3D world we experience is the computed image.

### The Logic

1. Black holes have entropy proportional to area, not volume.
2. No region can contain more information than a black hole of the same size.
3. Therefore, information content of any region is bounded by surface area.
4. Therefore, the bulk degrees of freedom are not independent-they can be reconstructed from boundary data.

If you accept Bekenstein-Hawking entropy, holography follows.

## 3.5 Black Holes and Horizons

Let's make sure we understand what a horizon is-and why every observer has one.

### The Event Horizon

A black hole is not a physical object in the usual sense-it's a region of spacetime. The **event horizon** is the boundary of that region. Once you cross it, you cannot escape.

The Schwarzschild radius of a black hole of mass $M$ is:

$$R_s = \frac{2GM}{c^2}$$

For the Sun, this is about 3 kilometers. For Earth, it's about 9 millimeters. Any mass compressed within its Schwarzschild radius becomes a black hole.

The horizon is not a physical surface. You could cross it without noticing anything special. But once you're inside, the geometry of spacetime is such that all paths-even light paths-lead inward.

Here's a way to think about it: near a black hole, space is falling inward like a waterfall. The event horizon is where the water falls faster than you can swim.

### Other Horizons

Black holes are not the only source of horizons.

**Cosmological horizons**: The universe is expanding. Beyond a certain distance-about 46 billion light-years-galaxies are receding faster than light. Light from those regions will never reach us.

**Acceleration horizons**: If you accelerate continuously, there is a region behind you from which light can never catch up. You have a **Rindler horizon**. This produces the **Unruh effect**: an accelerating observer perceives the vacuum as a warm bath of particles.

In each case, the horizon is a boundary that limits what the observer can access. It is the edge of their observable universe.

### Every Observer Has a Screen

Here's the key insight: every observer has a horizon, and therefore every observer has a screen.

For an observer in our universe:
- There's a cosmic horizon at roughly 46 billion light-years
- If they're near a black hole, there's an event horizon
- If they're accelerating, there's a Rindler horizon

The horizon is approximately spherical. The area of this sphere bounds the amount of information the observer can access.

This is a deep shift in perspective. Instead of thinking about space as a fixed container, we think about each observer's horizon as their fundamental interface with reality.

## 3.6 Why a Sphere?

The screen is always (approximately) spherical. This is not an arbitrary choice-it follows from causality.

Light travels at the same speed in all directions. If you stand at a point and wait, the light that can reach you from a time $t$ ago forms a sphere of radius $ct$ around you.

Your past light cone-the set of events that could have influenced you-has spherical cross-sections. Your future light cone also has spherical cross-sections.

The sphere is a consequence of the geometry of causality.

### The Cosmic Microwave Background

The cosmic microwave background (CMB) illustrates this beautifully.

The CMB is light from about 380,000 years after the Big Bang, when the universe cooled enough for atoms to form and light to travel freely. This light appears as a sphere around us-the **last scattering surface**.

We're at the center of this sphere, but so is everyone else. Every observer in the universe sees themselves at the center of their own CMB sphere.

The CMB sphere is our screen. Every cosmological observation is, ultimately, a measurement of patterns on this 2D surface.

## 3.7 The Geometry of the 2-Sphere

The mathematical object describing the screen is the 2-sphere, $S^2$.

$$S^2 = \{(x, y, z) \in \mathbb{R}^3 : x^2 + y^2 + z^2 = 1\}$$

We can parameterize it with spherical coordinates $(\theta, \phi)$:
- $\theta$ is the polar angle, from 0 at the North Pole to $\pi$ at the South Pole
- $\phi$ is the azimuthal angle, from 0 to $2\pi$ around the equator

The metric is:

$$ds^2 = d\theta^2 + \sin^2\theta \, d\phi^2$$

### Spherical Harmonics

Any function on the sphere can be expanded in **spherical harmonics**, $Y_\ell^m(\theta, \phi)$. These are the natural modes of vibration of the sphere.

The CMB temperature variations are analyzed by expanding in spherical harmonics. The **power spectrum**-how much power at each angular scale $\ell$-tells us about the early universe.

### Finite Resolution

If the screen has a smallest length scale-a pixel size at the Planck length-then there is a maximum $\ell$:

$$\ell_{max} \sim \frac{R}{\ell_P}$$

The total number of independent modes is roughly $\ell_{max}^2 \sim R^2/\ell_P^2$-proportional to area in Planck units, exactly what Bekenstein-Hawking says.

The finite resolution of the screen means our experience of a continuous world is an approximation. At the smallest scales, space is pixelated.

## 3.8 Patches and Overlaps

You cannot see the whole screen. Some parts are hidden by your horizon or by instrumental limits. You only access a **patch**-a portion of the sphere.

Another observer, at a different location or with different instruments, accesses a different patch. Where patches overlap, observers can compare notes.

If the screen is a sphere $S^2$ and observer $i$ sees patch $P_i$, then two observers can compare data on the overlap $P_i \cap P_j$. That overlap is the seed of consistency.

### A Concrete Example

Consider two astronomers on opposite sides of Earth. During the night, they see different parts of the sky. But some stars are visible to both-stars near the horizon for each observer.

These shared stars provide a link. The astronomers can calibrate by comparing their observations of the overlap region. Once they agree on the overlap, they can combine their observations into a consistent map of the whole sky.

### Coordinate Charts and Atlases

A sphere cannot be covered by a single smooth coordinate system. If you try to put latitude-longitude coordinates on a sphere, you run into problems at the poles.

Mathematicians handle this by using multiple overlapping coordinate charts, called an **atlas**. Each chart covers part of the sphere. Where charts overlap, there are transition functions that tell you how to convert coordinates.

This is exactly analogous to our observer patches. Each observer has a local description. Where observers overlap, they must agree on how to translate between their descriptions.

Physics is the art of finding descriptions that work in many charts and have consistent translations between them.

## 3.9 What Is an Observer?

We've talked about "observers" and their "patches." But what exactly IS an observer in this model?

### Not External Watchers

In classical physics, observers are implicitly outside the system-disembodied measurers who don't affect what they measure. This won't work here. Observers must be part of the system they observe.

### Observers as Patterns in the Data

An observer is a special kind of pattern in the horizon data-a subsystem with three key properties:

**1. Bounded access**: The observer can only interact with a finite patch P of the screen. This patch defines what the observer can measure, know, and act upon. The boundary of the patch is the observer's horizon.

**2. Stable records**: The observer contains internal correlations that persist over time-memory. When you measure something and remember the result, your brain has become correlated with the measured system. These correlations are the "records" that define measurement outcomes.

**3. Self-modeling**: An observer can build compressed representations of its environment. Your brain doesn't store raw sensory data; it builds a model of the world.

### The Vortex Analogy

Think of observers as stable vortices in a fluid.

The fluid is the quantum state on the horizon-constantly evolving, highly correlated. A vortex isn't separate from the fluid; it's a pattern within the fluid. It persists over time. It has a definite location. It interacts with other patterns.

An observer is like that. It's not a ghostly presence watching from outside. It's a stable, self-reinforcing pattern within the data on the screen. The pattern has access to a local region (its patch), maintains internal structure (its records), and can interact with nearby patterns (other observers, measured systems).

### Movement and Time

Do observers "move around" on the sphere?

Not in a simple sense. Different patches represent different observers, or the same observer at different moments. "Movement" is actually a sequence of overlapping patches with consistent marginals.

What creates the sense of time? The internal structure of the quantum state provides a natural flow: the **modular flow** from quantum statistical mechanics. For a thermal state, modular flow generates time evolution. The thermal time principle (Connes and Rovelli) shows this is the origin of experienced time.

### Why This Matters

This definition of observers resolves several puzzles:

**No external reference frame**: Observers are internal to the system, so there's no need for an external "God's-eye view."

**Measurement is physical**: When an observer measures something, correlations form between subsystems within the horizon data. There's nothing mysterious about "collapse"-it's just the establishment of records.

**Consistency follows from structure**: Two observers whose patches overlap must agree because they're both patterns in the same underlying state. The state is self-consistent, so the observers are consistent.

### Reality from Computation

Here's a concrete way to think about the screen and its observers.

Imagine the screen as a **gauge-invariant quantum system** on the 2-sphere, something like a quantum cellular automaton but with important structure. Triangulate the sphere into tiny cells. At each edge of the triangulation sits a finite-dimensional quantum system (a qudit). At each vertex, a gauge constraint (Gauss's law) restricts which configurations are physical. Not all states survive; only those satisfying the constraint at every vertex.

**Observer patches** are subsystems defined by boundary-gauge-invariant algebras. Each patch is like a computational thread, a connected region where an observer can ask questions and get answers. The algebra $\mathcal{A}(R)$ defines what that observer can measure: the operators that commute with the boundary gauge transformations.

**Overlap consistency** is automatic. Where two patches intersect, they access the same gauge-invariant observables. Both observers are reading the same underlying data, just from different angles. The gauge redundancy at boundaries is what makes gluing non-trivial and gives rise to the "edge modes" that carry geometric information.

**The dynamics** comes from MaxEnt: among all states consistent with the constraints, nature selects the maximum entropy state. This is like a Gibbs state $\rho \propto e^{-H}$ where $H$ is a sum of local terms. The system thermalizes at the UV scale, and the macroscopic physics emerges from this equilibrium.

**The 4D bulk isn't on the sphere.** It emerges from the entanglement structure between patches. When you look around and see three-dimensional space, you're experiencing a compressed encoding of how your patch is entangled with others. Distance in the bulk is entanglement on the boundary.

*The screen is the computation. Observer patches are the threads. Reality is what they agree on.*

This computational picture is made rigorous through **quantum link models**, a class of lattice gauge theories with finite-dimensional Hilbert spaces that realize exactly the structure we need. The technical details are in the paper; the intuition is that the sphere is running a quantum computation, and we are processes within it.

## 3.10 Entanglement Creates Depth

We have said the 3D world emerges from 2D boundary data. But how? What creates the feeling of depth?

The answer is **entanglement**.

### Quantum Entanglement

Quantum entanglement is correlation with no classical analogue. The simplest example is the Bell state:

$$|\Psi\rangle = \frac{1}{\sqrt{2}}(|00\rangle + |11\rangle)$$

Two qubits, each of which can be 0 or 1. If you measure the first qubit, you get 0 or 1 with equal probability. But once you measure the first and get (say) 0, the second is guaranteed to also be 0.

The key point: each qubit alone looks random. But together, the pair is perfectly correlated. The whole knows more than the parts.

### Entanglement Entropy

If you have a large quantum system and divide it into two parts, A and B, you can ask: how entangled are they?

The measure is **entanglement entropy**:

$$S_A = -\text{Tr}(\rho_A \log \rho_A)$$

where $\rho_A$ is what you get when you trace out system B.

### The Area Law

In many quantum systems-particularly ground states of local Hamiltonians-entanglement entropy follows an **area law**:

$$S_A \propto \text{Area}(\partial A)$$

The entanglement between region A and its complement is proportional to the boundary between them, not the volume.

This is deeply connected to holography. If entanglement scales with boundary area, then bulk degrees of freedom are not all independent.

### Ryu-Takayanagi: Distance from Entanglement

In 2006, Shinsei Ryu and Tadashi Takayanagi made the connection precise. In the AdS/CFT correspondence:

$$S_A = \frac{\text{Area}(\gamma_A)}{4 G_N}$$

The entanglement entropy of a boundary region A equals the area of the minimal surface $\gamma_A$ in the bulk.

This is stunning. **Entanglement on the boundary creates geometry in the bulk.** If two regions of the boundary are highly entangled, the minimal surface connecting them is short-they are close in the bulk. If weakly entangled, they are far apart.

**Distance in the emergent space is entanglement on the screen.**

### Tensor Networks

This can be visualized with **tensor networks**. Imagine the screen as a 2D grid of nodes, each representing a qubit. Nodes are connected by bonds representing entanglement.

Stack layers of such grids, with each layer connected to the one below. The pattern of bonds forms a network. You can define a "distance" based on how many bonds you must cross.

For certain tensor networks (like MERA), this distance matches the geometry of curved spacetime. The network encodes an emergent spatial dimension.

The bulk is decoded from the boundary. Depth-the sense of being "inside" 3D space-is how we render entanglement patterns into a coherent world.

### Testable Consequences

The holographic principle makes sharp, testable predictions:

**1. Area law vs. volume law**: If the holographic principle is correct, entanglement entropy in gravitational systems must scale with area, not volume. Tensor network models that respect the holographic bound produce area-law scaling. Generic quantum field theory states produce volume-law scaling. This is a discriminating test.

**2. The Bekenstein bound is saturated by black holes**: No system can have entropy exceeding $S = A/(4\ell_P^2)$. Black holes saturate this bound-they are maximally entropic for their size. Any violation identifies a measurement contradiction with the model.

**3. Information is finite**: The observable universe contains at most $\sim 10^{122}$ bits. This is enormous but finite. Any evidence of truly infinite information content would contradict holography.

These predictions have been tested in every context where we can check:
- Black hole thermodynamics confirms area-entropy
- AdS/CFT calculations match both sides of the duality
- No violation of Bekenstein bounds has ever been observed

The holographic principle started as a proposal. It is now one of the most tested ideas in theoretical physics.

## 3.11 The Reverse Engineering

Let's trace the reverse engineering explicitly.

**The intuitive picture**: Information scales with volume. Space is the fundamental container.

**The hint**: Black hole entropy scales with area. The Bekenstein bound limits information by surface area.

**The lesson**: The fundamental degrees of freedom are 2D, not 3D. The boundary is primary; the bulk is emergent.

**The first-principles reframing**:

1. Each observer has a horizon-a spherical screen bounding their accessible information
2. The screen carries the fundamental data, limited by $S \leq A/(4\ell_P^2)$
3. Entanglement patterns on the screen create the geometry of the emergent 3D bulk
4. Different observers have different screens, but consistency on overlaps makes the emergent 3D world shared and stable

The holographic principle is not a philosophical preference. It's what the hints force us to conclude.

## 3.12 Pixel Limits

Let's put numbers on this.

The Planck length is $\ell_P \approx 1.6 \times 10^{-35}$ meters-about $10^{20}$ times smaller than a proton. The Planck area is $\ell_P^2 \approx 2.6 \times 10^{-70}$ m².

**The observable universe**: Radius $R \approx 4.4 \times 10^{26}$ m. Horizon area $A \approx 2.4 \times 10^{54}$ m². Number of bits: $N \approx 10^{122}$.

This is a truly enormous number-but it is finite. The observable universe contains a finite amount of information.

**A solar-mass black hole**: Schwarzschild radius $R_s \approx 3$ km. Number of bits: $N \approx 10^{77}$.

This is still huge, but much smaller than the observable universe. Yet it's far more than the entropy of the Sun as a normal star (about $10^{58}$). Collapse increases entropy because the horizon has vastly more microstates than ordinary matter.

The finite resolution of the screen means continuous space is an approximation. At the smallest scales, reality is digital.

## 3.13 Where We Go Next

We have established that:
- Information lives on horizons, not in volumes
- Horizons are spherical, a consequence of causality
- The amount of information is finite, bounded by area
- Entanglement patterns on the screen create emergent 3D geometry

But we haven't yet explained dynamics. The screen we've described is static-it encodes information. What makes things happen? What creates the arrow of time?

The answer involves entropy again, but now entropy's role in dynamics. The Second Law says entropy increases. But why? And what does this have to do with the screen?

In the next chapter, we explore the edge of the screen-the boundary conditions that govern what can happen. We will see how entropy growth is not just a statistical tendency but a geometric constraint, built into the structure of horizons themselves.

The reverse engineering continues in Chapter 4: Entropy on the Edge.
