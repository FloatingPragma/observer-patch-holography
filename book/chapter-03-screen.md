# Chapter 3: The Screen and the Sphere

## 3.1 The Volume Hint

The ordinary intuition says more space should hold more data.

A bigger hard drive stores more files, a bigger warehouse holds more boxes, and a bigger brain, presumably, holds more memories. The amount of stuff you can fit into a container should scale with its volume.

This is the **intuitive picture**: information content scales with volume.

$$\text{Information} \propto V$$

If you have a box and you divide it in half, each half should hold half the information. If you double the size of a room, you should be able to fit twice as many things in it.

The symbol $V$ means volume. The proportionality sign says that, in this
ordinary intuition, information capacity should grow in direct proportion to
the amount of three-dimensional space available.

This seems so obvious that nobody questioned it for most of physics history.

And it's wrong.

The universe gave us a spectacular, unexpected hint that information does not
work this way at all. The hint came from the strangest objects in the cosmos:
black holes.

## 3.2 The Teacup Problem: The Hint

In 1972, John Wheeler put a simple thought experiment to his Princeton graduate student Jacob Bekenstein.

Imagine a cup of hot tea. The tea has entropy: it is hot and messy, with many
microscopic arrangements of molecules that produce the same macroscopic state.

Lower the cup into a black hole. This is, among other things, a terrible
thing to do to a cup of tea.

The tea crosses the event horizon and vanishes. No one outside can ever see it again. If the tea is gone, so is its entropy. The total entropy of the observable universe has decreased.

But wait. The Second Law of Thermodynamics says total entropy never decreases. The Second Law is the rule that gives time its direction.

If a black hole can erase entropy, the Second Law is wrong.

### Bekenstein's Bold Response

Bekenstein proposed that black holes must have entropy. When the tea falls in,
the entropy does not disappear. It shows up as an increase in the black hole's
own entropy.

But where could a black hole's entropy hide?

Black holes are supposed to be simple. In general relativity, a black hole is
fully described by just three numbers: its mass, its electric charge, and its
spin. Wheeler called this the "no-hair theorem": black holes have no
distinguishing features.

So where are the microstates? Where is the internal structure that entropy requires?

Bekenstein looked at the only thing that changes when you throw stuff in: the
size of the event horizon. He made an educated guess, constrained by
dimensional analysis and theoretical consistency: the entropy is proportional
to the **area** of the horizon:

$$S \propto A$$

The volume does not appear. The area does.

### Hawking Confirms It

Stephen Hawking was skeptical. He set out to prove Bekenstein wrong by showing black holes have no temperature.

He studied quantum fields near a black hole horizon. What he found shocked him.

The vacuum of quantum field theory seethes with virtual particle pairs that pop
into existence and annihilate. Near a horizon, one particle can fall in while
the other escapes. To a distant observer, the black hole emits radiation,
**Hawking radiation**.

Hawking calculated the temperature:

$$T_H = \frac{\hbar c^3}{8\pi G M k_B}$$

In Hawking's temperature formula, $T_H$ is the black-hole temperature, $M$ is
the black-hole mass, $G$ is Newton's gravitational constant, $c$ is the speed
of light, $\hbar$ is Planck's constant divided by $2\pi$, and $k_B$ is
Boltzmann's constant. Larger black holes are colder because $M$ sits in the
denominator.

Once a black hole has temperature, it must have entropy. From thermodynamics, Hawking derived:

$$S_{BH} = \frac{A}{4 \ell_P^2}$$

where $\ell_P = \sqrt{\hbar G/c^3} \approx 1.6 \times 10^{-35}$ m is the Planck length.

The entropy of a black hole is proportional to its surface area, measured in Planck units.

### The Surprising Conclusion

**The hint**: Information scales with area, not volume.

**The lesson**: The intuitive picture that information content scales with the
size of a container fails for gravitational capacity. Black-hole entropy and
related bounds push strongly toward a boundary-sensitive description.

**The first-principles reframing**: The 3D world we experience can be an emergent bulk reconstructed from boundary data.

## 3.3 Why Entropy Points to the Boundary

Entropy counts how many microscopic arrangements fit one macroscopic description. Chapter 4 develops that idea carefully. The screen chapter needs one narrower lesson.

For black holes, the entropy is set by horizon area:

$$S_{BH} = \frac{A}{4 \ell_P^2}.$$

The surprise is that area supplies the natural counting measure for the most
extreme gravitating objects. Any observer-centered account of accessible
information must therefore take boundaries seriously.

Bekenstein sharpened the point further. Pack enough energy into a region and the region becomes a black hole. The black hole then supplies the maximum entropy compatible with that size. Area becomes the natural ceiling for accessible information in gravitational settings.

## 3.4 From Area Scaling to Holography

The jump from thermodynamics to geometry begins when the largest possible
entropy in a region is controlled by its boundary. A boundary-first description
then becomes the natural bookkeeping choice. The bulk may be the world we
experience, while the independent data is organized more economically on the
boundary.

This is the holographic idea in its simplest form. A two-dimensional surface can encode a three-dimensional description, just as a hologram stores depth information on a film.

Chapter 8 returns to holography in full. For the present chapter, the conclusion is simpler. The horizon is the right place to organize the data available to an observer.

## 3.5 Black Holes and Horizons

A horizon is the boundary of what one observer can ever check.

### The Event Horizon

A black hole is a region of spacetime rather than a physical object in the usual sense. The **event horizon** is the boundary of that region. Once you cross it, you cannot escape.

The Schwarzschild radius of a black hole of mass $M$ is:

$$R_s = \frac{2GM}{c^2}$$

For the Sun, this is about 3 kilometers. For Earth, it's about 9 millimeters. Any mass compressed within its Schwarzschild radius becomes a black hole.

$R_s$ is the Schwarzschild radius. The formula says that the critical radius
grows linearly with mass $M$. Compress the same mass inside that radius and
the escape velocity at the boundary reaches light speed.

A horizon is a causal boundary. You could cross it without noticing anything special. Once you're inside, the geometry of spacetime is such that all paths, even light paths, lead inward.

Near a black hole, space is falling inward like a waterfall. The event horizon is where the water falls faster than you can swim.

### Other Horizons

Black holes are not the only source of horizons.

**Cosmological horizons**: The universe is expanding, and cosmology distinguishes the observable-universe scale from the future event horizon. Regions from which light cannot reach us make observer access finite.

**Acceleration horizons**: If you accelerate continuously, there is a region behind you from which light can never catch up. You have a **Rindler horizon**. This produces the **Unruh effect**: an accelerating observer perceives the vacuum as a warm bath of particles.

In each case, the horizon is a boundary that limits what the observer can access. It is the edge of their observable universe.

### Every Observer Has a Screen

Finite observer access naturally suggests an effective screen picture.

The word "screen" needs care. OPH does not mean that every observer owns a
separate physical sphere, and it does not require a literal ball with data
painted on its surface. In the idealized description there is a shared
observer-facing screen net, often drawn with the symmetric chart $S^2$. A finite
observer has access only to a bounded local cut of that net. That cut is the
patch's support chart. The observer patch itself is the operational package
carried there: its accessible algebra and state, interfaces, readable records,
mismatch comparisons, repair moves, and checkpoint data.

For an observer in our universe, the accessible boundary can take several
forms. There is an observer-dependent cosmological horizon scale. Near a black
hole there is an event horizon. Under sustained acceleration there is a
Rindler horizon.

In the simplest symmetric situations, the relevant causal boundary is approximately spherical. The area of this sphere bounds the amount of information the observer can access.

This is a deep shift in perspective. Space is not a fixed container. Each observer's horizon is a fundamental interface with reality.

## 3.6 Why a Sphere?

The order matters here. OPH is not allowed to assume a light cone and then
discover a sphere; that would smuggle in the answer. Instead the framework
asks the patches themselves to certify their own geometry: do the repaired
pieces knit into one closed surface with no holes and no twists in its
stitching? When they do, the surface is a sphere and the sphere is round.
Without that certificate, repair consensus alone can just as well settle on a
torus or another support. Whether the twelve-port machinery delivers the
certificate is work in progress.

Only after the sphere has been certified does its symmetry group hand back
Lorentz kinematics. Light cones come out of that construction; they do not go
into it.

### The Cosmic Microwave Background

The cosmic microwave background (CMB) illustrates this beautifully.

In the standard FLRW reconstruction, the CMB is light from the era when the
cosmic plasma cooled enough for atoms to form and photons to travel freely.
This light appears as a sphere around us: the **last scattering surface**.

We're at the center of this sphere, but so is everyone else. Every observer in the universe sees themselves at the center of their own CMB sphere.

The CMB sphere is a useful cosmological proxy for thinking about an observer-centered screen picture. It is one especially vivid example of how observer-accessible information can be organized on an apparent 2D sky.

## 3.7 The Geometry of the 2-Sphere

On the certified global-support branch, the mathematical object describing the
support screen is the 2-sphere, $S^2$.

$$S^2 = \{(x, y, z) \in \mathbb{R}^3 : x^2 + y^2 + z^2 = 1\}$$

We can parameterize it with spherical coordinates $(\theta,\phi)$. The angle
$\theta$ runs from the North Pole to the South Pole, and $\phi$ runs around the
equator from $0$ to $2\pi$.

The metric is:

$$ds^2 = d\theta^2 + \sin^2\theta \, d\phi^2$$

$S^2$ means the two-dimensional surface of a unit sphere, not the solid ball
inside it. The set notation says: take all points $(x,y,z)$ in ordinary
three-dimensional real space whose distance from the origin is 1. The metric
$ds^2$ then tells you how to measure tiny distances along that curved surface.
The $\sin^2\theta$ factor records a familiar fact about globes: circles of
latitude shrink toward the poles, so a step in $\phi$ covers less ground the
closer you stand to either pole.

### Spherical Harmonics

Any function on the sphere can be expanded in **spherical harmonics**, $Y_\ell^m(\theta, \phi)$. These are the natural modes of vibration of the sphere.

The CMB temperature variations are analyzed by expanding in spherical
harmonics. The **power spectrum**, meaning how much power sits at each angular
scale $\ell$, records the hot dense side of the standard cosmological
reconstruction.

### Finite Resolution

At finite screen resolution, the smallest cell length sets a maximum angular
mode $\ell$:

$$\ell_{max} \sim \frac{R}{\ell_P}$$

Here $R$ is the radius of the region under discussion. The total number of
independent modes is roughly $\ell_{max}^2 \sim R^2/\ell_P^2$, proportional to area in Planck units, in line with the area scaling suggested by Bekenstein-Hawking.

Our experience of a continuous world is the large-scale limit of this finite
screen description.

### The Finite Machine Under the Chart

Three objects are easy to confuse. The **abstract observer patch** is the
operational package: what it can read, which state it carries, which records it
can consult, which boundary data it exposes, and which repairs it can perform.
The **support patch** is the cap, collar, or other geometric region that
displays those operations on a sphere. The **carrier patch** is the physical or
digital machine that realizes them. The observer is identified with this
visible operational package. A hidden circuit or material change is physically
silent only when it preserves the complete exposed structure: port incidence
and orientation, readouts, records, dynamics, repairs, refinement behavior,
and checkpoint continuation. Presentation labels can wash out. Architecture
cannot.

The reference carrier is called **echosahedral**, echo crossed with
icosahedral, a carrier that echoes its reads back through icosahedrally
arranged ports. It has a finite internal algebra and state, twelve overlap
ports, a readout at each port, an
observer-readable central record register, a finite mismatch score, a finite
menu of update and repair moves, and checkpoint data. In the icosahedral
screen sieve, the ports sit at the vertices of a regular icosahedron and
opposite ports form six axes. The twenty
triangular faces and thirty edges organize incidence, edge-sector, and collar
bookkeeping. Recurrent loops inside a patch or a small group of ports can
supply local memory and winding-sensitive behavior. With physical phase data
and a coupling law, they can also expose phase-lock-sensitive observables.
They are local subchannels, not a claim that the whole universe is a torus or
that every federation has proved phase locking.

There is a precise sense in which this arrangement, twelve ports at the
corners of an icosahedron with each paired against the one directly opposite,
is the cheapest one the bookkeeping allows, and the symmetry group of the
arrangement falls out for free. The details live in the technical papers. The
arrangement theorem does not by itself supply physical currents; attaching
those currents to the ports is work in progress.

The body is also not an observer by geometry alone. It becomes an elementary
observer only when its readback, records, feedback, predictions, and checkpoint
continuation pass the operational tests. A connected subfederation can satisfy
the same contract and function as a larger observer.

At one interface, the machine performs a short cycle:

```text
read local state
      ↓
expose a boundary packet
      ↓
compare it with the neighboring packet
      ↓
repair a checkable mismatch
      ↓
write the accepted result
      ↓
checkpoint and read again
```

A federation routes many such ports into shared interfaces. Accepted repairs
lower a declared mismatch score. Different compatible repair orders reach the
same public normal form only when the protected dependencies are complete,
conflicting updates have one coherent combined payload, the local-diamond test
passes, and the repair menu is complete.
Nested icosahedral subdivisions provide one convenient finite regulator for
the spherical support chart. They are not the twelve-port carrier federation,
and the map connecting repaired carrier data to the certified smooth sphere is
work in progress.

## 3.8 Patches and Overlaps

You cannot see the whole screen net. Some parts are hidden by your horizon or by
instrumental limits, and some parts are simply not in your operational domain.
The support chart displays that access domain as a connected region
$P_O\subset S^2$.

Another observer, at a different location or with different instruments,
has a different support region. Where support regions overlap, observers can
compare notes.

If the screen net is charted by $S^2$ and observer $i$ has support region $P_i$,
then two observers can compare data on the overlap $P_i \cap P_j$. That overlap
is the seed of consistency. Observer $i$ sees the local algebra, state, and
records displayed on $P_i$. What $i$ and $j$ can make public is the part of
those descriptions that agrees on $P_i \cap P_j$.

![Two observer support regions on the S^2 chart share a lens-shaped overlap where their descriptions can be compared.](../assets/book_diagrams/s2-screen.svg){width=74%}

### A Concrete Example

Consider two astronomers on opposite sides of Earth. During the night, they see different parts of the sky. But some stars are visible to both, the ones near the horizon for each observer.

These shared stars provide a link. The astronomers can calibrate by comparing their observations of the overlap region. Once they agree on the overlap, they can combine their observations into a consistent map of the whole sky.

### Coordinate Charts and Atlases

A sphere cannot be covered by a single smooth coordinate system. If you try to put latitude-longitude coordinates on a sphere, you run into problems at the poles.

Mathematicians handle this by using multiple overlapping coordinate charts, called an **atlas**. Each chart covers part of the sphere. Where charts overlap, there are transition functions that tell you how to convert coordinates.

This is analogous to our observer patches. Each observer has a local
description. Where observers overlap, they must agree on how to translate
between their descriptions. The atlas is the mathematical way to say what the
informal picture says physically: one shared net can be covered by many local
charts, and no finite chart has to contain the whole net.

Physics is the art of finding descriptions that work in many charts and have consistent translations between them.

## 3.9 What Is an Observer?

We've talked about "observers" and their "patches." But what exactly *is* an observer in this model?

### Not External Watchers

In classical physics, observers are implicitly outside the system: disembodied measurers who don't affect what they measure. This won't work here. Observers must be part of the system they observe.

### Observers as Patterns in the Data

An observer is an operational pattern in a patch or a connected patch
federation. It has bounded access to a local algebra, stable records it can
reread, an internal readback or state estimate, and future behavior that
depends on those records. It exposes predictions at its boundary and carries
enough checkpoint information to continue after a controlled repair. A human
observer can instantiate this role, as can a detector or a software process.
The definition names a self-reading structure rather than a species.

### The Vortex Analogy

Think of observers as stable vortices in a fluid.

The fluid is the finite carrier state exposed through the horizon chart:
constantly changing and highly correlated. A vortex is a pattern within the
fluid. It persists over time. It has a definite location. It interacts with
other patterns.

An observer is like that: a stable, self-reinforcing pattern within the
observer-visible data. The pattern has bounded access, maintains internal
structure and records, and can interact with nearby patterns through exposed
interfaces.

### Movement and Time

Do observers "move around" on the sphere?

Not in a simple sense. The support region can change from one record slice to
the next. Records, external interfaces, and checkpoint continuation carry the
operational identity while the spherical chart displays that change as motion.

What creates the sense of time? The internal structure of a quantum state
provides a natural algebraic ordering through **modular flow**. On the certified
geometric branch, an independent clock instrument calibrates that ordering as
physical duration. The thermal time principle supplies the guide; the state by
itself does not supply a wristwatch.

Chapter 11 explains the extra geometry, event correspondence, and calibration
needed before that parameter becomes the clock read by an observer.

### Why This Matters

This definition of observers resolves several puzzles:

**No external reference frame**: Observers are internal to the system, so there's no need for an external "God's-eye view."

**Measurement is physical**: When an observer measures something, correlations form between subsystems within the horizon data and stable records are created. That record formation captures the main physical content behind textbook collapse language.

**Consistency can follow from a verified structure**: Two observers may be restrictions of
one underlying state or begin from independently presented local data. In both
cases, translation, comparison, and certified repair can make their shared
descriptions agree on overlaps.

### Reality from Computation

A concrete screen picture looks like this.

Represent the screen chart by a **gauge-invariant quantum system** on the
2-sphere, something like a quantum cellular automaton with additional
structure. Triangulate the chart into tiny cells. At each edge of the triangulation sits a
finite-dimensional quantum system (a qudit). At each vertex, a gauge constraint
(Gauss's law) restricts which configurations are physical. Not all states
survive; only those satisfying the constraint at every vertex.

In plainer language, a qudit is a small quantum register with finitely many
possible readouts, like a qubit with more than two options. A gauge constraint
is a local accounting rule: many internal descriptions may be allowed, but only
the combinations that leave the shared physical quantities unchanged count as
physical.

**In this quantum-link regulator**, an abstract observer patch is
represented by a subsystem with a boundary-gauge-invariant algebra. Its support
region is like a computational thread through the chart, a connected place
where an observer can ask questions and get answers. The algebra
$\mathcal{A}(R)$ defines what that observer can measure: the operators that
commute with the boundary gauge transformations.

**Overlap consistency** becomes a physical interface test. Two patches expose
their gauge-invariant boundary packets, translate them into a common frame, and
compare the shared observables. A mismatch triggers one of the allowed local
repair moves and writes the accepted result to the record layer. Gauge
redundancy at the boundary makes the gluing non-trivial and gives rise to the
edge modes that carry geometric information.

**State selection and repair play different roles.** Maximum entropy selects
the realized state within a declared constraint family. At fixed cutoff this is
modeled by a Gibbs-like state $\rho \propto e^{-H}$, where $H$ is a sum of local
terms and $\rho$ is a quantum probability table, unpacked properly next
chapter. The patch update and repair instruments perform the synchronization
work when exposed records disagree. Evidence that both layers come from one
source and remain compatible under refinement is needed to keep state selection
and repair aligned across scales.

**The 4D bulk is not painted on the sphere.** On the successful geometric
branch, event spacetime is reconstructed from refined overlap records and their
entanglement structure. When you look around and see three-dimensional space,
you are seeing the effective geometry encoded by those relations. In the
constructions emphasized later, bulk distance is read from boundary
entanglement structure.

*The patch federation does the work. The screen is the chart. Reality is what observer patches agree on.*

In this book I sometimes call the chart readout a **folded screen**. The phrase
has a narrow meaning. Local patches compare overlaps, repair mismatches, and
settle to a quotient normal form. The spherical screen chart displays that
settled state as caps, collars, cuts, edge sectors, and boundary records.
Folding is the screen-facing presentation of repair.

This chart does real work. Caps and collars on the sphere identify the local
questions an observer can ask, and overlaps between caps identify the data two
observers can compare. On the certified smooth branch, the conformal
symmetries of the same sphere become the Lorentz symmetries of the shared
spacetime description, so sky directions and cap sides transform together. The
finite patch federation supplies the machine underneath that chart.

### A Physical or Digital Realization

The carrier can be built from quantum registers, optical channels, software
state machines, or another bounded material. Such implementations count as the
same echosahedral carrier only when they preserve the full exposed structure
listed earlier, from port incidence to checkpoint continuation, within their
operating tolerance. More general observer patches can have other bounded interface counts, and
their downstream symmetry and geometry need not match this branch.

A bench device makes the self-reading loop visible through its wiring,
calibration, raw readouts, controls, repair log, and repeatable checks. Those
records connect the physical carrier to the abstract observer-patch protocol.

## 3.10 Entanglement Creates Depth

The screen gives a boundary. Three-dimensional depth appears when entanglement
starts arranging the data into an interior.

When parts of a quantum state are strongly correlated, they behave as one connected structure. In holographic settings this relation becomes quantitative: boundary entanglement constrains bulk geometry. The Ryu-Takayanagi formula and related results make that statement precise in the regimes where they apply.

One lesson is enough here. Depth is read off from correlation structure.
Strongly linked regions count as nearby in the emergent bulk. Weakly linked
regions count as distant.

Chapter 9 develops this in detail. In the present chapter, entanglement does one job. It explains why a screen can support an interior world, not a flat catalog of data.

## 3.11 The Sphere Ladder

The screen notation can be read as a small ladder of logical roles.

$S^0$ is the first distinction, the seed of a readout. Something can be marked
as this, with silence as the alternative. $S^1$ is recurrence: a loop in which a record can
return to itself and be checked again. $S^2$ is the screen, the public archive
where finite observers expose overlap data. The final rung is the
three-dimensional observer-frame space, whose canonical kinematic chart is
hyperbolic space $H^3$. The ladder names operational roles rather than a global
$S^3$ topology.

This ladder is a teaching map, with the full particle taxonomy carried by the
later chapters. The photon belongs to the unbroken electromagnetic sector, gluons to color gauge
transport, the graviton to emergent geometry and diffeomorphism structure,
$W$, $Z$, and $H$ to the electroweak and Higgs sectors, and hadrons to QCD
composites. The ladder explains how OPH moves from seed, to loop, to screen, to
observer-frame geometry. Event spacetime, physical gauge currents, and particle
fields enter through their separate later constructions.

## 3.12 The Reverse Engineering

The reverse-engineering trail is short. The intuitive picture says
information scales with volume and space is the container. The hint is that
black-hole entropy scales with area, with gravitational entropy bounds pushing
toward boundary-limited information. The lesson is a boundary-first
description. In the symmetric construction used here, each observer has an
effective horizon naturally charted by a
spherical screen. The finite patch data exposed through that chart is limited
by $S\leq A/(4\ell_P^2)$. Its entanglement patterns create the geometry of the
emergent three-dimensional bulk, and overlap consistency makes that bulk shared
and stable across observers.

The holographic principle is the explanatory reading of these converging
thermodynamic and geometric facts.

## 3.13 Pixel Limits

The numbers are small enough to state directly.

The Planck length is $\ell_P \approx 1.6 \times 10^{-35}$ meters, about $10^{20}$ times smaller than a proton. The Planck area is $\ell_P^2 \approx 2.6 \times 10^{-70}$ m².

**The de Sitter horizon**: Add up the capacity of the horizon and you get
about $10^{122}$ bits, give or take a convention. The number is enormous, and
it is finite. The observable universe contains a finite amount of information.
It is a number worth remembering.

**A solar-mass black hole**: Schwarzschild radius $R_s \approx 3$ km. Number of bits: $N \approx 10^{77}$.

This is huge, but much smaller than the observable universe. Yet it's far more than the entropy of the Sun as a normal star (about $10^{58}$). Collapse increases entropy because the horizon has vastly more microstates than ordinary matter.

One reminder before leaving the numbers. Continuous space is an effective
large-scale approximation, and the spherical screen is the chart on which the
repaired patch data is displayed.

## 3.14 Where We Go Next

We have assembled four linked statements, each carrying its own load.
Gravitational
entropy bounds and holographic arguments push away from naive volume counting
and toward horizon-sensitive information organization. Separately, the
declared incidence, mesh, and cross-ratio receipts select a spherical support
on the certified branch. Lorentzian causal structure comes after that support
reconstruction, not before it. Area bounds constrain accessible information,
and entanglement patterns are the candidate data from which the effective
three-dimensional geometry must be recovered.

This chapter has separated the screen geometry from the synchronization cycle.
What gives the surviving records a direction? What creates the arrow of time?

The answer involves entropy again, this time in dynamics. The Second Law says entropy increases. But why? And what does this have to do with the screen?

In the next chapter, we explore the edge of the screen, the boundary conditions that govern what can happen. Entropy growth appears as a geometric constraint built into the structure of horizons themselves, alongside its statistical reading.

Chapter 4 turns the screen from a storage surface into a thermodynamic one.
