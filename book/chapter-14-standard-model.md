# Chapter 14: The Standard Model from Consistency

## 14.1 The Intuitive Picture: Particles and Forces Are Fundamental

The intuitive picture is straightforward. The universe is made of particles.
Forces act between them. The Standard Model is the final inventory of what
exists.

In this picture, an electron is a tiny object with definite properties, and
fields are invisible fluids that fill space. You learn the Standard Model as a
catalog: quarks, leptons, gauge bosons, the Higgs. That is the whole picture.

This view works for calculations. It also hides what is actually strange about
our best theory of matter.

## 14.2 The Surprising Hint: The Standard Model Is Not Fundamental

The Standard Model is extremely successful, and it carries deep warnings. Its
vacuum energy and loop integrals blow up in the ultraviolet. Its couplings run
with scale. Its anomaly cancellations are delicate. Its chirality is startling.
These clues point to an emergent effective description. The Standard Model is
not the foundation.

## 14.3 The Quantum Revolution

To understand what the Standard Model really says, we need to start with
quantum mechanics itself. Quantum mechanics is deeply, irreducibly weird.

### Planck's Desperate Act

In December 1900, Max Planck presented a formula to the German Physical
Society. He called it "an act of desperation."

The problem was blackbody radiation. When you heat an object, it glows. At low
temperatures, it glows red. Hotter, it glows white. The question was: how much
light at each wavelength?

Classical physics gave a disastrous answer. The Rayleigh-Jeans formula
predicted infinite energy at short wavelengths. Ovens should emit deadly gamma
rays. This was the "ultraviolet catastrophe."

Planck found a formula that fit the data extremely well. To derive it, he had
to assume something absurd: energy comes in discrete packets. Light of
frequency $f$ carries energy in multiples of $hf$, where $h$ is a tiny
constant.

$$E = nhf, \quad n = 0, 1, 2, 3, \ldots$$

Planck did not believe this was real physics. He thought it was a mathematical
trick. It took Einstein to show it was genuine.

### Einstein's Light Quanta

In 1905, Einstein explained the photoelectric effect. When light hits metal,
electrons pop out. The energy of those electrons depends only on the light's
frequency, not its intensity. Brighter light produces more electrons, not
faster ones.

Einstein's explanation: light really does come in packets. A photon of
frequency $f$ carries energy $hf$. One photon kicks out one electron. The
photon's frequency determines the electron's energy.

This was radical. For two centuries, physicists had proven that light was a wave. Young's double-slit experiment showed interference patterns. Maxwell's equations described electromagnetic waves. Einstein was saying light was particles?

Both were true. Light is neither purely wave nor purely particle. It's something new that exhibits both behaviors depending on how you probe it.

### Bohr's Atom

In 1913, Niels Bohr proposed a model of the hydrogen atom. Electrons orbit the
nucleus, but only in specific orbits. When an electron jumps between orbits, it
emits or absorbs a photon.

The model was frankly bizarre. Why should only certain orbits be allowed? Bohr
had no answer. He declared that angular momentum must be quantized:

$$L = n\hbar, \quad n = 1, 2, 3, \ldots$$

The model worked brilliantly for hydrogen. It explained the Balmer series, the
specific wavelengths of light that hydrogen emits. It failed for everything
else. Helium was a mess. The model was obviously incomplete.

### de Broglie's Audacity

In 1924, Louis de Broglie made a wild proposal in his PhD thesis. If light
waves can behave like particles, maybe particles can behave like waves.

He proposed that every particle has an associated wavelength:

$$\lambda = \frac{h}{p}$$

where $p$ is momentum. For everyday objects, this wavelength is absurdly tiny.
A baseball's de Broglie wavelength is about $10^{-34}$ meters. For electrons,
it is comparable to atomic sizes.

In 1927, Davisson and Germer proved de Broglie right. They bounced electrons off a nickel crystal and saw interference patterns. Electrons really do behave like waves.

### Schrödinger's Equation

Erwin Schrödinger took de Broglie's idea and ran with it. If electrons are
waves, what is waving?

Schrödinger proposed that electrons are described by a wave function
$\psi(x,t)$. The equation governing this wave is:

$$i\hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m}\nabla^2\psi + V\psi$$

This is the Schrödinger equation, and it works spectacularly well. It predicts
atomic spectra, chemical bonds, and semiconductor behavior. It is the
foundation of quantum chemistry and materials science.

What is $\psi$? Schrödinger initially thought it described a smeared-out
electron, spread across space like a cloud. Max Born had a different
interpretation: $\psi$ squared gives the probability of finding the electron at
each location.

$$P(x) = |\psi(x)|^2$$

Operationally, the wave function does not assign a classical trajectory. It gives the probabilities for different measurement outcomes.

The early formulas introduce the basic quantum dictionary. In Planck's
$E=nhf$, $E$ is energy, $n$ is a whole-number quantum count, $h$ is Planck's
constant, and $f$ is frequency. In Bohr's $L=n\hbar$, $L$ is angular momentum
and $\hbar=h/(2\pi)$. In de Broglie's $\lambda=h/p$, $\lambda$ is wavelength
and $p$ is momentum. In Schrödinger's equation, $\psi$ is the wave function,
$m$ is mass, $V$ is potential energy, and $\nabla^2$ measures spatial
curvature of the wave. Born's rule, $P(x)=|\psi(x)|^2$, turns the wave
function into a probability density for detection at position $x$.

That dictionary was assembled by many people under pressure from experiment.
Planck's blackbody curve, Einstein's photons, Bohr's spectral lines, de
Broglie's matter waves, Schrödinger's wave mechanics, Heisenberg's matrices,
Born's probability rule, Dirac's relativistic equation, and Feynman's diagrams
are different steps in one long reconstruction. The Standard Model inherits
that whole history.

### Heisenberg's Uncertainty

Werner Heisenberg approached quantum mechanics differently. He focused on observables: things you can actually measure.

In June 1925, suffering from hay fever on the island of Helgoland, Heisenberg developed matrix mechanics. Observable quantities became matrices. When he tried to calculate, he discovered something strange: the order of multiplication matters.

Position times momentum is not the same as momentum times position:

$$XP - PX = i\hbar$$

This commutation relation is the mathematical heart of quantum mechanics. It implies the uncertainty principle:

$$\Delta x \cdot \Delta p \geq \frac{\hbar}{2}$$

You cannot simultaneously know both position and momentum with arbitrary precision. This is a fundamental feature of reality. There is no state that has both precise position and precise momentum.

### The Copenhagen Interpretation

Bohr and Heisenberg developed what became the "Copenhagen interpretation." The wave function doesn't describe objective reality. It describes our knowledge. When we measure, the wave function "collapses" to a definite value.

This interpretation was never universally accepted. Einstein famously objected: "God does not play dice." The mathematics works. Quantum mechanics makes predictions, and those predictions are confirmed to extraordinary precision.

The core lesson is operational. Quantum theory gives probabilities for measurement outcomes with extraordinary accuracy. What those probabilities mean ontologically depends on the interpretation.

## 14.4 From Particles to Fields

Quantum mechanics describes particles. But particles can be created and destroyed. An electron and positron can annihilate into photons. A photon can create an electron-positron pair. How do you write a wave function for a variable number of particles?

You don't. You need quantum field theory.

### Dirac's Equation

In 1928, Paul Dirac sought a relativistic version of Schrödinger's equation. He found something deeper.

The Dirac equation describes spin-1/2 particles like electrons:

$$i\hbar \gamma^\mu \partial_\mu \psi - mc\psi = 0$$

The equation had a problem: it predicted states with negative energy. An electron could fall into these states, releasing infinite energy.

The matrices $\gamma^\mu$ are Dirac gamma matrices. They package spin and
relativity into one algebraic object. The derivative $\partial_\mu$ measures
change in the spacetime direction $\mu$. The field $\psi$ is a spinor
field, not a single nonrelativistic wave, and $mc$ carries the particle
mass scale. Dirac's compact line says that spin, antimatter, and special
relativity belong together.

Dirac's solution was audacious. The negative energy states are filled. The vacuum is a sea of negative-energy electrons. What we call a "positron" is a hole in this sea.

This prediction was confirmed in 1932 when Carl Anderson photographed positron tracks in a cloud chamber. Antimatter exists.

### Second Quantization

The Dirac sea was a stepping stone. The modern view is cleaner: fields are the fundamental objects, and particles are excitations of fields.

Consider a violin string. The string can vibrate in different modes. Each mode has a definite frequency. When you pluck the string, you excite various modes.

Quantum fields work similarly. The electromagnetic field can be decomposed into modes. Each mode is a quantum harmonic oscillator. Exciting a mode means adding photons.

The vacuum isn't empty. It's the ground state of all fields. Every mode is in its lowest energy state. But even the ground state has fluctuations. These zero-point fluctuations are real and measurable.

### Feynman Diagrams

Richard Feynman developed a beautiful pictorial language for particle physics. Draw space horizontally and time vertically. Particles are lines. Interactions are vertices where lines meet.

An electron emitting a photon:

```
    e- ---•--- e-
          |
          γ
```

The power of Feynman diagrams is that each diagram corresponds to a mathematical expression. You can calculate by drawing pictures.

To find the probability of a process, you draw all possible diagrams and add them up. This is perturbation theory. It works when interactions are weak.

### Renormalization

There's a catch. When you calculate loop diagrams, you get infinities.

Consider an electron. It's surrounded by a cloud of virtual photons. These photons affect the electron's mass and charge. When you calculate this effect, you get infinity.

The solution is renormalization. You absorb the infinities into the definition of mass and charge. The "bare" parameters are infinite, but the physical parameters are finite.

This sounds like cheating, but it works with astonishing precision. Quantum electrodynamics (QED) predicts the electron's magnetic moment to 12 decimal places. The prediction agrees with experiment to extraordinary precision.

Renormalization works for some theories (called "renormalizable") but not others. The Standard Model is renormalizable. Perturbative Einstein gravity is not. This is one reason gravity remains outside the Standard Model.

### Running Couplings

A strange consequence of renormalization: coupling constants change with energy.

The fine-structure constant measures the strength of electromagnetism. OPH keeps
three readings separate. The undressed source/root pixel gives
$136.9948351646\ldots$. Combining that root with the finite-screen unified
gauge width evaluated at the CODATA-derived comparison pixel gives the
mixed-provenance diagnostic $137.0359595008\ldots$. The number measured in the
laboratory is the long-distance Thomson endpoint, $137.035999177(21)$. A
source-only endpoint requires the missing same-scheme hadronic spectral backend.

That low-energy number sits inside the same particle sector as the weak bosons.
Once the electroweak transport family is read from the selected local fixed
point, electromagnetism is the unbroken channel left after the weak and
hypercharge sectors mix. The fine-structure constant therefore belongs to the
same transport family that yields the $W$ and $Z$ benchmarks.

The strong force coupling runs the opposite way. At low energies, it's strong (hence the name). At high energies, it weakens. This is "asymptotic freedom," discovered by Gross, Wilczek, and Politzer in 1973.

Running couplings mean the "constants" of physics aren't constant. They depend on the scale at which you probe.

## 14.5 The Standard Model Zoo

The Standard Model organizes all known particles into a coherent model.

### Fermions: The Matter Particles

Matter is made of fermions: particles with spin 1/2. They obey the Pauli exclusion principle. No two identical fermions can occupy the same quantum state. This gives atoms structure, gives us the periodic table, and keeps you from falling through the floor.

**Quarks** come in six flavors. Up, charm, and top carry charge $+2/3$. Down,
strange, and bottom carry charge $-1/3$.

Quarks are never found alone. They're always bound into hadrons by the strong force. Protons are (uud), neutrons are (udd).

**Leptons** also come in six types. The electron, muon, and tau carry charge
$-1$. Their three neutrinos are neutral.

The electron is stable. The muon and tau decay quickly.

### Three Generations

The fermions come in a strange pattern: three copies. The up and down quarks, plus the electron and its neutrino, form the first generation. The charm and strange quarks, plus the muon and its neutrino, form the second. The top and bottom, plus the tau and its neutrino, form the third.

Historically, the Standard Model by itself does not explain why there are three generations. OPH later argues that, on its realized one-Higgs branch, the minimal viable count is three. The charged members of the second and third generations are heavier copies of the first, while the neutrino sector has its own mixing pattern. Almost all ordinary matter uses only first-generation particles.

### Bosons: The Force Carriers

Forces are mediated by bosons: particles with integer spin.

**Photon** (spin 1): The observed quantum of the electromagnetic field. In the
ordinary Maxwell vacuum it has a massless pole, travels at the invariant null
speed, and couples to electric charge.

**W and Z bosons** (spin 1): Carry the weak force. W has charge plus or minus 1. Z is neutral. Both are massive: about 80-90 GeV. The weak force is weak at low energies because its carriers are heavy.

**Gluons** (spin 1): Carry the strong interaction in perturbative descriptions.
There are eight color components. The pure Yang--Mills quadratic action has no
hard mass term, but confined QCD has no isolated free-gluon particle pole.

The Yang-Mills mass gap is a statement about the spectrum of the compact
nonabelian gauge theory, separate from assigning a hard mass to the gluon. In
OPH, the gap is accounted for by repair dynamics on the declared
controlled compact-gauge construction: exact local repair gives a positive
Euclidean relaxation generator, and its first nonzero repair eigenvalue is the
first nonzero Yang-Mills energy.

**Higgs boson** (spin 0): The source of mass for W, Z, and fermions. Discovered at CERN in 2012. Mass about 125 GeV.

**Graviton** (spin 2): The hypothetical carrier of gravity. Not part of the Standard Model. Never directly detected.

### The Gauge Groups

The Standard Model is organized by symmetry. One usually writes the gauge group as:

$$G_{SM} = SU(3)_C \times SU(2)_L \times U(1)_Y$$

The notation names three continuous accounting systems. $SU(3)$ is a
three-component special-unitary symmetry, $SU(2)$ is its two-component cousin,
and $U(1)$ is the circle-like symmetry behind a single conserved charge. The
subscripts say which physical bookkeeping each factor carries.

$G_{SM}$ means "the Standard Model gauge group." The subscript $C$ means color.
The subscript $L$ means left-handed weak isospin. The subscript $Y$ means
hypercharge, the charge that mixes with weak isospin to produce ordinary
electric charge after symmetry breaking. The product sign means the three
symmetry systems are present together.

**SU(3)_C** is the color group. Quarks carry color charge: red, green, or blue. Gluons carry color-anticolor combinations. The strong force binds quarks into colorless combinations.

**SU(2)_L** is the weak isospin group. It acts only on left-handed particles.
The weak force therefore violates parity.

**U(1)_Y** is the hypercharge group. It combines with SU(2)_L to give electromagnetism after symmetry breaking.

The subscripts matter. L means "left-handed." The weak force distinguishes left from right. This is one of nature's deepest asymmetries.

## 14.6 Chirality: Nature's Handedness

Nature treats left and right differently. This is one of the deepest
asymmetries in the Standard Model.

### What Is Chirality?

A relativistic fermion has a left-handed face and a right-handed face. For
massless particles, that split lines up with helicity, with spin either
tracking the motion or leaning against it. For massive particles the relation
is subtler, but the left-right split remains built into the theory.

Helicity is the easy visual version: compare the direction of a particle's spin
with its direction of travel. Chirality is the deeper field-theory label. For
massless particles they coincide; for massive particles they do not have to.

### The Weak Force Discriminates

The charged weak interaction carried by the $W$ boson couples only to
left-handed fermions. A right-handed electron sits out those charged-current
processes.

This was discovered through parity violation experiments in 1956-1957. Chien-Shiung Wu studied the beta decay of cobalt-60. If parity were conserved, electrons should emerge equally in both directions along the spin axis. They didn't. More electrons came out opposite to the spin.

Lee and Yang had predicted this. Wu proved it. Parity violation earned Lee and Yang the Nobel Prize. Wu, who did the experiment, was not included.

### Why Chirality Matters

Chirality matters everywhere. It is essential to weak parity violation and to
anomaly-cancellation constraints, and it sharply restricts which fermion mass
terms can appear without extra structure.

## 14.7 Anomaly Cancellation: Why the Charges Are What They Are

Consider the electric charges of quarks and leptons. At first glance they look
arbitrary: the up quark at $+2/3$, the down quark at $-1/3$, the electron at
$-1$, the neutrino at $0$. The real explanation is anomaly cancellation.

### What Is an Anomaly?

A quantum theory can look symmetrical in its classical dress and tear at
the seams once quantization is done. That failure is an anomaly. If it hits a
gauge symmetry, the theory stops being self-consistent.

### The Cancellation

The Standard Model survives because one generation of quarks and leptons
cancels every dangerous hypercharge contribution at once. Color, weak charge,
the cubic hypercharge sum, and the gravitational sum all close together.

The famous charges do not float freely. Thirds of an electron
charge are not decorative details. They are the values that let the structure
hold.

### Connection to OPH

The same issue appears in geometric dress. Glue observer patches around
a loop and return to the starting point. If some leftover phase remains, the
gluing tears. Field theory calls that failure an anomaly. The screen picture
calls it bad loop bookkeeping. Either way the cure is the same: the charge
assignments must make the loop close cleanly.

The Standard Model's hypercharges look so crisp for that reason. Up to overall
normalization, they are the solution that lets the gluing hold together.

## 14.8 The Higgs Mechanism

The Standard Model has a puzzle. A pure gauge kinetic action has massless
quadratic modes, yet the W and Z are massive. Gauge redundancy can coexist
with their mass because the Higgs field changes the physical phase.

### Spontaneous Symmetry Breaking

Consider the Higgs potential:

$$V(\phi) = -\mu^2|\phi|^2 + \lambda|\phi|^4$$

This is symmetric under rotations in field space. But the minimum isn't at zero. It's in a circular valley at radius v = mu/sqrt(lambda).

$\phi$ is the Higgs field. $\mu$ and $\lambda$ are parameters of the
potential. The negative quadratic term pushes the field away from zero, while
the positive quartic term keeps the energy from running off to infinity. The
nonzero radius of the valley is the vacuum expectation value that feeds masses
to the weak bosons and fermions.

The field "falls" to some point in this valley. The symmetry is broken spontaneously. The equations are symmetric; the ground state is not.

That settled nonzero value is called the vacuum expectation value: the
background value of the Higgs field that other particles move through.

### Eating Goldstone Bosons

When a continuous symmetry is spontaneously broken, massless particles appear: Goldstone bosons. They correspond to motion along the valley.

In a gauge theory, something special happens. The gauge bosons "eat" the Goldstone bosons and become massive. This is the Higgs mechanism.

For the electroweak group SU(2) x U(1), three Goldstone bosons get eaten. The
W+, W-, and Z become massive. One combination of generators remains unbroken.
On the ordinary vacuum with the Maxwell kinetic term, that electromagnetic
combination has the familiar massless transverse pole.

### Fermion Masses

Fermions also get mass from the Higgs. The Yukawa couplings connect left-handed and right-handed fermions through the Higgs field:

$$\mathcal{L}_{Yukawa} = y_e \bar{L} \phi e_R + y_u \bar{Q} \tilde{\phi} u_R + y_d \bar{Q} \phi d_R + \text{h.c.}$$

This line is a compact part of the Lagrangian, the formula that says which
field interactions are allowed. The $y$ values are Yukawa couplings. They set
how strongly each fermion talks to the Higgs field, and therefore how much mass
that fermion gets after symmetry breaking.

The barred fields are conjugate fields. $L$ is a left-handed lepton doublet,
$Q$ is a left-handed quark doublet, and $e_R$, $u_R$, and $d_R$ are
right-handed charged-lepton, up-type-quark, and down-type-quark singlets.
$\tilde{\phi}$ is the Higgs doublet arranged with the conjugate weak charge.
"h.c." means Hermitian conjugate, the companion term required to make the
Lagrangian real.

When the Higgs gets a vacuum expectation value, these terms become mass terms. The masses are proportional to the Yukawa couplings.

Why do the Yukawa couplings have the values they do? Why is the top quark so much heavier than the electron? The Standard Model leaves this unexplained.

## 14.9 From Overlaps to Gauge Structure

The OPH connection is direct.

### Gauge as Gluing Redundancy

In the standard presentation, gauge symmetry is a postulate. You write down a Lagrangian that's invariant under local transformations.

A local transformation is a change of internal description that can vary from
point to point. Gauge symmetry says such changes must not alter physical
predictions.

Gauge symmetry emerges from the redundancy in how observers glue their patches
together.

Different observers describe the same overlap region using different frames. The transformation between frames is a gauge transformation. The freedom that leaves overlap observables invariant forms the gauge group.

This is "gauge-as-gluing." Gauge symmetry isn't fundamental. It's the grammar of how patches fit together.

The same receipt-certified branch also carries the Yang-Mills form. Compact-gauge reconstruction,
local holonomy data, the four-dimensional continuum scaling chart, the
reflection-positive ordinary vacuum sector, and the absence of additional
gauge-invariant relevant dimension-four pure-gauge operators beyond the
curvature-squared invariant produce the Euclidean action. Repair collars supply
the spectral part: leaving the repaired vacuum sector costs a positive amount
of Euclidean relaxation energy, and that cost is the Yang-Mills mass gap on the
controlled compact-gauge construction.

Once the edge transports behave like a compact gauge group, the long-distance
field theory has the usual Yang-Mills action. That action is the cost of
bending the gauge connection. The mass gap is the first nonzero energy cost for
leaving the repaired vacuum sector. The OPH claim keeps the branch conditions
visible because a slogan about gluing is weaker than a field-theory continuum
statement.

### Edge-Center Completion

When you have a boundary between patches, there are degrees of freedom that live on the edge. These edge modes carry "charges" that label how the two sides connect.

Technically, the Hilbert space decomposes:

$$\mathcal{H}_{collar} = \bigoplus_\alpha (\mathcal{H}_{left}^\alpha \otimes \mathcal{H}_{right}^\alpha)$$

The direct-sum symbol means the boundary data split into sectors labeled by
$\alpha$. The tensor-product symbol joins the left and right sides inside one
sector. The formula is just the precise way to say that an edge carries a label
both neighboring patches must respect.

The labels alpha are the edge charges visible in this one-collar algebra. In
the bosonic gauge picture they are the seed carriers for the tensor category
from which the boundary gauge group is reconstructed. They are not assumed to
list every charge generated when seeds are fused. A fixed collar can therefore
have finitely many visible blocks while the tensor-generated category has
infinitely many simple charge types.

One tensor category uses one common strict transport choice at each cutoff.
Only visible blocks with trivial closed-loop action under that same choice are
retained as seeds. Blocks that become transportable only after incompatible
strictifications belong to alternative fixed-stage categories; their union is
not silently closed under tensor product.

The letter $\mathcal H$ names a Hilbert space, the quantum state space for a
piece of the system. "Collar" means the thin overlap zone near a boundary. The
superscript $\alpha$ says that each left and right Hilbert space belongs to one
shared edge-charge sector. The formula is the bookkeeping form that makes
boundary agreement possible, not an extra postulate about particles.

### Fusion Rules Define the Group

When you concatenate collars, edge charges fuse. The fusion rules:

$$\alpha \otimes \beta = \bigoplus_\gamma N_{\alpha\beta}^\gamma \, \gamma$$

define a tensor category. The Tannaka-Krein reconstruction theorem says,
roughly, that if the fixed-cutoff seed carriers generate a rigid tensor
category, remain transportable across patches, persist coherently under a
certified refinement system, and carry compatible finite-dimensional
forgetful fibers, then the resulting refinement-limit category recovers the
compact symmetry group behind them. The fusion table is central, and it is
used together with the refinement transport and carrier data. The group is
read off from how charges behave in that full persistent sector package.

For intuition, treat the fusion rules as a multiplication table for charges.
If you know how every charge combines with every other charge, you have enough
information to recover the symmetry that those charges are representing.

The labels $\alpha$, $\beta$, and $\gamma$ are charge sectors. The tensor
symbol $\otimes$ means "combine these sectors." The integers
$N_{\alpha\beta}^{\gamma}$ count how many times sector $\gamma$ appears when
$\alpha$ and $\beta$ fuse. A tensor category is the organized collection of
these sectors, their fusions, their duals, and their consistency rules.
It is a bookkeeping machine for charges: which charges exist, how they combine,
which charge is the mirror of which, and which combinations count as the same
operation in different orders.

On the receipt-certified branch, some compact gauge group is reconstructed from
the persistent sector data rather than guessed in advance.

![Tannaka-Krein reconstruction reads a compact gauge group from the way edge sectors fuse and represent one another.](../assets/book_diagrams/tannaka-krein.svg){width=82%}

The refinement step needs its own receipt; ordinary state coarse-graining does
not provide it. At each finer cutoff the receipt must certify that every coarse
record has a fine extension, exhibit a surjective map of compact boundary
actions, give a center-compatible embedding of the finite collar blocks with
explicit multiplicities, and provide compatible finite tensor realizations and
forgetful fibers. Pulling a representation back along the surjective group
map preserves its carrier and cannot split a coarse simple. These pullbacks
compose along a cofinal refinement chain, whose endpoint is the directed
colimit. Without that receipt, only the fixed-stage categories exist and the
refinement-limit gauge group remains conditional.

For example, if a U(1) collar sees charge one and the common strict transport
retains it, tensor closure contains charge
two and, in fact, every integer charge. The charge-two projector belongs to the
two-fold tensor-word or concatenated-collar realization; it need not appear in
the original one-collar center. Its forgetful fiber is still one-dimensional
and remains the same one-dimensional carrier at the next certified refinement.

### The Standard Model Factors

Why does the realized group have the form SU(3) x SU(2) x U(1) up to finite quotient?

A quotient means that some formally different group elements act the same on
all physical states and are counted once. It is like discovering that two labels
on a wiring diagram name the same actual connection. The Standard Model quotient
removes that duplicate counting across color, weak isospin, and hypercharge.

On that receipt-certified branch, the reconstruction step by itself gives some compact gauge group from the
transportable sector category. The Standard Model choice enters through minimal admissible
selection: among the admissible one-Higgs
low-energy sector packages, choose the minimal package that satisfies the
transport, anomaly, refinement, Yukawa, CP, and weak-sector conditions.
Equivalently, the logic factors into classification and selection. Overlap
gluing first tests whether triangle or higher associator defects can be
strictified. It then tests the ordinary holonomy across the allowed strict
representatives and requires at least one to act trivially around every closed
loop. To form one tensor category, one such representative is then fixed for
the whole stage, and only sectors on which that same representative has
trivial loop action are retained. Only those sectors form the
transportable category from which tensor reconstruction reads a compact group.
The selector acts on realized matter packages only after that.

The two tests detect different mismatches. On the central branch, the degree-2
class records the projective triangle defect, and its vanishing permits a strict
edge 1-cocycle. On the noncentral branch, the full crossed-module orbit can admit
several strict representatives related by higher-gauge edge changes, so it need
not determine one ordinary 1-cocycle class. The associator-only test asks whether
at least one strict representative exists. In either case an allowed strict
cocycle may still have nontrivial holonomy around a noncontractible loop, and
strict transport requires that at least one allowed representative have trivial
represented holonomy.
"Zero obstruction" in the gauge-reconstruction branch therefore means both
associator strictifiability and the existence of a trivial-holonomy strict
representative for an individual sector; tensor reconstruction additionally
fixes one common representative for all seeds in the chosen category.

Once you ask for the smallest matter sector that can carry color, weak
interactions, chirality, and ordinary charge, the answer is forced into a
color triplet, a weak doublet, and one abelian charge direction. The weak
factor has to behave like $SU(2)$ because weak doublets come in the right
two-dimensional pseudoreal form. The color factor has to behave like $SU(3)$
because color triplets need a genuinely complex three-dimensional action. Once
those two are in place, the remaining commuting charge direction is $U(1)$,
and the sixth-integer hypercharge pattern sharpens the result to the Standard
Model quotient.

The representation words only say how a particle multiplet transforms. A weak
doublet is a two-entry object rotated by the weak symmetry. A color triplet is a
three-entry object rotated by the color symmetry. "Pseudoreal" and "complex"
distinguish whether the mirror representation is effectively the same object or
a genuinely different one.

The same low-energy sector also fixes the counting. The minimal coupled carrier makes the quark doublet a color triplet and therefore fixes $N_c=3$. On the same one-Higgs quark branch, intrinsic CKM CP capability requires at least three generations, weak-sector ultraviolet consistency keeps the count finite, and the smallest viable answer is $N_g=3$. The Witten anomaly then remains as a consistency check on the resulting triplet-doublet structure. This anomaly is a global $SU(2)$ obstruction: the theory is consistent only when the number of left-handed weak doublets is even.

## 14.10 Hypercharge from Gluing Consistency

Given the gauge group, what determines the matter content?

### The Anomaly Condition Again

Loop-coherent gluing requires a trivial central obstruction class and at least
one allowed strict edge transport with trivial represented holonomy around every
closed overlap loop. This combined transportability condition is not a Standard
Model selector. In a chiral effective field theory, the same consistency burden
reappears as anomaly
cancelation, but the full bridge between the two is a separate step.

Given one generation of chiral fermions with
$SU(3)\times SU(2)\times U(1)$ charges, and requiring Yukawa couplings to a
Higgs doublet, the hypercharge ratios are determined. A standard normalization
then fixes the absolute lattice.

### The Derivation

Start with Yukawa invariance. The Higgs coupling has to be neutral under
hypercharge, so the charges of the left-handed doublets, right-handed
singlets, and Higgs must add up in the allowed way:

$$Y_u = Y_Q + Y_H, \quad Y_d = Y_Q - Y_H, \quad Y_e = Y_L - Y_H$$

Add the anomaly cancellation conditions. The first line says that the weak
doublets cannot leave a mixed weak-hypercharge anomaly:

$$N_c Y_Q + Y_L = 0 \quad (SU(2)^2 U(1))$$

The second line is the mixed gravitational condition. It says the chiral
hypercharge assignment must remain consistent when the fermions couple to
gravity:

$$2N_c Y_Q - N_c Y_u - N_c Y_d + 2Y_L - Y_e = 0 \quad (\text{gravitational})$$

Solving those constraints first fixes the lepton and Higgs charges in terms of
the quark-doublet charge:

$$Y_L = -N_c Y_Q, \quad Y_H = N_c Y_Q$$

The right-handed singlet charges then follow from the Yukawa relations:

$$Y_u = (N_c+1)Y_Q, \quad Y_d = -(N_c-1)Y_Q, \quad Y_e = -2N_c Y_Q$$

With $N_c=3$ and standard normalization, the familiar lattice appears:

$$\boxed{Y_Q = \frac{1}{6}, \quad Y_L = -\frac{1}{2}, \quad Y_u = \frac{2}{3}, \quad Y_d = -\frac{1}{3}, \quad Y_e = -1, \quad Y_H = \frac{1}{2}}$$

These are exact rationals, the Standard Model hypercharges, with the ratios
fixed by anomaly freedom together with Yukawa invariance and the absolute
values fixed by standard normalization. There is nothing to tune. The
sixth-integer lattice is exactly the one compatible with the physical quotient
\((SU(3)\times SU(2)\times U(1))/\mathbb Z_6\).

The $Y$ symbols are hypercharges. $Q$ labels the left-handed quark doublet,
$L$ the left-handed lepton doublet, $H$ the Higgs doublet, and $u$, $d$, and
$e$ the up-type quark, down-type quark, and charged lepton singlet sectors.
$N_c$ is the number of colors. The boxed line is the familiar charge lattice
written before electroweak mixing turns hypercharge and weak isospin into
ordinary electric charge.

That is what makes the derivation satisfying. The equations are doing one job:
they explain why the charges come out in the strange pattern we observe.
Quarks carry third-integer charges because the weak interaction, the Higgs
couplings, and anomaly cancellation all have to coexist in one self-consistent
chiral theory.

## 14.11 The Number of Colors: Why N_c = 3

In the full argument, the color count is fixed directly by the same coupled
carrier that emits the $SU(3)$ factor. The global $SU(2)$ anomaly is an
important check on the realized structure, while the coupled carrier determines
the count.

### The Coupled Color Carrier

The weak sector needs a pseudoreal doublet. The color sector needs a genuinely
complex nonabelian role. The smallest common carrier that supports both on one
block is

$$\mathbb C^3 \otimes \mathbb C^2.$$

That fixes the quark doublet to be a color triplet:

$$\boxed{N_c = 3}$$

This is the decisive structural step. The color count is emitted by the same
minimal coupled carrier that produces the $SU(3)$ factor, not by a later oddness
argument.

### The Witten Check

The global $SU(2)$ anomaly must cancel on the realized branch. Each
generation contributes $N_c$ quark doublets and one lepton doublet, so the
number of left-handed $SU(2)$ doublets per generation is

$$N_c + 1.$$

With $N_c=3$, this becomes

$$N_c + 1 = 4,$$

which is even. So Witten's anomaly is satisfied generation by generation. In
this derivation it confirms the realized triplet-doublet structure. It does not
select the color count.

## 14.12 Why Three Generations?

Anomaly cancellation works generation by generation. Each generation independently satisfies the conditions. So why three?

### CKM CP Capability Requires Three

The CKM matrix describes how quarks mix under the weak force. In general, it is
a unitary $N_g\times N_g$ matrix. The number of physical CP-violating phases
is:

CP means charge-parity reversal: swap particles with antiparticles and mirror
space. A CP-violating phase is a built-in complex phase that lets those mirrored
processes differ. In ordinary language, it is one source of particle-antiparticle
rate differences in weak interactions.

$$\text{(CP phases)} = \frac{(N_g - 1)(N_g - 2)}{2}$$

For $N_g=1$ or $N_g=2$, the formula gives zero phases. For $N_g=3$, it gives
one phase. The third generation is the first case with intrinsic CKM CP
capability.

So the realized quark branch requires at least three generations:

$$N_g \ge 3$$

### Weak-Sector UV Completability Limits

Too many generations spoil asymptotic freedom. The $SU(2)$ beta function
coefficient is:

Asymptotic freedom means an interaction gets weaker at shorter distances or
higher energies. The beta function is the bookkeeping rule for how a coupling
changes with scale.

$$b_{SU(2)} = \frac{22}{3} - \frac{1}{3}N_g(N_c + 1) - \frac{1}{6}$$

The final $-1/6$ is the contribution of one Higgs doublet. For
$b_{SU(2)} > 0$ (asymptotic freedom):

$$N_g(N_c + 1) < \frac{43}{2}.$$

With $N_c=3$, this becomes

$$4 N_g < \frac{43}{2} \implies N_g \le 5$$

Combining the lower and upper bounds gives the viable window:

$$3 \le N_g \le 5.$$

### The Minimal Viable Window

CKM CP capability and weak-sector UV completability define the viable window.
Here UV completability means that the theory can keep making sense at shorter
distances and higher energies, with no immediate breakdown when the resolution
is increased.

A minimal admissible realization principle then picks the smallest viable
realization. "Minimal admissible" means the smallest option that satisfies the
listed consistency tests:

$$\boxed{N_g = 3}$$

The one-Higgs slot also has a clean local geometric carrier. On the selected
electroweak branch, the weak screen chart can be modeled locally as
$\mathbb{CP}^1$ with the minimal positive Hopf line bundle $\mathcal O(1)$.
Borel-Weil then gives

$$H^0(\mathbb{CP}^1,\mathcal O(1))\cong\mathbb C^2.$$

This is the first nontrivial holomorphic section space, so it supplies the
weak doublet carrier. OPH fixes the hypercharge convention with the neutral
component condition $Q(\phi^0)=T_3+Y=0$, giving $Y=+1/2$. A nonzero section
direction is a point of the same projective line, but that ray cannot determine
the unbroken electromagnetic group. Hypercharge multiplies the whole doublet by
a common phase, which projectivization erases. For the nonzero lower-component
vacuum vector

$$\phi_0=\frac{v}{\sqrt 2}\binom{0}{1},\qquad v\ne0,$$

the electroweak action is

$$e^{i\alpha T_3}e^{i\beta Y}\phi_0
=e^{i(\beta-\alpha)/2}\phi_0.$$

Independent $T_3$ and hypercharge phases therefore fix the ray
$[\phi_0]$. Its projective stabilizer is a two-torus,
$U(1)_{T_3}\times U(1)_Y$, modulo the inherited finite center. The vector itself
is fixed only when $\beta=\alpha$, locally, leaving the electromagnetic
$U(1)_Q$ generated by $Q=T_3+Y$. The projective line explains the scalar
carrier; the nonzero vacuum vector explains why electromagnetism remains
unbroken. This construction does not explain the Higgs mass, the quartic, or
the weak scale; those belong to the OPH hierarchy and Higgs/top quantitative
branch.

![The generation-count window starts at three for intrinsic CP capability and closes above five from weak-sector ultraviolet consistency.](../assets/book_diagrams/generation-count.svg){width=84%}

Refinement stability explains why extra unfixed Yukawa structure is disfavored. Among the allowed options, the smallest viable one wins. With $N_c=3$ and $N_g=3$, each generation carries four left-handed weak doublets, an even number, so the Witten anomaly is automatically satisfied on the realized branch.

## 14.13 Why Chirality?

Why does nature distinguish left from right?

### Mass Terms Are Relevant

A Dirac mass term connects left and right chiralities:

$$m\bar{\psi}\psi = m(\bar{\psi}_L\psi_R + \bar{\psi}_R\psi_L)$$

If both chiralities exist in conjugate representations, this term is allowed. Under the renormalization group, it's a "relevant" deformation. It grows at low energies.

### Refinement Stability

Relevant operators that aren't forbidden by symmetry or constraints get turned
on under refinement. They can't be kept at zero without fine tuning.

If a mass term is allowed, it will generically appear. The fermion will become massive. At low energies, it will decouple.

To keep fermions light without fine tuning, the mass term must be forbidden. The cleanest way: make the fermion chiral. If only one chirality exists, there's no partner to couple to. No mass term is possible.

The Standard Model fermions are chiral for that reason. Chirality protects their masses from running to the cutoff scale.

## 14.14 What Particles Are in This Model

Before discussing which particles the model predicts, we need to be clear about what a "particle" even means in our approach. The answer is both more precise and more radical than the intuitive picture shows.

In the conventional view, particles are fundamental objects, tiny balls of
stuff that move through space. Fields fill the gaps, and particles are what
detectors click on. This picture is useful for calculations, but it gets the
ontology backwards. Particles are patterns first. They are not primitives.

Think about what an observer actually sees. Each observer has a patch on the
holographic screen and a collection of allowed questions. When the answers
settle into a stable excitation that survives local time evolution, keeps its
identity across overlaps, and transforms in a repeatable way under the emergent
symmetries, the theory has found a candidate carrier pattern. Calling that
pattern a quantum particle requires the physical-Hilbert-space, spectrum, and
pole or asymptotic-state receipts stated below.

There is a finer record-level question underneath that sentence. A stable
sector can be transported across a chart, yet two detector clicks on opposite
sides of the chart boundary have to be linked as the same continuing track. OPH
treats that as a separate stitching problem. The chart, clock,
interface, charge transport, and refinement must leave one track clearly
preferred. If they do not, the theory should say "ambiguous" instead of forcing
a particle history.

In ordinary language, a particle is a recurring role in the screen data. A
quantum particle claim needs more: a physical Hilbert-space representation, a
positive-energy spectrum, and the appropriate stable pole or asymptotic state.
An electron has such a massive quantum role in the Standard Model. A photon has
the corresponding massless spin-one role only after the electromagnetic
quantum-pole receipt passes.

That picture has teeth. The model does not place representations on the stage
and then ask whether they fit. It reads candidate charge and carrier roles from
the way the algebra net closes on itself; actions and physical spectra decide
which of those roles propagate as particles.

### The Particle Structure In One Picture

The particle picture can be told as one continuous line. The framework first
reconstructs compact gauge data from sectors that can be transported without a
loop mismatch. The small-economy selection then picks the realized Standard
Model quotient, the hypercharge lattice, and the color-generation count. The
same structure identifies electromagnetic, color, and metric carrier roles.
Explicit action and phase receipts give their classical transverse or TT
modes; quantum particle poles require a stronger receipt. From there the
pixel fixed point organizes the electroweak sector, the Higgs/top relation, the
conditional running-quark support surface, the charged-lepton story, and a
target-informed weighted-cycle neutrino comparison candidate rejected by the
correlated NuFIT 6.1 profile. Its source kernel and charged basis are open, so
physical mixing, mass ordering, Majorana phases, and absolute masses are
unformed.
Hadrons add the strong-binding problem on top of that particle-level picture.

The sphere ladder from Chapter 3 is useful here only as a logic map. It says
seed, loop, screen, bulk. It does not say photon, gluon, graviton, hadron.
Those role labels come from the recovered Lorentz and gauge structure. The
unbroken electromagnetic direction, color directions, and metric tensor mode
become photon, gluon, and graviton particle labels only after their respective
quantum receipts pass. $W$ and $Z$ are massive weak carriers on the declared quantitative branch, the
Higgs is the scalar electroweak excitation, and hadrons are QCD composites.

### How the Concrete Particle Entries Arise

Stable patterns on the screen matter because they land on the particle entries a
physicist actually cares about. First comes the structural side. Chapter 15
supplies Lorentz kinematics, so stable excitations sort themselves by the usual
labels of mass, spin, and helicity. The realized gauge quotient, hypercharge
lattice, and generation-color counting supply the particle-side structure. Together
they decide which charged excitations can exist and how they transform.

Then comes the fixed-point closure. The same screen cell is read twice. From the
outside it is a small displacement from golden-ratio balance. From the inside it
is the electromagnetic observation step available to observers in the encoded
world. On the endpoint branch, matching those two descriptions displays
$\alpha^{-1}(0)=137.035999177(21)$ and
$P\simeq1.6309682094$.

The route is a chain of translations. Golden-ratio entropy balance fixes the
screen side. The boundary width sets the electromagnetic detuning. The source
map then carries that value through unification, running gauge couplings, the
electroweak anchor, and the electromagnetic channel used at low energy.
The value is forced because the same local pixel has to satisfy the outside
geometry equation and the inside electromagnetic endpoint equation at the same
time.

The fine-structure constant belongs here beside the weak sector. It is the
local electromagnetic width of the observer-supporting pixel. From there the
same construction continues into the weak sector, the Higgs/top relation, the
conditional running-quark support sector. The weighted-cycle neutrino entry is
a target-informed comparison candidate rejected by the correlated NuFIT 6.1
profile, with no physical prediction status. Hadrons come later, because protons
and mesons are bound states of quarks and gluons. Their masses live in the
strong-binding problem, not in the bare quark table.

For that reason, a laboratory observer does not measure the undressed source
diagnostic as the fine-structure constant. A low-energy electromagnetic
measurement sees the dressed Thomson current after lepton vacuum polarization,
confined-quark/hadron spectral transport, and finite endpoint matching. The
mixed root-plus-comparison-pixel diagnostic is $137.0359595008\ldots$, just
below the public endpoint. The CODATA-minus-diagnostic gap names the numerical
burden on the same-scheme hadronic endpoint transport; it does not calculate
that transport from source data.

The screen cell wants to sit at the golden-ratio balance point, the exact
self-similar equilibrium of the local entropy hierarchy. A universe with
observers cannot remain perfectly silent. It needs a
small displacement so that records can form, photons can carry information, and
measurements can leave durable traces. The fine-structure lane reads the
Thomson-limit value as that small displacement.

## 14.15 What the Electromagnetic Branch Supplies

When two observer patches share a charged region, they may use different local
descriptions without changing the shared data. The recovered charge
bookkeeping closes on an unbroken $U(1)$ factor. That result identifies the
electromagnetic symmetry and connection role. A group label alone does not yet
provide a propagating field or a photon.

The next step is an explicit Maxwell branch. If the low-energy action contains
the usual positive $F^2$ kinetic term, and the selected vacuum has no Higgs,
Stueckelberg, medium, or nonlocal mass operator, gauge reduction leaves two
transverse classical waves. Their quadratic Green function has a pole at
$\omega^2=c_*^2|\mathbf k|^2$. This is a precise massless classical
carrier-mode statement.

A photon is the quantum interpretation of that mode. To earn that word from
the OPH construction, one still needs a positive-energy quantization, a
physical Hilbert space, a positive-residue two-point-function pole, and a stable
asymptotic state. The present group and Maxwell-equation derivation does not by
itself supply those objects.

## 14.16 What the Gravitational Branch Supplies

Chapter 15 explains how modular screen geometry leads to a classical Einstein
branch. On a flat background, the additional pure Einstein--Hilbert action can
be linearized and gauge-reduced. The result has two transverse-traceless
classical wave modes, conventionally called the plus and cross polarizations,
with the same invariant null speed $c_*$.

That calculation is not yet a graviton construction. A classical Einstein
equation does not quantize the metric, define a graviton Hilbert space, or prove
a positive-residue particle pole. Diffeomorphism invariance also does not rule
out every massive or additional mode: extra fields, higher-derivative terms,
bimetric sectors, and Stueckelberg descriptions can preserve a gauge
redundancy while changing the physical spectrum.

## 14.17 Why This Matters: Comparison to String Theory

String theory provides a useful contrast. After the worldsheet theory is
quantized, its physical spectrum can contain a massless spin-two state. The
state, its norm, and its pole belong to the same quantum construction.

The present OPH core reaches a different checkpoint. It reconstructs a gauge
classification, an electromagnetic action branch, and a classical Einstein
branch. Those results support classical transverse and spin-two waves under
their stated action and phase assumptions. A later string or field-theory lift
may supply the missing quantization and particle pole, but it must do so
explicitly. The particle interpretation is not automatic from overlap
redundancy alone.

## 14.18 Why Composite Masses Are Different

Consider the proton. Its mass is 938.272 MeV, measured to extraordinary
precision. Can OPH compute it from the same quadratic carrier analysis?

No. Even the carrier analysis is conditional on an action, background, phase,
and, for a particle, a quantum pole receipt. The proton is harder still. It is
a bound state, and bound
states ask for the full nonperturbative drama of quarks, gluons, and
confinement.

That difference matters. Some results in the framework are structural and
sharp. Others depend on solving the strong-coupling machinery in detail. The
electroweak sector sits close enough to the local fixed-point readout that masses and
couplings can be pinned down cleanly. Hadrons sit deeper in the strong-coupling problem.

A promising route into that jungle uses edge entanglement. It does not
weight charge sectors arbitrarily. It assigns each one a local geometric cost
set by the gauge group itself. Read those costs carefully enough and the
effective gauge couplings can be inferred from the vacuum.

In simple test cases such as $\mathbb Z_5$ and $S_3$, that weighting pattern
shows up with tight numerical accuracy. Even the golden-ratio fingerprint of
$\mathbb Z_5$ appears where the group geometry says it should. Entanglement
geometry leaves visible marks on the coupling structure.

The same golden-ratio motif returns on the fixed-point side. Perfect
self-similar balance would sit exactly at $\phi$. A lived universe with durable
records sits nearby, carrying the slight detuning that makes structure and
history possible. Reliable extraction of gauge couplings from entanglement
therefore sharpens the quantitative picture without breaking the narrative
spine of the chapter.

## 14.19 Gauge Unification and the Proton

One of the great puzzles of particle physics is why the three gauge couplings (for the strong, weak, and electromagnetic forces) have such different strengths at low energies, yet seem to converge when extrapolated to high energies.

In the 1970s, physicists noticed a numerical tease. If you run the couplings
upward using the renormalization group equations, they almost meet at a single
point around $10^{16}$ GeV. This suggested that all three forces might unify at
high energies, the dream of Grand Unified Theories.

The snag was immediate. With just the Standard Model particle content, the
three couplings do not quite meet. They miss each other. In the 1990s,
physicists discovered that adding supersymmetric partners fixes this: with
MSSM-like particle content, the couplings unify beautifully, predicting
$\alpha_s(M_Z) \approx 0.117$, close to the measured value of
$0.1177 \pm 0.0009$.

OPH separates two ideas that are often fused together. Couplings can display
unification-like running without the Standard Model being embedded in a larger
simple group. A heat kernel is a standard way of weighting group
representations with a diffusion-like smoothing parameter. In the edge-mode
construction, that weighting reproduces MSSM-like one-loop running: entropy
weights a representation by one copy of its dimension because one side of the
entanglement cut is traced over, while loop corrections see both indices of the
representation block. A second factor of the dimension returns in the running.
That is what lets the beta-function shifts land near the familiar unification
benchmark.

At the unification-scale heat-kernel parameter $t_U \approx 1.64$, this gives:

$$\Delta b_{\text{edge}} \approx (2.49,\ 4.38,\ 3.97)$$

compared to the MSSM target $(2.50,\ 4.17,\ 4.00)$. The agreement is within 5%
for all three coefficients in this edge-mode picture. What emerges here is
unification-like running behavior, not an MSSM spectrum hidden inside OPH.

MSSM means the Minimal Supersymmetric Standard Model, a popular extension of the
Standard Model. OPH adds no MSSM particle spectrum here. It compares the
running pattern of the couplings.

The sharper structural prediction concerns *how* any unification-like closure would happen.

### Why Protons Don't Decay

Traditional Grand Unified Theories achieve unification by embedding the Standard Model gauge group into a larger simple group like SU(5) or SO(10). This embedding has a dramatic consequence: it introduces new gauge bosons called X and Y bosons that can turn quarks into leptons. Protons should decay, with minimal SU(5) predicting lifetimes around $10^{31}$ years.

But Super-Kamiokande has been watching for proton decay since 1996. The experimental limit is
$\tau_p > 10^{34}$ years, a thousand times longer than predicted. The simplest GUTs are dead.

OPH takes a different path. The gauge group is not embedded in anything
larger. On the bosonic branch carrying the explicit compact-gauge refinement
receipt, Tannaka-Krein
reconstruction builds the gauge group from the persistent charge-sector data,
yielding the product structure:

$$G_{\mathrm{phys}} = \mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1) / \mathbb{Z}_6$$

There is no larger group. No X and Y bosons. No leptoquark generators. Any
coupling closure happens geometrically, with all three couplings sharing a
common "diffusion time" on the edge, not algebraically through group
embedding.

The prediction is stark: **gauge-mediated proton decay is forbidden**.

This is one of the cleanest experimental forks in the road. A simple-group
unification scheme predicts new gauge bosons that mediate proton decay. The OPH
route predicts that those bosons never exist. The difference is a detector
signal in those channels versus a structural null result.

The claim is unusually valuable. Many high-energy ideas differ
mainly in elegance or ultraviolet taste. Proton decay is harsher. Either the
detector sees the relevant channel or it keeps not seeing it. OPH lands on the
null-decay side for structural reasons.

This is a unique experimental signature. Standard SUSY GUTs predict both
precision unification and proton decay. OPH separates those questions:
the full connected gauge group has only the product-group adjoint content and
no mixed leptoquark generators, so gauge-mediated proton decay is forbidden,
while the edge-mode construction can display MSSM-like unification-style
running without simple-group embedding. Hyper-Kamiokande null results in those
channels, combined with precision measurements that favor unified couplings,
support geometric edge-sector running over algebraic unification.

## 14.20 What the Model Explains

The framework accounts for several concrete facts.

**The integers.** Why three colors? Why three generations? Why those specific hypercharges? These are consequences of consistency requirements, not free parameters. Anomaly cancellation and Yukawa invariance fix the hypercharge lattice, the minimal coupled carrier fixes the color triplet, and CKM CP capability together with weak-sector ultraviolet consistency fixes the generation count.

**The carrier modes.** On the declared Maxwell vacuum, the quadratic
electromagnetic action has two transverse massless classical modes. On the
declared pure-Einstein flat branch, the quadratic metric action has two
transverse-traceless massless classical modes. These statements concern the
displayed actions. They do not follow from symmetry alone and do not yet
construct photon or graviton particles. The experimental photon-mass and
gravitational-wave dispersion bounds become tests of OPH only after the
corresponding quantum spectral receipt is supplied.

**The particle structure.** Section 14.14 gives the concrete structure. The
framework fixes the carrier roles and states the conditional classical-mode
receipts. It then organizes the fine-structure
fixed point, the weak-boson benchmark pair, a Higgs/top relation, a theorem-grade
two-modulus quark source-spread obstruction with target-audit mass textures, and a target-informed
weighted-cycle neutrino candidate rejected by the correlated NuFIT 6.1 profile.
The charged-lepton basis is open and has nearly degenerate singular values, so
the physical PMNS matrix is unformed. Neutrino mass ordering, Majorana phases,
and absolute masses are unresolved. Charged leptons, direct-top comparison,
hadrons, and strong CP sit around that same map with different levels of support.

The reason these numbers belong in one chapter is that the framework organizes
them with one local fixed-point structure. The same pixel ratio feeds the
electroweak scale, the low-energy electromagnetic endpoint, and the effective
gravitational coupling. The point does not require every intermediate symbol.
OPH ties
electroweak relations, the Higgs/top quantitative relation,
electromagnetism at low energy, and Newton's constant into one common
structure.

The hierarchy relation is part of that same local story. The weak scale is emitted
by the pixel transmutation law:

$$
v/E_\star=P_\star^{-1/2}\exp[-2\pi/(4\alpha_U(P_\star))].
$$

Here $v$ is the weak scale, $E_\star$ is the branch energy scale used for the
display, and $\alpha_U(P_\star)$ is the unified coupling read from the pixel
branch. The source branch gives a unique value for that coupling. The selected
source-to-Higgs coarse-graining square has
zero defect, so the usual bare-Higgs-plus-cutoff split is a regulator
coordinate description on that branch. The SI gravity decimal is a separate
clock-scale readout.

There is also a global capacity-side register story. OPH treats reality as a
continual repair of disagreements between local observers. In the hierarchy
story, a local physics pixel has to synchronize with the cosmic
record-capacity screen. The screen side supplies twelve irreducible curvature
ports. Reversible write/check orientation turns those ports into a 24-slot
oriented repair register. The weak hierarchy uses this 24-fold normalization to
read a local piece of the global screen depth. This is the book-level picture
behind the local/global hierarchy link: particle physics samples the local
screen readout of the same closure that fixes the whole horizon record.

The 24 slots have a simple group-theory origin. The visible product branch has
eight color channels, three weak channels, and one hypercharge channel.
Reversible repair reads each channel in two orientations, write and verify.
Twelve one-way channels therefore become a 24-slot oriented repair register.
The SU(5) adjoint lands on the same single-orientation number for a different
support; its X/Y mixed gauge channels are absent from the OPH product branch.

The exact bridge capacity is about $3.53235\times10^{122}$ on the public
endpoint branch. The rounded $3.31\times10^{122}$ value is the de Sitter
capacity-scale display. The screen-sieve theorem explains why the weak-scale
exponent is divided by twelve: the finite spherical screen exposes twelve
icosahedral sampling ports, the refinement-stable orbit of its unit fivefold
curvature defects.

At the formal level, this bridge value is also a fixed point in log-capacity
coordinates. The physical point is simpler: the exact bridge value belongs to
the repair rhythm, while the rounded value is the de Sitter capacity-scale
display.

The readback side has a matching finite statement. At a finite cutoff, the
closed observer record selects one public capacity atom, so the delivery
resolution is a single value. No menu of incompatible scales remains. In the
refinement limit, that single atom becomes the horizon-to-cell factor used by
the repair tick.

**Charge quantization.** All color-singlet particles have integer electric charge. No fractional charges like $\pm 1/3$ can exist outside hadrons. This follows from the global structure of the gauge group.

**Gauge-mediated proton decay.** Gauge-mediated proton decay is forbidden. The gauge group is a product group with no embedding in a larger simple group, so no leptoquark generators exist. Published experimental limits ($\tau_p > 10^{34}$ years) are consistent with this prediction.

**Why hadrons are harder.** Quark masses are short-distance parameters.
Hadrons are bound states. Their masses come from the nonperturbative dynamics
of confined quarks and gluons. The OPH hadron story therefore has to pass
through strong binding, the same hard physics that makes most of the proton's
mass come from confinement rather than from the bare quark masses.

## 14.21 The Big Picture

The Standard Model looks like the answer to a very specific question: What is
the simplest admissible low-energy sector package that can be transported by
OPH's gluing rules, reconstructed as compact gauge data, and survive under
refinement?

The reconstruction fixes symmetry and geometric roles, not quantum particles.
On an additional Maxwell-action branch, gauge reduction leaves two transverse
classical null modes. On an additional pure Einstein--Hilbert branch about a
suitable flat background, constraint reduction leaves two
transverse-traceless classical null modes. A photon or graviton claim needs
more: a positive-energy physical Hilbert space, a positive-residue pole in the
physical spectrum, and the relevant stability or asymptotic-state conditions.
Gauge or diffeomorphism redundancy by itself does not supply those receipts,
nor does it exclude Higgs, Stückelberg, medium, higher-derivative, or
extra-field effects outside the declared branch.

The quark and lepton charges are fixed on the minimally selected
gauge-consistency branch. Three colors and three generations follow there from
the combined demands of anomaly cancellation,
chirality, the minimal coupled color carrier, CKM CP capability, and
ultraviolet consistency.

The result is concrete. The book points to a specific gauge
structure, charge pattern, color count, and generation count. It also reaches
conditional action-level carrier modes, the $W$ and $Z$ benchmark pair, a
Higgs/top quantitative relation, and a quark source-spread non-identifiability
theorem on a chosen physical basis. MAR gives every positive rescaling of the
two quark profiles the same structural score, so it does not select their
Yukawa eigenvalues. Numeric quark rows remain withheld, and the
mixed-convention GeV matrices are target-audit mass textures rather than
physical Yukawas. The weighted-cycle neutrino point belongs to the comparison ledger:
NuFIT 6.1 rejects it. The charged basis is open and has nearly degenerate
singular values, leaving physical mixing, mass ordering, and absolute masses
unresolved. Strongly coupled bound states add the QCD problem on top of that
particle-level structure.

Particles emerge from the screen as stable patterns that transform under emergent symmetries. The natural sequel is spacetime itself. If the particle inventory is fixed by consistency, can geometry be fixed the same way?

That's the question of **Chapter 15: Relativity from Modular Time**.
