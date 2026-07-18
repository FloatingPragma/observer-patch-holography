# Chapter 13: The de Sitter Patch

## 13.1 The Intuitive Picture: The Universe Is Static or Decelerating

Start with the old cosmological picture.

The universe is either static, with things staying roughly as they are, or
decelerating, with gravity pulling everything together and slowing expansion.
This is the natural expectation from Newton through Einstein.

Einstein himself added a "cosmological constant" to his equations in 1917 to create a static universe, a universe that neither expanded nor contracted. When Hubble discovered the universe is expanding, Einstein dropped the constant, calling it his "greatest blunder."

Even after accepting expansion, the expectation was deceleration. Gravity
attracts. The mutual pull of all the matter in the universe should slow the
expansion, like a ball thrown upward gradually slowing. The expansion might
stop or even reverse.

Supernova data broke that picture.

## 13.2 The Surprising Hint: The Universe Is Accelerating

### The 1998 Supernova Observations

In January 1998, two teams of astronomers independently announced results that overturned our understanding of the cosmos.

Saul Perlmutter led the Supernova Cosmology Project. Brian Schmidt and Adam
Riess led the High-Z Supernova Search Team. Both groups had spent years hunting
Type Ia supernovae, the "standard candles" of cosmology.

Everyone expected to find that expansion is slowing. The data showed the opposite.

Distant supernovae were fainter than expected, farther away than a decelerating
universe would predict. The universe is **speeding up**.

Something is pushing the cosmos apart. Something is fighting gravity and winning. The teams called it "dark energy."

The supernova result also rested on historical work. Henrietta Leavitt's study of
Cepheid variable stars gave astronomy a way to climb the cosmic distance
ladder. Edwin Hubble's expansion law made the universe dynamical. Walter Baade,
Allan Sandage, Vera Rubin, Kent Ford, and many others sharpened the large-scale
picture long before the 1998 teams found acceleration. Cosmology is a relay
race. The de Sitter clue entered this book through a century of measurements,
calibrations, and arguments about what the sky was actually saying.

### The Cosmological Constant Returns

In the simplest late-time fit, a positive cosmological constant, Lambda greater
than zero, creates a repulsive large-scale tendency that grows with distance.
In the conventional late-time reconstruction, matter once dominated the
large-scale stress budget. As the effective expansion diluted matter, the
dark-energy component took over.

The standard fit places the transition to accelerated expansion about 5 billion
years back along our reconstructed branch and summarizes the cosmic budget as
about 68% dark energy.

If that component is a true cosmological constant, the future approaches a de
Sitter phase. The DESI-era data keep dynamical dark energy on the table, so this
chapter treats de Sitter as the standard branch used for the OPH capacity
calculation.

## 13.3 The First-Principles Reframing: De Sitter Is the Natural Screen

The deeper question is why late-time cosmology points toward a de Sitter patch.

### The Static Patch

What does one observer actually experience in de Sitter space?

As you look outward, galaxies recede faster and faster. At a critical distance $r_H = c/H$, the recession velocity equals the speed of light. Beyond this radius, light can never reach you.

Here $H$ is the Hubble expansion rate for the de Sitter phase. The formula says
that expansion itself creates a distance beyond which signals cannot overcome
the stretching of space.

This defines your **cosmological horizon**, the boundary of your causal access.

![A de Sitter static patch gives one observer a finite horizon screen whose capacity scales with area.](../assets/book_diagrams/desitter-static-patch.svg){width=78%}

Inside the horizon, you can use static coordinates. This region, the
**static patch**, is all of de Sitter space that you can ever access.

### De Sitter Fits OPH

**The de Sitter horizon is the natural holographic screen.**

The fit is tight. Observers have finite patches, and the static patch
is bounded by a horizon. The patch boundary is an $S^2$, exactly the geometry
the framework wants. The entropy is finite through the Gibbons-Hawking area
law. No observer sees beyond the horizon, so there is no God's-eye view.
Observer equivalence is built in because de Sitter is maximally symmetric.
Time is patch-dependent because there is no preferred global clock.

The static patch is the natural arena for physics from an observer's perspective.

## 13.4 The Gibbons-Hawking Temperature

In 1977, Gary Gibbons and Stephen Hawking proved that the cosmological horizon
radiates like a black body:

$$T_{dS} = \frac{\hbar H}{2\pi k_B}$$

For our universe, this is about 10^{-30} Kelvin, far below direct detection.
In conventional inflationary cosmology, horizon-scale quantum fluctuations are
stretched and later seed structure formation. OPH assigns that job to finite
observer-screen synchronization: screen geometry and source release supply
candidate large-scale coherence, and ordinary Boltzmann transfer carries finite
source data to the sky.

The symbols echo earlier horizon physics. $T_{dS}$ is the de Sitter
temperature. $\hbar$ is Planck's constant divided by $2\pi$. $H$ is the Hubble
rate for the de Sitter phase. $k_B$ is Boltzmann's constant, which converts
energy units into temperature units. The denominator $2\pi$ is the same
circle factor that appears in Unruh and Hawking horizon temperatures.

Empty space is not glowing brightly around us. An observer confined to one
static patch sees the horizon as a thermal environment. Part of the quantum
state is inaccessible beyond the horizon, and that loss of access has the same
thermodynamic signature that horizons have elsewhere in gravitational physics.

### Why This Temperature? The Unruh Connection

The Gibbons-Hawking and Unruh formulas are closely related, but the identification has to be stated carefully.

A geodesic observer at the center of the static patch has zero proper acceleration, while a generic observer held at fixed radius has a radius-dependent proper acceleration. Near a horizon, the local de Sitter temperature reduces to the corresponding Unruh form:

$$T_U = \frac{\hbar a}{2\pi c k_B}$$

So the de Sitter and Unruh temperatures are locally linked, but they should not be identified by assigning every static-patch observer the same acceleration $a = cH$.

This has an important implication for OPH: **de Sitter horizons satisfy the same thermodynamic relations as Rindler horizons**. This is standard Gibbons-Hawking thermodynamics.

### Finite Entropy

If the horizon has temperature, it must have entropy:

$$S_{dS} = \frac{A}{4\ell_P^2} = \frac{\pi c^5}{G\hbar H^2}$$

This is the entropy associated with one de Sitter static patch, the logarithm
of the effective number of states accessible within that patch.

Here $A$ is the horizon area, $\ell_P$ is the Planck length, $c$ is the speed
of light, $G$ is Newton's gravitational constant, and $H$ again sets the
de Sitter expansion rate. The first expression is the area law. The second is
the same law after writing the horizon radius in terms of $H$.

For the late-time horizon of our universe, $R_{dS} \approx 1.66 \times 10^{26}$ m. The bare radius-squared count is

$$N_{\text{patch}} = \left(\frac{R_{dS}}{\ell_P}\right)^2 \approx 1.05 \times 10^{122}.$$

The entropy capacity includes the area factor:

$$N_{\text{scr}} = S_{dS} = \pi N_{\text{patch}} \approx 3.31 \times 10^{122},$$

or about $4.77 \times 10^{122}$ bits.

That is the practical meaning of the formula. It is a capacity statement. The patch does not contain an infinite amount of information hidden in a smooth continuum. It contains a finite number of distinguishable states, and the area of the horizon tells you how large that state space can be.

This finite entropy has major implications. An observer's accessible patch has a
finite information capacity. The smooth continuum starts to look like an
effective description laid over a screen with a hard budget.

There is one naming trap. The cap-normal theorem in the papers also uses the
phrase "de Sitter," but there it names a bookkeeping space for oriented round
caps and circles. It is not the cosmological de Sitter spacetime whose static
patch and entropy are discussed in this chapter.

### Why This Matters for Gravity

Jacobson's derivation of Einstein's equations requires horizons with a
temperature proportional to surface gravity, an entropy proportional to area,
and a first law tying heat to entropy. De Sitter thermodynamics supplies that
structure. In OPH it becomes the natural thermodynamic backdrop for the
gravity relation recovered from observer-patch consistency.

## 13.5 The Problem of Time in De Sitter

In Anti-de Sitter space, there's a boundary at spatial infinity that provides a universal time reference.

De Sitter has no spatial boundary. The only boundary is the horizon-and the horizon is observer-dependent.

### Horizon Complementarity

Leonard Susskind and collaborators proposed **de Sitter complementarity**: operationally, physics should be described patch by patch, without privileging a single global observer description.

Alice describes physics in her patch using her Hilbert space. Bob describes physics in his patch using his Hilbert space. Where their patches overlap, their descriptions must be consistent. In the complementarity reading adopted here, patch-relative descriptions are primary.

A Hilbert space here is the quantum state space for the degrees of freedom
accessible inside one observer's horizon, not a private mental space.

This fits naturally with OPH. Reality is a collection of consistent patches. You can't step outside and view the universe from nowhere.

## 13.6 Static Patch Holography

Where should we put the holographic screen in de Sitter?

A natural candidate: on the cosmological horizon.

For an observer at $r = 0$, the horizon is a sphere at $r = c/H$. This sphere has area $4\pi c^2/H^2$ and the entropy capacity above.

The three-dimensional bulk inside the horizon is treated holographically as data organized on the two-dimensional horizon.

When an object falls toward the horizon, it gets redshifted and appears to freeze onto the surface, its information smeared across the screen.

The horizon is the natural screen for cosmology. It is the last place where
an observer can trade signals with the rest of the patch. If physics is
organized around what observers can compare, then the cosmological horizon is
exactly where that comparison structure has to live.

### Why This Is Not dS/CFT

When physicists say "de Sitter holography is unsolved," they typically mean: no AdS/CFT-like duality with a clean boundary CFT at infinity is available. The classic dS/CFT proposal puts a Euclidean CFT at future infinity. This leads to notorious problems: potential non-unitarity, complex weights, and no clear operational access for any observer.

OPH takes a different path. It uses static-patch holography with
positive Lambda. The boundary is the observer's horizon, not future infinity.
The construction asks for local algebras and overlap consistency, without one
global CFT. Each observer has a horizon screen, and observer-dependence is
part of the setup.

This is a different target. The "unsolved problem" of dS holography is about finding a global boundary theory at infinity. OPH asks how local observer patches, each bounded by a horizon, yield consistent physics.

### The Lambda-Capacity Relation

A key point: the cosmological constant belongs to global capacity rather than
local patch data. Null modular probes reconstruct the stress tensor only up to
a term proportional to the metric itself, so $\Lambda g_{ab}$ enters as the one
global scale the local construction cannot erase.

The symbol $\Lambda$ is the cosmological constant, the part of Einstein's
equation that acts like a uniform large-scale tendency for space to accelerate.
It belongs to the proposed global capacity branch.

On that branch the dimensionless Lambda-capacity relation is posed as a
**global** self-reading constraint: the outside total horizon capacity must
equal the inside observer-accessible public record. With the selected OPH scale
bridge, the target is displayed as

$$N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}}),\qquad
\Lambda_{\mathrm{CRC}}=\frac{3\pi}{G N_{\mathrm{CRC}}}.$$

The direction of inference runs from measurement to working value. Through the
scale bridge, the measured $\Lambda$ locates the working capacity
$N_{\mathrm{CRC}}$,
and that located value enters the program as an input, never as one of its
outputs. The
dimensionless capacity relation is the global size parameter carried by every
consistent patch if the readback closure is derived.

The global closure has a sharp target. Let $F(N)$ denote the active horizon
capacity read back by stable observers inside an OPH universe supplied with
capacity $N$. A source-derived construction of $F$ is absent. The required
cosmic record-closure equation is

$$
N_{\mathrm{CRC}}=F(N_{\mathrm{CRC}}).
$$

The proposed finite-count representation uses terminal observer normal forms
that close on themselves at capacity $N$, normalized by the full screen
Hilbert-space size:

$$
\Pi(N)=|\Omega^{\mathrm{sc}}_N|e^{-N},
\qquad
N_\star=\operatorname*{argmax}_{N\ \mathrm{admissible}}
\left(\log|\Omega^{\mathrm{sc}}_N|-N\right).
$$

$\Omega^{\mathrm{sc}}_N$ is the set of self-consistent observer-supported
screen records at capacity $N$. The absolute-value bars count how many such
records there are. The exponential factor $e^{-N}$ removes the trivial growth
from increasing the raw screen capacity. The selected capacity $N_\star$ is the
admissible value where this balance is largest.

The same proposal can be written as a fixed-point map. Define
$\ell(N)=\log|\Omega^{\mathrm{sc}}_N|-N$. The declared formal map
$T_\eta(N)=N+\eta\ell'(N)$ moves $N$ in the direction indicated by the slope of
that balance function. Under the stated derivative-sign condition this formal
map has a unique stable fixed point. OPH has not derived the required count or
shown that this formal map is the physical observer readback. Identifying its
fixed point with the de Sitter entropy capacity requires three open premises,
including the counting statement.

The image is the important part. The small positive Lambda is proposed as the
size label of a finite self-reading horizon. It is the number every observer
patch would inherit if the universe reconstructs its own boundary from inside.
In the book’s chain of consistency requirements this is the capacity self-reading
step: a world that reads itself has to read back its own size, without deficit
or slack.

This is the informal $N_{\mathrm{CRC}}$ story. From the outside the datum is
total horizon capacity. From the inside the datum is total observer-accessible
public record. Closure says that the universe must be able to reconstruct its
own boundary. Observers exist inside it, and infer geometry, horizons, entropy,
$\Lambda$, history, and records from the information available inside the
universe.

The capacity relation belongs to global screen size: a larger
self-consistent record capacity corresponds to a smaller positive cosmological
constant once the selected scale bridge is used for the SI display.

### Many Observers, One Lambda

The philosophical stance of OPH, no objective camera angle and only perspectives that must agree on overlaps, maps naturally onto de Sitter static-patch intuition. Each timelike observer has a horizon and a patch. There is no operational access to a single global description.

On that same proposed branch, Lambda is the global quantity that
**can** be shared across overlaps. It is a capacity constraint that all
consistent overlapping descriptions inherit. Different observers see
different patches, and they all see the same Lambda encoded in the finite
size of their horizons.

### The Cosmology Picture

The cosmology picture is easiest to state in plain language. When the
entropy-maximizing state is rotationally symmetric for an observer, the
large-scale stress tensor looks like a perfect fluid. When the same isotropy
holds across observers, the spatial slices have constant curvature. Combined
with the gravity relation from the earlier chapters and a positive
cosmological constant, that gives the familiar FLRW geometry used in
cosmology.

### What "Early Universe" Means Here

In OPH, the Big Bang names an effective boundary in a reconstruction made from
records inside finite observer patches, rather than a God's-eye first instant
at which spacetime objectively begins. The CMB, light element abundances,
galaxy surveys, supernovae, and lensing data are real records. The global
story built from them is a powerful model of our branch, with no view from
outside the universe smuggled into the account.

So the "early universe" in this book means the hot dense side of the
reconstructed record: the side where the effective plasma is opaque, radiation
dominates the bookkeeping, and the CMB surface is not transparent. It names a
record-side limit rather than an absolute origin of being. Likewise, "cosmic
time" is an effective clock that appears when many local records synchronize,
an internal clock rather than an external stopwatch measuring spacetime from
outside.

Inflation is where this distinction becomes operational. Standard cosmology
introduces an inflaton field to explain flatness, horizon coherence, and the
nearly scale-invariant pattern in the CMB. OPH assigns those jobs to finite
observer-screen synchronization.
Flatness belongs to the selected spatial branch, horizon coherence to shared
boundary repair, the near scale-invariant pattern to screen geometry and
release, and the acoustic peaks to ordinary Boltzmann evolution once finite
source data are fixed. The CMB and growth curves are empirical tests of that
source story.

The screen part of this proposal has a precise conditional form. The scalar
field is one third of the logarithmic volume change on a uniform-density
release surface, with the average and dipole removed. A source-side maximum
entropy law fixes its amplitude from the expected release energy. The screen
precision then gives an exact angular spectrum through a ratio of gamma
functions.

The selected dimensionless pixel ratio is denoted by $P_*$. If a finite source
supplies the required infinitesimal edge-center reserve density, the scalar
tilt is

$$n_s=1-\frac{P_*}{48}=0.9660214956\ldots$$

A single spherical snapshot loses depth information. Even a perfect sky map can
compare points only across distances up to the sphere's diameter. Many different
three-dimensional patterns can therefore cast the same pattern onto one screen.

OPH has two ways to recover the missing depth. The first follows the screen as
it is refined and checks that this change really corresponds to zooming through
physical wavelengths. When that check passes, the scaling law fixes the shape
of the three-dimensional spectrum and the screen fixes its amplitude. The
second combines correlations from many radii, much like medical tomography
combines many views to reconstruct an interior.

The mathematics of these two routes is complete. A finite source carrying the
reserve density, release law, stress, clock, freezeout surface, physical scale,
and one of the radial certificates is work in progress. Temperature and
polarization spectra also require ordinary Boltzmann transfer and a complete
likelihood. The decimal above is a conditional source value. It has no status
as a fitted CMB result or physical prediction.

## 13.7 Scrambling and Chaos

De Sitter space is a **fast scrambler**-perhaps the fastest possible.

Information sent toward the horizon gets thermalized, mixed with all the other quantum information. The scrambling time is:

$$t_{scrambling} \sim \frac{1}{H}\ln S \sim \frac{280}{H}$$

For our universe, this is about 4 trillion years. Black holes are the standard saturators of the chaos bound in holographic settings, and de Sitter is often discussed as a fast-scrambling horizon with analogous scaling.

$t_{scrambling}$ is the time needed for initially localized information to
become thoroughly mixed across the horizon degrees of freedom. The symbol
$\sim$ means "scales like," not exact equality. $S$ is the de Sitter entropy.
The number 280 comes from the logarithm of the huge entropy associated with
our late-time horizon.

The smooth, empty appearance of the de Sitter vacuum can be read as highly scrambled information in this perspective.

## 13.8 The Swampland and Anthropic Selection

String theory has difficulty producing stable de Sitter vacua.

Swampland arguments suggest that stable de Sitter vacua may be impossible in consistent quantum gravity. If true, our universe would be slowly rolling down a potential hill.

Even if de Sitter vacua exist, why is Lambda so small (10^{-122} in Planck units)?

The **anthropic principle** offers an answer: if Lambda were much larger,
galaxies couldn't form. If it were negative, the universe would recollapse. In
the simplest fit, we find ourselves in a universe with small positive Lambda
because that is where observers can exist.

## 13.9 Reverse Engineering Summary

Historical cosmology expected expansion to slow under gravity. The sky disagreed.
The supernova data and the standard late-time fit point toward de Sitter
behavior, and de Sitter fits the observer-first picture with almost suspicious
neatness. Each observer has a static patch, a horizon, a temperature, an
entropy budget, and finite accessible information. The cosmological horizon is
the natural screen in this reading.

## 13.10 The Dark Sector: Static Law and Physical Gaps

In 1933, Fritz Zwicky looked at galaxies in the Coma Cluster and found a
problem. They were moving as if the cluster contained far more mass than the
telescopes could see. Four decades later, Vera Rubin and Kent Ford found the
same kind of problem inside spiral galaxies. The outer stars were orbiting too
fast. If the visible stars and gas were all the mass there was, those outskirts
should have behaved differently.

The modern evidence is broader than rotation curves. Galaxy clusters bend light
more strongly than their visible matter allows. The Bullet Cluster separates
hot gas from the gravitational lensing map. The cosmic microwave background
carries the imprint of extra gravitational weight on the hot dense side of the
standard reconstruction. This is what physicists call the dark sector: the part
of the cosmic budget inferred from gravity rather than from light.

The standard answer is dark matter: new particles that clump around galaxies,
pull on stars and light, and mostly ignore electromagnetism. That answer works
well in large-scale cosmology, but the particles themselves have not shown up
in detectors. MOND takes the opposite route and changes the low-acceleration
law of gravity. It captures a surprising amount of galaxy phenomenology, while
clusters and precision cosmology remain difficult terrain.

OPH reads the same evidence through the observer screen. A galaxy is luminous
matter, dust, gas, and dark gravitational behavior inside a finite de Sitter
patch. The question becomes: what does gravity read when the boundary
bookkeeping is nearly, but not perfectly, additive?

### The Modular Leftover

Chapter 11 used the branch-qualified entanglement-equilibrium route to recover
the local Einstein relation. In the ideal case, the geometric readout,
null-stress bridge, small-ball area identity, and tensor upgrade close cleanly
enough that local gravity takes the familiar form.

The recovery analysis has two versions. In the clean case, the boundary lines up
with the way the state divides and the leakage across it is exactly zero. In the
broader case, the leakage falls off exponentially with the thickness of the
boundary layer, provided the boundary loses its correlations fast enough as you
move inward. Ordinary fast decay of two-point correlations does not by itself
establish that. A thicker layer suppresses the leak, and a larger boundary
raises the prefactor.

A nonzero leakage can therefore serve as a finite diagnostic of imperfect
recovery. That leakage is a single number with no direction, so it cannot be
inserted into Einstein's equation as a source of gravity. A gravitational source of this kind has to
be a proper conserved tensor, reconstructed with directions built in and passing
the usual normalization and coupling tests. The finite quantities a simulation
can measure are stand-ins for that reconstruction, and they do not amount to a
full dark source on their own.

### The Repair-Charge Medium

The proposed dark-sector completion treats the scalar repair register like a
rotor. The integer $n$ counts signed repair occupation. The angle $\theta$
tracks its conjugate phase. Together they form the number-and-phase pair needed
for dynamics. Finite OPH supplies the scalar register, while this pair and its
action are proposed additions.

Given that action, changing the phase gives an exact balance law for repair
charge, and changing the geometry gives a stress tensor. In a dilute,
homogeneous phase, charge conservation gives $\rho_R\propto a^{-3}$ and
$w_R\simeq0$. The medium therefore dilutes like pressureless matter. Its cosmic
abundance and perturbations are work in progress.

### Why Galaxies Flatten

In the spherical deep-gradient regime, a cubic repair-phase link energy gives

$$a_R=\sqrt{a_0g_b}$$

where $g_b$ is the Newtonian acceleration from baryons. For a galaxy this gives
the flat-rotation-curve behavior astronomers actually see.

$g_{\text{obs}}$ is inferred from the rotation curve, and $g_b$ is the
Newtonian acceleration of stars, gas, and dust. The same asymptote gives the
baryonic Tully-Fisher scaling

$$V^4 = GM_ba_0.$$

These are conditional consequences of the proposed action. The OPH derivation
of the canonical pair, source constants, complete constitutive law,
relativistic limit, lensing, clusters, Solar-System response, abundance, and
perturbations is work in progress. Cosmology is outside the falsification
program.

The same action lets coherent material enter the repair channel through a
fixed dimensionless efficiency. That number is not a force law. A device force
also needs a dimensional calibration and either an external repair gradient,
net repair charge, or momentum carried away by the repair field. A closed,
repair-neutral device cannot support itself through an internal coherence
contrast.

---

The arena is a finite static patch bounded by a holographic horizon. What populates this arena? What are the particles and forces we observe, and why do they have the peculiar properties they do?

The next chapter treats the Standard Model of particle physics as an effective
structure. The path has three steps. Gluing observer patches together sorts out
which charge sectors can be carried around cleanly. Refining that picture lets a
reconstruction read a compact symmetry group from the surviving sectors.
Choosing the smallest single-Higgs matter content that works lands on the
realized Standard Model.

This is **Chapter 14: The Standard Model from Consistency**.
