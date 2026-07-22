# Chapter 7: The Recovery Rule

## 7.1 The Intuitive Picture: Information Can Be Copied Freely or Lost Forever

Start with the ordinary picture of information loss.

Information can be freely copied or irreversibly destroyed. When you write a
letter, you can make as many copies as you like. When you burn a book, the
information is gone forever. These are two distinct fates: duplication or
annihilation.

This is the commonsense view embedded in our everyday experience. You can photocopy a document infinitely. You can record a conversation and play it back endlessly. Information is cheap to replicate. When the Library of Alexandria burned, the information vanished into the void; a crashed hard drive and a fading memory go the same way. Destruction is final.

Classical physics supported this intuition. The state of a system is a point in phase space. You can, in principle, measure it exactly and write down as many copies as you wish. And entropy increases, meaning organized information degrades into random noise. The past becomes inaccessible as the universe forgets.

Physics broke that picture from both directions.

## 7.2 The Surprising Hint: No-Cloning, Yet Information Cannot Be Destroyed

### The No-Cloning Theorem

The first shock came from quantum mechanics. In 1982, William Wootters and Wojciech Zurek proved the **no-cloning theorem**: there is no quantum operation that can copy an unknown quantum state.

If you have a qubit in state |psi> and want to create |psi>|psi>, you cannot. The linearity of quantum mechanics forbids it.

The obstruction is a fundamental law rather than a shortfall in our tools.
Quantum information cannot be copied. You cannot make a backup of a quantum
state. You cannot read it out and write it elsewhere without disturbing the
original.

This seems catastrophic for building reliable systems. Classical computers work precisely because we can make redundant copies. If one bit flips, the backup catches it. How can you protect information you cannot copy?

### The Black Hole Information Paradox

The second shock came from black holes, and it pointed in the opposite direction.

In 1974, Stephen Hawking made a disturbing discovery. Black holes emit faint
radiation due to quantum effects near the event horizon. This **Hawking
radiation** has a precise temperature:

$$T = \frac{\hbar c^3}{8\pi G M k_B}$$

For a solar-mass black hole, this is about 60 nanokelvin, undetectably cold. But for small black holes, the temperature can be significant. The radiation carries energy away. Black holes evaporate.

$T$ is the Hawking temperature. $M$ is the black-hole mass. The constants
$\hbar$, $c$, $G$, and $k_B$ are Planck's constant divided by $2\pi$, the
speed of light, Newton's gravitational constant, and Boltzmann's constant.
Because $M$ is in the denominator, smaller black holes are hotter.

The problem was severe. Hawking's calculation showed the radiation is thermal: random, uncorrelated noise carrying no information about what fell in. If you throw a book into a black hole and wait for evaporation, all you get out is random static.

If this is true, information is destroyed. A pure quantum state (the book) becomes a mixed thermal state (the radiation). This violates **unitarity**, the foundational principle that quantum evolution preserves information.

Hawking was willing to accept this. Most other physicists were not.

### A Holographic Resolution Perspective

After decades of debate, the broad holographic lesson is that black-hole evaporation need not destroy information. In semiclassical holographic models, the Hawking radiation is not truly random: it carries subtle correlations, so information that looked lost can instead be encoded in the radiation.

This lesson was sharpened by the Page-curve and island calculations developed in the 2010s. In semiclassical holographic models, they support encoded-information viewpoints and show how information that seemed lost to the black hole interior can instead be carried by correlations among the outgoing radiation particles.

Information cannot be copied (no-cloning), yet information cannot be destroyed (unitarity). These twin constraints require a specific structure: **quantum error correction**.

## 7.3 The First-Principles Reframing: Error Correction Structure Preserves Information

The deeper question is why nature forbids copying and yet refuses to lose
information.

### The Library of Alexandria Revisited

In 48 BC, Julius Caesar's troops set fire to the Egyptian fleet in Alexandria's
harbor. The flames spread to warehouses, then to buildings, and according to
legend consumed the Great Library, the ancient world's greatest repository of
knowledge. Hundreds of thousands of scrolls burned, and ash
drifted over the Mediterranean.

We intuitively understand this loss is permanent. Once a book is burned, the information is destroyed. Entropy increases, smoke disperses, and time ensures we cannot run the movie backward.

But is the information *really* gone?

This question haunted Ludwig Boltzmann in the 1870s. His colleague Josef
Loschmidt pointed out something troubling: the fundamental laws of physics are
reversible. Newton's equations run equally well forward or backward. If you
knew the exact position and momentum of every molecule of smoke and ash,
including every atom that had been paper and ink, you could in principle
reverse their trajectories and reconstruct the scrolls.

The information isn't destroyed. It's scrambled. Hidden in correlations among billions of particles, diluted into the environment until no practical measurement could extract it. But mathematically, physically, it remains there.

### The Universe's Error Correction

**The universe appears to use error-correcting structure that preserves information even when it appears lost.**

In quantum mechanics, this requirement is non-negotiable. Closed-system quantum
evolution is **unitary**. If information were genuinely destroyed in that
setting, the standard quantum-mechanical evolution law would fail.

So the universe must preserve information, even when it looks scrambled beyond
recognition. Some mechanism, a "Save Game" feature, must in principle allow the
smoke to remember what the scroll said.

But how can information be preserved if it cannot be copied? Perfect copies were never required. It is enough to encode the information **redundantly**, in a way that survives local errors.

## 7.4 Claude Shannon's Discovery

The recovery thread begins in 1948, in an office at Bell Telephone Laboratories, at the desk of Claude Shannon.

His colleagues worried about static on phone lines and how to squeeze more calls onto a cable. Shannon was asking what information *is*: whether it can be measured, and how a message can survive a channel that tries to destroy it.

Shannon had spent World War II working on cryptography, trying to make messages secure from eavesdroppers. He then attacked the opposite problem: how to make messages survive noise that corrupts them randomly.

His 1948 paper, "A Mathematical Theory of Communication," is one of the most
influential scientific works of the twentieth century. It founded information
theory, and buried in its pages was the recovery idea that matters here.

### The Noisy Channel

Imagine you're sending a message through a bad phone line. You say "yes," but static might make it sound like "mess" or "ness." How can you guarantee your message gets through?

Shannon's answer was that you can't eliminate noise, but you can beat it with **redundancy**.

A simple example is repetition coding. Send a single bit three times. A zero
becomes `000`. A one becomes `111`.

Suppose noise flips one bit. You receive "010." Majority vote says the original was "0," with two zeros versus one one. The information survives.

This seems obvious, although Shannon proved something surprising: every noisy
channel has a **capacity**, a maximum rate at which information can be sent
reliably. Below that capacity, codes exist whose error rate can be made
arbitrarily small.

The trick is clever encoding. Spread information across many symbols in subtle patterns. The receiver can reconstruct the original even when individual symbols are corrupted, because the patterns survive even when specific symbols don't.

### The Cost of Reliability

Redundancy isn't free. Extra symbols mean slower transmission. Extra bits mean more storage. And there's a fundamental cost: Landauer's principle says erasing a bit requires at least kT ln 2 of energy, about 3 times 10 to the negative 21 joules at room temperature.

The universe has finite resources. Recovery must be efficient, local, bounded. You can't store infinite backups of infinite data.

This constraint shapes reality. The area law says a boundary can only carry so many bits. If information capacity is bounded by area, then recovery must respect geometry. Distant regions can't share unlimited redundancy.

**Spacetime can be read through a Shannon-code analogy.** Gravity then acts
like an error corrector, keeping the global account consistent even when local
observations are noisy.

## 7.5 The Mathematics of Redundancy

The mathematics starts with ordinary information.

### Shannon Entropy

Shannon defined the information content of a random variable X with outcomes {x} and probabilities {p(x)}:

$$H(X) = -\sum_x p(x) \log p(x)$$

This measures uncertainty: how many yes/no questions you'd need to ask, on average, to learn the outcome.

$X$ is the random variable, $x$ labels one possible outcome, and $p(x)$ is the
probability of that outcome. The sum adds the uncertainty contribution from
each possible outcome.

Examples make the meaning concrete. A fair coin has $H=1$ bit, one yes-or-no
question. A heavily loaded coin at 99% heads has about $0.08$ bits, because
there is very little uncertainty left. A certain outcome has $H=0$.

### Mutual Information: Shared Predictive Content

The mutual information between X and Y measures how much knowing one tells you about the other:

$$I(X:Y) = H(X) - H(X|Y) = H(X) + H(Y) - H(X,Y)$$

If X and Y are independent, I(X:Y) = 0: knowing one tells you nothing about the other. If they're perfectly correlated, mutual information equals entropy: knowing one determines the other.

$H(X|Y)$ means the uncertainty left about $X$ after $Y$ is known. $H(X,Y)$ is
the joint entropy of the pair. Mutual information is the amount of uncertainty
that disappears when one variable predicts the other.

### Conditional Mutual Information: The Recovery Metric

Recovery enters through conditional mutual information, CMI for short, which measures correlation between X and Y *given* knowledge of Z:

$$I(X:Y|Z) = H(X|Z) + H(Y|Z) - H(X,Y|Z)$$

If I(X:Y|Z) = 0, then X and Y are **independent given Z**. Once you know Z, learning Y tells you nothing new about X.

The vertical bar again means "given." CMI asks how
much extra connection remains between $X$ and $Y$ after the mediator $Z$ is
supplied.

This is the mathematical definition of "Z screens X from Y." All information that Y has about X is contained in Z.

Small CMI means approximate given-data independence, and approximate given-data independence enables recovery.

## 7.6 Markov Chains and Screening

We say X goes to Y goes to Z forms a **Markov chain** if X and Z are independent given Y:

$$p(x,z|y) = p(x|y) \cdot p(z|y)$$

This is equivalent to I(X:Z|Y) = 0.

The equation says that once $Y$ is known, the joint probability for $X$ and
$Z$ factors into two separate probabilities. In ordinary language, $Y$ carries
all the information that connects the two ends.

### The Screening Property

When $X$ leads to $Y$ leads to $Z$, we say $Y$ screens off $X$ from $Z$. Once
you know $Y$, $X$ adds nothing new about $Z$. All $X$-$Z$ correlation is
mediated through $Y$. The middle system carries everything relevant.

With $Y$ in hand, $X$ can be discarded without losing what $X$ could have told
you about $Z$.

### Physical Examples

Consider three locations along a copper wire: A, B, C, with B in the middle. In thermal equilibrium, B's temperature screens A from C. Heat from A reaches C only through B. If you know B's temperature precisely, knowing A's temperature adds nothing to your prediction of C's.

This is **locality**. Effects propagate through space. Distant regions communicate only through intermediates.

Your skin is a Markov blanket. It screens your internal organs from the external world. Everything the world knows about your liver, it knows through your skin (and other body surfaces). Everything your liver knows about the world, it knows through your skin.

An observer's patch works the same way. It carries all accessible information
about what lies beyond. In the ideal recovery limit, the patch is a sufficient
summary of what can be recovered from outside.

## 7.7 Quantum Recovery: The Petz Map

### From Classical to Quantum

Everything we've discussed has quantum analogs.

For a quantum state described by density matrix rho, the von Neumann entropy is:

$$S(\rho) = -\text{Tr}(\rho \log \rho) = -\sum_i \lambda_i \log \lambda_i$$

where the lambdas are the eigenvalues of rho.

This is Shannon entropy with quantum bookkeeping. The eigenvalues
$\lambda_i$ are the weights of the independent quantum alternatives after the
density matrix is diagonalized. The trace expression is the coordinate-free way
to compute the same uncertainty without choosing a favorite basis.

The quantum CMI is the same formula as before, now built from quantum entropies:

$$I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC)$$

### Strong Subadditivity: The Miracle Theorem

In 1973, Elliott Lieb and Mary Beth Ruskai proved one of the most important theorems in quantum information:

**Strong Subadditivity**: For any quantum state, I(A:C|B) is greater than or equal to 0.

CMI is never negative.

The statement sounds obvious. The proof took years of hard functional analysis, and it is the foundation of quantum recovery.

Strong subadditivity says B can only help, never hurt: at worst it is useless, and it never creates confusion that was absent before.

### The Petz Map: Physical Recovery

In 1986, the Hungarian mathematician Dénes Petz asked a natural question: if I(A:C|B) = 0 exactly, can we physically reconstruct the state?

Petz constructed an explicit procedure, later called the **Petz recovery map**:

$$R_{B \to BC}(\sigma) = \rho_{BC}^{1/2} (\rho_B^{-1/2} \sigma \rho_B^{-1/2} \otimes I_C) \rho_{BC}^{1/2}$$

This formula is shown because it is the repair operation made explicit. It is
the quantum version of saying: take the state available on $B$, rebalance
it using the reference correlations, then rebuild the missing $C$ side.

In the formula, $\sigma$ is the state you actually have on $B$, while
$\rho_{BC}$ is the reference correlation pattern telling the map how $B$ and
$C$ fit together. The square roots and inverse square roots are matrix
operations that rebalance the known state before rebuilding the missing side.

The formula's details are secondary for the main story. This is a physical
operation, in principle something you could implement on a quantum device. Given
only B's state sigma, the Petz map outputs a comparison state on BC that
reproduces the reference correlations in the exact Markov case.

Think of it like calibrating a distorted photograph. The original image (BC) got scrambled into a noisy version (B alone). The Petz map knows what the original "should" look like (from the reference state rho_BC) and applies the inverse distortion.

### Approximate Recovery: The Fawzi-Renner Theorem

Perfect recovery requires I(A:C|B) = 0 exactly. But in physics, nothing is exact. What if CMI is merely small?

In 2015, Omar Fawzi and Renato Renner proved the theorem that turns recovery into a working tool:

**Theorem**: For any state rho_ABC with I(A:C|B) less than or equal to epsilon, there exists a recovery map R such that:

$$\|\rho_{ABC} - (\mathbb{I}_A \otimes R_{B \to BC})(\rho_{AB})\|_1 \leq 2\sqrt{2\epsilon}$$

Small CMI implies approximate recoverability. The smaller I(A:C|B), the better the recovery.

The norm bars measure how distinguishable the original state and the recovered
state are. $\mathbb I_A$ means that subsystem $A$ is left alone while the
recovery map acts on $B$. The small number $\epsilon$ is the allowed CMI error.
The theorem turns a small information leak into a concrete reconstruction
bound.

This is the mathematical heart of the recovery rule: **redundancy implies reconstruction**.
Where you draw the cut matters as much as how much leaks across it.

## 7.8 Example Calculations

The recovery rule becomes clearer in a toy state.

### A Bell Pair Plus Extra Qubit

Let A and B be entangled in a Bell state, and let C be an independent qubit.

Since C is independent, knowing B tells you everything B could possibly tell you about C, which is nothing. So I(A:C|B) = 0 exactly. B screens A from C perfectly.

Recovery is trivial here: C has nothing to do with A, so "recovering" C from B just means C can be anything.

### The GHZ State: Maximum Correlation

The GHZ state is different:

$$|\text{GHZ}\rangle = \frac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$$

Compute $I(A:C|B)$ for this state.

For a pure state |psi> of ABC, we have S(ABC) = 0 (pure states have zero entropy).

The reduced state on AB is:
$$\rho_{AB} = \frac{1}{2}(|00\rangle\langle00| + |11\rangle\langle11|)$$

This is a classical mixture, not entangled. Its entropy S(AB) = 1 bit.

Similarly, S(BC) = 1 bit and S(B) = 1 bit.

So:
$$I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC) = 1 + 1 - 1 - 0 = 1$$

The GHZ state has nonzero, genuinely tripartite CMI. B does not screen A from C
at all. The correlation between A and C is genuinely tripartite: you need all
three systems to see it.

You cannot recover C from B alone. The GHZ state is non-Markov.

## 7.9 The Fourth Axiom: Local Markov/Recoverability

We can state the recovery rule as a physical principle. One piece of
vocabulary first: a **collar** is a thin buffer strip wrapped around a region,
like a moat around a castle.

**Axiom 4 (Local Markov/Recoverability)**: For any three support regions $P_A$,
$P_B$, and $P_C$ on the screen chart, where $P_B$ topologically separates
$P_A$ from $P_C$:

$$I(A:C|B) \leq \varepsilon(B)$$

$A$, $B$, and $C$ label the subsystems carried by those support regions. $B$
is the separator, typically a collar, and the small quantity $\varepsilon(B)$
is the allowed leakage past it. Exact Markov screening sets it to zero. In the
finite thermal screen model, $\varepsilon(B)$ is a boundary-size prefactor
times a factor that decays exponentially with the collar's thickness.

The framework's repair tool is the Fawzi-Renner recovery map from earlier in
this chapter. As the collar thickens, the repair error shrinks toward zero,
and the recovered state approaches the exact zero-leak form that the ideal
splitting-and-recombining identities assume.

### Two Routes to Zero Leakage

Sometimes the cut lines up perfectly. When the energy stored at the boundary
depends only on a conserved label shared by both sides, the natural way to cut
the system matches the way the state itself divides, and the leakage across
the boundary is exactly zero.

Otherwise, thickening the buffer beats the leakage. The concrete picture is
the one in the figure below: a cap $A$ on the screen, a collar $B$ wrapped
around it, and the exterior $D$, everything outside the moat. For a
well-behaved thermal state whose boundary correlations mix uniformly as one
moves inward, the leakage falls off exponentially with the collar's thickness,

$$I(A:D\,|\,B)\ \le\ c\,|\partial B|\,e^{-\delta/\xi},$$

where $\delta$ is the thickness of the collar, $\xi$ is a fixed decay length,
and $|\partial B|$ counts the cells along the collar's boundary. A thicker
collar drives the leakage down, while a longer boundary raises the prefactor
in front. The collar has to thicken faster than the boundary grows.

### Screening Through the Separator

If region B sits between regions A and C, then B approximately screens A from C. The correlations between A and C are almost entirely mediated through B.

The “almost” is quantified by $\varepsilon(B)$. Collar thickness and boundary
size pull in opposite directions: thickness improves screening, while a longer
boundary supplies more places for a residual correlation to cross.

![A cap, a thin collar, and the exterior form the A-B-D split used in recoverability arguments.](../assets/book_diagrams/collar-tripartition.svg){width=70%}

### Constructive Gluing (Tree Covers)

In the finite-dimensional code-subspace setting, the arena of the
error-correcting codes below, Axiom 4 yields a clean constructive result for
patch families arranged like a tree. Each new patch touches the already-glued
union only through a single separator $B$, the induced $A$-$B$-$C$ split is a
genuine tensor product at each step, and recovery maps glue the patches into
a global state.

The reconstruction error per step is bounded by

$$\|\rho_{ABC} - (\mathrm{id}\otimes\mathcal R)(\rho_{AB})\|_1 \le 2\sqrt{\ln 2\; I(A:C|B)}$$

(CMI in bits), and errors accumulate at most additively (capped by 2).

Loopy covers add a global closure condition. If several overlaps wrap around
and return to the starting point, their transformations must close cleanly on
the full loop. Pairwise gluing handles the local joins; loop closure removes a
global defect. In a chiral effective field theory, the same consistency burden
reappears as anomaly cancellation.

This matches holographic expectations. In AdS/CFT, entanglement between
boundary regions scales with the area of the minimal surface connecting them.
Here, recoverability bounds scale with separator size rather than bulk volume.

### Why This Matters

The recovery rule does a surprising amount of work. If the interior of a region can be
recovered from its boundary, bulk physics is encoded in boundary physics. If
$I(A:C|B)$ is small, then $A$ and $C$ behave independently once $B$ is known,
which is exactly the operational face of locality. Ground states of local
Hamiltonians often show area-law entanglement, and recoverability helps explain
why correlations can remain organized in a boundary-sensitive way. Classical
facts become the records that survive because they are redundantly encoded.

## 7.10 The Black Hole Information Paradox Reframed

The recovery perspective sharpens one of physics' most famous puzzles.

### Hawking's Puzzle, Revisited

The black-hole puzzle is an extreme version of the recovery question. A book
falls past the horizon. Long afterward the black hole has almost evaporated.
If the outgoing radiation is only heat, the book's information seems gone. If
quantum evolution preserves information, the radiation must carry correlations
that a coarse thermal description misses.

### The Page Curve

In 1993, Don Page proposed a resolution. If information is preserved, the entropy of Hawking radiation should follow a specific curve.

Early on, radiation entropy increases. Each photon emitted is uncorrelated with previous photons.

At the **Page time**, roughly when the black hole has lost half its mass,
radiation entropy should start *decreasing*. Later photons become correlated
with earlier ones. The radiation starts "remembering" what fell in.

Page's curve is the shape unitarity demands: entropy rises until Page time,
falls after it, and returns to zero for a final pure state.

The Page curve long stood as a consistency requirement for unitarity, not as a direct calculation.

### The Recovery Perspective

The recovery rule gives a natural OPH reading of holographic interior encoding.

Label the systems in the usual way. $A$ is the information thrown into the
black hole, Alice's diary. $B$ is the early Hawking radiation. $C$ is the late
Hawking radiation.

Initially, B is small. The collected radiation is too small to decode the diary information carried jointly by the full evaporation state.

As time passes, B grows. More radiation is emitted, and the correlations needed for decoding become increasingly accessible in the radiation subsystem.

At Page time, B becomes large enough to screen A from C. The CMI I(A:C|B) is
then expected to drop, and later radiation becomes recoverable from earlier
radiation once the separator grows large enough to do its screening work.

### Islands

In 2019, several groups (Penington; Almheiri, Engelhardt, Marolf, and Maxfield) made this precise using a concept called "islands."

When computing entropy in theories with gravity, you should include contributions from **island regions** inside the black hole.

Before Page time, no island contributes. Radiation entropy tracks the naive Hawking calculation and keeps increasing.

After Page time, an island appears. The black-hole interior, called the
**island**, is encoded in the radiation. Including the island contribution,
radiation entropy decreases.

The island formula reproduces the Page curve in semiclassical holographic
models and makes the encoding picture vivid. Alice's diary may be physically
inside the black hole, yet its information does not need to live in an
autonomous interior tensor factor. The black-hole lesson is that recovery and
encoding belong to the basic architecture.

## 7.11 Spacetime as Error Correction

The black hole resolution points to a deeper truth: spacetime may have the structure of a quantum error-correcting code.

### Quantum Error Correction

In quantum computing, you can't copy quantum information (no-cloning theorem). So how do you protect qubits from noise?

**Quantum error correction** spreads information across many physical qubits in
entangled configurations. If some qubits are corrupted, the others can
reconstruct the original.

The simplest example is the three-qubit code. Logical $|0\rangle$ becomes
$|000\rangle$, and logical $|1\rangle$ becomes $|111\rangle$.

If one qubit flips, majority vote recovers the original. This is just classical repetition. Quantum codes are more sophisticated, protecting against both bit-flips and phase errors.

### The HaPPY Code

In 2015, Fernando Pastawski, Beni Yoshida, Daniel Harlow, and John Preskill
built a toy model of holography using error correction, the **HaPPY code**.
The acronym is assembled from their surnames, which is roughly what passes
for wordplay in quantum gravity.

They constructed a tensor network in which the bulk is the logical
information and the boundary is made of the physical qubits.

Information in the bulk is redundantly encoded in the boundary. Erase part of
the boundary and bulk information survives in a form recoverable from the
remaining boundary.

This is exactly the recovery rule: I(Bulk : Erased | Remaining) is approximately 0.

The "gravity" in the HaPPY code emerges from the code structure. Regions of the bulk are closer when they share more boundary support. Distance becomes a property of information, not something fundamental.

---

## 7.12 The Indestructible Past

The recovery rule has a startling implication: in this recoverability picture,
nothing is simply deleted from the full quantum description.

If the universe is unitary and holographic encoding is stable, information is
redistributed into increasingly nonlocal correlations of the full quantum state.

The Library of Alexandria? The scrolls burned, and the information scrambled into smoke, heat, and light. That radiation spread across the cosmos at light speed. It is diluted across an unimaginably vast region of space. In principle, with a computer the size of the observable universe, you could run the Petz map and watch the smoke reconstitute into Sophocles. You would have to wait a while for the answer, and then check whether you still remembered the question.

Paleontology and astronomy use weak versions of this. Fossils preserve information about creatures from millions of years ago. Astronomy records light that has traveled for billions of years before reaching our telescopes. The cosmic microwave background is one vivid example of very old information preserved in radiation.

The recovery rule says this is not accident or luck. In a unitary encoding
picture, the past is carried forward in increasingly scrambled form.

### The Structural Constraint

Practical recovery is impossible. The computation required to
recover the Library of Alexandria would exceed any conceivable technology.
Chaos amplifies tiny errors. A single misplaced bit in trillions grows into
garbage.

The distinction matters. The past is recoverable in principle but
inaccessible in practice. That gives us both unitarity and the lived arrow of
time. The past is encrypted with a key we will never find. Patterns, it turns
out, are harder to destroy than the things that carry them.

## 7.13 Reverse Engineering Summary

Information can remain recoverable without being freely copied. No-cloning
blocks duplication. Recovery survives because the information is encoded across
extended correlations, and the mathematics underwrites the claim: strong
subadditivity guarantees that CMI cannot go negative, and Fawzi-Renner and
Petz show that when the missing correlation is small enough, there is a map
that rebuilds what looked lost.

The physics mirrors the mathematics. Black-hole evaporation is read through
the Page curve. Entanglement wedges, the boundary regions that carry
everything about a given piece of the bulk, reconstruct interior data from
edge data. Quantum error correction works in the lab, which means encoded
recovery has operational teeth. That is how a noisy world carries history, how
observers can agree on a past they never saw, and why black holes do not
behave like cosmic shredders. It is also why spacetime starts to look like a
code: a structure whose geometry and stability are tied to the same redundancy
that protects information.

Shannon started with a practical problem, sending messages over noisy phone
lines. His solution, redundancy, reappears as one of the strongest organizing
analogies for spacetime and holographic encoding.

---

We have the Screen. We have the Algebra. We have the Consistency Rules. We have Recovery.

Where does space come from? Where does time come from? How does the abstract
structure of quantum information become the geometry we inhabit?

The next chapters turn recovery into geometry. We'll see how boundaries encode interiors, how entanglement draws the map, and how the consistency conditions we've developed start to look suspiciously like gravity.
