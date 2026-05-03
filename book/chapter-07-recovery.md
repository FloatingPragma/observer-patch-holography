# Chapter 7: The Recovery Rule

## 7.1 The Intuitive Picture: Information Can Be Copied Freely or Lost Forever

Start with the ordinary picture of information loss.

Information can be freely copied or irreversibly destroyed. When you write a
letter, you can make as many copies as you like. When you burn a book, the
information is gone forever. These are two distinct fates: duplication or
annihilation.

This is the commonsense view embedded in our everyday experience. You can photocopy a document infinitely. You can record a conversation and play it back endlessly. Information is cheap to replicate. Conversely, when the Library of Alexandria burned, when a hard drive crashes, when memories fade with age, the information vanishes into the void. Destruction is final.

Classical physics supported this intuition. The state of a system is a point in phase space. You can, in principle, measure it exactly and write down as many copies as you wish. And entropy increases, meaning organized information degrades into random noise. The past becomes inaccessible as the universe forgets.

Physics broke that picture from both directions.

## 7.2 The Surprising Hint: No-Cloning, Yet Information Cannot Be Destroyed

### The No-Cloning Theorem

The first shock came from quantum mechanics. In 1982, William Wootters and Wojciech Zurek proved the **no-cloning theorem**: there is no quantum operation that can copy an unknown quantum state.

If you have a qubit in state |psi> and want to create |psi>|psi>, you cannot. The linearity of quantum mechanics forbids it.

This does not reflect any shortfall in our tools. It is a fundamental law.
Quantum information cannot be copied. You cannot make a backup of a quantum
state. You cannot read it out and write it elsewhere without disturbing the
original.

This seems catastrophic for building reliable systems. Classical computers work precisely because we can make redundant copies. If one bit flips, the backup catches it. How can you protect information you cannot copy?

### The Black Hole Information Paradox

The second shock came from black holes-and pointed in the opposite direction.

In 1974, Stephen Hawking made a disturbing discovery. Black holes aren't quite black-they emit faint radiation due to quantum effects near the event horizon. This **Hawking radiation** has a precise temperature:

$$T = \frac{\hbar c^3}{8\pi G M k_B}$$

For a solar-mass black hole, this is about 60 nanokelvin-undetectably cold. But for small black holes, the temperature can be significant. The radiation carries energy away. Black holes evaporate.

Here's the problem. Hawking's calculation showed the radiation is thermal-random, uncorrelated noise carrying no information about what fell in. If you throw a book into a black hole and wait for evaporation, all you get out is random static.

If this is true, information is destroyed. A pure quantum state (the book) becomes a mixed thermal state (the radiation). This violates **unitarity**-the foundational principle that quantum evolution preserves information.

Hawking was willing to accept this. Most other physicists were not.

### A Holographic Resolution Perspective

After decades of debate, the broad holographic lesson is that black-hole evaporation need not destroy information. In semiclassical holographic models, the Hawking radiation is not truly random: it carries subtle correlations, so information that looked lost can instead be encoded in the radiation.

This lesson was sharpened by the Page-curve and island calculations developed in the 2010s. In semiclassical holographic models, they support encoded-information viewpoints and show how information that seemed lost to the black hole interior can instead be carried by correlations among the outgoing radiation particles.

Information cannot be copied (no-cloning), yet information cannot be destroyed (unitarity). These twin constraints require a specific structure: **quantum error correction**.

## 7.3 The First-Principles Reframing: Error Correction Structure Preserves Information

The deeper question is why nature forbids copying and yet refuses to lose
information.

### The Library of Alexandria Revisited

In 48 BC, Julius Caesar's troops set fire to the Egyptian fleet in Alexandria's harbor. The flames spread to warehouses, then to buildings, and according to legend, consumed the Great Library-the ancient world's greatest repository of knowledge. Hundreds of thousands of scrolls burned. Sophocles' lost plays, Aristotle's missing books, Euclid's unfinished theorems-gone. Ash drifted over the Mediterranean.

We intuitively understand this loss is permanent. Once a book is burned, the information is destroyed. Entropy increases, smoke disperses, and time ensures we cannot run the movie backward.

But is the information *really* gone?

This question haunted Ludwig Boltzmann in the 1870s. His colleague Josef Loschmidt pointed out something troubling: the fundamental laws of physics are reversible. Newton's equations run equally well forward or backward. If you knew the exact position and momentum of every molecule of smoke and ash-every atom that had been paper and ink-you could, in principle, reverse their trajectories and reconstruct the scrolls.

The information isn't destroyed. It's scrambled. Hidden in correlations among billions of particles, diluted into the environment until no practical measurement could extract it. But mathematically, physically, it's still there.

### The Universe's Error Correction

**The universe is built with error-correcting structure that preserves information even when it appears lost.**

In quantum mechanics, this requirement is non-negotiable. Quantum evolution is **unitary**-reversible by definition. If information could genuinely be destroyed, unitarity would break. If unitarity breaks, probabilities don't sum to 1. Physics collapses into nonsense.

So the universe must preserve information, even when it looks scrambled beyond recognition. There must be a mechanism-a "Save Game" feature-that allows, in principle, the smoke to remember what the scroll said.

But how can information be preserved if it cannot be copied? The answer: you don't need to copy information perfectly to protect it. You need to encode it **redundantly** in a way that survives local errors.

## 7.4 Claude Shannon's Discovery

The recovery thread begins in 1948, in a cramped office at Bell Telephone Laboratories in Murray Hill, New Jersey.

Claude Shannon was not like other engineers. While his colleagues worried about practical problems-how to reduce static on phone lines, how to compress calls onto cables-Shannon was thinking about something deeper. What *is* information? Can it be measured? How do you send a message reliably when the channel tries to destroy it?

Shannon had spent World War II working on cryptography, trying to make messages secure from eavesdroppers. He then attacked the opposite problem: how to make messages survive noise that corrupts them randomly.

His 1948 paper, "A Mathematical Theory of Communication," is one of the most
influential scientific works of the twentieth century. It founded information
theory, and buried in its pages was the recovery idea that matters here.

### The Noisy Channel

Imagine you're sending a message through a bad phone line. You say "yes," but static might make it sound like "mess" or "ness." How can you guarantee your message gets through?

Shannon's answer: you can't eliminate noise, but you can beat it with **redundancy**.

Here is the simplest example. Instead of sending a single bit, send it three
times. A zero becomes `000`. A one becomes `111`.

Suppose noise flips one bit. You receive "010." Majority vote says the original was "0"-two zeros versus one one. The information survives.

This seems obvious, but Shannon proved something surprising: every noisy channel has a **capacity**-a maximum rate at which you can send information reliably. If you send slower than capacity, there exist codes whose error rate can be made arbitrarily small.

The trick is clever encoding. Spread information across many symbols in subtle patterns. The receiver can reconstruct the original even when individual symbols are corrupted, because the patterns survive even when specific symbols don't.

### The Cost of Reliability

Redundancy isn't free. Extra symbols mean slower transmission. Extra bits mean more storage. And there's a fundamental cost: Landauer's principle says erasing a bit requires at least kT ln 2 of energy-about 3 times 10 to the negative 21 joules at room temperature.

The universe has finite resources. Recovery must be efficient, local, bounded. You can't store infinite backups of infinite data.

This constraint shapes reality. The area law says a boundary can only carry so many bits. If information capacity is bounded by area, then recovery must respect geometry. Distant regions can't share unlimited redundancy.

**Spacetime itself behaves like a Shannon code.** Gravity acts like an error corrector, keeping the global account consistent even when local observations are noisy.

## 7.5 The Mathematics of Redundancy

Let's build up the mathematics step by step.

### Shannon Entropy

Shannon defined the information content of a random variable X with outcomes {x} and probabilities {p(x)}:

$$H(X) = -\sum_x p(x) \log p(x)$$

This measures uncertainty-how many yes/no questions you'd need to ask, on average, to learn the outcome.

Examples make the meaning concrete. A fair coin has $H=1$ bit, one yes-or-no
question. A heavily loaded coin at 99% heads has about $0.08$ bits, because
there is very little uncertainty left. A certain outcome has $H=0$.

### Mutual Information: The Key Quantity

The mutual information between X and Y measures how much knowing one tells you about the other:

$$I(X:Y) = H(X) - H(X|Y) = H(X) + H(Y) - H(X,Y)$$

If X and Y are independent, I(X:Y) = 0-knowing one tells you nothing about the other. If they're perfectly correlated, mutual information equals entropy-knowing one determines the other.

### Conditional Mutual Information: The Recovery Metric

Here's where recovery comes in. The conditional mutual information measures correlation between X and Y *given* knowledge of Z:

$$I(X:Y|Z) = H(X|Z) + H(Y|Z) - H(X,Y|Z)$$

If I(X:Y|Z) = 0, then X and Y are **conditionally independent given Z**. Once you know Z, learning Y tells you nothing new about X.

This is the mathematical definition of "Z screens X from Y." All information that Y has about X is contained in Z.

Small conditional mutual information means approximate conditional independence-and approximate conditional independence enables recovery.

## 7.6 Markov Chains and Screening

We say X goes to Y goes to Z forms a **Markov chain** if X and Z are conditionally independent given Y:

$$p(x,z|y) = p(x|y) \cdot p(z|y)$$

This is equivalent to I(X:Z|Y) = 0.

### The Screening Property

When $X$ leads to $Y$ leads to $Z$, we say $Y$ screens off $X$ from $Z$. Once
you know $Y$, $X$ adds nothing new about $Z$. All $X$-$Z$ correlation is
mediated through $Y$. The middle system carries everything relevant.

This matters. It means you can throw away X and still have full access to anything X could have told you about Z-as long as you keep Y.

### Physical Examples

Consider three locations along a copper wire: A, B, C, with B in the middle. In thermal equilibrium, B's temperature screens A from C. Heat from A reaches C only through B. If you know B's temperature precisely, knowing A's temperature adds nothing to your prediction of C's.

This is **locality**. Effects propagate through space. Distant regions communicate only through intermediates.

Your skin is a Markov blanket. It screens your internal organs from the external world. Everything the world knows about your liver, it knows through your skin (and other body surfaces). Everything your liver knows about the world, it knows through your skin.

An observer's patch works the same way. It carries all accessible information about what lies beyond. In the ideal recovery limit, the patch isn't just a window-it is a sufficient summary.

## 7.7 Quantum Recovery: The Petz Map

### From Classical to Quantum

Everything we've discussed has quantum analogs.

For a quantum state described by density matrix rho, the von Neumann entropy is:

$$S(\rho) = -\text{Tr}(\rho \log \rho) = -\sum_i \lambda_i \log \lambda_i$$

where the lambdas are the eigenvalues of rho.

The quantum conditional mutual information is:

$$I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC)$$

### Strong Subadditivity: The Miracle Theorem

In 1973, Elliott Lieb and Mary Beth Ruskai proved one of the most important theorems in quantum information:

**Strong Subadditivity**: For any quantum state, I(A:C|B) is greater than or equal to 0.

Conditional mutual information is never negative.

This sounds obvious but it's not. The proof took years and required sophisticated functional analysis. And it's the foundation of quantum recovery.

Strong subadditivity says B can only help, never hurt. If you want to learn about correlations between A and C, knowing B cannot make things worse. In the worst case, B is useless. But B can never create confusion that didn't exist before.

### The Petz Map: Physical Recovery

In 1986, Hungarian mathematician Denes Petz asked a natural question: if I(A:C|B) = 0 exactly, can we physically reconstruct the state?

The answer is yes, and Petz constructed the explicit procedure later called the **Petz recovery map**:

$$R_{B \to BC}(\sigma) = \rho_{BC}^{1/2} (\rho_B^{-1/2} \sigma \rho_B^{-1/2} \otimes I_C) \rho_{BC}^{1/2}$$

In the formula, $\sigma$ is the state you actually have on $B$, while
$\rho_{BC}$ is the reference correlation pattern telling the map how $B$ and
$C$ fit together. The square roots and inverse square roots are matrix
operations that rebalance the known state before rebuilding the missing side.

Don't worry about the formula's details. This is a physical operation, something you could implement with a quantum computer. Given only B's state sigma, the Petz map outputs a state on BC that correctly reproduces all correlations with A.

Think of it like calibrating a distorted photograph. The original image (BC) got scrambled into a noisy version (B alone). The Petz map knows what the original "should" look like (from the reference state rho_BC) and applies the inverse distortion.

### Approximate Recovery: The Fawzi-Renner Theorem

Perfect recovery requires I(A:C|B) = 0 exactly. But in physics, nothing is exact. What if conditional mutual information is merely small?

In 2015, Omar Fawzi and Renato Renner proved a powerhouse theorem:

**Theorem**: For any state rho_ABC with I(A:C|B) less than or equal to epsilon, there exists a recovery map R such that:

$$\|\rho_{ABC} - (\mathbb{I}_A \otimes R_{B \to BC})(\rho_{AB})\|_1 \leq 2\sqrt{2\epsilon}$$

Small conditional mutual information implies approximate recoverability. The smaller I(A:C|B), the better the recovery.

This is the mathematical heart of the recovery rule: **redundancy implies reconstruction**.

## 7.8 Example Calculations

Let's see the recovery rule in action.

### A Bell Pair Plus Extra Qubit

Let A and B be entangled in a Bell state, and let C be an independent qubit.

Since C is independent, knowing B tells you everything B could possibly tell you about C-which is nothing. So I(A:C|B) = 0 exactly. B screens A from C perfectly.

Recovery is trivial here: C has nothing to do with A, so "recovering" C from B just means C can be anything.

### The GHZ State: Maximum Correlation

The GHZ state is different:

$$|\text{GHZ}\rangle = \frac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$$

Let's compute I(A:C|B).

For a pure state |psi> of ABC, we have S(ABC) = 0 (pure states have zero entropy).

The reduced state on AB is:
$$\rho_{AB} = \frac{1}{2}(|00\rangle\langle00| + |11\rangle\langle11|)$$

This is a classical mixture, not entangled. Its entropy S(AB) = 1 bit.

Similarly, S(BC) = 1 bit and S(B) = 1 bit.

So:
$$I(A:C|B) = S(AB) + S(BC) - S(B) - S(ABC) = 1 + 1 - 1 - 0 = 1$$

The GHZ state has nonzero, genuinely tripartite conditional mutual information. B doesn't screen A from C at all. The correlation between A and C is genuinely tripartite-you need all three systems to see it.

This means you can't recover C from B alone. The GHZ state is non-Markov.

## 7.9 The Fourth Axiom: Local Markov/Recoverability

We can state the recovery rule as a physical principle.

**Axiom 4 (Local Markov/Recoverability)**: For any three patches P_A, P_B, P_C on the screen, where P_B topologically separates P_A from P_C:

$$I(A:C|B) \leq \varepsilon(B)$$

Here $\varepsilon(B)$ measures how much correlation can bypass the separator.
Its form follows the geometry of the separator itself, with natural candidates
including boundary-size scaling or exponential decay with separation.

### Screening Through the Separator

If region B sits between regions A and C, then B approximately screens A from C. The correlations between A and C are almost entirely mediated through B.

The "almost" is quantified by ε(B). Larger separators allow more "leakage"-more correlation that bypasses the screen.

### Constructive Gluing (Tree Covers)

In the finite-dimensional code-subspace setting, Axiom 4 yields a clean
constructive result for tree-ordered covers. Each new patch overlaps the
existing glued union only on a single separator $B$. The induced $A$-$B$-$C$
split is a genuine tensor product at each step, and recovery maps glue the
patches into a global state.

The reconstruction error per step is bounded by

$$\|\rho_{ABC} - (\mathrm{id}\otimes\mathcal R)(\rho_{AB})\|_1 \le 2\sqrt{\ln 2\; I(A:C|B)}$$

(CMI in bits), and errors accumulate at most additively (capped by 2).

Loopy covers ask for one more check. If several overlaps wrap around and return
to the starting point, the gluing has to close cleanly on the full loop, not
just pair by pair. If it fails, the reconstruction accumulates a genuine
global defect. In a chiral effective field theory, the same consistency burden
reappears as anomaly cancellation, although the precise bridge is a later EFT
step.

This matches holographic expectations. In AdS/CFT, entanglement between
boundary regions scales with the area of the minimal surface connecting them.
Here, one natural scaling candidate ties recoverability bounds to separator
size and not to bulk volume.

### Why This Matters

The recovery rule has dramatic consequences. If the interior of a region can be
recovered from its boundary, bulk physics is encoded in boundary physics. If
$I(A:C|B)$ is small, then $A$ and $C$ behave independently once $B$ is known,
which is exactly the operational face of locality. Ground states of local
Hamiltonians tend toward area-law entanglement because recovery keeps the
correlations under control. Classical facts become the records that survive
because they are redundantly encoded.

## 7.10 The Black Hole Information Paradox Resolved

The recovery rule resolves one of physics' most famous puzzles.

### Hawking's Calculation

In 1974, Stephen Hawking made a disturbing discovery. Black holes aren't quite black-they emit faint radiation due to quantum effects near the event horizon.

Here's the problem. Hawking's calculation showed the radiation is thermal-random, uncorrelated noise carrying no information about what fell in. If you throw a book into a black hole and wait for evaporation, all you get out is random static.

If this is true, information is destroyed. A pure quantum state (the book) becomes a mixed thermal state (the radiation). This violates unitarity-the foundational principle that quantum evolution preserves information.

### The Page Curve

In 1993, Don Page proposed a resolution. If information is preserved, the entropy of Hawking radiation should follow a specific curve.

Early on, radiation entropy increases. Each photon emitted is uncorrelated with previous photons.

But at the **Page time**-roughly when the black hole has lost half its mass-something changes. Radiation entropy should start *decreasing*. Later photons become correlated with earlier ones. The radiation starts "remembering" what fell in.

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

At Page time, B becomes large enough to screen A from C effectively in the heuristic picture. The conditional mutual information I(A:C|B) is then expected to drop.

This motivates an encoded-information picture. Later radiation becomes
recoverable from earlier radiation once the separator grows large enough to do
its screening work.

### Islands: The Mathematical Proof

In 2019, several groups (Penington; Almheiri, Engelhardt, Marolf, and Maxfield) made this precise using a concept called "islands."

When computing entropy in theories with gravity, you should include contributions from **island regions** inside the black hole.

Before Page time, no island contributes. Radiation entropy equals naive Hawking calculation-increasing.

After Page time, an island appears. The interior of the black hole-the **island**-is encoded in the radiation. Including the island contribution, radiation entropy decreases.

The island formula reproduces the Page curve in semiclassical holographic
models and makes the encoding picture vivid. Alice's diary may be physically
inside the black hole, yet its information does not need to live in an
autonomous interior tensor factor. The black-hole lesson of this chapter is
that recovery and encoding belong to the basic architecture.

## 7.11 Spacetime as Error Correction

The black hole resolution points to a deeper truth: spacetime may have the structure of a quantum error-correcting code.

### Quantum Error Correction

In quantum computing, you can't copy quantum information (no-cloning theorem). So how do you protect qubits from noise?

The answer is **quantum error correction**: spread information across many physical qubits in entangled configurations. If some qubits are corrupted, the others can reconstruct the original.

The simplest example is the three-qubit code. Logical $|0\rangle$ becomes
$|000\rangle$, and logical $|1\rangle$ becomes $|111\rangle$.

If one qubit flips, majority vote recovers the original. This is just classical repetition. Quantum codes are more sophisticated, protecting against both bit-flips and phase errors.

### The HaPPY Code

In 2015, Patrick Hayden, Sepehr Nezami, Fernando Pastawski, John Preskill, and Beni Yoshida built a toy model of holography using error correction-the **HaPPY code**.

They constructed a tensor network in which the bulk is the logical
information and the boundary is made of the physical qubits.

Information in the bulk is redundantly encoded in the boundary. Erase part of the boundary and bulk information survives-you can recover it from the remaining boundary.

This is exactly the recovery rule: I(Bulk : Erased | Remaining) is approximately 0.

The "gravity" in the HaPPY code emerges from the code structure. Regions of the bulk are closer when they share more boundary support. Distance becomes a property of information, not something fundamental.

## 7.12 What Recovery Implies

Recovery sits on a strong foundation. No-cloning blocks naive copying. Strong
subadditivity guarantees that conditional mutual information cannot go
negative. Fawzi-Renner and Petz show that when the missing correlation is
small enough, there is a map that rebuilds what looked lost.

The physics mirrors the mathematics. Ordinary quantum evolution keeps
information in play. Black-hole evaporation is read through the Page curve.
Entanglement wedges reconstruct bulk data from boundary data. Quantum error
correction works in the lab, which means the core logic of encoded recovery is
not speculative hand-waving. The world keeps telling us the same thing:
information does not need to sit in one place to survive.

---

## 7.13 The Indestructible Past

The recovery rule has a startling implication: in this recoverability picture, nothing is ever truly lost.

If the universe is unitary and holographic encoding is robust, information is not simply destroyed; it is redistributed into increasingly nonlocal correlations of the full quantum state.

The Library of Alexandria? The scrolls burned, but the information scrambled into smoke, heat, and light. That radiation spread across the cosmos at light speed. It is diluted across an unimaginably vast region of space-but it's still there. In principle, with a computer the size of the observable universe, you could run the Petz map and watch the smoke reconstitute into Sophocles.

Paleontology and astronomy use weak versions of this. Fossils preserve information about creatures from millions of years ago. Astronomy records light that has traveled for billions of years before reaching our telescopes. The cosmic microwave background is one vivid example of very old information preserved in radiation.

The recovery rule says this is not accident or luck. It's structural: the past is encoded in the present in increasingly scrambled form.

### The Structural Constraint

Of course, practical recovery is impossible. The computation required to
recover the Library of Alexandria would exceed any conceivable technology.
Chaos amplifies tiny errors. A single misplaced bit in trillions grows into
garbage.

This distinction matters enormously. The past is recoverable in principle but
inaccessible in practice. That gives us both unitarity and the lived arrow of
time. The past is not erased. It is encrypted with a key we will never find.

## 7.14 Reverse Engineering Summary

Information does not have to be freely copied to remain recoverable.
No-cloning blocks duplication, but recovery survives because the information is
encoded across extended correlations. That is how a noisy world can still
carry history. It is why observers can agree on a past they never saw. It is
why black holes do not behave like cosmic shredders. And it is why spacetime
starts to look like a code, a structure whose geometry and stability are tied
to the same redundancy that protects information.

Shannon started with a practical problem-sending messages over noisy phone lines. His solution, redundancy, turns out to be built into spacetime itself. The universe is the ultimate error-correcting code.

---

We have the Screen. We have the Algebra. We have the Consistency Rules. We have Recovery.

But where does space come from? Where does time come from? How does the abstract structure of quantum information become the geometry we navigate?

The next chapters turn recovery into geometry. We'll see how boundaries encode interiors, how entanglement draws the map, and how the consistency conditions we've developed start to look suspiciously like gravity.
