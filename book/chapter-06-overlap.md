# Chapter 6: Overlap and Agreement

## 6.1 The Intuitive Picture: Local Causes Explain Correlations

Start with the classical picture.

When two distant events are correlated, there must be a common cause in their
shared past. Correlations come from shared history, hidden connections, or
pre-existing properties. If I flip two coins and they always match, either the
coins were manufactured together with matching weights, or someone is signaling
between them. There's no spooky magic at a distance.

This is the worldview of classical physics and common sense. Einstein himself held it dear. Objects have definite properties whether or not we measure them. Measurements reveal pre-existing facts. If two particles are correlated when measured far apart, they must have carried that correlation with them from the start, like matched gloves packed in separate suitcases.

The technical term for this intuition is **local realism**. "Local" means that
nothing influences distant events faster than light. "Realism" means that
properties exist independently of observation.

Local realism is so natural that questioning it seems absurd. Of course the moon exists when nobody's looking. Of course a particle has a definite spin before you measure it. Of course distant correlations require either a shared cause or a connecting signal.

Bell's theorem broke that picture.

## 6.2 The Surprising Hint: Bell's Theorem and Nonlocal Correlations

### Einstein's Challenge: The EPR Paper

To understand why quantum consistency is hard, we need to visit 1935.

Albert Einstein was fifty-six years old and deeply troubled. He had helped create quantum mechanics (his 1905 paper on the photoelectric effect was one of the founding documents), but he never accepted its implications. "God does not play dice," he famously declared.

In May 1935, Einstein, Boris Podolsky, and Nathan Rosen published what became known as the EPR paper. Its title was dry: "Can Quantum-Mechanical Description of Physical Reality Be Considered Complete?" The New York Times was less restrained, announcing the work under the headline "Einstein Attacks Quantum Theory."

EPR constructed a thought experiment. Take two particles created together and let them fly apart. Quantum mechanics says they can be *entangled*, with correlations that have no classical analog. Measure a property of one particle, and you instantly know the corresponding property of the other, even if they're light-years apart.

The puzzle is sharp. According to quantum mechanics, the particles don't have definite values until measured. But if I measure particle A and find it has spin-up, I instantly know particle B has spin-down-without ever touching particle B. Did my measurement somehow affect particle B instantaneously? Einstein called this "spooky action at a distance" and found it absurd.

EPR concluded that quantum mechanics must be incomplete. The particles must have had definite values all along-values we just didn't know. There must be "hidden variables" underneath the quantum description.

Most physicists shrugged and went back to calculating. Niels Bohr wrote an impenetrable response. The debate seemed philosophical, not scientific.

For nearly thirty years, everyone assumed it couldn't be settled by experiment. Then along came John Bell.

### Bell's Breakthrough

John Stewart Bell was an Irish physicist working at CERN in the 1960s. He was quiet, precise, and privately dissatisfied with the foundations of quantum mechanics. In his spare time, between designing particle accelerators, he worked on a problem everyone else had abandoned.

In 1964, Bell published a short paper that changed everything. He proved that
the question was empirical and designed an experiment that could distinguish
quantum mechanics from the relevant class of local hidden-variable theories.

The key was correlation. When two observers measure entangled particles, their results are correlated. Bell showed that local hidden-variable theories set a ceiling on how correlated the results can be. This ceiling is called the Bell inequality:

$$|S| \leq 2$$

The quantity S combines correlations from four different measurement settings. Local hidden-variable models of Einstein's "reasonable" picture cannot produce correlations stronger than 2.

Quantum mechanics predicts something stronger:

$$S = 2\sqrt{2} \approx 2.83$$

That's a 41% violation, large enough to test.

This $S$ is the Bell-CHSH correlation score, not entropy. The vertical bars
mean absolute value. A local hidden-variable model can arrange many patterns,
but it cannot push this score above 2. Quantum mechanics can reach
$2\sqrt{2}$ for the right entangled state and measurement angles.

### What Makes This So Strange

Let me be concrete. Alice and Bob each receive one particle from an entangled pair. They're far apart: different continents, different planets, it doesn't matter. Each chooses randomly whether to measure their particle along angle A1 or A2 (for Alice) or B1 or B2 (for Bob).

In the hidden variable picture, each particle carries a tiny instruction manual: "If measured at angle A1, give result +1. If measured at B2, give result -1." And so on. The instruction manual was written when the particles were created. The particles are like correlated coins-maybe both were programmed to give the same answers.

Bell's genius was realizing you could test this. Run the experiment thousands of times. Calculate the correlations. If the world is described by a local hidden-variable model of this Bell type, S <= 2, however the instruction manuals were written.

Quantum mechanics can beat that ceiling. When Alice and Bob choose the right measurement angles, entanglement produces correlations of 2 times the square root of 2.

### The Experiments

For two decades after Bell's paper, experimentalists raced to test it. The challenges were enormous. You needed to create entangled pairs reliably, separate them, measure them independently, and collect enough data to beat statistical noise.

Alain Aspect performed the definitive early tests at Orsay, outside Paris, in 1981-82. His team used pairs of entangled photons created by exciting calcium atoms. They measured polarizations and found S approximately equal to 2.70, well above 2 and consistent with quantum predictions.

But there were loopholes. What if the particles somehow communicated with each other? (The locality loophole.) What if only certain particles got detected? (Detection loophole.) What if the measurement choices weren't truly random? (Freedom-of-choice loophole.)

Over the following decades, experimenters closed the major loopholes one by one. The 2015 "loophole-free" Bell tests by teams in Delft, Vienna, and Colorado closed the locality and detection loopholes simultaneously, while using fast random setting choices that strongly constrained freedom-of-choice concerns. The particles were separated by large distances, the measurements were completed before any signal could travel between them, and the detection efficiency was high enough to rule out selection effects.

**Suitable entangled Bell experiments repeatedly violate Bell inequalities.**

This means at least one ingredient in the classical Bell-assumption set must
fail. The pressure falls on locality, realism, or the assumption that
measurement settings are not secretly pre-coordinated with hidden variables.

Many physicists read the Bell results as strong pressure against naive local
realism, although the exact interpretive lesson is contested. Accepting a
deeper nonlocal structure or preferred causal bookkeeping carries a high
explanatory cost for many physicists.

Quantum correlations exceed what any local hidden variable theory permits. The intuitive picture of pre-existing properties carried from a common past is experimentally contradicted.

## 6.3 The First-Principles Reframing: Consistency and Nonlocal Correlations

The reverse-engineering question is simple: why does nature behave this way? What principle would make such nonclassical correlations structurally natural?

### Objectivity Is Agreement

Stand on a street corner in New York City and look across the road. A bright
red Ferrari is parked there, gleaming and
expensive, the kind of car that makes people stop and stare. A second
observer, Bob, is standing fifty feet down the block. He sees the side profile
and the license plate. A third observer, Charlie, is looking out of a
second-story window and sees the roof of the car.

We take for granted that there's a single, objective "real" Ferrari sitting there. But ask a dangerous question: *How do we know the car is real?*

The only evidence any of you has is your own private sensory data, your patch.
You have the view from the corner. Bob has the view from the sidewalk. Charlie
has the view from above.

If Bob walked up to you and said, "That's a nice blue elephant," you would have a problem. If Charlie yelled down, "No, it's a green helicopter," the world would dissolve into chaos.

**Objectivity is simply the process of checking for agreement.**

If all three of you agree on the overlap of your visual fields and report "Red
Car," you conclude the car is real. The "object" emerges from the intersection
of your views. Reality is the consensus reached by a network of observers.

### Why Classical Consistency Is Easy

In classical physics, checking consistency is much simpler on basic overlap structures than in the quantum case.

The state of a classical system is a point in phase space, a list of all positions and momenta. If Alice knows the full state, so does Bob. They're reading from the same book.

When information is partial, we use probability distributions. Let rho_A be Alice's distribution, rho_B be Bob's. If they both measure observable O, their expected values must agree:

$$\langle O \rangle_A = \int O(s)\rho_A(s)ds = \int O(s)\rho_B(s)ds = \langle O \rangle_B$$

For tree-like overlap structures, agreeing marginals can be glued into a joint
distribution. If Alice's distribution over variable X matches Bob's marginal
over X, and Bob's distribution over variable Y matches Carol's marginal over Y,
there is a joint distribution P(X,Y,Z) that reproduces all the marginals.

The angle brackets mean expectation value, the average result predicted for
observable $O$. The variable $s$ labels a classical state, and $\rho_A(s)$ and
$\rho_B(s)$ are Alice's and Bob's probability distributions over those states.
The integrals add the observable's value over all possible states, weighted by
the probabilities each observer assigns.

In general overlap graphs, the classical marginal problem can fail and is computationally hard; agreement on pairwise overlaps is not always sufficient.

### Why Quantum Consistency Is Hard

Quantum mechanics is different.

Given reduced density matrices that are pairwise consistent on overlaps, does a global state exist that produces them all?

Unlike the classical case, the answer can be **NO**. This is the Quantum Marginal Problem (QMP).

Why can't you just glue quantum marginals together? The answer involves one of
quantum mechanics' sharpest constraints: **entanglement is monogamous**.

If particles A and B are maximally entangled, then A cannot also be maximally entangled with C. You can't share maximal quantum correlation with more than one partner.

One standard qubit monogamy relation is the Coffman-Kundu-Wootters inequality:

$$\tau_{A:B} + \tau_{A:C} \leq \tau_{A:BC}$$

In this qubit setting, A's pairwise entanglement budget with B and C cannot exceed its total entanglement with BC together.

Think of it like attention. If you're having a deeply intimate conversation with one person, you can't simultaneously have an equally deep conversation with someone else. Quantum correlations work the same way.

### The Consistency Filter

**Bell-violating correlations are treated here as a structural feature that may help an observer-consistency framework remain viable.**

Imagine the space of all possible local states-all assignments of density matrices to patches. This space is enormous. Most assignments are inconsistent; different patches disagree on overlaps.

Apply the overlap consistency condition. Any assignment where patches disagree gets filtered out. The consistent assignments form a tiny subset.

**Reality is the collection of local states that survives the consistency filter.**

The hardness of the Quantum Marginal Problem tells us the filter is doing real work. The constraints are genuinely restrictive. This helps explain why overlap consistency is a nontrivial structural requirement, with real content. It also suggests that the allowed state-space is highly structured.

Overlap conditions favor allowing correlations that exceed classical bounds.
Bell-violating correlations can then be read as part of the quantum structure
available to an observer-consistency framework, without appealing to a large
hidden-variable overhead.

In a universe built on observer agreement, the nonlocal correlations that so troubled Einstein are not inexplicable. They become part of the consistency structure, with no unexplained add-on.

## 6.4 Defining the Overlap

What does Bell's theorem have to do with observer patches?

Everything.

Bell showed that when two observers access the same entangled system, their
correlations can exceed classical bounds. They can't *communicate* faster than
light, each observer's local statistics look completely random, but when they
*compare notes*, patterns emerge that no classical account can explain.

This comparison is overlap. When Alice and Bob's patches both include information about an entangled system, their descriptions must be compatible in a very specific way.

Recall our setup. Alice has patch P_A with algebra A(P_A). Bob has patch P_B with algebra A(P_B). If their patches overlap, they share a region:

$$R = P_A \cap P_B$$

This region R is the "Looking Glass." It contains observables common to both. For reality to be consistent, **Alice and Bob must agree on the state of the Looking Glass.**

$P_A$ and $P_B$ are Alice's and Bob's patches. The symbol $\cap$ means
intersection. $R$ is the part both patches contain.

In a simple finite-dimensional toy model, Alice describes her patch with density matrix rho_A and Bob describes his with rho_B. Then consistency on the overlap can be pictured as equality of the reduced descriptions on region R:

$$\text{Tr}_{A \setminus R}(\rho_A) = \text{Tr}_{B \setminus R}(\rho_B)$$

This is only the toy-model picture. More generally, the right statement is that the two restricted states agree on the shared overlap algebra.

The set-minus symbol $\setminus$ means "remove this part." Alice traces out
the part of her patch outside the overlap, and Bob does the same. If the two
reduced density matrices match, their descriptions agree on the shared region.

![Overlap consistency means both patches assign the same state to the shared region.](../assets/book_diagrams/overlap-consistency.svg){width=76%}

### The Mathematical Translation

Here is the same equation in ordinary language.

A density matrix is quantum mechanics' way of describing partial knowledge. If you know a system is definitely in state |psi>, you use a pure state. If you only know the system is in state |psi_1> with probability p_1 or state |psi_2> with probability p_2, you use a density matrix:

$$\rho = p_1|\psi_1\rangle\langle\psi_1| + p_2|\psi_2\rangle\langle\psi_2|$$

The "trace" operation (Tr) is how you marginalize-how you focus on one part of a system while ignoring the rest. If Alice has access to particles A and B but Bob only has access to B, then "Tr_A" traces out particle A, leaving just the description of B.

In this density matrix, $p_1$ and $p_2$ are probabilities. The kets
$|\psi_1\rangle$ and $|\psi_2\rangle$ are possible pure states. The paired bra
and ket $|\psi\rangle\langle\psi|$ is the projector onto that state. A density
matrix is therefore a weighted quantum mixture of possible state assignments.

The consistency condition says: when Alice traces out everything Bob can't see, and Bob traces out everything Alice can't see, they'd better end up with the same description of the overlap.

### Overlap Is a Protocol

In 1665 Christiaan Huygens watched two pendulum clocks mounted on one support
settle into a shared rhythm. He called the effect an “odd sympathy.” The clocks
had not exchanged philosophical arguments. Their common beam carried the
coupling.

An observer overlap also needs a physical coupling. Two astronomers looking at
the same star from different continents need a shared reference frame,
synchronized clocks, and calibration conventions for their instruments. The
overlap becomes useful only when they agree on the translation between their
frames. Agreement always includes some shared dictionary.

Physics uses standardized units, coordinate systems, and calibration procedures because they are the protocols that make overlap possible.

The finite patch architecture turns that dictionary into a port protocol. A
patch does not hand its entire private state to its neighbor. It exposes one
boundary packet, receives the packet routed to the same interface, and compares
the two after the declared translation. A checkable mismatch selects an
allowed local repair. The patch then writes the result into a rereadable record
and exposes the boundary again. On the icosahedral screen-sieve branch, twelve
ports cover the icosahedral directions and $A_5$ symmetry keeps any one port
from being privileged by its label. Another material can realize the same
protocol only when the complete port geometry, translation, readout, dynamics,
repair, refinement, and checkpoint behavior agree.

This is operational synchronization: exposed packets settle to compatible
public records when the confluence conditions hold. Physical phase locking is
a narrower claim. It requires actual phase variables, a coupling law, and a
phase-sensitive test. Recurrent subchannels in an echosahedral carrier can
host such a test, but the consensus theorem does not manufacture it.

### Overlap Has a Cost

Sharing observations isn't free. You need energy to send signals and memory to store them. Every message takes time. That means overlap is always limited. You only share a slice of your full experience.

An observer has finite capacity. If you want to make your patch more consistent with others, you spend resources exchanging data. Agreement is work.

This cost will become important later when we discuss how classical reality emerges. The facts that get widely shared are the ones that can be copied cheaply, and quantum mechanics places fundamental limits on what can be copied.

## 6.5 The Quantum Marginal Problem Is QMA-Complete

In 2006, computer scientist Yi-Kai Liu proved that deciding whether quantum marginals are compatible is QMA-complete.

QMA is the quantum analog of NP. Just as NP captures problems where solutions are easy to verify but hard to find, QMA captures problems where a quantum computer could verify a quantum proof, but finding the proof is impossibly hard.

Being QMA-complete means the Quantum Marginal Problem is as hard as any problem in the class. If you could solve QMP efficiently, you could solve any QMA problem efficiently.

If the complexity labels are unfamiliar, keep the operational lesson: local
quantum consistency is serious bookkeeping. In the worst case it is as
hard as the hardest verification problems quantum computers are expected to
handle.

### Why the Hardness Matters

In classical physics, local data often determine global data on simple overlap structures, and compatibility checking is much easier than in the quantum case.

In quantum physics, local data constrain but don't determine global data. In
the worst case, checking consistency is computationally hard: there is no known
efficient general algorithm to decide whether quantum marginals are compatible.

This shows that quantum mechanics hides global structure in a deeply complex
way. You cannot easily deduce the whole from the parts.

## 6.6 A Concrete Counterexample: Three Qubits

One quantum case looks consistent locally but cannot be glued together.

Consider three qubits $A$, $B$, and $C$. Suppose every pair were maximally
entangled: $A$ with $B$, $B$ with $C$, and $A$ with $C$. Each pair looks
harmless on its own.

Each pair being maximally entangled seems fine. The reduced state of any single qubit is maximally mixed-equal probability of spin-up or spin-down. That's consistent.

Try to find a state |psi>_ABC that produces all three Bell pairs. You can't.

The obstruction is simple. For any pure state of three parties, there's a constraint:

$$S(\rho_A) = S(\rho_{BC})$$

The entropy of A equals the entropy of BC. This is a consequence of entanglement structure.

Here $\rho_A$ is the reduced density matrix of subsystem $A$, and
$\rho_{BC}$ is the reduced density matrix of the joint subsystem made from
$B$ and $C$. For a pure total state, the entropy of one side of a split equals
the entropy of the other side.

If AB is maximally entangled, then rho_A is maximally mixed: S(rho_A) = 1 bit.

So S(rho_BC) = 1 bit.

But if BC is maximally entangled, then rho_BC is pure, so S(rho_BC) = 0.

Contradiction. The marginals are individually valid but globally incompatible. Monogamy strikes again.

### GHZ and W: Two Ways to Share

There are different ways to distribute entanglement among three particles.

The **GHZ state**:
$$|\text{GHZ}\rangle = \frac{1}{\sqrt{2}}(|000\rangle + |111\rangle)$$

Look at any pair-say, qubits A and B. Trace out C. The reduced state shows no entanglement at all. AB looks completely classical. But when all three particles are measured together, perfect correlations emerge. It's an all-or-nothing state.

The **W state**:
$$|W\rangle = \frac{1}{\sqrt{3}}(|001\rangle + |010\rangle + |100\rangle)$$

In the W state, every pair has some entanglement, but none is maximal. The entanglement is spread around, diluted.

The GHZ state is named after Greenberger, Horne, and Zeilinger. The W state is
named for the shape of its entanglement pattern, not for a person. The
normalizing factors $1/\sqrt{2}$ and $1/\sqrt{3}$ make total probability 1.
The three slots in each ket are the three qubits. These two states show that
"three-party entanglement" is not one thing. It has distinct species.

Quantum agreement is a budget. Spend it on one overlap and you have less for another.

## 6.7 The Kochen-Specker Theorem

There's an even more direct demonstration that quantum mechanics resists classical consistency.

In 1967, Simon Kochen and Ernst Specker proved a theorem with a hard consequence: in a Hilbert space of dimension 3 or higher, there is no single noncontextual assignment of pre-existing values to all quantum observables.

### What Does This Mean?

Imagine trying to create a "cheat sheet" for a quantum system, a list saying "if you measure observable A, you'll get value a; if you measure observable B, you'll get value b; ..." and so on for every possible measurement.

Kochen-Specker says: no such cheat sheet exists.

The sharp lesson is narrower and more precise: there is no single noncontextual cheat sheet assigning pre-existing values to all observables at once. Any viable hidden-variable picture must therefore be contextual, and the measurement setting cannot be treated as irrelevant bookkeeping.

Kochen-Specker closes one of the most attractive escape routes from quantum
strangeness. A system cannot carry a full secret list of answers, with
measurement merely revealing whichever answer belongs to the question an
experimenter happened to ask. The question itself belongs to the physics.

The practical consequence is easy to state. Quantum systems do not come with a sealed answer key that every possible experiment simply reads out. A measuring setup selects a compatible family of questions, and the consistency conditions live inside that family. Change the family, and the bookkeeping changes with it. OPH leans on exactly that point. What observers can stably compare depends on the overlap algebra they actually share.

### The Peres-Mermin Magic Square

The Mermin-Peres square gives a vivid example. Arrange nine observables for two qubits in a 3x3 grid. Each row and each column contains three observables that can be measured together (they commute).

The product of observables in each row is +I (the identity).
The product of observables in each column is +I.
Except the last column, whose product is -I.

Try to assign definite values (+1 or -1) to each observable such that the product rules hold.

The product of all row products = (+1)(+1)(+1) = +1.
The product of all column products = (+1)(+1)(-1) = -1.

But each observable appears once in a row and once in a column. So the product of row products should equal the product of column products.

+1 does not equal -1. Contradiction.

No single noncontextual value assignment exists that satisfies these
constraints. Any viable account must therefore treat value assignment as
context-dependent. This is **contextuality** in the sense highlighted by the
theorem.

That is exactly why this theorem fits OPH so well. Observer overlap happens
between actual measurement contexts, actual accessible observables, and actual
records, not between perfect God's-eye inventories. Context belongs to the
structure.

That gives the theorem a direct narrative role in the book. The world does not stay coherent because every observer secretly samples one master spreadsheet. It stays coherent because local contexts can be glued where they genuinely overlap. Kochen-Specker tells you why the stronger fantasy fails.

## 6.8 Wigner's Friend: Consistency Between Nested Observers

The consistency challenge becomes sharper when observers themselves become part
of the system.

In 1961, Eugene Wigner proposed a thought experiment that puts the observer inside the measurement problem.

Wigner's friend is in a sealed laboratory, measuring a quantum system. From the friend's perspective, the measurement has produced a definite outcome record-say, spin-up. In standard textbook language, the friend would update the system to the corresponding outcome state.

But Wigner is outside the lab. He describes the entire lab-including his friend-using quantum mechanics. From Wigner's perspective, the lab is in a superposition: (friend sees spin-up and atom is spin-up) + (friend sees spin-down and atom is spin-down).

Who's right?

From the friend's view: the measurement record is definite.
From Wigner's view: the isolated lab can be described by a superposed quantum state until he interacts with it.

Both descriptions are internally consistent. The problem arises at the overlap, when Wigner opens the door and compares notes with his friend.

This is the nested-observer version of the whole book. One observer can carry a
finished record while another treats that record as part of a larger
superposed description, provided they have not compared notes. The tension
appears when communication begins and the two descriptions have to settle into
one common account.

Wigner's friend matters beyond foundations theater because it is the simplest
clean model of nested access. One observer inhabits a definite record, another
treats that record as part of a larger quantum state, and both
descriptions remain admissible until communication forces a common restriction.
The whole OPH project keeps asking how such restrictions line up without
contradiction.

At that moment, their descriptions must agree. The consistency condition forces a resolution. Before the door opens, they can maintain different descriptions. After it opens, they share an overlap, and quantum mechanics demands their states match on that overlap.

This is observer-relativity with constraints. The "facts" depend on who's
asking, and the overlap conditions determine which facts can coexist.

Frauchiger and Renner's 2018 no-go theorem, and the Wigner's-friend experiments that followed it, show that even sophisticated extensions of quantum mechanics struggle to maintain consistency when observers observe observers. The consistency conditions are doing real work.

## 6.9 Quantum Darwinism: How Overlaps Build Objectivity

If quantum mechanics is so resistant to consistency, how does the classical world emerge? How do we get the stable, objective facts that everyone agrees on?

The answer involves a concept called **quantum Darwinism**, developed by Wojciech Zurek.

The mechanism is environmental copying. A quantum system interacts with air molecules, photons, and everything around it. Some information about the system gets copied into the environment. Quantum mechanics forbids perfect copying of arbitrary unknown states, so the useful information is redundantly encoded in stable records.

Consider Schrödinger's cat. If the cat is alive, air molecules bounce off it in a certain way. Light reflects off it in a certain way. Heat radiates from it in a certain way. Each of these environmental fragments carries partial information about the cat's state.

When you look at the cat, you are reading information from environmental
fragments. Many observers can read different fragments and agree.

The information that gets redundantly copied is the information that becomes
"objective." It is the information that survives across multiple overlaps.
Quantum superpositions do not get copied this way. Only certain pointer states
remain stable under environmental interaction.

### The Birth of Classical Facts

A classical fact is quantum information that has been copied redundantly into
the environment, made available through multiple independent channels, and made
stable against small perturbations.

The red Ferrari is classical because trillions of photons have bounced off it, carrying correlated information to many observers. In the decoherence / quantum-Darwinism picture, the cat's environmentally stable pointer-state records become effectively classical for observers, while interference between alternatives becomes inaccessible in practice.

Classical objectivity is quantum redundancy. The facts everyone agrees on are the facts that got copied everywhere.

## 6.10 Reality as a Sheaf

The big picture is compact.

We've been building toward a radical view of reality. We do not begin by requiring a single, global "state of the universe." A useful mathematical analogy is a sheaf-like gluing picture in which local data are stitched together by consistency conditions.

A sheaf is the mathematics of local descriptions that may or may not glue into
one global description. Maps are the friendly example. If every local map agrees
on its overlaps with neighboring maps, you can assemble an atlas. If the overlap
rules conflict around a loop, no single map exists.

### The Internet Analogy

Think of the internet. There's no single file called "The Internet" stored somewhere. There are billions of computers, each with its own memory. They communicate via protocols. When my computer sends a packet to yours, we "agree" on the content. The "internet" is the emergent consistency of all these local interactions.

Reality can be treated as a collection of local states, one for each observer,
constrained to agree on overlaps. No single observer occupies a God's-eye view
of one global state.

When a global state exists, that is useful. But we do not require one. Local states satisfying consistency conditions are enough for physics.

### Living Without a Global Wavefunction

This is philosophically adjacent to relational quantum mechanics, QBism, and
the Copenhagen refusal to assign one wavefunction to the universe as a whole.
The distinctive move here is to make the overlap conditions explicit and let
them carry real mathematical weight.

The framework adds a precise mathematical model. Its consistency conditions
can be written explicitly, and the global gluing problem is nontrivial.

### Transitivity and Networks

With many observers, each pair of overlapping patches must agree on their intersection. This forms a web of constraints.

If Alice and Bob agree on their overlap (AB), and Bob and Carol agree on their overlap (BC), then Bob can mediate indirect compatibility between them on simple tree-like covers. Local pairwise consistency can help enforce global structure there, but loops or more general covers can produce frustration unless higher-order constraints are satisfied.

Loops add another condition. Go from Alice to Bob to Carol and back to Alice;
you should return with the same state on shared overlaps. Otherwise the system
has **frustration**: the local assignments cannot all hold simultaneously.

This is analogous to gauge theory and geometry. Move a vector around a loop; if it comes back rotated, there is nontrivial holonomy. Holonomy is the leftover mismatch that appears after a full circuit. Likewise, a loop that does not close cleanly signals an obstruction to global gluing, with no simple globally consistent assignment.

## 6.11 Formal Statement

The consistency condition can be stated directly.

### Setup

We have a screen $S^2$, a collection of patches $\{P_i\}$, an algebra
$A(P_i)$ of observables for each patch, and a state $\omega_i$ on each patch.

### The Condition

For any two patches P_i and P_j with non-empty overlap:

$$\omega_i|_{\mathcal{A}(P_i \cap P_j)} = \omega_j|_{\mathcal{A}(P_i \cap P_j)}$$

The restrictions to the overlap algebra must be the same state.

$\omega_i$ and $\omega_j$ are the states assigned by observers $i$ and $j$.
$\mathcal A(P_i\cap P_j)$ is the algebra of observables available on their
shared patch. The vertical restriction bar means "look only at this shared
algebra."

In plainer English: for any observable O that both Alice and Bob can measure:

$$\omega_i(O) = \omega_j(O)$$

They must assign the same expectation value.

$O$ is any observable in the overlap. The equation says that the two observers
may keep different private descriptions elsewhere. On the overlap, they must
make the same predictions for questions both can actually ask.

### The Patch Graph

The patches form a graph whose nodes are observers and whose edges connect
overlapping patches.

The topology of this graph determines what kind of global structure can emerge. Loops in the graph create constraints. On tree-like covers, local consistency is much easier to promote toward global structure, but in the quantum setting that depends on the relevant compatibility hypotheses. If there are loops, you need to check that going around each loop is consistent.

## 6.12 What Overlap Predicts

If overlap consistency is real, the world should look Bell-like, Markov-like,
and Darwinian all at once. Bell experiments should obey the quantum bound and
violate the classical one. Structured states separated by a good intermediary
region should show small CMI. A genuine global
state, when it exists, should automatically induce matching reduced states on
every overlap. Public facts should appear only where information has been
copied redundantly into the environment.

That is exactly the pattern physics presents. Bell stays within the Tsirelson
limit. Overlaps induced from explicit global states match by construction.
Redundancy is what turns a private quantum state into something many observers
can check. No experiment has revealed a stable public world without
environmental redundancy, and no experiment has pushed Bell beyond the quantum
bound.

## 6.13 Reverse Engineering Summary

Objectivity becomes agreement. Correlations do not need to be
explained by hidden classical instruction sets carried from a common past.
Bell's theorem makes that fantasy too small for the world we observe. Here,
the stronger quantum structure belongs to the consistency conditions linking
patches.

That matters because the gluing problem is hard. Quantum marginals can fail to
fit together even when pairwise overlaps look harmless, and the combinatorics
become vicious as more observers are involved. Bell-violating entanglement
becomes part of the efficient bookkeeping that allows many local perspectives
to remain compatible without an impossible burden of pre-coordination.

---

We have the Screen. We have the Algebra. We have the Consistency Rules.

But what if the web gets torn? What if I measure something here, and you measure something there, and we lose the connection? What if information seems to disappear into a black hole or leak out through quantum noise?

That tension is sharper here. Overlap is contextual, entanglement is limited,
and global compatibility can be brutally hard to check. A world built from such
ingredients should feel fragile. The next chapter explains why it does not.

The recovery chapter examines mechanisms that reconstruct missing information
and preserve the web of consistency when individual links appear broken.
