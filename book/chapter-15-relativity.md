# Chapter 16: Matter, Motion, and Classical Physics

## 16.1 The Intuitive Picture: Matter Is Stuff, Motion Is Force

Before we get technical, let's state the common-sense picture most of us grew
up with.

**The intuitive picture**: Matter is made of tiny objects moving around in
space. Each object has a position and velocity. Forces push them, pull them,
and bend their paths. Energy is a kind of fuel that keeps the motion going.

In this view, the world is a stage (space), time ticks forward, and matter is
the cast. Classical physics is the script: Newton's laws, conservation of
energy, and the principle of least action.

This picture works spectacularly well at everyday scales. So why not take it
as fundamental?

## 16.2 The Surprising Hint: The Classical World Is Not Fundamental

Quantum physics breaks the intuitive picture in three ways:

1. **Particles do not have definite paths.** In the double-slit experiment,
   each particle explores multiple paths at once. The classical trajectory
   appears only after interference and measurement.

2. **Fields are more fundamental than particles.** The same electron can be
   created or destroyed. "Particle" is not a permanent object but a
   long-lived excitation of a field.

3. **Energy is more than fuel.** It is a generator of time evolution, and it
   is tied to symmetry. In relativity, energy and momentum are components of a
   single object, and mass is energy at rest.

The hint is clear: the classical picture is an emergent approximation.
The question is not "why does classical physics work?" but "what makes it
work so well?"

## 16.3 The First-Principles Reframing: Matter as Stable Patterns

In OPH, **matter is a stable pattern in the screen data**.

Think of the screen as a high-resolution, quantum information canvas. Most
patterns are noisy and ephemeral. Some are stable: they survive overlap
consistency, persist under modular time, and can be tracked across patches.
Those stable patterns are what we call **particles**.
That's the conceptual picture. The detailed particle story is richer: the framework already fixes the gauge structure and several prominent particles, while the more delicate masses are still being filled in.

The key reframing is:

- **Matter is not a primitive substance.**
- **Particles are not tiny billiard balls.**
- **Matter is the set of robust, localized excitations of the net of
  algebras on the screen.**

A useful analogy is a ripple in a pond. The water is the substrate, but the
ripple is a pattern that moves and interacts. The ripple is not a separate
thing; it is a stable excitation. Particles play the same role in the
emergent EFT.

## 16.4 What Is a Particle?

In ordinary physics, a particle is defined by symmetry. Wigner showed that
"particle types" are **irreducible representations of the Poincare group**,
classified by mass together with spin or helicity. This is more than a
definition; it explains why particles have sharp representation labels.

In OPH, this appears after two steps:

1. **Lorentz kinematics is recovered on the extracted prime geometric subnet** from
   geometric modular flow on the screen, with the scaling-limit geometric cap-pair
   realization and ordered cut-pair rigidity scaffold still open beneath full
   UV internalization (as described in Chapter 15).
2. **Localized excitations organize into representations** of this emergent
   symmetry in the EFT regime.

So particles are **the representation theory of emergent symmetries**.
A "mass" is the representation label that tells you how the excitation
responds to time translations. A "spin" is the label for how it responds to
rotations.

Once Lorentz kinematics is in place, energy and momentum form a single
four-vector and mass is the invariant. This gives the familiar relation

$$E^2 = p^2 + m^2$$

(in units where $c=1$). This is why energy and mass are so tightly linked in
classical physics: they are two faces of the same symmetry.

This is why particles are universal: they are bookkeeping devices for symmetry
classes, not fundamental objects.

## 16.5 What Is Energy?

Energy is more than a number. It is the generator of time evolution.

In OPH, time is **modular flow**. The generator of modular flow is
the modular Hamiltonian:

$$K = -\log \rho.$$

In the EFT regime, this connects to the ordinary Hamiltonian and the stress
energy tensor. On null surfaces, the modular generator becomes an integral of
$T_{kk}$ (null energy density), and in a UV CFT it reduces to the standard
stress-tensor charge on small caps.

So energy has a clean meaning:

- **Energy is the charge that generates time translations**, and in the
  emergent EFT it is encoded in the stress tensor $T_{ab}$.

Conservation of energy then follows from symmetry: if the emergent action is
invariant under time shifts, Noether's theorem gives a conserved energy.

## 16.6 Motion and Forces: Why Things Move the Way They Do

Classical motion can be described in two equivalent ways:

- **Force laws**: $F = ma$.
- **Variational laws**: trajectories extremize an action.

Both are effective descriptions. In OPH, motion is a property of
stable patterns moving under modular flow, observed consistently across
patches. Forces describe how those patterns interact within the emergent EFT.

**Locality and consistency constrain motion.** Overlaps
force observers to agree on what happened. The Markov structure enforces local
relations between neighboring regions. These requirements leave very little
freedom in the form of effective equations of motion.

The broad shape of the low-energy laws is set by the same consistency
structure that gives us gauge symmetry in Chapter 14. The exact
coupling-by-coupling story is part of the continuing program.

## 16.7 Why the Principle of Least Action Appears

The principle of least action can sound mystical, but it is a direct
consequence of quantum interference.

In quantum mechanics, the probability amplitude for a particle to go from
A to B is a sum over all possible paths:

$$\mathcal A \sim \sum_{\text{paths}} e^{i S/\hbar}.$$

Here the action is

$$S = \int L(q, \dot q, t)\,dt,$$

where $L$ is the Lagrangian.

When the action $S$ is large compared to $\hbar$, phases oscillate rapidly and
cancel out. Only paths where $S$ is stationary survive. This yields the
Euler-Lagrange equations:

$$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot q}\right) = \frac{\partial L}{\partial q}.$$

So the "least action" rule is not a separate law. It is the classical limit
of quantum consistency. In OPH, the effective action is part of the EFT
bridge. Once that bridge is in place, the stationary-action description
follows in the usual semiclassical limit.

Historically it is called "least" action, but what really survives is
**stationary** action: small variations do not change the path to first order.

## 16.8 The Classical Limit: Why the World Looks Deterministic

Classical physics is an **emergent approximation** that appears when:

- The action is large compared to $\hbar$ (stationary phase dominates).
- The system is strongly entangled with its environment (decoherence).
- Observers coarse-grain over microscopic details (MaxEnt selection).

### Why Decoherence Matters for Consistency

Decoherence is crucial. In this framework it is part of what makes a stable shared classical description possible, not a lucky accident bolted on afterward.

Here's why. The overlap condition demands that observers agree on shared observables. But quantum mechanics permits states that are superpositions, "both A and B." If macroscopic interference remained broadly accessible at everyday scales, different observers sampling different environmental fragments would fail to recover a single robust public record.

Decoherence solves this by rapidly entangling macroscopic objects with their environments. This entanglement has a specific structure: it correlates the object's state with environmental "records" that can be accessed by multiple observers independently.

The key insight from **quantum Darwinism** (Chapter 6) is that only certain states, pointer states, get their information redundantly copied into the environment. These are the states that many observers can access and agree upon. Coherent superpositions do not get copied into stable public records; once the system is entangled with the environment, the interference becomes locally inaccessible.

**Classical facts are quantum states that pass the consistency filter.** A "classical" property is one that:
1. Gets redundantly encoded in the environment
2. Can be accessed through multiple independent channels
3. Produces agreement when different observers check

The pointer basis, the set of states that decohere into classical alternatives, is not arbitrary. It is constrained by the system-environment coupling and by which observables can be stably shared across patches. States that can't be consistently shared across patches don't survive as "real" in the intersubjective sense.

So classical physics is the **stable, compressible limit** of the deeper quantum structure: the patterns that survive the consistency filter. The world looks deterministic because only the consistent patterns-the ones that all observers can agree on-rise to the level of "facts."

### Why Classical Physics Isn't Fundamental

This resolves an old puzzle: why does the quantum world give rise to classical physics at all?

In the standard picture, classical physics is an approximation that breaks down at small scales. But our model inverts this: classical physics is what emerges when consistency constraints are satisfied. The classical world isn't the fundamental reality poorly approximating quantum mechanics-it's the consistent core that multiple observers can share.

The quantum world is larger but less shareable. Superpositions exist, but they can't be consistently communicated. When you try to share quantum information broadly, decoherence kicks in, and you're left with classical correlations.

**Classical physics is the public face of quantum reality.** It's not fundamental; it's the stable consistency regime that many observers can share.

## 16.9 Reverse Engineering Summary

| Intuitive Picture | Surprising Hint | First-Principles Reframing |
|---|---|---|
| Matter is fundamental stuff moving through space | Quantum interference and creation/annihilation show particles are not permanent objects | Matter is a stable excitation pattern in the screen net |
| Energy is just a fuel | Energy is a generator tied to symmetry | Energy is the charge of emergent time translations (stress tensor) |
| Motion follows force laws | Trajectories emerge from interference | Least action is the classical limit of quantum consistency |
| Classical physics is the fundamental layer | Quantum mechanics underlies classical physics | Classical physics emerges because only consistent patterns survive the overlap filter |

**The key reverse engineering insight**: classical physics is not the starting
point. It is what you get when quantum information on the screen organizes
into stable patterns, when modular time becomes geometric, and when overlap
consistency enforces locality. Particles, energy, and motion are the emergent
vocabulary of that stable regime.

**Why classical physics emerges**: The overlap condition demands that observers agree on shared observables. Decoherence, the rapid entanglement of macroscopic objects with their environments, is part of what makes stable shared records possible. Only pointer states that get redundantly copied into the environment can be consistently shared across patches. Classical facts are quantum states that pass the consistency filter. The deterministic, objective world of everyday experience is the public face of a quantum reality too fragile to be broadly shared.

We've seen that spacetime, particles, and classical physics all emerge from the screen through consistency requirements. But why these particular laws? Why these constants? Could the universe have been different?

The next chapter explores a radical idea: physical laws themselves are evolutionary survivors. Just as life evolves through natural selection, laws evolve through cosmic selection.

This is **Chapter 17: Darwin's Laws**.
