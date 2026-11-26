# Free-end current, one-quantum arrivals and reference phase

Number the n two-load dimers from the free end toward the fixed pin. In the previously derived Ising representation define two Hermitian Majoranas per dimer by prefixing its Z or Y by X on every earlier dimer. The physical tip current is the first Y, with empty prefix. Include the fixed-end pin as in the [native mode derivation](../../row_08_annihilation/daughter_modes/README.md). The Majorana generator couples the Z-family to the Y-family through a rectangular bidiagonal matrix with entries 2h and ±2J. The Y-family has no zero mode. The tridiagonal square on its n coordinates has eigenvectors proportional to sin(jk) with k = πα/(n+1), α = 1,…,n, up to an alternating sign convention. Its positive square-root eigenvalues are

$$
\epsilon_\alpha=2\sqrt{(J-h)^2+4Jh\sin^2(k/2)},\qquad
w_\alpha=2\sin^2(k)/(n+1).
$$

Here w is the squared overlap with the first Y coordinate. Since the Majoranas pairwise anticommute and square to identity, the Heisenberg-evolved tip satisfies

$$
\tfrac12\{Y(t),Y\}=C(t) I,\qquad C(t)=\sum_\alpha w_\alpha\cos(\epsilon_\alpha t).
$$

A first sign s projects onto Pₛ = (I+sY)/2. Sandwiching the anticommutator by Pₛ gives PₛY(t)Pₛ = s C(t) Pₛ. Hence the conditional probability of sign b at time t is (1+bsC(t))/2 on every nonempty first-sign branch. The product of signs has mean C(t) for *any* incoming state, including a state entangled with another wing. For n=3 the three weights are one quarter, one half and one quarter. This also explains the unchanged line weights on the native sham. The acquisition computes independent dense Schrödinger evolutions from the collision and does not use this identity to generate the counts. [Matrix tests](test_clock.py) · [Checked algebraic identities](FreeEnd.lean)

Every Y coordinate is a linear combination of creators and annihilators of the positive-frequency normal modes. It changes occupation parity. Thus QY(t)Q=0 for the sharp one-mode projector Q under native free waiting. Reading Y nonselectively sends a density operator to (ρ+YρY)/2. If QρQ=ρ, the second summand has opposite parity, and the probability of repeating the herald is exactly one half. A general fixed terminal effect need not be parity odd. Its one-mode matrix elements can carry only differences of one-mode energies. Those energies lie in [2(J−h),2(J+h)] for J>h. Thus every such difference is at most 4h = 2; on the actual three-dimer receiver the measured width is 33.21463444−31.80075564 = **1.41387880**. The lower edge of the infinite uniform wing dispersion is 31.5. These scales refer to excitation energy relative to the ordered material ground. The incoming one-quantum state can have coherence *between* modes.

Suppose an independent ancilla, initially uncorrelated with the wing, evolves freely and is coupled only at the fixed terminal read. If its state is stationary, tracing any joint effect against it gives a fixed wing effect. The same difference-frequency bound follows. A ticking ancilla can supply its own frequencies and their sidebands. For a target vacuum/one transition operator D, a sharp one-quantum target has zero expectation of D; an independent reference cannot turn this into relative vacuum/one coherence by a product read. An initial joint routing, coherent preparation or target-reference interaction is a different operational resource. This distinguishes an independent external clock, a two-time material response and a phase physically stored in the arriving excitation. The proof leaves open other particle sectors, heralds which preserve cross-number coherence, and relational instruments that intervene before waiting.
