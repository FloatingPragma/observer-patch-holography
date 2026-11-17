# Coulomb variance in the full ring lower bound

The Hamiltonian, references, clock and preparations are exactly those
of the frozen sixteen-CAR four-cell source. This package changes an
operator inequality. It changes no physical coefficient or mass.

Write the exact electric energy as

```
sum_x(E+q_x)^2/4=(E+qbar)^2+B,
B=sum_x(q_x−qbar)^2/4>=0.
```

Vector diamagnetism bounds the first term below by the scalar kinetic
energy of the full Fock-wavefunction norm. The remaining pointwise
matrix is therefore `W_N(theta)=H_matter,N(theta)+B_N`. Every Fock state
and both energy signs are retained in W_N.

For occupation n_i at cell c_i, put
`f_xi=1_(c_i<=x)−(4−c_i)/4`. The exact expression is

```
B=(1/4)sum_x[sum_i f_xi n_i+3−2x]^2.
```

Using n_i squared=n_i gives79 nonzero occupation monomials of degree
at most2. Squaring with the exact union rule for occupation supports
gives794 monomials of degree at most4. Both identities are verified on
all65536 words with exact integer arithmetic, using common denominators
64 and4096. The sector norm of B is its largest retained diagonal value.

At a fixed holonomy, let v be the unique full free N-particle ground
Slater state, with energy e and free Fock gap Delta. Its occupation
moments obey `expect(prod_(i in A)n_i)=det(C_A)`, where C is the complete
occupied one-particle correlation projector. This supplies
`b=<v,Bv>` and `s2=<v,B^2 v>` using determinants of order at most4.
It is an exact identity for this comparison vector, not a Gaussian
assumption about the interacting vacuum.

Split the entire fixed-N Fock space into v and its full orthogonal
complement. On the complement H_matter is at least e+Delta, and B is
positive. The off-diagonal norm is exactly sqrt(s2−b squared). Therefore

```
W_N >= e + max(0,
 [b+Delta−sqrt((b−Delta)^2+4(s2−b^2))]/2).
```

The zero alternative uses W_N>=H_matter. Degenerate free Fermi edges
use that unrefined bound. All such controls and all non-improvements
are retained. No orthogonal Fock state is deleted by the two-block
inequality.

The one-body eigensystem is checked against a1e−10 source comparison
envelope. For a separated Fermi edge, the admitted projector error is
`2 delta/(gap−2 delta)`. The finite-Slater angle bound and the exact
sector norms of B and B squared bound their expectation errors, with
an additional1e−6 moment-arithmetic envelope. The lower b, lower Delta
and upper off-diagonal norm are used in the quadratic root. This
arithmetic scope is explicit; the expectation calculation is not
described as interval diagonalization.

Crucially, the actual matrix W_N has the same holonomy derivative norm
bound as H_matter, since B is independent of theta. The1024-point lower
bounds thus feed the same global Lipschitz/Barta construction as the
frozen source. No derivative of the Slater frame is needed. Every
verification function a=0, ±1/4, ±1/2, ±1, ±2 is retained. N9 follows
from the exact particle-hole symmetry with theta reversal; all other
number sectors retain their complete-source bounds.

The resulting infinite-rotor lower edges are approximately

```
H_7,H_9 >= −89.733931646261,
H_8 >= −93.007054049945.
```

For the frozen zero-mass K4 preparation, the neutral energy-excess bound
becomes0.262866 and the complete charged resolvent gap bound becomes
3.01025685. Reusing its actual stored residuals gives a full-rotor
prepared Green error below3.753e−5. The electron-input case uses only
the source4|m| perturbation bound.

These inequalities tighten the same source. They do not establish a
neutral eigenvalue count or ground-state fidelity. A small spectral
residual alone does not exclude an unobserved lower neutral state.
The exact-vacuum-propagator acceptance fields remain false. The
fractional Wilson member's independently reported additive-mass
obstruction concerns that member; neither this wall source nor its
stronger energy inequality transfers a Wilson chiral-protection claim.
