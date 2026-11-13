# Positive transfer on the same pure-gauge source

## Full link space

Use exactly the source, Haar measure and Gauss range of
`openwave_pure_su2_source`. The electric heat kernel on one link has the
Peter--Weyl expansion

```
k_h(U V^-1)=sum_(n>=0) (n+1) exp[-h alpha n(n+2)/4] chi_n(U V^-1).
```

The series and its smooth derivatives converge for every positive `h`.
Every representation coefficient is strictly positive. The resulting
convolution is the injective positive self-adjoint electric semigroup;
products give all links. Define

```
T_h=exp(-hV/2) exp(-hE) exp(-hV/2).
```

The magnetic multipliers are bounded, positive and boundedly invertible at
fixed graph and step. Hence `T_h` is compact, injective and positive. Its
coordinate kernel is strictly positive. Gauge invariance of both factors
makes the physical range reducing. No finite group replaces SU(2), no
representation coefficient is dropped, and no sign-indefinite logarithm
is used as a surrogate transfer.

For any finite family of positive-time cylinder histories, reflect the
operator order and take adjoints at the reflection plane. The reflected
quadratic form is the squared norm of the sum of their half-history
vectors. A bond plane supplies a middle factor `T_h`, whose positive
square root gives the same norm identity. Thus arbitrary complex linear
combinations have nonnegative reflected form. The unique normalized
positive transfer vacuum supplies the stationary vacuum histories.
This proves time-reflection positivity on the actual full gauge space.

## Common clock

The electric domain and bounded magnetic perturbation satisfy the
Lie--Trotter product hypotheses, including on the reducing Gauss range:

```
strong-limit_(n->infinity) (T_(t/n))^n=exp[-t(E+V)], t>=0.
```

The finite Peter--Weyl polynomial spaces are dense electric form cores.
Their increasing Galerkin compressions converge in closed forms at fixed
graph. This statement takes the representation cutoff to infinity at
fixed spatial graph. It asserts no spatial-volume or lattice-spacing
limit and no physical rescaling of the source time parameter.

Time reflection positivity and this common clock are finite-regulator
instances of the Yang--Mills reconstruction requirements. They do not
supply all Osterwalder--Schrader axioms. In particular the heat-kernel
temporal action is explicitly declared; it is not identified with an
isotropic Wilson plaquette measure at finite step. The construction of
positive gauge transfer matrices is also described by
[Luescher](https://doi.org/10.1007/BF01614090). The argument here uses the
displayed electric heat kernel and bounded magnetic multiplication.

## Complete executed controls

The exact one-plaquette Gauss space is spanned by every character `chi_n`.
For the numerical control alone compress to `0<=n<=12`. Both compressed
electric and magnetic operators are positive. For positive bounded
`E_N,V_N`, Duhamel and contractivity give the conservative estimate

```
||(exp(-hV_N/2)exp(-hE_N)exp(-hV_N/2))^m
       -exp[-mh(E_N+V_N)]|| <= mh^2 ||[E_N,V_N]||/2.
```

Split the symmetric product into two Lie products of half duration and
use the integral commutator bound `s^2||[E_N,V_N]||/2` on each, followed
by telescoping. This also gives the displayed weaker coefficient.
The producer evaluates all transfer eigenvalues, complete 13-outcome
plaquette probabilities and a five-history reflected Gram matrix. Its
vacuum is the actual vacuum of each tested transfer. It does not relabel
different transfer vacua as an identical Hamiltonian ground state.
The finite numerical comparison checks the explicit construction; the
full-space proof supplies reflection positivity and the strong clock limit.

The self-adjoint operator `diag(1,-1/2)` gives a negative bond-reflection
form on its second vector. This control excludes the inference that
self-adjointness alone implies reflection positivity.
