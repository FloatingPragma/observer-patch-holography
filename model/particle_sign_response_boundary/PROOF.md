# What local Wilson weight does not control

The existing [local spectral theorem](../particle_local_wilson_density/PROOF.md)
bounds the state-averaged local weight of small Wilson eigenvalues and the
local Hilbert–Schmidt error of a sign approximation. It does not bound the
electric derivative of the sign function. The actual native gauge family
below exhibits the distinction, with every Dirac and internal column
retained. This is an obstruction to one proposed inference, not a failure
theorem for the configured mixed state or for an interacting record limit.

## 1. A gauge-coordinate family in the declared source

Use the same periodic one-cell coordinate source, its fifteen-dimensional
representation of S(U(3)×U(2)), source Clifford matrices, and spatial
beta=gamma_tuple[3]. The native weak representation consists of four
doublets and seven singlets. The code reads its exact exterior basis and
weak Casimir from the frozen source packet, rather than assigning a new
internal representation.

For `0<delta<1/2`, put `c=(1-delta)/2`, `s=sqrt(1-c²)` and let theta
be the angle with cosine c and positive sine s. Choose the actual weak
group elements on the three coordinate links as

```
U1(u)=cos(theta+u) I+i sin(theta+u) sigma1,
U2(u)=cos(theta-u) I+i sin(theta-u) sigma2,
U3(u)=I.
```

Their color blocks are identity. These are legitimate values of the
existing quantum link coordinates, with unchanged action coefficients.
A coordinate value is not a normalizable quantum gauge-field state.

On each weak doublet write `alpha_i=i beta gamma_i`,
`A_i=alpha_i tensor sigma_i`. A1 and A2 are commuting Hermitian
involutions; each anticommutes with beta. Define

```
P_-=(I-A1 A2)/2,       P_+=(I+A1 A2)/2.
```

Each has dimension four on an eight-dimensional Dirac/doublet block.
Across all four native doublets their ranks are 16 and 16. Keep also the
rank-28 weak-singlet projector P_s. The complete source kernel reduces as

```
H_W(u)|P_- = B(u) beta + D(u) A1,
H_W(u)|P_+ = B(u) beta + F(u) A1,
H_W(u)|P_s = -beta,
B(u)=1-2c cos(u),  D(u)=2c sin(u),  F(u)=2s cos(u).       (1)
```

At u=0 this gives the complete squared spectrum

```
H_W(0)² = delta² P_- + (3+2delta) P_+ + P_s,
ranks                 16                 16      28.   (2)
```

In particular, the small-mode projector is exactly P_-. None of its
multiplicity is a particle census or a rest-mass assertion.

## 2. The electric derivative is inversely proportional to the gap

The sign on the first block is

```
sign(H_W(u))|P_-
 = [B(u) beta+D(u) A1]/sqrt(B(u)²+D(u)²).
```

Since `B(0)=delta`, `B'(0)=0` and `D'(0)=1-delta`, differentiation is
exact and gives

```
d sign(H_W(u))/du |0 = (1-delta) A1 P_- / delta.         (3)
```

The derivative of the large block vanishes at zero; the singlet block
does not vary. Thus (3) is the full sixty-dimensional derivative.
The code checks both the literal native-link kernel and its complete
spectral expression, as well as the two derivative identities obtained
from `sign(H)²=I` and `[H,sign(H)]=0`.

The electric metric is the source metric `-Tr_M15(T_a T_b)`. The native
anti-Hermitian generator i sigma_i acts on four weak doublets and has
squared norm 8. The two varied links are independent, so the u-direction
has squared norm 16. Its unit electric vector field is `X=(1/4)d/du`.
Consequently

```
||X sign(H_W)|| = (1-delta)/(4delta),
||X sign(H_W)||_HS² = (1-delta)²/delta².                (4)
```

These formulas use all sixty source columns. They diverge while the
kernel and its first coordinate derivative remain bounded.

In a local coordinate chart away from singular kernels, electric momentum
is `p_X=-iX`. On smooth compactly supported vectors in that chart,

```
[p_X,M_sign] psi = -i M_(X sign) psi.                  (5)
```

Every invertible witness has an open neighborhood on which the derivative
is large. Hence no bounded global multiplication-commutator estimate can
be inferred from the boundedness of sign itself. Equation(5) states a
domain-qualified first-derivative identity. It does not assert that sign
preserves the entire electric Laplacian operator domain, nor that the
full second-order commutator is defined on every source state.

## 3. Small original weight does not remove the derivative requirement

The sum of squared Hilbert–Schmidt derivatives over an orthonormal
electric frame is gauge invariant: a gauge transformation conjugates the
kernel and orthogonally rotates the Lie-algebra frame. It is at least
the directional cost in (4). This is a scalar source diagnostic and can
be evaluated after a density twirl; no gauge choice removes the large
cost.

Mix a rare neighborhood of (1) with original weight w and a gapped
neighborhood of the identity with weight 1-w. Keep the preparation flag.
At the central coordinate symbols the rare contribution to the low-mode
trace is 16w, whereas its unit-electric derivative cost is

```
w (1-delta)²/delta².                                  (6)
```

Taking w=delta sends the first quantity to zero and the second to
infinity. This demonstrates why an unweighted small-mode count supplies
no upper bound on an inverse-gap-weighted derivative.

There is also a normalizable source interpretation of this counterfamily.
For each fixed positive delta, choose smooth normalized coordinate
wave packets supported inside the corresponding invertible neighborhoods.
The small-mode ranks and a strict lower bound for (4) persist there by
continuity. A mixed density, with the preparation flag and conjugate
reference retained, is gauge invariant after the full site-gauge twirl.
The gauge-invariant summed derivative cost remains large. Each smooth
packet has finite electric energy. No bound uniform in delta is asserted.
These are alternative preparations of the declared coordinate source;
they are not claimed to be the existing link-polynomial density family.

Finite source-word approximations suffice for bounded versions of this
statement. The faithful native matrix words are dense in the coordinate
Hilbert space. Approximate the smooth packet in norm, retain its original
mixture weight, and use a bounded clipped derivative-cost effect or the
finite-pulse experiment below. All complete bounded record probabilities
then converge. The finite word vector has finite Casimir energy and admits
the canonical finite conjugate purification. Approximation cannot be used
to infer convergence of an unbounded derivative moment; the clipped effect
is essential for that comparison. The package supplies no practical word
degree or pulse-compilation error for these alternative packets.

## 4. Complete finite-pulse coordinate controls

Take `delta<=1/4` and `u=arcsin(delta)`, so `0<u<=2delta`. Prepare the
auxiliary Dirac/M15 mode tracially, `I_60/60`. Retain the three block
labels P_-, P_+, P_s and both sign labels at each of two readouts, before
and after the specified gauge pulse. There are twelve original records.
For block b of rank r_b and correlation C_b their probabilities are

```
p(b,a,a') = (r_b/60)(1+a a' C_b)/4,   a,a' in {-1,+1},
C_- = B/sqrt(B²+D²),
C_+ = [delta B+2s F]/[sqrt(3+2delta)sqrt(B²+F²)],
C_s = 1.                                               (7)
```

Every record is nonnegative and their sum is one. The singlet switches
are retained zero labels. The program does not condition on the small
block. Equation(7) is an exact coordinate-symbol control; a point
coordinate is not relabeled as a normalizable preparation.
The fixed Pauli axes, left-translation pulse and block menus use a chosen
gauge frame. They are not asserted to commute with the full Gauss action.
A physical Gauss-preserving version needs a specified compensating
apparatus/reference or a covariant observable construction. This example
does not supply that apparatus; the gauge-invariant summed derivative
diagnostic of section 3 is independent of this extra pulse requirement.

The finite pulse gives a useful lower bound. Here
`delta<=B<=delta+delta²`, `D=(1-delta)delta`, so `D/B>=3/5`.
It follows that `C_-<=5/sqrt34<7/8`, and the conditional small-block
switch probability is at least 1/16. After retaining its original rank
weight and a rare preparation flag of weight w=delta, the original
rare-switch probability is at least `delta/60`. The unit-electric pulse
length is `4u<=8delta`; dividing that switch probability by the squared
pulse length gives a lower bound `1/(3840delta)`. This divergence is a
susceptibility statement. The original unconditioned switch probability
still tends to zero, and all probabilities remain bounded by one.

For each fixed delta, the complete pulse experiment can be approximated
on normalizable finite word states and a sufficiently large finite code
in the full coordinate carrier, without imposing an additional Gauss
restriction on its apparatus.
Left translations of the links preserve word degree. The two sign menus
are bounded coordinate instruments, whose complete finite compressions
converge strongly as the code grows. A strict probability margin allows
their total error to be chosen smaller than, for example, delta/120,
without discarding any flag or detector outcome. This existence argument
does not claim that the displayed coordinate law was executed on such a
large finite code, supply its resource modulus or discharge the gauge
constraint on the pulse apparatus.

## 5. The sufficient spectral input and the remaining task

For an invertible finite Hermitian H with eigenvectors a and eigenvalues
lambda_a, the Fréchet derivative of sign in a Hermitian direction V is
entrywise

```
[D sign_H(V)]_ab = 0                         (same-sign eigenvalues),
                 2 V_ab/(|lambda_a|+|lambda_b|) (opposite signs).
```

In particular its squared Hilbert–Schmidt norm is exactly

```
sum_(lambda_a lambda_b<0)
        4 |V_ab|²/(|lambda_a|+|lambda_b|)².              (8)
```

The quotient follows by differentiating the spectral functional calculus,
or by the involution and commutation tangent equations used in the code.
A state-weighted bound on these inverse-gap matrix elements is a genuine
additional premise. A uniform gap is sufficient but stronger than needed.
Spatially local response requires the corresponding local and electric
estimates, with volume-independent constants and the actual quantum state.

A smooth regulator offers a different declared route. For
`f_mu(x)=x/sqrt(x²+mu²)`, divided differences have absolute value at most
`1/mu`, since `0<=f_mu'(x)<=1/mu`. Thus the finite Hilbert–Schmidt bound
`||D f_mu,H(V)||_HS<=||V||_HS/mu` holds without a gap. This statement alone
is not an operator-norm bound or a full interacting time-evolution result.
Taking mu to zero again needs a controlled local comparison on the same
state. A polynomial action already has bounded coordinate derivatives at
each fixed degree, with constants that must be tracked under refinement.

The counterfamily does not invalidate the proved local spectral-weight
bound, the polynomial thermodynamic record limit, or any particular
configured preparation. A local Duhamel proof might avoid differentiating
sign directly. The exact remaining obligation is a volume-uniform local
record-response estimate for the full electric/link/matter dynamics,
including state and reference errors. Neither a particle mass nor failure
of that obligation is concluded here.
