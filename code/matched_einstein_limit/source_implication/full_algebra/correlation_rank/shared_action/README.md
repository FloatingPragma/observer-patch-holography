# A regional action meter must resolve a three-register correlation

The selected 80-cell preparation has three accessible registers: the left
reversal flag and two head coins. There is a two-sided family through this
preparation with exactly the same **three pairwise reduced densities** and
exactly the same four-station radar volume. At the reference its full entropy
increases by **0.006092217 to 0.006093092 nats per unit Hilbert–Schmidt
variation**. The independently implemented Pauli measurement that multiplies
both coin outcomes by the recorded reversal sign detects this change.
The example tests a proposed regional work or stress meter on one stationary
finite survey. [Source, ordinary gain, and reachable preparations](../../PROOF.md)
· [Native analyzer settings](../README.md)

## A signed preparation with identical pairwise records

Write the density in the native reversal blocks as
$X= X_e\oplus X_o$. The nine operators
$U_{ij}=\sigma_i\otimes\sigma_j/2$ are orthonormal and have zero partial
traces. On their span let $a=A_e-A_o$ and $v=V_e-V_o$ be the differences of
the two sectors' area and volume gradients. These are certified derivatives
of the independently timed twelve directed round trips with the six-unit
deadline, evaluated at the exact source. The same interval calculation
that certified the product-term count supplies both matrices.
[Source coefficient intervals](../results.json)

Inside the nine-dimensional affine slice $(X_e+D,X_o-D)$,
every trial that reads any two registers has exactly the same probability
law. Forgetting the tag gives the constant density $X_e+X_o$. Each coin
conditioned on either tag retains its marginal since both partial traces
of $D$ vanish. Thus the reversal weight, parked exchange energy and the
mean of every sum of at-most-two-register work terms are also constant.
This covers arbitrary independent repetitions on any of the three pairs.
A multi-step interaction accessing all three registers can create a
three-register effect even when each individual gate acts on a pair.

Set $w=a-v(a\cdot v)/\|v\|^2$. The certificate gives
$0.001347052\leq\|v\|\leq0.001347077$; the implicit function theorem
therefore gives an exactly volume-preserving two-sided curve in this slice
whose initial unit tangent is $(-w,+w)/(\sqrt{2}\|w\|)$. Both source blocks
and the local prepared-state neighborhood are interior at the reference.
Along this tangent the independently calibrated finite stationarity theorem
gives $dS=\eta\|w\|/\sqrt2$, positive within the reported enclosure. This
is the largest absolute entropy slope of any unit tangent in this blind
fixed-volume slice. For stable outward bounds the certificate evaluates
$\|w\|^2=\sum_{i<j}(a_i v_j-a_j v_i)^2/\|v\|^2$.

The directly available three-register observable
$Z_{\rm reversal}\,Y_{\rm left}Z_{\rm right}/2$ has signed slope in
$[1.106899536,1.107432379]$ along this *same* tangent. This setting is
chosen from the established native Pauli menu after analysis. It reads a
hidden correlation; the calculation supplies no physical unit of stress
for this single observable. As an independent numerical control the test
diagonalizes the recorded source center and evaluates $-\log X$ on the
tangent, reproducing the signed slope without the area law. The rigorous
enclosures above use the exact-source radar and ordinary-gain certificates.
[Ordinary calibration](../../certificate.json) · [Executable checks](test_shared_action.py)

## Allowing the same bilinear coupling in independently rotated bases

A more general meter can condition its two-coin interrogation on the tag.
Suppose its two sector work operators have Pauli correlation coefficients
$c_b O_b J P_b^T$, with an arbitrary common real $3\times3$ matrix $J$,
independent real gains $c_b$, and independent orthogonal $O_b,P_b$.
We allow arbitrary additional tagged one-coin and scalar terms. Even this
class misses the entropy derivative on some unit fixed-volume source tangent
by more than $8.16\times10^{-11}$ nats. The bound is conservative. It
excludes exact equality for this class at the specified source.

Indeed, matching on the joint fixed-volume tangents would require both
$C_b=A_b-tV_b$ to have proportional singular values at one common real
volume price $t$. Overall signs and the gain can be absorbed into the work
operator. For a nonzero matrix $C$ use the two invariants
$F=\|\operatorname{cof}C\|_F/\|C\|_F^2$ and
$G=|\det C|/\|C\|_F^3$. Acting on a unit-Frobenius matrix they have
Lipschitz constants 2 and 3 respectively: the second and third exterior
powers are compressions of the tensor powers to normalized exterior bases;
telescoping the tensor powers and applying the triangle inequality gives
the bounds. Both invariants are unchanged by left/right rotations and any
nonzero gain.

Put $h=\max(|F_e-F_o|/2,|G_e-G_o|/3)$. If a pair of common-shape
approximants has errors $d_e,d_o$ in Frobenius norm, normalizing each
block changes its error to at most $2d_b/\|C_b\|_F$. The triangle
inequality thus gives
$\sqrt{d_e^2+d_o^2}\ge h\min_b\|C_b\|_F/(2\sqrt{2})$.
For a zero approximant the error in its block is already its full norm;
$h\le 2$ suffices. Directed interval arithmetic verifies
$h\min_b\|C_b\|_F>10^{-11}$ on all four dyadic boxes used for
$A-tV$, $|t|\le1$, and all 159 used for $sA-V$, $|s|\le1$. For
$t=1/s$, dividing the latter matrix by $s$ only increases this bound.
Taking the infimum over $t$ projects out the single volume constraint.
Multiplying by the positive ordinary gain produces the stated entropy
error. A control with differently scaled, rotated copies of an identical
coupling correctly fails the exclusion.

The comparison constrains the *work operator at the reference*. A metric
variation of a Hamiltonian can introduce additional operators through
geometry-dependent bases, couplings or multi-step native dynamics. A
regional source mechanism remains possible if its independently derived
stress includes reversal-sensitive joint terms outside the tested class.
The calculation establishes a first variation on the chosen radar
tetrahedron. Ten timelike frames, shrinking proper balls, Ward conservation
and a nonlinear Einstein solution require further physical inputs.

Run `python certify.py` for the interval certificate and
`pytest -q test_shared_action.py` for the source-spectrum and shared-shape
controls. [Stored enclosure](results.json)
