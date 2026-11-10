# Scalar transition poles on the full physical local ground

## Source and exact scalar identities

The source is the thirty-mode native fermion Fock space tensor the
five-state quantum scalar. Its local Hamiltonian is `H=N+Y`, where
N is the actual scalar excitation projector and Y is the sum of all
three native Yukawa channels. Every coefficient is one. The physical
domain has all twelve gauge generators zero and no compensating
boundary charge. The clock is `exp(-itH)`, with hbar one.

The [complete local ground theorem](../local_reference/PROOF.md) gives
the global physical ground projector P_g of rank two, raw L=1,
B3=0 and physical spin one-half. The reference is `rho=P_g/2`.
Its selection is the declared preparation law; no arbitrary-input
ground-filter probability is supplied.

The actual full normal-ordered source polynomials satisfy

`N^2=N`, `{N,Y}=Y`, `Y^2=H^2-H`, `[N,Y^2]=0`.

The first two identities follow from the five-state scalar: Y joins
only the scalar blank and excited parts. The code verifies them on
the full thirty-mode polynomial algebra, with the actual Yukawa
coefficients and signs. It also verifies their exact reduced matrix
versions. N and Y are even, commute with all twelve gauge generators,
carry physical spin zero, and preserve raw L and B3. These statements
apply to the total Y. Separate responses for the three individual
Yukawa channels are not inferred.

The actual scalar current is `C_N=i[H,N]`. It is self-adjoint in the
physical metric and satisfies

`C_N^2=Y^2`, `i[H,Y]=-C_N`.

The receipt includes the complete rational energy, metric, scalar,
Yukawa and current matrices. The source includes both magnetic
partners of the full spin-one-half maximal-quark block and has
dimension thirty-two. N and total Y preserve that complete block;
each full operator image is checked through its invariant inclusion.

## Exact two-level response

Let e be the global physical ground energy and put

`xi=e(e-1)`, `omega=1-2e=sqrt(1+4xi)>0`.

The selected ground obeys `Y^2 P_g=xi P_g`. Since N and Y commute
with Y squared, both probe images lie in the eigenspace of
`H^2-H` with eigenvalue xi. Its only energies are e and 1-e.
The complete source spectrum contains each of these with rank two,
so the ground and its positive-energy partner exhaust this probe
response. The other seven energy outcomes have exactly zero weight.

Expectation of `{H,N}=N+H` on the stationary reference gives

`n=<N>=e/(2e-1)`.

The code independently checks this mean using the exact polynomial
ground projector and all physical matrix traces. The scalar variance
and the total-Yukawa mean are

`v=<N^2>-<N>^2=xi/(1+4xi)`,
`<Y>=e-n=-2xi/omega`.

Thus the unscaled elastic and inelastic weights are

| Operator | Elastic weight | Inelastic weight |
|---|---|---|
| N | `n^2` | `v` |
| Y | `4xi^2/(1+4xi)` | `v` |
| Scalar current C_N | `0` | `xi` |

The Y and N inelastic matrix elements have opposite signs because
`Y=H-N`. The current identity supplies the factor omega between
the current and scalar inelastic matrix elements. The stationary
current mean is exactly zero.

With the convention

`chi_AB(t)=-i theta(t) Tr(rho [A(t),B])`,

and Fourier-Laplace transform `integral_0^infinity exp(i z t) chi(t) dt`,
the exact raw retarded responses for `Im(z)>0` are

`chi_NN(z)=chi_YY(z)=2 omega v/(z^2-omega^2)`,
`chi_NY(z)=-chi_NN(z)`.

The poles are at plus and minus omega, with opposite residues v.
The elastic energy record at zero excitation contributes no
zero-frequency commutator pole. Numerically omega is about 23.4612
in the declared source units. This is a finite-clock scalar transition;
it is not identified with a relativistic rest mass.

## Complete spectral records

Writing `x=E(E-1)`, the squarefree source minimal polynomial is

`p(E)=(E-1)(x-105)(x^3-313x^2+31708x-1041468)`.

The full matrix annihilation identity is checked. Its six cubic-derived
roots have rank two each; the two roots of x=105 also have rank two,
and E=1 has rank sixteen. These ranks sum to thirty-two. The source
ground projector is evaluated from

`P_e=[p(H)/(H-e)]/p'(e)`

using the exact polynomial quotient, reduced modulo the ground energy
factor. Its trace is exactly two. Exact matrix traces check the means
of N, Y, Y squared and the current. Rational root isolation and
outward dyadic rounding supply every energy and probability interval.
No numerical eigenbasis enters the certificate.

Two recorded settings, each with probability one-half, use

`K_N=N/2`, `K_Y=Y/16`.

The joint elastic/inelastic probabilities are respectively

`(n^2/8, v/8)` and `(4xi^2/[512(1+4xi)], v/512)`.

Their sums are `n/8` and `xi/512`. The spin-degenerate elastic
record and the positive-energy partner each retain rank two. Every
other one of the nine energy outcomes is explicitly retained with
probability zero. The complement of the complete source block is
also retained, with zero probability for this preparation. Failure
has no energy postselection. The full alphabet has twenty-two joint
records: six positive-probability records and sixteen exact zeros.

## Full-carrier even apparatus and failures

The number probe has positive failure operator

`F_N=I-N+(sqrt(3)/2)N`.

For the Yukawa probe define the full positive finite-matrix square root

`F_Y=sqrt(I-Y^2/256)`.

This is defined on the entire native carrier. The global quark
character bound C_q at most 54, together with the complete lepton
squared-transfer bound six, gives

`||Y||^2 <= (sqrt(108)+sqrt(6))^2=114+36sqrt(2)<165<256`.

The quark bound includes every color representation and occupation;
no physical-singlet restriction is used for this instrument norm.
Both settings therefore satisfy the exact full-carrier identity

`K^dagger K+F^dagger F=I`.

K and F are commuting self-adjoint functions of the same even neutral
operator, so the pointer dilation

`U=[[F,-K],[K,F]]`

is unitary on both sides and preserves gauge charge, raw charges,
physical spin and fermion parity. The setting and outcome registers
are even, neutral finite apparatus registers. The odd source parity
is retained. The controller is an intervention, rather than part
of the idle Hamiltonian.

The scalar-number failure is applied as its full polynomial. The
Yukawa failure acts on the selected ground as
`F_Y P_g=sqrt(1-xi/256)P_g`, a consequence of the full source identity
`Y^2P_g=xiP_g`. It is not defined by compression to that ground.
The joint failure probabilities are `1/2-n/8` and `1/2-xi/512`.
Together with all successful and complementary outcomes they sum
to one exactly.

The full Yukawa square root has a deterministic uniform approximation.
For `Z=Y^2/256`, `0<=Z<rI`, `r=165/256`, set

`F_m=I-sum_(k=1)^m a_k Z^k`,

where `a_1=1/2` and `a_(k+1)=a_k(2k-1)/(2k+2)`. The positive
decreasing coefficients give

`||F_m-F_Y|| <= a_(m+1) r^(m+1)/(1-r)`.

The receipt gives exact rational errors for finite polynomial orders.
Since both F_m and F_Y are positive contractions, replacing F_Y in
the pointer matrix has unitary defect at most twice that error.
A finite polynomial is not claimed to be the exact instrument.
The exact instrument uses the unique positive square root of a
finite rational matrix; the uniform series supplies an explicit
evaluation algorithm with a full-carrier error bound.

## Verification and scope

Hostile tests reject a false scalar projector, a changed scalar
penalty, invalid square-root orders and missing custody inputs.
An independent generalized eigensystem compares every weight,
including the exact zeros, under the declared diagnostic tolerance
`2*10^(-10)`. That comparison does not widen a certified interval.

The theorem is a complete scalar response on the selected physical
local ground with the actual native clock. It does not supply a
many-cell stationary vacuum, momentum-dependent propagation, an
individual-Yukawa-channel response, or a measured particle
identification. All declared inputs, the inherited ground receipts
and the repository-local import closure are pinned.
