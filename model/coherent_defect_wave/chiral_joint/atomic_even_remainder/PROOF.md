# An exact even effective Hamiltonian and full-complement quartic bound

## Actual source and physical reference

The source is the native periodic two-cell-per-axis cube with all
thirty fermion modes and five quantum scalar states at each of eight
cells, and the full 119-state Peter-Weyl link on each of twenty-four
oriented edges. The onsite Yukawa coefficients, scalar penalty and
link penalty are one. The real hopping coefficient lambda is explicit:

`H_lambda=H0+lambda T`.

The [atomic reference theorem](../atomic_reference_control/PROOF.md)
proves the complete physical atomic ground projector

`P=P0=tensor_x P_ground,x tensor all_links_blank`,

of rank 256 and energy `E_at=8E0`. It covers nonblank links and their
charged endpoint sectors. Let `Q=I_physical-P`, with the identity on
the entire Gauss-constrained carrier. The theorem supplies a positive
rational lower bound g for the full atomic gap and a complete-carrier
hopping norm bound `||T||<=tau=360`. It also proves `PTP=0`.
This certificate recomputes the unrestricted local bound, lattice
isolation and full transporter norm controls used by that theorem.

Write `v=|lambda|tau` and assume `v<=g/5`. With the atomic energy
subtracted, the complete physical block matrix is

`H_lambda-E_at=[[0,B^dagger],[B,D]]`,

`B=lambda QTP`, `D=D0+Delta`,

`D0=Q(H0-E_at)Q>=gI`, `Delta=lambda QTQ`.

Thus `||B||<=v`, `||Delta||<=v` and `D>=gamma I`, where
`gamma=g-v>0`. All of Q is retained in these operators.

## Exact parity on the native carrier

Let A be the four cube cells whose coordinate sum is even and define

`U_A=(-1)^(sum_(x in A) N_fermion,x)`.

Every native local Yukawa word has even fermion degree. The scalar
and electric terms are also even. All local gauge-current words
have even fermion degree, so U_A preserves the physical Gauss domain.
The code checks these statements on the actual source words.

Every one of the twenty-four actual hopping terms joins opposite
sublattices and contains one creation and one annihilation operator.
Its Wilson link is bosonic. Consequently

`U_A H0 U_A=H0`, `U_A T U_A=-T`,

`U_A H_lambda U_A=H_(-lambda)`.

Every local physical ground vector has odd fermion parity. Four such
cells in A give `U_A P=P`. This scalar action on the entire reference
band is necessary for an even effective matrix, rather than merely
an even spectrum. U_A therefore preserves Q; write its restriction
there as U_Q. In particular,

`U_Q D(lambda) U_Q=D(-lambda)`, `U_Q B(lambda)=-B(lambda)`.

The larger spectral source does not obey this same hypothesis.
At side four its actual displacement-two coefficient is nonzero,
and its endpoints have the same coordinate parity. The certificate
records this source counterexample and the test rejects transferring
the bipartite argument to that lattice.

## Canonical graph symmetry

The atomic reference theorem constructs the unique invariant graph
in the norm ball `||X||<=r=2v/gamma<=1/2` by

`X=-D^(-1)B+D^(-1)X B^dagger X`.

Its contraction factor is at most `4v^2/gamma^2<=1/4`.
The full-source parity maps a solution at lambda to `U_Q X(lambda)`
at minus lambda, in the same norm ball. Uniqueness gives

`X(-lambda)=U_Q X(lambda)`.

The canonical isometry and effective Hamiltonian are

`J_lambda=(I;X)(I+X^dagger X)^(-1/2)`,
`K_lambda=J_lambda^dagger(H_lambda-E_at)J_lambda`.

Since `X(-lambda)^dagger X(-lambda)=X(lambda)^dagger X(lambda)`,
the normalization agrees at both signs. It follows that

`J_(-lambda)=U_A J_lambda`, `K_(-lambda)=K_lambda`.

This is an operator equality in the fixed canonical reference
coordinates. The same parity shows that the complete self-energy
`B^dagger D^(-1) B` is even in lambda.

## A controlled fourth-order bound

The graph equation first gives the exact norm estimate

`||X+D^(-1)B|| <= (v/gamma)r^2`.

Put `L=B^dagger X` and `M=I+X^dagger X`. Graph invariance gives
`K=M^(1/2) L M^(-1/2)`. The spectrum of M lies in `[1,1+r^2]`,
so

`||K-L|| <= [sqrt(1+r^2)-1/sqrt(1+r^2)] ||L||`

`         <= r^2 ||L|| <= v r^3`.

Combining the graph and normalization errors yields

`||K+B^dagger D^(-1)B||`

` <= v^2 r^2/gamma + v r^3 <=12v^4/gamma^3`.

The two contributions are bounded separately by `4v^4/gamma^3`
and `8v^4/gamma^3`. This estimate retains the exact D with all
charged and neutral complement states.

For the second-order term let

`K2=-PTQ D0^(-1) QTP`.

Set `C=D0^(-1/2) Delta D0^(-1/2)`, so C is self-adjoint and
`||C||<=v/g<1`. The exact finite-operator resolvent identity is

`[(D0+Delta)^(-1)+(D0-Delta)^(-1)]/2-D0^(-1)`

` =D0^(-1/2) C^2(I-C^2)^(-1) D0^(-1/2)`.

It is positive and has norm at most

`v^2/[g(g^2-v^2)]`.

The even self-energy equals the average of its values at the two
signs. Multiplying the preceding bound by the two B factors therefore
gives

`||B^dagger D^(-1)B-B^dagger D0^(-1)B||`

` <=v^4/[g(g^2-v^2)]`.

No formal series or finite complement truncation is used. Combining
the exact resolvent estimate with the canonical graph correction gives

`||K_lambda-lambda^2 K2||`

` <=R4=v^4[1/(g(g^2-v^2))+12/(g-v)^3]`.

For `u=v/g` in `[0,1/5]`, the dimensionless coefficient

`1/(1-u^2)+12/(1-u)^3`

has positive derivative and endpoint value `1175/48<25`. Hence

`||K_lambda-lambda^2 K2|| <=25v^4/g^3`.

The receipt evaluates the exact sharp and uniform bounds at zero,
both admitted endpoints and a smaller nonzero coefficient. Negative
coupling uses the same bound. The original full hopping coefficient
one is outside this small-coupling inference.

## Stable gap-transfer API

`uniform_quartic_remainder(gap,tau,coupling)` returns the exact rational
bound `25(|lambda|tau)^4/g^3`, rejecting nonrational inputs and values
outside `v<=g/5`. Its hypotheses include the proved scalar reference
parity and complete physical atomic gap; supplying numbers alone does
not establish those source hypotheses.

If a consumer separately proves a positive second-order spectral gap
delta2, `choose_quartic_coupling(gap,tau,delta2)` chooses the largest
coefficient of the form `2^(-k)` satisfying

`|lambda|tau<=g/5`,

`lambda^2<=g^3 delta2/(100 tau^4)`.

It therefore guarantees

`R4<=lambda^2 delta2/4`,

`lambda^2 delta2-2R4>=lambda^2 delta2/2>0`.

Weyl's finite self-adjoint eigenvalue bound can then transfer that
separately established gap. The parity theorem alone does not
identify a unique interacting ground or a second-order spin spectrum.

## Verification and domain

Exact rational finite-block tests verify the symmetric resolvent
identity, its positive norm bound and the cancellation with the
actual parity pattern. Hostile tests reject a same-sublattice source
term, a nonscalar parity action on the reference, an excessive
hopping coefficient and invalid gap inputs. A separate generalized
eigensystem diagnostic checks canonical graph evenness and the
quartic bound in a noncommuting finite block, with absolute tolerance
`2*10^(-12)`. That diagnostic supplies no certified source input.

All inherited atomic bounds and source parity checks are recomputed,
and declared inputs plus their repository-local import closure are
hashed with missing files rejected. The theorem applies to the
specified two-cell-per-axis Peter-Weyl source, its complete physical
complement and its declared hopping interval. No unit-hopping,
larger-lattice or thermodynamic statement is transferred.
