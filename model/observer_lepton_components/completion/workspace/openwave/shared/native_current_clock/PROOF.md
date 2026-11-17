# Current queries and their complete walk records

## Source operation and query complement

The carrier and its CAR ordering map are exactly those of `spin_gauge_join`.
On one seam the native source operation U is a signed involution, with
whole-event holds at the electric-window boundary. Let Q be the diagonal
realized number current. Its spectrum is{-1,0,1}, U Q U=-Q, and the moving
support Q squared is invariant under U. The Hermitian operator J=i U Q
satisfies J squared=Q squared. Set P=I-J squared. Then J P=P J=0 and P
is a projection. The two-query-bit operation

    C = [[J,P],[P,-J]]

is Hermitian and C squared=I. On a moving column it performs the actual
source hop with phase i*j*CARsign*(-1)^query. On a hold or overflow it
flips the query bit without moving matter, memory or links. These columns
specify the operation on every matrix unit and every untouched reference.
There is no normalization or deletion of the query-one output.

The finite source slice retains every27 electric assignment and every six
memory values, not just Gauss-zero or a selected character. It is closed
under each of the three seams. All other source sectors are in an explicit
orthogonal invalid flag. This flag is the identity under the full walk and
is a failure when read; it is not accepted as an electron state.

## SELECT, reflection and the source spectrum

The recovered response frame gives equal seam lengths ell and unit axes n_e.
The existing port-spin attachment assigns the Hermitian involution sigma.n_e.
Define SELECT as the direct sum of (sigma.n_e) C_e on the three valid
seam addresses, and identity on the fourth padded address. On the invalid
source flag SELECT is identity. The spectator's independent spin factor
tensors identity throughout. SELECT is a Hermitian involution.

Let T embed the signal space by adjoining uniform3 in the seam register and
query0. On the admitted source, R=2TT*−I; on the invalid flag R=I. Then

    W=SELECT R,   W*W=I,
    T*SELECT T = A = (1/3) sum_e (sigma.n_e) J_e = D/(3 ell).

The source supplies the actual J_e, including current sign, CAR sign,
active memory and electric-window modifications. The choice of coherent
reflection and spin-current query is declared observer architecture. This
does not select D, D squared or arccos(A) as physical energy. The same finite
source also has its established scalar Dirichlet description, with the
physical attachment hypotheses stated in the Lean modules.

For an eigenvector v of A with eigenvalue lambda, write u=Tv. When
beta=sqrt(1-lambda squared)>0, define w=(SELECT u-lambda u)/beta. The
vectors u,w are orthonormal and

    W u=lambda u+beta w,   W w=−beta u+lambda w.

Thus (u−i w)/sqrt2 has walk phase lambda+i beta. At lambda=+1 or−1,
u itself is an eigenstate. The endpoint cases remain in the algorithm.
In particular T*W^n T equals the Chebyshev polynomial T_n(A) for every
nonnegative integer n. This is a consequence of the complete walk, including
its complement, not an evolution installed from a fitted eigenvalue list.

## Complete reference observation

For a normalized initial state psi, controlled W^n and a retained reference
phase z in{1,i} give full terminal arrays

    p_plus(x)=|(psi(x)+z (W^n psi)(x))/2| squared,
    p_minus(x)=|(psi(x)−z (W^n psi)(x))/2| squared.

Every event, query, spin, source and invalid coordinate x is kept. Their
sum is one, and the plus-minus difference is Re(z <psi,W^n psi>). Query-one
and event-padding coordinates are terminal outcomes. The query register
is coherently retained during the walk, so committing it after each tick
would define another instrument. The reference prepares and controls an
experiment; it is not a self-starting physical electron clock.

The native M6 source words execute uniform preparation, the correct-sign
eight-address reflection, the individual spin gates and their full event-
controlled multiplex. The exact projector expansion on the two event bits
compiles the latter, including its padding identity. Their source-word
errors are numerical diagnostics. The separate interval audit certifies
every emitted primitive and the exact geometric targets, with corrected
compiled-piece half-diamond bound1.815e-7. Its geometric bound exceeds the
original conditional1e-12 allowance; the original receipt is preserved and
the corrected bound is separately named. Integer current evaluation,
coherent source/tag/link arithmetic,
oracle inversion, control and final pointer chart remain declared. A
full M6 decomposition of those larger oracle services is not claimed.

## Electric continuity with coherent selector terms

Let E_f be an actual integer link coordinate and let Q_f be its realized
number current on the source. On a source hop, Delta E_f=−q Q_f; holds and
overflow have Q_f=0. Within each query block, moving and hold source domains
are orthogonal, and the spin gate acts on a separate factor. Therefore

    SELECT* E_f SELECT − E_f
      = −q (|f><f|_event tensor I_query,spin tensor Q_f) = −q D_f.

R commutes with E_f, since it changes only event/query coordinates. Hence
the full walk identity is

    W* E_f W − E_f = −q R* D_f R.                         (1)

The right-hand side includes all coherent selector cross terms. It is not
obtained by replacing a coherent seam superposition by classical weights.
It can be read from the realized source-current flag immediately before
SELECT, after applying R, using fresh copies for incompatible observations.
For a stationary density rho with finite electric expectation, (1) implies

    q Tr(rho R*D_f R)=0

on every link. In particular every finite-window stationary q=1 state has
zero mean realized charge transport on each seam. The current used here is
the actual transport diagonal; it is not J_f=iU_fQ_f, a coherent-current
quadrature. Nor is it the existing oriented two-path area response or a
coherent loop quadrature. Those observables can remain nonzero. This result
does not assert that all circulation or magnetic observables vanish on
all OPH candidates.

All Gauss sectors are preserved by SELECT, and R commutes with the Gauss
operators. The numerical stationary diagnostic restricts only its prepared
signal eigenstates to Gauss zero and uniform memory character0. The full
query and walk certificates continue to retain all other sectors. Every
eigenvalue in that fixed diagnostic subspace is reported.

The exact constraint identifies what this electric-only construction still
needs for a stationary sourced magnetic field: a source-derived divergence-
free field response that can balance the matter contribution in (1). Merely
substituting the Pauli square for this missing dynamics would not establish
it. `PauliMinimalCoupling.sigmaDot_mul_self` is the exact algebraic precursor;
its physical commutation hypothesis still needs an independently measured
field and the same clock. No universal Maxwell-only angular-momentum gate
is introduced for the already conserved CAR field spin.
