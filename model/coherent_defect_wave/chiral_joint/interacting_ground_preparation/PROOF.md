# A finite algebraic preparation and its complete probability bounds

## Declared preparation service

This construction uses the same full side-two Hamiltonian as
[collective_spin_response](../collective_spin_response/PROOF.md),
`H_lambda=H0+2^-34 T`. The ground gap, graph isometry and trial-state
bounds refer to this Hamiltonian and its physical Gauss subspace. There
are two full source copies. All scalar, link and fermion occupations,
both source copies, spin references and filter registers are retained.

The controller admits exact finite algebraic state synthesis, finite
matrix inversion of the specified source, and controlled applications of
the resulting unitary. These are declared architecture services. No
efficient native gate compilation or numerical execution on the full
source carrier is claimed. The construction does not require the exact
interacting ground vector or its absolute energy as input.

## Actual local ground and an explicit trial controller

The [403-class neutral construction](../neutral_local_spectrum/PROOF.md)
gives an invariant intertwiner from 403 quotient coordinates to all
113948 basis states in the actual physical seed component. Its normalized
seed is `psi_up=((2048,1)-(1024,2))/sqrt(2)`, its full-H minimal
polynomial p_H has degree 49, and its least visible energy E has positive
seed weight. It is the least root of

`[E(E-1)]^3-313[E(E-1)]^2+31708 E(E-1)-1041468=0`.

The [global physical local theorem](../local_reference/PROOF.md) proves
that the same simple root in the magnetic-up sector is the actual
physical local ground, and that its complete magnetic multiplicity is
two. The certificate checks the exact factor, simple-root property,
matching energy isolators, dimensions, positive overlap and source pins.

An explicit finite algebraic ground-vector recipe is

`q_E(z)=p_H(z)/(z-E)`,
`g_up=q_E(H_local) psi_up/||q_E(H_local) psi_up||`,
`g_down=S_minus g_up`.

The denominator is nonzero because the seed's ground weight is positive.
Spectral decomposition of the finite invariant cyclic subspace shows
that q_E annihilates every other energy. The actual local spin
commutators and spin-half seed representation fix the normalized down
partner. The complete native source inclusion of this vector is inherited
from the 403-class intertwiner; no bare one-fermion state replaces it.

At each site retain two gauge-neutral spin-half fermionic reference modes
of raw L=-1 and zero reference Hamiltonian, with their full four-state
Fock space. The local source/reference singlet

`s=(g_up ref_down-g_down ref_up)/sqrt(2)`

and the full local vacuum both have trivial Gauss charges, total spin
zero, total raw L zero and even total parity. The unitary

`I-|vac><vac|-|s><s|+|vac><s|+|s><vac|`

swaps them and acts as identity on the complete complement. Apply it at
all eight sites with links blank. The resulting joint state Psi is the
product of eight singlets. In the source ground-doublet coordinates,
`Psi=(1/16) sum_j |j> tensor J|j>`, where J is the tensor product of
spin-half time reversal on the references. These reference systems are
retained rather than traced out by the controller.

Let u be the normalized native-coordinate cube singlet: the unique cube
ground at the specified least cubic root, transformed by the actual
bond-axis unitary U and embedded through the local ground doublets.
The cube spectral projector gives u by finite polynomial algebra and a
positive normalization. Its phase has no effect on the target
`eta=u tensor J u`. The exact overlap is

`<Psi,eta>=1/16`.

Set `w=(Psi-eta)/sqrt(15/8)`. The complete Householder unitary
`R=I-2|w><w|` obeys `R Psi=eta`. It preserves Gauss charges, total raw
L and total parity because both vectors have the same labels. It also
preserves joint translations and proper cubic rotations: Psi is invariant;
the unique cube ground transforms in a one-dimensional source symmetry
representation, and its time-reversed reference has the conjugate
character. Source translations of u have momentum zero. Only the
rotations of the actual cube construction enter this argument; a
continuous rotated spin symmetry of H_lambda is unnecessary.

The source factor is the pure trial u. The full-source graph and ground
projector bound gives a phase choice with

`||g_lambda-u||<=eta_g`,
`p=|<g_lambda,u>|^2 >= p0=(1-eta_g^2/2)^2>1/2`.

At the stated coupling eta_g is about 0.0001693193, so p0 exceeds
0.99999997. The trial uses exact local and cube algebra, not the unknown
full interacting eigenvector. Prepare two copies by the same controller.

## Finite energy agreement without a ground-energy input

On the complete source carrier, the triangle inequality gives

`||H_lambda||<=B=8*57+24+lambda*360`.

The actual local norm bound 57 includes all native Yukawa terms and the
scalar penalty; all twenty-four link penalties are present. With t=1/B,
define the algebraic Cayley unitary

`C=(I-i H_lambda/B)(I+i H_lambda/B)^(-1)`,
`W=C tensor C^dagger`.

The inverse exists on the entire carrier because H_lambda is Hermitian.
For two source energies a,b in [-B,B], their Cayley values obey

`|C(a)-C(b)| = (2|a-b|/B)/sqrt((1+a^2/B^2)(1+b^2/B^2)) >= |a-b|/B`.

Let Delta be the certified full physical ground gap. For every mixed
ground/excited pair, the corresponding eigenvalue z of W satisfies
`|1-z|>=Delta/B`. For an explicit positive integer M use

`F_M(W)=(1/M) sum_(k=0)^(M-1) W^k`.

Every eigenvalue of this operator has modulus at most one, and a mixed
ground/excited amplitude is bounded by `2B/(M Delta)`. The certificate
chooses `M=ceil(2B/(epsilon Delta))` with epsilon=1/1024 and stores the
exact resulting integer and the tighter actual leakage bound.

A complete implementation prepares an M-state ancilla uniformly,
applies controlled W^k, reverses the uniform-state preparation, and
measures every ancilla basis state. Ancilla zero realizes F_M. All
nonzero outcomes are retained, either individually or under a failure
flag. An exact algebraic uniform preparation exists by the same
Householder construction on the finite ancilla. There is no invocation
of a ground-energy spectral projector in this filter.

Write the trial as `u=sqrt(p) g+sqrt(1-p) e`, with e orthogonal to g.
The two-copy input has mutually orthogonal ground/ground, mixed, and
excited/excited sectors. F_M preserves those sectors, acts identically
on ground/ground, has norm at most epsilon on each mixed sector, and
has norm at most one on excited/excited. Therefore the success
probability q satisfies `q>=p^2`, and the conditional ground-pair
fidelity satisfies

`F >= p^2/[p^2+(1-p)^2+2 epsilon^2 p(1-p)]`.

This expression increases with p on [0,1], since its reciprocal is
`1+((1-p)/p)^2+2 epsilon^2 (1-p)/p`. Substitution of p0 and the
actual finite-M leakage gives the receipt's rational lower bound. The
fidelity of either retained source marginal is at least the pair
fidelity. Its trace distance from the exact ground density is at most
`sqrt(1-F)`, enclosed by an outward rational square-root bound.

## Readout and retained failures

On preparation success, the collective-spin apparatus acts on one source
copy. The second copy and all spin references are retained. For every
joint detector measurement, its conditional output distribution differs
in total variation by at most the stated trace-distance bound from the
exact-ground distribution. Unitarily evolving the prepared source for
any real time preserves that bound because the target ground is
stationary. Each line/outside/success/failure marginal is consequently
enclosed by the exact-ground interval expanded by this one error. Both
line and complete outside weights have positive lower bounds.

The prepared source is not asserted to be exactly stationary. Its
approximate-ground readout has its own conditional intervals. With
preparation success q and failure 1-q, the full law obeys

`(1-q)+q*(p_line+p_outside+p_pump_fail)=1`.

The receipt gives `q>=p0^2` and retains all failure outputs. The original
maximally mixed input/filter experiment in collective_spin_response
keeps its exact success 1/D; this distinct constructive preparation
does not change that law. Neither experiment proves a relativistic
mass, an empirical particle assignment or a thermodynamic limit.
