# The actual memory complement is dissipative

The five-state memory in `native_repair_mass` advances with probability
p=1/12 and otherwise resets. Old memory is a retained output record. It is
a classical-quantum instrument; replacing it by a coherent sterile sector
would change the source. This audit reconstructs its complete ideal
180-dimensional CQ superoperator, with all five memory labels and all36
Nambu/family density coordinates. The existing native implementation has
its own transfer/error receipt; this calculation does not erase that error
or claim a new literal native execution.

The stationary distribution is πj=(1−p)p^j/(1−p^5). In memory coordinates
(π,e1−e0,...,e4−e0), T splits into1 and a4-dimensional block whose exact
characteristic polynomial is

    x^4 + p x^3 + p² x² + p³ x + p^4.

Its four roots are p times the nontrivial fifth roots of unity. The gap
to the unit circle is11/12. This is a relaxation gap. It is not a positive
sterile mass.

Let M(h)=diag(exp(h Lj))(T tensor I36), with
Lj=−i[X tensor Pj,.]. After the same memory-coordinate change, write its
blocks as A(h),B(h),C(h),D(h), where A has dimension36 and D has dimension144.
The complete projected resolvent is exactly

    [zI−A(h)−B(h)(zI−D(h))^−1 C(h)]^−1.

The producer checks this identity against inversion of the full180-state
matrix at three nonreal z values for h=1/8,1/64,1/512. It also retains every
pole and projected spectral residue norm. The maximum direct Schur error
is below5e−15. Spectral reconstruction becomes ill-conditioned as h shrinks
(eigenbasis condition about1.12e6 at h=1/512); these eigenvalue listings are
floating diagnostics, not interval certificates.

At h=0 the marginal pole is36-fold. Expanding M and eliminating its full
complement gives the slow generator

    K(h)=K1+h K2+O(h²),
    K1=M1,PP=−i[X tensor Y,.],
    K2=M2,PP+M1,PQ(I−D0)^−1 M1,QP−K1²/2,

where Mk is the coefficient of h^k in the transformed M. All terms are
computed from the source projectors and exact reset probabilities.

The anti-self-adjoint part of K2 is a Hamiltonian commutator. Its Hermitian
coefficient has the form I2 tensor iA3, where A3 is real antisymmetric.
This follows already from [X tensor Pj,X tensor Pk]=I2 tensor[Pj,Pk].
Its expectation in every real eigenvector of Y is zero. Therefore the
first correction in h to the coherent mass eigenvalues vanishes. The
symmetric part of K2 remains: retaining only its Hamiltonian component
would discard damping from the actual recorded source.

Numerically each tested complete channel has36 slow poles and144 fast
poles. Only two poles have modulus within1e−10 of one. This numerical
observation is not promoted to an exact peripheral-spectrum theorem.
None of the channels is unitary. The leading positive Hamiltonian
coefficients remain approximately(0.00206415,0.0636598,0.934276), not a new
light spectrum. Even at the largest tested step the displayed coherent
approximation changes the smallest by only about1e−7 in source units.

Consequently this complement does not supply a three-state unitary neutral
mass/current frame. No PMNS matrix is inferred by removing damped poles.
It also does not modify the separate432-state charged register: the source
attachment needed for that claim is absent. A new coherent sterile service
must declare its Hamiltonian and currents explicitly. The sibling
`neutral_chain` package does that as a separate construction; it does not
rename the recorded reset chain as a seesaw.
