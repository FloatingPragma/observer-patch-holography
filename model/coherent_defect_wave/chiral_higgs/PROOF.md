# A quantum Higgs interaction for the native chiral module

The matter space at each spatial site is the full exterior algebra of
`C2 tensor R15`, with spin-major mode order `spin*15+internal`.
The [local native encoding](../chiral_source/PROOF.md) realizes that space
and its twelve currents with actual charged native clouds. Its selected
Weyl kinetic symbol is `h_L(p)=-sigma.p`, as in the
[spectral chiral construction](../chiral_regulator/PROOF.md).
The Yukawa tensors are exactly those of the
[owned current calculation](../physical_attachment/PROOF.md), with the
retained exterior-subset basis. No measured particle masses enter.

A Hamiltonian, scalar memory, coefficients and reference are specified
construction data. The native idle repair fixes the even charged cloud and
does not select these additional operators. The interaction below is a
finite quantum law, with a separate conditional constant-background
calculation. Their references cannot be interchanged by replacing a field
operator with its expectation value.

## 1. Local Weyl pairs and exact gauge invariance

For the antisymmetric spin tensor
`epsilon=[[0,1],[-1,0]]`, define

`A_channel,k = (1/2) c^T [epsilon tensor C_channel,k] c`.

Each retained `C_channel,k` is a symmetric fifteen-by-fifteen matrix;
therefore `Delta=epsilon tensor C` is antisymmetric, as required for CAR
pairing. There are three independent channels: up, down and charged lepton.
The up channel uses H, and the down and charged-lepton channels use its
complex conjugate. All fields are left Weyl fields. A singlet conjugate
label is not a supplied right Weyl field with the same gauge representation.

For every one of the twelve Hermitian matter currents T and Higgs generator t,
the retained tensors obey

`T^T C_k + C_k T + sum_j t_jk C_j=0`,

using `t -> -t^T` in a conjugate-Higgs channel. The code checks all 72
component identities after adjoining the two spin components. Since
`[dGamma(T),c]=-Tc`, these are full-Fock operator identities, rather than
checks limited to a one-particle matrix.

Every complex two-by-two matrix S satisfies
`S^T epsilon S=(det S)epsilon`. Thus the pair is a Lorentz spin singlet
under the left-handed `SL(2,C)` field transformation. Executed controls
include nonunitary boosts as well as a rotation. This proves the algebraic
vertex transformation when H is a scalar. It supplies neither the scalar's
relativistic dynamics nor an interacting continuum theorem.

## 2. A finite quantum doublet memory

Use five scalar states: vacuum, two particle states and two antiparticle
states. This is the total occupation at most one subspace of two particle
and two antiparticle oscillators. The projection onto this total number
cutoff commutes with the internal gauge action. Define

`phi_i=|vac><particle_i|+|antiparticle_i><vac|`,

`Q_H(t)=diag(0,t,-t^T)`, `N_H=diag(0,1,1,1,1)`.

Direct matrix multiplication gives `[Q_H,phi_i]=-sum_j t_ij phi_j` on the
entire five-dimensional memory. The finite ladder does not satisfy exact
canonical bosonic commutation relations at its cutoff. Its covariance is
exact despite that cutoff.

The declared local interaction is

`H_Y=sum_channel,k [y_channel A_channel,k Phi_channel,k + h.c.]`,

where Phi is phi for up and `phi^dagger` for down and charged lepton.
The pair tensor identity and scalar commutator cancel exactly, so
`[Q_matter+Q_H,H_Y]=0` for all twelve currents. The coefficients may be any
complex numbers. The interaction is even in the fermion operators and
lifts through the local source encoding without adding support outside
that site. A spatial scalar/gauge kinetic law is a separate declared part
of a lattice construction. The scalar memory itself does not establish a
Klein-Gordon limit or a physical Higgs phase.

The Higgs hypercharge is one half. If a reference density rho commutes with
the hypercharge generator, cyclicity of trace gives
`Tr(rho phi_i)=0`. A gauge-invariant finite state therefore cannot be
replaced by a nonzero charged expectation without changing its state or
observation convention. A gauge-fixed background is a conditional response
setting; its physical realization must instead be expressed through an
appropriate gauge-invariant state and observables.

There is also an exact variance control. On
`chi=(|vac>+|particle_0>)/sqrt2`, the finite field has
`<phi_0>=1/2` and
`<phi_0^dagger phi_0>-|<phi_0>|^2=3/4`.
A nonzero mean alone does not turn the operator coupling into a classical
mass term. The fluctuating part can entangle matter and scalar records.

## 3. Executed quantum interaction and two different references

Set `y_e=1`, `y_u=y_d=0`, and scalar energy `H_H=N_H` in source units.
The active matter factor consists of the six CAR modes for
`e^c,L_1,L_2`, with both spin components. All 64 occupations are retained.
The scalar gives a full 320-dimensional Hamiltonian. The remaining 24 quark
CAR modes are uncoupled spectators in this onsite calculation; they are
not erased from the full construction. The full local matter/scalar dimension
is `5*2^30`, and every displayed active-block multiplicity is multiplied by
`2^24` when all spectator occupations are counted. The prepared spectator
reference is the empty quark state; a spatial kinetic extension changes the
Hamiltonian and must calculate its own reference. Exact source matrices verify all
twelve total-current commutators, both Higgs components and the entire
particle/antiparticle cutoff.

Ordering states by scalar vacuum versus scalar occupation writes the same H as

`H=[[0,W^dagger],[W,I256]]`.

The exact eigenvalues of the 64-dimensional `W^dagger W` are

| Value r | Multiplicity |
| --- | ---: |
| 0 | 6 |
| 1 | 20 |
| 2 | 22 |
| 3 | 8 |
| 4 | 2 |
| 5 | 4 |
| 6 | 2 |

For r greater than zero, H has energies
`lambda_±(r)=(1±sqrt(1+4r))/2`, each with the displayed multiplicity.
The remaining energies are zero with multiplicity six and one with
multiplicity 198. These count all 320 states. The spectral projectors of
`W^dagger W` are the exact polynomials
`P_r=product_(k!=r)(W^dagger W-kI)/(r-k)`, and the receipt checks their
idempotence, ranks, eigenvalue identities and completeness.

The full ground energy is minus two, with multiplicity two, and the gap
above its entire ground space is `(5-sqrt21)/2`. For a vector v in P6,
`(v,-Wv/3)` is a ground vector; its squared norm is `5/3` times the squared norm of v.
This constructs the exact rank-two ground projector and the normalized
reference rho. All twelve gauge generators annihilate its range, and the
range has odd fermion parity. `H+2I` is positive on the entire active source.
This finite, degenerate onsite ground is not asserted to be the Lorentz
vacuum of an interacting field theory, and its onsite excitation gap is
not a relativistic mass.

Lepton charge, assigning minus one to e^c and plus one to each L component,
is conserved by this interaction. Restricting to zero lepton charge and
zero gauge generators gives a five-dimensional sector. The code solves the
actual constraint kernels, without assuming that the full source ground
belongs to them. Its characteristic polynomial is
`(E-1)^3(E^2-E-4)`, up to the harmless determinant sign convention. Its unique
ground energy is `(1-sqrt17)/2`, strictly above minus two. That reference is
therefore a distinct sector ground. Charged apparatus operations can leave
that sector, so its energy origin cannot be used to assert positivity of
all charged excitations without the apparatus and compensating charges.

The initially empty matter/scalar state is not stationary. Its actual
coordinate support closes on five source basis states. In that basis H is
a star with center diagonal zero, four leaf diagonals one, and four real
couplings minus one. On the two-dimensional cyclic subspace of the empty
state, `(H-I/2)^2=17I/4`. At `t=pi/sqrt17` the complete 320-outcome basis
readout has empty probability `1/17`, four matter-pair-plus-scalar outcomes
of probability `4/17` each, and all other outcomes zero. The no-conversion
alternative is retained. This is an executed quantum conversion, without
calling it asymptotic particle scattering.

## 4. Complete charged energy records on the full quantum reference

A fermionic reservoir of the e^c spin-zero mode has initial empty/occupied
settings, each recorded with probability one half. The full fermionic mode
SWAP gives the four source Kraus operators

`c^dagger/sqrt2`, `-n/sqrt2`, `(I-n)/sqrt2`, `c/sqrt2`.

Their squared operators sum to the identity on all 320 states. The
reservoir carries the same hypercharge as the exchanged mode; it is a
weak and color singlet. Thus its exchange preserves the total charge when
the reservoir is included. The source instrument is driven; no claim of
energy conservation without its external work supply is made.

The exact energy PVM has fourteen distinct outcomes, including all zero
probabilities on a prepared state. For r greater than zero, the isometry
proportional to `(v,Wv/(lambda-1))` maps P_r to the corresponding H energy
space. The zero energy space is P0 in the scalar-vacuum block; the remaining
scalar subspace has energy one. These spaces are orthogonal and exhaust
all 320 dimensions. The code evaluates their weights on each of the four
actual post-preparation densities and retains all 56 joint records.

The reference is the normalized full ground projector. A complete energy
PVM on `I320/320` prepares it with original success probability `2/320=1/160`;
all fourteen energy outcomes have probabilities equal to their full
multiplicities divided by 320. This supplies a complete finite preparation
service, without selecting it by idle dynamics or evolving the empty even
state into an odd reference. Every recorded
excitation energy is `lambda+2`, so it is nonnegative. This is a complete
finite quantum energy response on the same source and reference, distinct
from the constant-Higgs spectrum below. A physical spatial dispersion,
interacting reference limit and operational energy standard are needed to
interpret any line as a measured rest mass.

## 5. Conditional background masses without adding mirror fields

For a declared neutral background `H=(0,v)`, use the actual symmetric matrix

`M=v(y_u C_up,1+y_d C_down,1+y_e C_lepton,1)`.

For real positive coefficients, exact multiplication gives
`M^dagger M=diag(v^2 y_u^2 [six],v^2 y_d^2 [six],v^2 y_e^2 [two],0 [one])`
in the retained basis order. The location of each diagonal entry is retained
in the code; the bracket notation denotes multiplicity. There are three
independent complex Yukawa coefficients for one family. No gauge-current
identity determines their ratios or selects v.

With `Delta=epsilon tensor M` and the canonical Nambu coordinates
`(c(p),c^dagger(-p))`, the exact response matrix is

`B(p)=[[h_L(p) tensor I15,-Delta*],[Delta,-h_L(-p)^T tensor I15]]`.

It is Hermitian, satisfies its particle-hole redundancy, and has
`B(p)^2=|p|^2 I60+I4 tensor M^dagger M` for the real coefficients tested.
For general complex symmetric M, Takagi singular values give the same
mass statement, with the appropriate conjugate mass-squared block. The
Nambu matrix duplicates operators; it does not introduce thirty independent
additional Weyl species. Pairing a charged left doublet component with a
conjugate singlet yields a physical Dirac response, with distinct weak
representations in its two chiral components. The resulting conditional
content is three up-color Dirac components, three down-color components,
one charged-lepton Dirac field and one massless neutral Weyl field.
The minimal fifteen-state module has no right-handed neutrino and supplies
no neutrino mass mechanism.

Analytic zero-momentum controls with declared `(y_u,y_d,y_e,v)=(1,2,3,1)`
verify masses one, two and three. These values test the executable matrix;
they are neither empirical inputs nor predicted Standard Model masses.
The finite antiperiodic spectral regulator excludes exact zero momentum;
rest parameters are recovered from its dispersion or the infinite-volume
limit. Its smallest finite-volume energy must not be substituted for the
rest mass. The quantum scalar ground above does not select this background
or reproduce it by taking an expectation. An interacting chiral continuum,
physical Higgs state, source-selected Yukawa coefficients and GeV scale
are separate physical requirements.
