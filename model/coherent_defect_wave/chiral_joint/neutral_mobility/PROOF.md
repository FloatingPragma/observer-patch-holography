# A physical odd composite with controlled arrival

The construction adds an explicit gauge-invariant composite hopping term
to the complete [119-state quantum-link Hamiltonian](../oriented_links/PROOF.md).
Every native fermion occupation, scalar state, quantum link and Yukawa
channel participates in the dynamics. A prepared physical odd excitation
has a strictly positive probability to reach a neighboring detector. The
calculation retains the no-arrival outcome and bounds the full-source
remainder. No ground state, mass or particle identity is inferred from the
preparation.

## 1. Construction and parameter ancestry

Let the scalar operators phi be the actual five-state quantum memories,
and let native indices 10,11 be the lepton doublet of the twelve actual
native currents. For each Weyl component s define

`J_x,s = phi_x,0 c_x,s,11 - phi_x,1 c_x,s,10`, `D_x,s=J_x,s/sqrt(2)`.

The [full-CAR composite construction](../higgs_composites/PROOF.md)
supplies these exact operators. Their gauge commutators vanish for all
twelve currents, and their spin commutators are

`[S_i,J_s] = -sum_t (sigma_i/2)_st J_t`.

The code independently recomputes the 24 gauge identities and six spin
identities in the complete normal-ordered canonical anticommutation
relation algebra, retaining every scalar matrix unit. These are operator
identities on all thirty native fermion modes, not a one-particle
approximation. D is odd. Its product across distinct sites is even and
commutes with every local Gauss generator.

On the two relevant fermion modes of one spin component and all five
scalar states, `J^dagger J` has eigenvalues 0,1,2 with multiplicities
14,4,2. The remaining native fermion modes are graded spectators for
this norm calculation. Thus `||D||=1` on the complete local carrier.
Its anticommutator with its adjoint is not the identity; no canonical
fermion algebra for the composite is assumed.

For an even periodic side N, let `q(x,i)=-1` on the antiperiodic fermion
seam and `+1` elsewhere. The additional Hamiltonian is

`H_N = (i lambda/2) sum_(x,i) q(x,i) D_x^dagger sigma_i D_(x+e_i)
        + adjoint`.

The complete idle generator is `H=H_oriented+H_N`. Lambda is a declared
finite exact real source coefficient. The default value is one, independent
of any mass or comparison data. The base source has spacing pi, all three
Yukawa coefficients one, scalar penalty one and link penalty one. Those
interactions remain present during the complete arrival experiment. The
clock is `exp(-it H)` in source units with hbar one.

The additional term is a permitted architecture choice. It is not derived
from the star-preserving idle source, a classical Higgs expectation, or a
bare projected Hamiltonian. Its local range is one cell. It acts as identity
on each full 119-state link because D is itself a local gauge singlet;
the separate base gauge transport operators remain fully quantum. The
standard operator estimate gives `||H_N|| <= 6 N^3 |lambda|`: each directed
edge has two unit-norm composite products before its adjoint. This bound
retains arbitrary occupations and scalar states.

## 2. Gauss covariance, parity and spatial symmetries

Every term in H_N has two odd fermion factors and is gauge invariant at
each endpoint separately. Its scalar factors are exact matrices. The
normal-ordered compiler includes all 768 nonzero elementary words at N=2
and verifies equality with its complete adjoint. Hence H is self-adjoint,
preserves total fermion parity and commutes with every local Gauss
generator on the full carrier. Its global lepton number is conserved.

For a proper cubic rotation R with spin lift S, the one-particle spin index
uses S, with creation and annihilation operators transforming as adjoints.
The antiperiodic seam sign comes from the fermion factor; the scalar is
spatially permuted and has spin zero. The same calculation
as for the base Hamiltonian gives `S sigma_i S^dagger=R_ji sigma_j`.
When an axis reverses, the hopping term becomes its adjoint, and the
purely imaginary coefficient changes sign. The antiperiodic seam phases
make the boundary terms agree. The code checks every complete geometric
word under all 24 proper rotations for N=2 and N=4. Three literal
nonblank, occupied full-Hamiltonian controls check pi rotations, and
three check translations. Uniform coefficients preserve these symmetries
without conjugating a matter charge. A two-pi rotation is fermion parity.
Continuous rotations and Lorentz covariance are outside this finite
symmetry theorem.

The base source conserves

`C_x=(-1)^N_matter,x product_(e incident to x) (-1)^N_excited,e`.

Since D_x is odd and does not act on links, `C_x D_x C_x=-D_x`.
An added hopping word flips C at its two endpoints. Consequently H_N
breaks this additional local conservation law while preserving its global
product, total fermion parity. The exact prepared-state control below gives
`||[H_N,C_0] psi||^2=12` at lambda one. This establishes the mechanism
that permits motion in the specified channel.

Literal Peter-Weyl multiplication on the 119-state sector selection would
also contain excitation-preserving nonvacuum Clebsch transitions, as the
base-link proof demonstrates. Such transitions can break this star parity.
They require their own full operator construction and norm bound. The
present result uses the specified composite hopping term alone.

## 3. A complete physical preparation

The all-blank source state Omega has no fermions or scalar excitations and
satisfies every Gauss constraint. It is a declared preparation input, not
a ground of the interacting source. At the origin,

`psi_s = D_0,s^dagger Omega`

`      = (c_0,s,11^dagger |scalar 1> - c_0,s,10^dagger |scalar 2>)/sqrt(2)`.

Both states are normalized gauge singlets, have odd source fermion parity,
and have lepton number one. To prepare them without violating total parity
or lepton number, retain a two-mode fermionic apparatus reference r with
gauge charge zero, lepton number minus one and spin one-half. This is an
apparatus memory, not an additional propagating matter family.

With the reference initially blank, define

`B^dagger=(D_up^dagger r_down^dagger - D_down^dagger r_up^dagger)/sqrt(2)`,

`K=B^dagger+B`.

The source-plus-reference state `Phi=B^dagger Omega_total` is normalized,
gauge invariant, total even, lepton neutral and a spin singlet. On the
span of Omega_total and Phi, K exchanges the two states. The second
creation is zero because the scalar memory excludes a second scalar
particle. Thus the exact gate `exp(-i pi K/2)` prepares `-i Phi`.
The code retains all four active local fermion modes, all five scalar
states and all four reference states in a 320-dimensional matrix check.
It verifies all twelve gauge commutators, all three spin commutators,
total parity and total lepton number. Graded spectator extension gives
the same full-native-domain result.

This gate is a declared preparation intervention before the idle clock
starts. No claim is made that uncontrolled H generates it. During idle
evolution the reference is retained and H acts on the source. The reduced
source input is the equal mixture of psi_up and psi_down; neither branch
is selected or discarded. No stationary-vacuum assumption enters the
calculation.

## 4. A complete local arrival detector

At the target cell y, let Omega_y denote the empty local matter/scalar
state, and define

`P_packet,y = sum_s D_y,s^dagger |Omega_y><Omega_y| D_y,s`.

The two vectors are orthonormal, so this is a rank-two projector on the
complete local matter/scalar carrier. It is gauge invariant, rotation
invariant under spin rotations, and supported in odd local matter parity.
It annihilates all other local occupations, including colored spectators;
those states participate in the complementary outcome.

The arrival effect is

`P_arrival=P_packet,y (I-C_y)/2`, `P_no_arrival=I-P_arrival`.

The factors commute. The second factor checks the retained local star
tag, which, on the packet subspace, requires even incident-link excitation
parity. This specifies a physical neutral packet channel and rejects a
change hidden in an unmatched incident flux record. Its spatial support
is one site and the six incident links. All fields outside this support
are retained without restriction. Both effects commute with all Gauss
generators and total parity, and sum exactly to identity on every full
source state.

An initialized binary flag records the outcome through the complete
unitary `P_arrival tensor X + P_no_arrival tensor I`. The reference and
source are retained. The two flag records are orthogonal and may be
copied to another observer's classical record. No arbitrary source
quantum state is copied, and no dynamical failure sector is removed.

At lambda zero, the prepared state has `C_y=+1` and the full base H
preserves it. Therefore the arrival probability of this specified
detector is exactly zero at every time. This statement concerns the
stated packet/flux channel, not all possible charged or string-dressed
motion of the base quantum source.

## 5. Exact full-source execution and remainder bound

Use N=2, origin `(0,0,0)` and target `(1,0,0)`, with all source coefficients
as specified in section 1. The antiperiodic nearest-neighbor terms on
the two edges between each pair add coherently. For psi_up the complete
source action gives

`P_arrival H psi_up = -i D_target,down^dagger Omega`,

so its squared norm is exactly one. The full H is used, including the
nonstationary scalar and Yukawa response on every other cell. Its sparse
executions have:

| Quantity | Exact value |
| --- | --- |
| Terms of H psi_up | 439 |
| Squared norm of H psi_up | 1675/8 |
| First-action terms with quark occupations | 348 |
| Terms of H squared psi_up | 64206 |
| Squared norm of H squared psi_up | 15948017/128 |

Quantum link excitations occur in the first action and are retained.
The executable checks all 96 cell/current constraints on the preparation
and on its complete first Hamiltonian image. The full normal-ordered
Gauss proof covers the entire evolution, including the second image.
The sparse support is an execution method, not a projected carrier.

For any finite self-adjoint H, the exact integral remainder gives

`||exp(-itH) psi - psi + it H psi|| <= t^2 ||H^2 psi||/2`.

Here `||H^2 psi||^2=15948017/128 <353^2`. Since
`P_arrival psi=0` and `||P_arrival H psi||=1`,

`(t-353 t^2/2)^2 <= p_arrival(t) <= (t+353 t^2/2)^2`

whenever the lower amplitude is nonnegative. At the exact source time
`t=1/353`, the complete probability intervals are

`1/498436 <= p_arrival <= 9/498436`,

`498427/498436 <= p_no_arrival <= 498435/498436`.

The norm error of the first-order state comparison is at most `1/706`.
Every intermediate full-source state in H squared psi contributes to this
bound. It does not substitute a compressed zero for an omitted residual.
A proper pi rotation about the arrival axis maps psi_up to psi_down,
preserves H and fixes the detector, so the same bounds hold for the
second reference branch and the unselected physical preparation Phi.

The positive lower bound proves finite physical motion on this declared
interacting source. It neither identifies a propagating asymptotic particle
nor proves a continuum field. Its coefficient lambda is a declared input,
and the finite source time has no conversion to measured seconds or GeV.
The base axial kinetic term retains its long path range. A stationary
reference, thermodynamic control, a relativistic dispersion, generation
multiplicity and an observed-particle dictionary require results on this
same Hamiltonian and apparatus domain.

The [receipt](results.json) stores exact rational/algebraic outputs and
the executable import ancestry. Its check command reconstructs the full
experiment and both retained outcome bounds. Frozen link or composite
spectra and grounds are not used.
