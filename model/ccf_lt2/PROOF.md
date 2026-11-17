# Common-source algebra and finite execution

## One source, with its ancestry preserved

The specified 114-factor snapshot is expanded with exact binary-rational
arithmetic. Its complete 22,273-monomial export must equal the source
commitment in [the specification](../SPEC.md). Write its normal-ordered polynomial as
`P=cI+P_0+P_1+P_01`, where each nonconstant local term is retained in its
own cell and P_01 is the remainder. The local block is `cI/2+P_0` and the
other is its translated copy. The audit verifies that equality rather than
assuming it. Replication sums local blocks and actual-edge copies of P_01.
This gives a literal finite-range definition and exact one-edge recovery.
It does not preserve the positive-square decomposition on arbitrary
higher-degree graphs. At each fixed finite graph the matter terms are bounded.

Each component is a full normal-ordered CAR polynomial. The native
Dirac addition retains all terms after its invertible canonical map. The
leptonic and quark current squares retain every term, including constants.
The assembled sum includes these components; subtracting them recovers
the specified matter polynomial. The native empty vacuum need not reduce the quark contact; no such
reduction is used in the full-source checks.

## Gauge, field and current identities

For a CAR word O let s_v be its signed native-third charge change. Its
dressed word is T=O product_e U_e^l_e. The builder checks B l=s, so

    [G_v,T] = ((B l)_v-s_v) T = 0.

Adjoint terms have the reverse shifts and conjugate coefficients. At a
hard wall, a positive/negative shift is the actual corresponding power of
the raising/lowering contraction. `[E,U]=U` remains exact. No cancelled
formal product U U† is used as an identity at the wall. Faces obey B C^T=0,
so every loop and its adjoint preserve Gauss. Electric squares do also.
The sum is Hermitian and Gauss preserving, at finite cutoff and on the
usual infinite-rotor electric domain on a finite graph. Bounded matter and
loop perturbations do not change that self-adjointness domain.

The shared phase derivative defines the matter current. Current terms
come from the Hamiltonian, not from a separately supplied observer
signal. A Hermitian wrong-link mutation tests the local charge identity
independently of Hermiticity. `source.phase_derivative` multiplies each
link exponent by `i l_e`; a second derivative multiplies by `-l_e l_f`.
It includes magnetic loops for the derivative of the complete source,
and exposes the matter current separately with `matter_only=True`.
Electric energy and source scalars have zero derivative. The exact checks
verify Hermiticity, Gauss, local-shard sums, mixed-derivative symmetry and
native charge continuity on each assembled graph.

## Literal complete-source columns

`apply_source_column` executes every actual generated term on a native
occupation bit mask and an electric-flux tuple. The mask represents all 96V
fermion modes. Creation/annihilation signs follow the complete CAR action;
no occupation cut or selected successful branch is applied. The field
Hamiltonian contributes its diagonal electric energy and both orientations
of every magnetic loop, including its scalar term.

The selected native-empty and single native-neutral-occupation preparations
are physical Gauss states with zero flux and reference. Every resulting
output is checked against the full local constraints. These are native
occupation columns, not asserted electron/neutrino mass eigenstates.
Their mean-energy difference includes interaction normal-order terms and
is not identified with a physical neutrino mass.

For a normalized basis state i, the computed coefficients give
`<H>=H_ii` and `variance=sum_f |H_fi|^2-H_ii^2` exactly. Executing columns
is not arbitrary-state time evolution or a bound on many-body simulation
cost. It does, however, evaluate the literal enlarged source including its
spectators and pair channels.

## Positive short-time transition on the complete matter-plus-field source

The numerical bound in this subsection is for `H_S` before geometry is
coupled. It is not a bound at the same time for the nested `H_total`.

On a hard-wall carrier of cap L every CAR word and link-shift product has
norm at most one. Hence an integer upper bound for the complete source is

    M = sum_words ceil(|Re coefficient|+|Im coefficient|)
        + ceil(a |edges| L^2/2 + 2 b |faces|).

The audit constructs it with exact arithmetic. For the Hermitian source,
the integral remainder and unitarity give

    ||exp(-i tau H)-I+i tau H|| <= tau^2 M^2/2.

Take an actual off-diagonal column entry z=H_fi !=0, with f orthogonal to i,
and set tau=|z|/(2M^2). Then

    |<f|exp(-i tau H)|i>| >= 3 tau |z|/4,
    P_f(tau) >= 9 tau^2 |z|^2/16 >0.

The target is a cross-cell native charged-pair/flux output actually present
in the FULL assembled column. Removing the Dirac addition removes this same
matrix element. That control does not assert zero all-time amplitude under
all possible indirect processes. The conservative time can be extremely
small: this proves an allowed unconditioned transition, not a practical
rate, on-shell annihilation, emitted photon or physical particle identity.

## Local energy allocation and geometry

Each Hermitian matter word is shared equally over the vertices touched by
its word and link support. Its adjoint has the same support. Pure scalars
are shared uniformly, a declared energy-density convention. Electric edge
energy is shared by its endpoints, and each magnetic face energy by its
vertices. Therefore every h_v is gauge invariant and Hermitian and
`sum_v h_v=H_S`. The source-column audit checks the latter identity exactly.

This local allocation is part of the construction. Redistributing scalar
energy at fixed local reference energies can change local response, even
if the total H_S is unchanged. It is not a coordinate relabeling.

The geometry modules act on the SAME source-column callback and retain
every output key. Their independent checks compare the complete callback
action with complete finite Kronecker matrices. Applying those modules to
the master source columns is a literal tensor-source execution without
enumerating the astronomical full tensor dimension. Hermitian, Gauss-
preserving source shards and commuting auxiliary factors retain these
identities in the joined law. Complete total energy is conserved by the
time-independent total Hamiltonian; subsystem energies need not be.

The finite lapse response and conditional static Green-function calculation
are not proofs of Einstein dynamics, common electromagnetic/matter light
cones or cosmological selection. In particular, uniform field weights and
edge-length Dirac weights on an irregular complex have not been proved to
approximate the same metric. The RER smooth-observation theorem requires the joint source and
observation premises stated in [the spacetime bridge](EINSTEIN_BRIDGE.md).

## Current of the local-capacity law

The exact matter-phase current is
`sum_v M_v N_v (partial h_v,matter/partial theta_e)`.
The second derivative has the same vertex-specific lapse factors. Fixed
references and auxiliary kinetic/potential terms have zero link-phase
derivative. `apply_total_phase_column` retains every nonzero output, with
`matter_only=False` including magnetic loop derivatives. The full check
compares every first and mixed second phase derivative with the respective
flux-change factors of the complete total column. It checks the native
charge continuity identity using the matter current. The two native columns
include a nonuniform Phi register and all source components.
