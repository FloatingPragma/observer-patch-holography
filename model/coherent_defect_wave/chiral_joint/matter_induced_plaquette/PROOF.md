# A contractible loop from the complete matter Hamiltonian

## One fixed source and clock

The source is the [literal Peter-Weyl construction](../peter_weyl_links/PROOF.md)
on an antiperiodic cube of side two and spacing pi. It retains all 240 source
CAR modes, eight five-state quantum scalars and 24 complete 119-state quantum
links. All three Yukawa coefficients and both excitation penalties equal
one. The same source Hamiltonian H and clock `exp(-itH)` serve both settings.
No classical scalar or link value is substituted. The matrix elements of the
compressed group coordinates include all nonvacuum Clebsch transitions.

A static reference has thirty CAR modes, transforming in the dual native
R15 gauge representation with an ordinary two-component spin index. Its
gauge currents are the second quantization of `-I_spin tensor T^T`; its
Hamiltonian is zero. Every reference occupation is retained. The complete
source-plus-reference carrier therefore has dimension

`2^270 * 119^24 * 5^8`.

The reference carries no fitted energy. Couplings, spacing, preparation,
observation time and detector settings are architecture inputs independent
of empirical masses or gauge couplings.

## Physical preparation and the actual square

The two recorded settings choose the native e^c block `{9}` of dimension
one or the colored block `{0,1,4}` of dimension three, each with probability
one-half. These blocks are specified by the actual twelve current matrices.
For either block B of dimension d, prepare

`psi_B = sum_(a in B) (|up,a>_source |down,a>_ref`
`                         - |down,a>_source |up,a>_ref)/sqrt(2d)`.

The source fermion is at cell `(0,0,0)`; all other source matter, scalars
and links are blank. The reference is retained. This is an exact gauge and
internal-spin singlet with even total fermion parity. Both source and
reference separately have odd parity. The selected spatial square is an
apparatus orientation; internal spin cancellation is not an assertion of
rotational invariance of that oriented apparatus.

The declared preparation controller swaps the orthonormal total vacuum and
psi_B and acts as the identity on their full complement. Explicitly,

`U_prep=I-|vac><vac|-|psi><psi|+|vac><psi|+|psi><vac|`.

Normalization and orthogonality give `U_prep^2=I`. Since both states are
gauge/internal-spin singlets of even total parity, the controller commutes
with those generators on the full carrier. This is a specified finite
controller, separate from the idle Hamiltonian H.

The square follows

`(0,0,0) -> (1,0,0) -> (1,1,0) -> (0,1,0) -> (0,0,0)`.

Its four stored edge indices are `(0,13,6,1)`. The first two matter hops
use adjoint coordinates; the final two use forward coordinates. All four
geometric edges have their anchor coordinate zero in the relevant axis,
so none crosses an antiperiodic seam. The loop has zero winding in every
periodic direction and is contractible.

Each actual kinetic coefficient attached to these stored edges is `i/4`.
Taking the two required adjoints gives coefficient product `1/256`.
The chronological spin action gives

`sigma_y sigma_x sigma_y sigma_x = -I2`.

These signs and orientations are checked from the actual path compiler and
native Pauli matrices.

## Full quantum loop operator

Let `V_e=sum_ab |a><b| tensor U_e,ab` for the selected native block, with
the literal compressed multiplication coordinates U of the source. Define

`W_B=Tr_internal(V_1 V_6 V_13^dagger V_0^dagger)`.

The trace is unnormalized. In coordinates its four factors are

`U_0,ab^dagger U_13,bc^dagger U_6,dc U_1,ad`,

summed over a,b,c,d in B. Their order across distinct link registers is
irrelevant. Every vacuum coordinate amplitude is the actual `1/sqrt(d)`;
the forward and adjoint branches create different flux orientations. Thus
`W_B |blank links>` has d^4 orthogonal terms of amplitude `1/d^2` and has
norm one. The code obtains these columns from the complete 119-state
matrices, rather than guessing labels or replacing the link by a number.

The full-domain endpoint Ward identities of the compressed coordinates
give `[G_x^A,W_B]=0` for every cell and each of the twelve actual currents.
`loop_ward_checks()` expands their commutators into four-coordinate words
and checks exact cancellation of every coefficient. These are 96 identities
on the entire four-link carrier, not only on its blank state. Independent
executions check all 96 Gauss equations on the resulting loop state.
Each V is a contraction on its complete internal/link carrier, so
`||W_B||<=d`; its specified blank-state norm is exactly one.

Define the normalized final state

`chi_B = psi_B tensor (W_B |blank links>)`,

with scalar memories blank. The matter/reference factor in this expression
does not include an additional copy of the links. Initial and final states
are orthogonal. Their source and reference gauge currents vanish exactly.
The colored preparation has six coordinates, its loop has 81 coordinates,
and its complete final state has 486 coordinates. All are retained.

## Exact full-H moments through degree five

For both settings, the exact matrix elements are

`<chi_B,H^n psi_B>=0`, for `n=0,1,2,3`,

`<chi_B,H^4 psi_B>=-1/(256d)`,

`<chi_B,H^5 psi_B>=-5/(128d)`.

Four distinct link registers differ between initial and final states.
Every actual kinetic path at side two contains exactly one edge; scalar
energies, link energies and Yukawa operators cannot create another changed
link. Fewer than four H factors therefore have zero matrix element. At
fourth order every relevant factor must be a hop on one of the four selected
edges, each used once. The implementation executes all 24 chronological
edge orderings using `peter_weyl_links.model.hopping_column`, including both
forward and adjoint branches and every internal index. Only the stated
cyclic order contributes to the selected oriented trace state.

The spin product supplies the minus sign. The reference contraction gives
the normalized internal trace `Tr/d`; the target loop trace state has norm
one. The resulting factors are `-1/256` for the singlet and `-1/768` for the
colored block. No color direction or reference outcome is selected away.

At fifth order, a relevant history requires at least four kinetic
factors. One Yukawa factor changes the matter number and scalar memory and
cannot return to the final state; two Yukawa factors leave too few hops.
Five kinetic factors cannot return the single fermion to its starting cell:
every actual side-two path flips the bipartite cell parity. The only possible
remaining factor is a diagonal source energy. Its five positions along the
four-hop sequence have values `0,1,2,3,4`, because that many distinct links
are excited. Their sum is ten. All five diagonal insertions are executed for
every edge ordering, giving the stated fifth moment. Actual full-H columns
also check the extracted diagonal action.

These exclusions are exact support statements for these matrix elements and
orders. No Yukawa term, occupied-link transition or scalar state is removed
from H. All of them contribute to the full remainder bound below.

## Complete matter elimination and a nonzero loop kernel

Let P project the entire source/reference matter factor onto the fixed
normalized psi_B while retaining every scalar and link state. Q=I-P on the
full carrier. Their unrestricted ranks are

`rank(P)=119^24*5^8`,

`rank(Q)=(2^270-1)*119^24*5^8`.

P commutes with the gauge currents because its matter/reference state is
an exact singlet. Hence the same decomposition restricts to each complete
Gauss sector. No physical Gauss-kernel dimension is inferred from these
unrestricted ranks. Every eliminated source and reference occupation is
retained in Q, including states with zero coupling to the preparation.

`A=PHP` contains only the scalar and link excitation energies. A hop moves
the sole source fermion to another cell; a Yukawa factor changes its matter
occupation by two. Both leave P. In particular, A has zero matrix element
between the blank and square-loop retained states. The complete electric
and scalar diagonal Hamiltonian leaves the prepared initial state fixed.

Write `D=QHQ`, `B_0=QHP`, `C_0=PHQ`. The inherited complete source word
bound is `||H||<=b=1896`. It uses only unit norms of the scalar/CAR factors
and contraction of each full compressed coordinate. The same bound holds
with the static reference and for all these compressions.

For real `z<-b`, exact elimination is

`H_eff(z)=A-C_0(D-zI)^(-1)B_0`,

with full reconstruction `(H-zI)(u,-(D-zI)^(-1)B_0 u)=(F(z)u,0)` and
`F(z)=H_eff(z)-zI`. No component of Q is projected out of the inverse.

In its norm-convergent large-|z| expansion, the two- and three-factor loop
coefficients vanish. The four-factor coefficient equals the full fourth
moment. Each intermediate state of the only contributing cycle has its
matter fermion away from the initial cell. It is therefore orthogonal to
the full fixed matter/reference-spin singlet P, not just to a selected color
coordinate. This proves Q-irreducibility of the leading history.

Put `R=1024*d*b^5` and `z=-R`. The leading loop kernel entry is the positive
real number

`ell=1/(256*d*R^3)`.

The complete remainder after that term has norm at most

`epsilon=b^5/(R^4*(1-b/R))`.

This follows by summing every higher power of D and using
`||B_0||,||C_0||,||D||<=b`. Since

`epsilon/ell = 1/(4*(1-b/R)) < 1`,

the exact loop matrix element has real part in
`[ell-epsilon,ell+epsilon]` and imaginary part of absolute value at most
epsilon. Its real part is strictly positive, so the entry is nonzero.
Reality of the full entry is not assumed. The receipt contains the exact
rational intervals and the resolvent-distance lower bound `R-b`.

This is an energy-dependent off-shell kernel on a space that retains a
matter/reference preparation. The leading term varies as `1/z^3` and tends
to zero at large negative z. A primitive magnetic interaction is not added.
An autonomous energy-independent pure-gauge generator is not derived.

## Complete records under the original source clock

Use the common time

`t=1/(64*3*b^3)`.

The selected detector has success projector `E_B=|chi_B><chi_B|` on the
complete source/reference carrier and failure `I-E_B`. Both are gauge
invariant and even in total fermion parity. A binary pointer records the
outcome; the controlled pointer flip `E_B tensor X+(I-E_B) tensor I` is
unitary. Every unnormalized success and failure state is retained. The
setting label also remains a recorded outcome.

The degree-five unitary Taylor remainder obeys

`||exp(-itH)-sum_(n=0)^5 (-itH)^n/n!|| <= delta=b^6*t^6/720`.

This is the integral remainder with a unitary factor, valid for the complete
bounded self-adjoint H. It includes all source and reference states. The
fourth-order amplitude is negative real, and the exact fifth-order amplitude
is purely imaginary. Write

`L=|M4|t^4/24`, `J=|M5|t^5/120`.

The exact amplitude modulus lies between `L-delta` and `L+J+delta`.
For this common time, `delta/L=d/4320`, which is strictly below one for both
settings. The original setting-and-success probability therefore lies in

`[(L-delta)^2/2, (L+J+delta)^2/2]`.

The original setting-and-failure interval is one-half minus that interval
with reversed endpoints. All four outcomes are retained and their complete
probability enclosure contains one. No success normalization or discarded
reference record enters this construction. The small time gives a certified
nonzero finite response; it is not a fitted propagation rate. Evolution uses
the original H, not the energy-dependent Feshbach operator.

## Scientific scope and replay

The construction supplies exact gauge-invariant contractible-loop responses
on the complete interacting Peter-Weyl source, including a colored native
channel and explicit references. The first closed-loop term is generated by
matter hops in H. The gauge-field excitation is read through
its actual native currents and quantum link states. It is not identified
with a physical photon, gluon, W or Z, and no Coulomb, confining or
relativistic phase is proved. No mass, ground state or continuum theorem
transfers from another link architecture.

The full-current loop identities, normalized physical states, all minimal
chronological histories, exact higher-order enclosures and complete records
are executable. Check mode validates every consumed Peter-Weyl receipt pin,
binds the transitive source closure and this package's code and proof, and
independently recomputes the exact receipt without writing. All amplitudes
and bounds use exact rational or algebraic values; replay has no numerical
tolerance or eigensolver dependency.
