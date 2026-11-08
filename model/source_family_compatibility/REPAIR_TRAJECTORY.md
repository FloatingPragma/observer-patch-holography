# The repair clock and a mobile odd-occupation sector

## 1. Which program the source theorem consumes

`model/native_program/interpreter.py` accepts a literal carrier, preparation
input, record-only repair or hold tape, rational law weights and finite FIFO
policy. Its admission step precedes the supplied source state. The default
local program has C weight `1/2`, repair probability `1`, exploration
probability `1/2` and FIFO capacity one. These readable constants determine
the corresponding program. They do not specify a spin Hamiltonian, a
positive-field ground preparation or a capacity-dependent change of action.
The interpreter's coherent-word service has a separate declared work
register and proper-code loading interface.

The protected-spin theorem consumes the full two-flavour CAR extension of
`defect_field_bridge/full_domain.py`, signed hopping, interaction and the
specified trajectory `c_E=J E^4`, `U_E=J E^8`. Its repair clocks have a declared
rate Gamma. The positive-field construction additionally supplies
`p=n^-2`, `B=J n^-5`, ground preparation and the clock `K_B/(Jp)`.
The generic native-program grammar does not automatically include this
two-flavour field extension. Its control semantics and that extension's
actual branch maps are kept distinct here.

There are two exact restrictions on a repair-only synthesis. Each committed
native repair branch is a weighted partial permutation. A diagonal source
input remains diagonal under any finite record-adapted composition, holds
and classical FIFO operations. The conclusion follows branch by branch and
does not require an unread channel to be unital. On the two-flavour
single-occupancy space P, the stronger existing identity is

    M_(coin,raw) P = delta_(raw,0) P/sqrt(2).                 (1)

Every such full-domain repair-only program consequently gives
`rho_SR tensor classical_record_law` on inputs supported in P, including
an untouched reference R. The record law is independent of the spin state.
Record-dependent seam selection, holds, finite stopping and lossless FIFO
export preserve this factorization. No such program prepares a particular
spin ground from arbitrary spin inputs or generates its nontrivial exchange
clock. A supplied initialization or a supplied coherent intervention changes
these hypotheses and has to be included in the program specification.

This is a statement about the executable classes just specified. It is not
an exclusion of a different explicit source-selection principle.

## 2. Equal complete capacity does not select the law mixture

`program_selection.py` executes two admitted one-attempt programs on H66.
They have the same carrier, uniform port/seam weights, fixed seam `(0,1)`,
one FIFO slot, state preparation and service counts. The C-law weights are
`c=1/4` and `c=3/4`. Both laws keep positive weight; every padded output
block is retained. For the uniform full H66 input, each program reaches
exactly ten distinct `(law,coin,raw)` records among its 25 allocated output
blocks. The stipulated terminal continuation is halt and exact full-record
reread. The existing compound-confusability evaluator gives ten isolated
record atoms and the same correctable capacity `ln(10)`.

Let a be occupation at ports `{0,4}`, b occupation at `{2,4}`, and prepare
`(|a,0>+|b,1>)/sqrt(2)` with a reference qubit. Direct execution gives

    <a,0|rho_out|b,1> = c/4 = 1/16 or 3/16.                (2)

The reference marginal is I/2 and each complete output trace is one.
All integer service counts, including program-description byte count, agree.
The complete capacity and those service counts therefore do not select this
coherent response parameter. Reading either program supplies its exact c;
the response of a fixed program is not ambiguous. This countercontrol does
not assert that either finite terminal program is the full cosmological
capacity fiber or that the bare capacity equation is OPH's only principle.

The stronger finite theorem
`RepairWordSchedule.selected_factorizes` removes schedule parameters when
its complete quotient-visible alphabet, full word-simplex feasibility and
counting-reference minimization premises hold. The two capacity controls
are not countermodels to those extra premises. Section 5 specializes this
theorem to the literal seam/law/coin alphabet; the alphabet's physical
completeness is not inferred from the capacity calculation.
`CapacityClosurePrinciple.closure_value_eq_of_unique` separately consumes
a specified source return map and its fixed-point uniqueness hypothesis.
The equal-capacity control constructs neither a complete cosmological
return map nor two models satisfying that uniqueness hypothesis.

## 3. Exact repair on one missing occupation

Use the same full-domain two-flavour branch function, without changing its
collision rule. On V sites define H_h to consist of exactly one empty site
and exactly one fermion at every other site. Its dimension is
`V 2^(V-1)` and its source occupation is `V-1`. Write P_x for the initial
hole-position projector. Each complete C or E repair preserves H_h: a seam
with two occupied endpoints holds, and a seam touching the hole either holds
or transports the neighbouring fermion with the actual CAR sign. The input
spin and an arbitrary reference are retained.

For a seam e and distinct holes x,y set `a=1_(x in e)`, `b=1_(y in e)`.
The complete unread channel Phi_e satisfies, for every spin matrix unit and
hence every operator A with `A=P_x A P_y`,

    Phi_(e,C)(A) = [1-(a+b)/2] A,
    Phi_(e,E)(A) = (1-a)(1-b) A.                            (3)

If neither hole touches e, both source columns are unchanged. If both
touch e, their raw currents differ for each coin. If exactly one touches e,
the C raw record agrees only for its hold coin, which contributes A/2.
The E record also encodes the odd endpoint total, so the two records never
agree. The retained contribution has no hop and CAR sign +1. This proves
(3) for every spin block, without imposing a spin-polarized input.

The diagonal hole probabilities have a different evolution, identical for
C and E: the hole traverses the seam with probability one half. Its carried
spin is transported with the literal source permutation. This supplies a
mobile source-odd readout relative to the filled reference, rather than an
identification with an even spin excitation. H_h is not asserted to be an
invariant ground band of the interacting Hubbard clock: generic hopping
can create doublons in this sector.

## 4. A complete hole-position detector preserving internal spin

Choose a deterministic sequence of E seams covering every vertex, with m
seams. Retain every coin and all sixteen raw slots at each step. Any two
different initial hole positions have disjoint full-record supports. To
prove this, take the first seam touching either initial hole. Before that
seam neither hole has moved. If it touches only one, the E parity records
differ; if it touches both, their currents differ for either coin.

For each nonzero fine history h there is thus a unique initial hole x(h).
Its Kraus map on H_h is

    K_h = 2^(-m/2) U_h P_(x(h)),                            (4)

where U_h is the signed spin transport along the recorded hole path. The
branch is injective, so its restriction extends to a signed permutation of
H_h. A record-controlled inverse transport, if supplied as a coherent
service, gives `U_h* K_h=2^(-m/2)P_(x(h))`. There are exactly `2^m` positive
histories for each initial hole and `V2^m` positive histories overall.
All `(32)^m` allocated histories remain. Sparse dictionary entries omitted
from storage denote zero maps; `scan_branch` supplies zero-default access
to every allocated slot.
The complete corrected channel is

    rho -> sum_x P_x rho P_x.                              (5)

Equation (5) preserves every internal-spin/reference matrix element within
each position block, with its original probability. It is a source-derived
position instrument on a common odd-occupation and spin register. The code
replays its signed branch maps; no inverse-transport native pulse tape is
claimed executed. Without that supplied correction the source executes the
same position measurement followed by its recorded spin transport.

For a coherent superposition of two hole positions, the source plus
classical record is block diagonal in those alternatives. Erasing or
compressing an exported classical record does not reconstruct the lost
source/reference coherence. Retaining a coherent premeasurement register
and undoing the entire dilation would be a different record interface.
The full isometric dilation retains the input information; the dephasing
statement concerns the committed classical-record algebra and its source
marginal, rather than global information loss in a closed quantum source.

## 5. The actual repair clock gives diffusion and dephasing

There is an exact native-attempt statement requiring no Poisson clock.
For the declared alphabet `(seam,law,coin)` with `4|E|` visible identities,
the counting-selected law of the finite theorem is IID uniform. With
`s_xy=|star(x) union star(y)|`, averaging (3) gives the one-attempt
decrement and k-attempt cross-position multiplier

    delta_xy = [deg(x)+deg(y)]/(4|E|)+s_xy/(2|E|),
    Phi_bar^k(P_x rho P_y) = (1-delta_xy)^k P_x rho P_y.     (6a)

All factors are exact rationals on a finite graph. This directly consumes
the finite theorem's factorization conclusion, while retaining its
full-alphabet, full-simplex and counting-reference hypotheses. Native
attempt count is not identified with a measured physical time.

Specify independent Poisson attempts at rate Gamma on every seam, with E
weight theta and C weight `1-theta`. The source theorem uses such a declared
Poisson clock; neither its physical rate nor theta is selected by the
capacity calculation in section 2. With no intervening coherent control,
the generator on a finite graph is exactly

    L_rep = Gamma sum_e [(1-theta)Phi_(e,C)+theta Phi_(e,E)-I].

Equation (3) diagonalizes every cross-position spin block:

    exp(t L_rep)(P_x rho P_y) = exp(-r_xy t) P_x rho P_y,
    r_xy = Gamma sum_e [(1-theta)(1_(x in e)+1_(y in e))/2
                        +theta 1_(e meets {x,y})].          (6)

No record is postselected. For a regular degree-d graph and x different
from y, the C contribution to r_xy is `Gamma d`. The E contribution is
`Gamma(2d-1)` for neighbours and `2Gamma d` otherwise. These rates hold on
all internal-spin coherences, including arbitrary reference coefficients.

The hole probabilities solve the continuous-time nearest-neighbour random
walk with rate Gamma/2 per edge. On the infinite cubic graph, the diagonal
Fourier multiplier and displacement second moment are

    lambda(k)=Gamma sum_(i=1)^3 (cos k_i-1),
    E |X_t-X_0|^2 = 3 Gamma t.                              (7)

The second formula follows from the generator on the squared displacement;
the first is its action on a plane wave. A scaled position `pX` has a finite
diffusive limit when `q Gamma p^2` stays finite for clock factor q. Its
cross-position coherence on that clock decays with `q r_xy`. This is an
exact open-system repair calculation, not a relativistic dispersion law.

In particular, the bare hole channel cannot approach identity on all such
coherences unless `q Gamma ->0` on a fixed-degree connected lattice. For
`(|a>+|b>)/sqrt(2)`, the norm-one test
`|a><b|+|b><a|` gives trace-norm distance from the unchanged input at least
`1-exp(-q r_xy t)`. This lower bound uses the original state weight. For the
spin ballistic clock `q=1/(Jp)`, fixed positive Gamma violates that necessary
condition. A repair rate sent to zero, a coherently retained record or a
different protected code would be an explicit modification of the selected
source family. The calculation does not exclude a different collective
matter observable or decide a Hubbard-plus-repair limit with nontrivial
coherent propagation between repairs.

## 6. A fair global-attempt schedule suppresses the disturbance

The number of global attempts and a rate per seam have different growth
laws. There is a positive complete-record comparison for the native global
schedule. Let V be the exact one-attempt dilation with uniform seam,
E weight theta, fair coin and all raw slots. Let V0 leave H_h unchanged and
write the same seam/law/coin distribution with raw zero. The latter is the
identity comparison instrument, rather than a postselected branch. Direct
calculation with the actual raw-zero holds gives

    V0*V = I - (1+theta)/(2|E|) sum_x deg(x) P_x,
    (V-V0)*(V-V0) = (1+theta)/|E| sum_x deg(x) P_x.          (8)

For C, one half of the hit coins hold with raw zero; for E, all hit coins
have nonzero raw code. Unhit seams give identity. These are the only terms
in V0*V; they have positive sign. Both dilations are isometries, so their
Gram difference is exactly (8). The identity tensored on any reference
does not change the operator norm.

Interleave m attempts with any common source unitaries preserving H_h,
and, if desired, common record-adapted channels that preserve this sector.
The source operations and their policy must be identical on each classical
record in the two comparisons. Telescoping the complete dilations, or the
corresponding reference-stable channel bounds, gives

    ||rho_native,source+archive - rho_identity,source+archive||_1
       <= min(2, 2m sqrt((1+theta)d_max/|E|)).              (9)

The archive includes every seam/law/coin/raw record. No conditioning on
raw zero occurs in the native calculation. Complete classical commitment,
FIFO export and subsequent common channels cannot increase (9). The same
argument also bounds the coherent dilations before commitment. Errors in
a supplied native pulse realization of either source propagation or repair
must be added separately; (8) is the exact branch-map comparison.

For example, take open cubes of side `n^2+1`, `p=n^-2`, and at most `n^2`
global attempts. They have `|E|=3n^2(n^2+1)^2`, maximum degree six, and
for every `0<=theta<=1`,

    (1+theta)d_max/|E| <= 4/n^6,
    complete source/archive trace error <= min(2,4/n).     (10)

This explicit fair schedule preserves one-hole coherence and spin over
that attempt budget in the growing-source limit. It is a workaround to
the fixed-per-seam-rate disturbance, using the same repair map. The
corresponding no-control hole diffusion also becomes small on this budget;
repair alone has not generated a propagating quantum field. Equation (10)
does not choose a physical duration per attempt or derive this cube/capacity
trajectory. Counting selection supplies the uniform schedule only under
the full-alphabet and word-simplex hypotheses stated in section 2.

The propagation premise in (9) must not be transferred to the full Hubbard
Hamiltonian on arbitrary H_h inputs. That Hamiltonian can produce doublons.
Its fully spin-aligned one-hole sector does stay invariant by number and
flavour conservation and the Pauli exclusion rule, so (9) applies to that
actual source restriction for every c,U. It supplies no interacting ground
sector: the aligned restriction has no double occupations and is a free
hole hopping problem. The general interacting low-energy and physical
particle attachments require separate source estimates.

## 7. Consequence for the common-source construction

The exact source offers both protected spin information and an odd-hole
position instrument. Its repair-only records leave the first invariant
and measure the second. The high-U full-ground theorem additionally sends
the odd-CAR energy gap to infinity on the spin ballistic clock. These are
separate constraints: (6) concerns the repair instrument on H_h, whereas
the energy exclusion concerns ground-subtracted Hubbard dynamics.

A construction retaining original odd CAR and the spin fields must give
one actual state and clock outside that excluded projection regime, and
control both its coherent propagation and the complete repair disturbance.
The native capacity theorem and readable default program do not supply
such a capacity-to-action trajectory. A fixed specified program can be
tested against equations (1)--(7); assigning the missing trajectory by
reading a numerical particle mass would change the source-selection task.
