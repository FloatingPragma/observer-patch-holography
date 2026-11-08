# Original charged response at the default coupling

The centered full CAR source at `U/c=8` has a quantitative obstruction to
retaining its entire original fermion response at finite ballistic energy.
The result uses the actual original and composite CAR operators on the
same ground state and clock. It supplies no many-body charged gap and
does not discard a possible finite-energy component.

## Exact local force on the full occupation space

Write `K=T+U sum_x n_up n_down-(U/2)N`, with each Weyl hopping link
`t_xy=(c/2)U_xy`, where `U_xy` is a unitary two-by-two matrix. The signed
Hubbard chart is included. For a site x of degree d, let

    r^2=sum_y |t_xy|^2=d c^2/4,
    d_s=r^-1 sum_(y,t) (t_xy)_(s,t) a_(y,t).

For r>0 the two d_s satisfy the CAR and anticommute with the two central
site modes. This follows from `sum_y t_xy t_xy*=r^2 I_2`. The full star
coupling is precisely

    T_x=r sum_s(a_(x,s)* d_s+d_s* a_(x,s)).            (1)

Terms on edges not incident to x commute with the odd local field
`b_(x,s)=(2n_(x,opposite)-1)a_(x,s)`. The contact term obeys
`[b_(x,s),K-T]=(U/2)a_(x,s)`. Thus its hopping force
`y=[b_(x,s),T]` acts on four CAR modes, even when the full source has
arbitrarily many sites. At r=0 this force is zero.

The exact four-mode matrices give

    spectrum(y* y)/r^2 = {0 [8], 1 [4], 5 [4]},
    spectrum({y*,y})/r^2 = {1 [10], 5 [4], 9 [2]}.    (2)

Multiplicities refer to the sixteen-dimensional Fock space. Consequently
`||y||^2=5r^2` and `0<={y*,y}<=9r^2 I`. The code also reconstructs the
literal three-site Weyl action for each Pauli-link axis, at `c=1,U=8`,
and verifies (2) for both original flavours. No singly occupied
projection or hidden-spin replacement of the repair map is used.

Use the invariant centered full ground and its retarded response
construction in `COMPOSITE_CAR.md`. Let `v=V0 e_(x,s)`,
`u=V1 e_(x,s)`, and let Y denote the retarded force vector, with the
addition and removal components retained. The same signed retarded
generator L satisfies

    ||v||=||u||=1,  <v,Y>=0,
    L u=alpha v+Y,  alpha=U/2,
    ||Y||^2=< {y*,y} > <=9r^2.                       (3)

The orthogonality in (3) is the invariant-ground identity `V0*Y=0`
from the composite construction; it is not an arbitrary-input operator
identity. The local force inequality itself holds on every full Fock
input. Both statements pass to the stated invariant thermodynamic ground
subsequence because their operators have finite support.

## A projection bound using the original probability

More generally, suppose normalized v,u, a vector Y orthogonal to v,
and a selfadjoint L satisfy `Lu=alpha v+Y`, with alpha>0 and
`||Y||<=beta`. Set

    P=1_[-E,E](L),  w=||Pv||^2,  E>=0.

Since `P v-w v` has squared norm `w(1-w)`,

    alpha w
       <= |<u,L P v>|+|<Y,Pv-wv>|
       <= E sqrt(w)+beta sqrt(w(1-w)).              (4)

The case w=0 is immediate. Otherwise
`alpha sqrt(w)<=E+beta sqrt(1-w)`. For E<alpha this yields

    w <= [(alpha E+beta sqrt(alpha^2+beta^2-E^2))
                   /(alpha^2+beta^2)]^2.             (5)

For E>=alpha use w<=1. The value at E=0 is
`beta^2/(alpha^2+beta^2)`. The executable evaluates an outward rational
enclosure through the equivalent monotone test

    D=(alpha^2+beta^2)w-beta^2-E^2 >=0,
    D^2>=4 E^2 beta^2(1-w),                          (6)

with exact arithmetic. It performs no division by preparation success.
The original signed measure on [-E,E] equals the sum of the original
positive addition and removal measures on [0,E], including both zero
energy contributions.

For the local source (3), `beta^2=9dc^2/4`. On a ballistic family
`c=v_physical/a`, `U=8c`, with fixed physical E and any degree at most six,

    limsup w <= 27/59,
    liminf original weight above E >=32/59.         (7)

These bounds are uniform in the finite volume. Particle-hole and internal
spin symmetry equate the two local charge measures, each of total mass
one half. Their separate low-window upper bound is 27/118 and their
separate escaped-weight lower bound is 16/59. Equation (7) alone does not
bound a spatially smeared packet: force cross correlations matter there.

## Uniform packets on the specified cubic ground average

`MOMENTUM_RESPONSE.md` constructs the arithmetic proper-cubic average of
the translation-invariant centered ground, keeping its same source action
and clock. Magnetic translations, hidden spin, particle-hole symmetry and
reality are preserved. This one averaged state is the reference for the
following packet statement. A symmetry-broken pure reference is not
identified with it.

The actual force Gram on this state has the exact range-two symbol

    G(k)=[g0+2g2 sum_j cos(2k_j)] I,
    g0=c^2(9/2-12qstar),  g2=c^2 q2.                 (8)

Here qstar and q2 are the stated central-neighbour-orbital and second
neighbour spin-plus-eta correlations. Their exact local spectra give
`g0<=27c^2/2`, `|g2|<=3c^2/4`, and therefore `0<=G(k)<=18c^2 I`.
The identity (8), including the cancellation of mixed-axis coefficients,
is proved and checked from the source in that owning file. A bound on
g0 alone would be insufficient to infer this operator estimate.

For every normalized lattice coefficient packet f, the source embeddings
give `v=V0 f`, `u=V1 f`, `Y=Y f`, with exactly the orthogonality in (3).
Equation (5) therefore holds uniformly with `beta^2=18c^2`. At
`c=v_physical/a`, `U=8c`, it gives, for every fixed E,

    limsup_(a->0) mu_(a,f_a)([-E,E]) <=9/17,
    liminf_(a->0) mu_(a,f_a)(outside [-E,E]) >=8/17.   (9)

The packets f_a may depend on a, including kinetic-valley preparations.
Their norms are one. The source state is the stated invariant ground
family; (9) does not supply a quantitative growing-box approximation to
that state. Arbitrary smeared addition and removal measures have combined
mass one, but they need not separately have mass one half. The separate
local factors in (7) are not transferred to them.

There is a direct clock consequence. Let the complete signed retarded
correlation be `phi_a(t)=integral exp(-itE) dmu_(a,f_a)(E)`, so
`phi_a(0)=1`. For any delta>0, average its real part over [-delta,delta].
The averaged deficit is

    1-(1/(2delta)) integral_(-delta)^delta Re phi_a(t) dt
       =integral [1-sinc(delta E)] dmu_(a,f_a)(E)
       >=(1/2) mu_(a,f_a)({|E|>2/delta}).             (10)

The integrand is nonnegative everywhere, and is at least one half on the
displayed exterior. Equations (9)-(10) imply that the liminf of this
deficit is at least 4/17 for every delta>0. If phi_a converged pointwise
to a function continuous at zero with value one, bounded convergence in
(10), followed by delta->0, would instead give deficit zero. Hence the
entire normalized, unfiltered original response has no such regular
limit on this ballistic default-coupling family.

This is an energy-tightness obstruction with complete probabilities. It
does not exclude an energy-filtered response with nonzero original
success weight, a renormalized asymptotic field, or a source-even
collective field. A filtered field must retain the complementary outcomes
and establish its own spectral and momentum convergence. The proof does
not assign a mass or infer a unique capacity trajectory.

## What the estimate does not imply

The three-dimensional response matrix

    L=[[0,4,0],[4,0,3],[0,3,0]],
    v=(1,0,0), u=(0,1,0), Y=(0,0,3)

satisfies the identities used in (4). Its original spectral weights are
8/25,9/25,8/25 at signed energies -5,0,5. It attains (5) at E=0.
Thus these response identities and force bounds permit a nonzero
zero-energy weight. This matrix is an implication countercontrol, not
a claimed realization of the full Hubbard source.

The exact finite eight-site gap is a result for that finite source.
Neither (7) nor finite-volume diagonalization gives a volume-uniform
many-body charged gap at `U/c=8`. The local theorem instead constrains
the original response weight which any continuum construction must
account for. Filtered responses, their complete complements and original
success probabilities require separate spectral control. Source-even
collective fields and other source families are outside this exclusion.
