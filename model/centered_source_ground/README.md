# The source-selected charge clock and its new ground

Requiring particle-hole invariance of the original bipartite Weyl/contact
source fixes a unique scalar charge-clock centering: `K=H-(U/2)N`.
This package constructs a **new** full ground for that explicitly changed
charged dynamics. It does not assume that the old ground stays a ground.
On the actual cube the old `N=4` ground fails that test, and the new
global ground is proved to lie entirely in `N=8`. An exact integer
certificate proves that it is unique and has hidden spin zero, using
the source spin algebra and a certified trial overlap. The explicitly
imported Lieb theorem supplies a separate ground-spin cross-check.
Its actual one-site state is faithful on the complete four-state CAR
algebra, with a certified least eigenvalue greater than `1/15000`.
The full finite gap is strictly larger than `1/128`; its separation from
every other charge sector exceeds `1/8` model energy unit at `c=1/8,U=1`.
The reference approximation has certified infidelity below `10^-6`.
These finite source energies do not identify relativistic masses.

The original action and earlier scientific results remain unchanged. The
subtraction is an explicitly recorded comparison dynamics, motivated by
the necessary CPT centering in
[continuum_cpt_constraint](../continuum_cpt_constraint/README.md).
The centered full-ground construction has its own state and domain.
It removes that signed first-moment obstruction; it does not prove an
interacting Lorentz theory, physical particle identities, or gravity.

## Deriving the centering from the source symmetry

On a finite open bipartite subset with `V` sites, write

```
H=dGamma(h)+U D,   D=sum_x n_xup n_xdown,   N=sum_x(n_xup+n_xdown).
```

In the physical Weyl frame, let `eta_x=(-1)^(x1+x2+x3)` and let
`J=eta (i sigma_y) conjugation` act on one-particle wave functions.
Spin time reversal commutes with the original hopping: its complex
conjugation changes the sign of `i`, while `i sigma_y` reverses each
conjugated Pauli matrix. Multiplication by `eta` reverses every nearest-
neighbor hopping. Hence `J h J^(-1)=-h`. The canonical particle-hole
automorphism `P(c(f))=c^*(Jf)` consequently preserves `dGamma(h)`;
the possible constant is `Tr(h)=0`. It exchanges onsite occupation with
the opposite-spin vacancy. Therefore

```
P(N)=2V-N,       P(D)=V-N+D,
P(H)=H+U(V-N),
P(H-kappa N)-(H-kappa N)=(U-2kappa)(V-N).             (1)
```

Thus `kappa=U/2` is necessary and sufficient for invariance on the full
Fock algebra. Restricting in advance to `N=V` would erase this test and
is not used. In the real scalar hidden-spin frame, the equivalent
particle-hole transform sends each `c_xsigma` to `eta_x c_xsigma^*`;
a common spin rotation relates the two choices. The code implements the
signed bit-complement unitary and checks (1) against every matrix entry
of the source's full 65,536-state cube Hamiltonian.

For each finite volume, its full centered-ground projector `P_K` is
defined by the lowest spectral value of `K`, without choosing a charge
or spin sector. Functional calculus gives `P(P_K)=P_K`. The normalized
state `P_K/Tr(P_K)` is also gauge and hidden-spin invariant. Equation
(1) gives `omega_K(N)=V`; this expectation alone would not prove sharp
half filling. The cube has the stronger charge-sector proof below.

The gauge-invariant ground has a positive GNS generator for `K`, with
ground energy subtracted only in the GNS representation. This generator
is obtained from the new state. Positivity is not inferred by subtracting
charge from the old state's generator.

## The actual cube changes from four particles to eight

The imported [all-number cube theorem](../interacting_vacuum/README.md)
proves the exact operator inequality `H>=-sqrt(3)/4` on every Fock state
at `c=1/8,U=1`. Its old ground has exactly `N=4` and therefore centered
energy `-2-sqrt(3)/4`, which is greater than `-4`. Any singly occupied
`N=8` product state has centered expectation `-4`: the contact term and
all hopping expectations vanish. Thus the old ground definitely ceases
to minimize the changed clock.

For every `N<=7`, the same all-number inequality implies

```
K|_N >= -sqrt(3)/4-N/2 >= -sqrt(3)/4-7/2 > -4.
```

The exact particle-hole unitary maps `N>=9` to `16-N<=7`, preserving
`K`. Therefore **every global centered ground vector lies in `N=8`**.
This uses the full old theorem as an operator bound, not its old state
as the new reference. No numerical eigenvalue determines the conclusion.

A perfect matching of the cube has four disjoint bonds, each with
`|t|=1/16`. Put the exact centered two-site singlet ground described
below on each bond, with the bond sign absorbed into its local orbital
phase. Each factor has exactly two particles. Every inter-dimer hopping
expectation is therefore zero, and the onsite and matching-bond energies
sum. This gives

```
-4-sqrt(3)/4 <= E_K <= -2-sqrt(17)/2,
inf spectrum(K|_(N !=8))-E_K
 >= sqrt(17)/2-sqrt(3)/4-3/2 > 1/8.                 (2)
```

The last rational bound follows from `sqrt(17)>103/25` and
`sqrt(3)<26/15`. These inequalities alone do not prove uniqueness or
a gap *inside* `N=8`.
An optional balanced-spin `N=8` Ritz calculation reports approximate
energies and residuals as diagnostics; (1)--(2) are its independent
admission bounds and do not depend on eigensolver success.

## An exact full finite gap and a certified reference approximation

The complete centered cube has a unique ground and a gap **strictly greater
than `1/128`** in the declared model energy units. The charge-sector bound
above is `>1/8`; the new calculation covers the remaining neutral sector.
The proof first treats the complete 4900-dimensional `N=8,S_z=0`
block, without assuming a singlet or using Lieb's theorem. It then
establishes spin zero and transfers the bound to every `N=8` spin sector.

`gap.py` obtains a numerical proposal, rounds its trial vector and Gram
factor to dyadic integers, and independently verifies the full matrix
identity over integers. With dyadic trial `v` and factor `Q`, it verifies

```
<v,Kv>/<v,v> < -1055/256,
K + 2105/512 I + vv* = QQ* + R,    ||R|| <= ||R||_infinity < 1/512.
```

Thus `K+1053/256 I+vv*` is strictly positive. A positive rank-one
correction can remove at most one nonpositive eigenvalue, while the trial
Rayleigh quotient lies below that threshold. Thus exactly one balanced
eigenvalue lies below `B=-1053/256`, its ground vector `g` is unique up to
phase, and every other balanced eigenvalue exceeds `B`.
No observed eigensolver residual or reported numerical
eigenvalue is used as an acceptance criterion. Sparse duplicate sums,
every integer product and every row reduction are checked before overflow;
all factor columns enter the proof. The large proposal arrays are
recomputed, not stored as unexplained evidence.

The same unprojected rational trial supplies a certified approximation to
this balanced ground. Let `rho=<v,Kv>/<v,v>` and
`s2=||(K-rho)v||^2/<v,v>`. Exact integer evaluation gives the variance
`s2`. The balanced spectral bound implies

```
1-|<g,v/||v||>|^2 <= s2/(B-rho)^2 < 1/1000000,
rho-s2/(B-rho) <= E_balanced <= rho.
```

The source spin algebra is checked independently on the complete Fock
space: `[S_+,S_-]=2S_z`, `[2S_z,S_+]=2S_+`, and the Hamiltonian commutes
with these generators. Exact frame conjugation matches all three spin
operators in the physical frame and the Hamiltonian to this hidden-frame
CAR construction.
The balanced addresses are independently reconstructed in the source's
combinations order; the rectangular `S_+` matrix from its 4900 states to
the 3136 states with `S_z=1` is checked against the full operator and the
Hamiltonian restriction against the same addresses.

Since `S^2` commutes with the balanced Hamiltonian, its unique ground `g`
has a definite spin. Direct integer evaluation of the same trial gives

```
<v,S^2 v>/<v,v> = ||S_+v||^2/||v||^2 < 1.
```

All spins at `N=8` are integers. Nonnegativity of `S^2` and the certified
overlap would therefore give a trial expectation greater than
`2(1-1/1000000)>1` if `g` had any spin `S>=1`. This contradiction proves
`S=0`. Every integer-spin multiplet contains `S_z=0`; hence the balanced
singlet is the unique ground of the complete `N=8` sector, and every
other energy in that sector exceeds `B`. Subtracting the trial upper
bound proves the neutral gap `>1/128`. The charge bound completes the
proof on all 65,536 states, independently of the imported Lieb theorem.

The replay certifies `-206203/50000 < E_ground < -82481/20000`, trial
infidelity below `10^-6`, half trace distance below `1/1000`, and expectation
error below `||A||/500` for every bounded observable `A`.
The [complete native preparation](../centered_native_preparation/README.md)
executes a signed amplitude tree and the full local inverse frame for this
`N=8` ground on all 65,536 field rows. Its complete preparation error is
less than `10^-7`, with every first-export failure retained. It uses supplied
blank registers and coherent source controls; autonomous ground selection
is not established.

For any explicitly declared Hermitian perturbation with norm below `1/256`,
minmax also preserves uniqueness and gives gap `>1/128-2||V||`. This finite
spectral stability does not establish protection under the repair instrument,
a topological invariant, a thermodynamic gap or a relativistic rest mass.

## An independent imported ground-spin cross-check

We also explicitly import Lieb's repulsive Hubbard ground-spin theorem;
the package checks its finite hypotheses, not a new proof of the theorem.
The primary statement is Theorem 1, printed page 8, of
[Lieb's 1993 review](https://arxiv.org/pdf/cond-mat/9311033), following
the [1989 original result](https://doi.org/10.1103/PhysRevLett.62.1201)
and its [erratum](https://doi.org/10.1103/PhysRevLett.62.1927.5).
For real hopping on a connected bipartite graph, positive onsite
repulsions and even half filling, it gives spin
`S=abs(|A|-|B|)/2` and no degeneracy beyond the spin multiplet.
The review's definition permits arbitrary real hopping signs.

The exact hidden-frame cube hopping is real, symmetric, spin independent,
connected and bipartite with four vertices in each part; every onsite
repulsion is one. `lieb_premises()` checks these source hypotheses.
The preceding full-Fock argument independently establishes half filling,
`N=8`. The review's centered interaction convention adds only the scalar
`UV/4=2` to our `K`, so it has the same eigenvectors. The theorem therefore
independently gives a unique hidden-spin-zero ground. It is not an
admission premise of the direct finite gap proof, but remains an explicit
import elsewhere in the package. Neither argument identifies physical
spin or establishes thermodynamic ground uniqueness.

## Charged response, currents and the complete repair

The unique centered ground is an `S=eta=0` singlet. Spin zero follows from
the direct finite certificate; uniqueness and the commuting eta algebra
force eta zero, since a nonzero eta multiplet would supply additional
degenerate charge sectors. Exact full-Fock CAR commutators then place every
linear creation/removal response, and its entire `K`-cyclic span, in
`S=eta=1/2`. The complete joint singlet isotypic space has dimension 588,
but the ground itself is unique. The charged joint isotypic space contains
1344 copies of the four-dimensional representation; no claim is made that
every copy is reachable by a linear probe. These are internal spin labels,
not an identification with physical Lorentz spin.

The actual particle-hole unitary preserves `K` and the ground. For every
Borel excitation-energy set, it relates the complete addition and removal
matrix measures, including off-diagonal probe coherences. In the physical
onsite frame, `C c C*=F c*` gives `A(B)=F B(B) F*`, where the removal
matrix is defined as `<g|c_i* P_B c_j|g>`; the retarded-response convention
uses its transpose. All local positive addition/removal measures agree.
All twelve oriented source bond currents reverse under the same unitary,
and their eight exact continuity equations use the same centered charge.
External Peierls links transform as `K(A)->K(-A)`; this does not fix their
physical electromagnetic normalization.

The full admitted cube repair has support only in its declared `N=6`
embedding. Every native raw map is zero on the complete `N=7` and `N=9`
sectors, and every outside map is the identity there. The check covers
22,880 input columns and all 130 C/E seam/outcome slots; arbitrary quantum
references and finite histories follow by the operator identities. These
charged responses are therefore repair-transparent, but this repair alone
produces no nontrivial particle record. Transparency is not a localized
topological invariant or a completed detector implementation.

The [complete onsite mass catalogue](../native_mass_channels/README.md)
has four channels. Their actual finite-ground order expectations are all
zero: the scalar staggered excess density is particle-hole odd, and the
three hidden-frame Neel orders form a spin vector in a singlet. This
allows nonzero fluctuations and the rigorously positive many-body gap.
It prevents interpreting that gap as proof of a selected classical mass
texture. The bounded-coupling continuum remains free and massless.

## Actual cube one-site faithfulness and its modular generator

This is a statement about the actual centered eight-site ground, not just
the dimer trial. It also follows for its normalized full-ground projector
without using the imported uniqueness result. Write `T=dGamma(h)` and
`D=sum_x n_xup n_xdown` in the sharp `N=8` sector, so `K=T+D-4`.
The source matrix has `A^2=3I` and zero trace. Its two spin copies give
`||T||=sqrt(3)/2` on this sector. Let `P0` project onto the 256 configurations
with exactly one particle on every site. Every nonzero hop creates a
doublon/hole pair, hence `P0 T P0=0`; the code checks all 3,072 nonzero
directed spin hops. The exact signed vertex-flip symmetries of the hopping
are transitive, so the invariant ground projector has the same doublon
expectation `d` at every site and `omega(D)=8d`.

Set `p=omega(1-P0)<=omega(D)` and
`delta=sqrt(17)/2-2>3/50`. Cauchy--Schwarz on the two off-diagonal blocks
and the norm bound on the remaining block give

```
|omega(T)| <= 2||T||sqrt(p(1-p))+||T||p
           <= 3||T||sqrt(p).
```

The dimer trial implies `omega(T)+omega(D)<=-delta`. Therefore
`delta<=|omega(T)|<=3||T||sqrt(omega(D))`, while also
`omega(D)<=||T||-delta`. These estimates prove

```
d >= delta^2/54 >1/15000,       d <sqrt(3)/16 <1/8.
```

Gauge, hidden-spin and particle-hole invariance now give the complete
one-site density and its finite modular generator:

```
rho_x=diag(d,1/2-d,1/2-d,d),
K_x=-log(1/2-d) I + log((1/2-d)/d) (N_x-1)^2,
lambda_min(rho_x)>1/15000,      ||K_x||<log(15000).
```

Thus the actual centered cube reference removes the old polarized
reference's zero-doublon obstruction at one site. Its entropy satisfies
`log(2)<S(rho_x)<3log(2)-(3/4)log(3)`, because the displayed entropy
increases with `d` on `(0,1/4)`. For every same-site density `sigma`, put
`Delta=sigma-rho_x` and `rho_epsilon=rho_x+epsilon Delta`. The full-CAR
entropy first law and a quantitative remainder are

```
dS(rho_epsilon)/d epsilon|_0 = Tr(Delta K_x),
0 <= D(rho_epsilon||rho_x)
  = epsilon Tr(Delta K_x)-[S(rho_epsilon)-S(rho_x)]
  <= 7500 epsilon^2 ||Delta||_HS^2/(1-epsilon),   0<=epsilon<1.
```

Indeed `rho_u>=(1-u)rho_x` bounds the logarithm Hessian by
`15000||Delta||_HS^2/(1-u)`, and integration of its Taylor remainder
gives the stated factor. This covers noncommuting directions and does
not discard any local state. It neither transfers the numerical bound
to an infinite ground nor identifies `K_x` with a geometric cap charge.

## Complete two-site ground, energies and local faithfulness

For hopping `t>0` and `U>=0`, let `s=sqrt(U^2+16t^2)`. In the singlet
and symmetric-doublon basis, the original two-site matrix is
`[[0,2t],[2t,U]]`. With

```
q=4t/(U+s),
|G>=(|singlet>-q|symmetric doublon>)/sqrt(1+q^2),
E_K=-(U+s)/2,
```

the normalized vector is the unique full-Fock centered ground. The
triplets have energy `-U`, the other two singlets have energies
`(-U+s)/2` and zero, the `N=1,3` energies are `-U/2 +/-t`, and `N=0,4`
have energy zero. Since `s>U` and `s>=4t`, all other levels are strictly
higher. The complete sixteen-state algebra is retained.

The one-site density, in vacuum/up/down/pair order, is

```
rho_site=diag(d,1/2-d,1/2-d,d),
d=(1-U/s)/4 >0.                                    (3)
```

It is faithful on the full local four-state CAR algebra for `t>0`.
Its actual modular generator is the diagonal `-log(rho_site)`. This
faithfulness belongs to the new reference; the old polarized cube's
zero doublon eigenvalue is not discarded or reinterpreted. No geometric
cap or physical temperature is obtained from (3).

For a fixed local spin, both addition and removal have the same two
positive excitation energies and weights:

```
g_- = s/2-t,       w_-=(1+4t/s)/4,
g_+ = s/2+t,       w_+=(1-4t/s)/4.                   (4)
```

Each charge branch has total weight `1/2`; the full signed response
has the four atoms `+/-g_-`, `+/-g_+`, with the corresponding weights.
At `U=0` the upper weight is zero and remains a retained null outcome.
At `U=1,t=1/16`, `g_+/-=(2sqrt(17)+/-1)/16`. These exact source energies
are symmetric under particle-hole exchange. The signed mean is zero
and the second moment is `t^2+U^2/4=65/256`. Their relativistic dispersion
and physical species are not established.

## A centered thermodynamic ground, without reusing the old vacuum

Use normalized full centered-ground projectors in growing open boxes.
Extend each by the tracial product state outside the box; the exterior
generator is zero. This makes the exterior particle-hole invariant
without affecting the interior ground inequality. Average over translates
as in [thermodynamic_ground](../thermodynamic_ground/README.md) and take
one weak-star subsequence. Boundary/volume estimates and the existing
finite-range dynamics convergence pass the local ground inequality to
the limit. The resulting state has a positive centered GNS generator.

All finite states are gauge and hidden-spin invariant. Translating the
particle-hole map changes `eta` by at most a common sign, absorbed into
the gauge symmetry, so particle-hole invariance also survives translation
averaging and the weak-star limit. Thus every site's density is exactly
one and its one-body spin density is `I/2`. This is a new selected
ground subsequence, not a proof that the old unshifted ground works.

The original and centered dynamics differ by a gauge phase. Consequently
the original locality estimates hold with identical norms. More strongly,
with the exact canonical composite `d=(2n_opposite-1)c`,

```
[K,c(f)]=-c(hf)-(U/2)d(f),       ||d(f)||=||f||.       (5)
```

For completeness, the onsite unitary `Z=I-2n_up n_down` obeys
`Z c_sigma Z=-d_sigma`. Products of these commuting even onsite unitaries
give the CAR for the whole `d` family; norm completion gives the stated
identity for every `ell^2` smearing. This norm premise is checked directly
on the complete two-site Fock space.

The full physical-smearing time bound is therefore
`v sum_j||partial_j f||+(Umax/2)||f||`. The same diagonal all-word
construction as [continuum_field_net](../continuum_field_net/README.md)
supplies a centered eight-valley continuum net with positive energy,
continuous translations and the proved propagation cone. Particle-hole
symmetry acts on all valleys (it pairs complementary valley phases),
so no valley is deleted to obtain it. Its averaged addition/removal
energy measures agree. This removes the previous CPT *moment* obstruction;
it does not derive a full spacetime CPT action or Lorentz covariance.

The centered moment identities are exactly

```
M1=h,             M2=h^2+U^2/4.                    (6)
```

They follow either from (5) and local CAR or by translating the original
quadratic-form identities with `r=1/2`. The finite first moment passes to
the centered continuum. The second moment has the upper bound
`||h_W f||^2+(U0^2/4)||f||^2`; equality need not pass through a weak limit.

## A uniform centered occupation estimate

There is a useful estimate stronger than the variational pocket bound.
On the positive kinetic band write `epsilon(k)=c s(k)>0`, where
`s(k)=sqrt(sum_j sin^2 k_j)`. For any measurable band subset with
`epsilon>=delta>0`, define maps to the centered GNS space by
`C f=c(P_+f)Omega` and `D f=d(P_+f)Omega`, using a conjugate input space
to accommodate the antilinear annihilator. The exact ground equation is

```
L C+C epsilon=-(U/2)D.
```

Integrating the derivative of `exp(-tL) C exp(-t epsilon)` gives

```
C=-(U/2) integral_0^infinity exp(-tL) D exp(-t epsilon) dt,
||C|| <= U/(2delta).                                (7)
```

The boundary at infinity vanishes because `L>=0` and `epsilon>=delta`;
`||D||<=1` follows from exact composite CAR. At fixed spacing, all
smeared fields lie in the derivation domain, so the integrated identity
is justified first on vectors and then in operator norm. Translation
invariance makes `C^*C` a Fourier multiplier. Applying (7) on arbitrary
measurable subsets and then on rational positive thresholds proves
the positive-band occupation bound `U^2/(4c^2 s(k)^2)` almost everywhere.
The negative-band creator has the identical bound for its hole density.
No pointwise field operator or pointwise GNS fiber is assumed.

The total kinetic occupation defect therefore obeys

```
w(k)<=U^2/(2c^2 s(k)^2),
D=integral w(k) <= min(2,11(U/c)^2/8).               (8)
```

For the integral bound, fold the Brillouin zone into the eight node
cells `[-pi/2,pi/2]^3`, use `|sin q_j|>=2|q_j|/pi`, and enclose the
cell by the radius `sqrt(3)pi/2` ball. With normalized Brillouin measure,
`integral 1/s(k)^2 <= pi sqrt(3)/2 <11/4`. The isolated nodes have
measure zero. This is specific to the three-dimensional integral.
For `c=v/a`, (8) is `O(a^2)` at bounded coupling. It still does not
uniformly control a normalized packet concentrated in a shrinking
momentum pocket or prove that a nonzero interaction survives.

There is an explicitly separate improvement for the original uncentered
ground. Its equation uses `U B` in place of `(U/2)d`, and the exact
anticommutator has expectation `r||f||^2`. Thus the two GNS maps formed
from `B` and its adjoint each have norm at most `sqrt(r)`. Repeating (7)
gives `w(k)<=2r(U/c)^2/s(k)^2`. The original bound `r<=1/2` yields
`D<=min(2,11(U/c)^2/4)`. This improves the original ground's microscopic
occupation rate without changing its clock, density or nonzero first
moment. It does not identify that original ground with the centered one.

## Replay and outstanding physical scope

```
.venv/bin/python -W error -m pytest -q model/centered_source_ground
.venv/bin/python -W error -m model.centered_source_ground --write
.venv/bin/python -W error -m model.centered_source_ground --check
```

The receipt binds the existing full-source ground bound, the new exact
particle-hole matrix check, imported Lieb theorem hypotheses, actual
cube one-site faithfulness and entropy estimates, the complete dimer
calculation, complete finite gap and trial-state error, full charged selection
and repair identities, and the centered-state construction. Floating proposals
are independently checked by exact integer arithmetic. The
[native preparation and charged readout](../centered_native_preparation/README.md)
retain the full centered carrier and every failure record. An effective
thermodynamic marginal convergence bound, surviving continuum interactions,
physical particle identities and masses, and the same-state geometric,
modular and Einstein requirements are not supplied by this finite result.
