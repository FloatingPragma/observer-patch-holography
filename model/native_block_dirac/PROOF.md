# A free Dirac readout from actual odd native spectral states

## Fixed native constituent and source-selected energy

Each constituent is the complete eight-cell native Hamiltonian
`h=H0+2^-34 T` and the old/dual/port local reference in the native abelian
join. Its reference has one shared vacuum, the complete original and dual
occupied Fock sectors with energy4n, and the two even, zero-energy gauge
ports. The full physical constituent has unique ground g with reference
vacuum and energy Eg. Its complete old, dual, port and neutral complements
have the proved positive lower bounds. No native scalar, link or matter
occupation is removed.
Explicitly the full physical constituent gap is
`min(delta_neutral,a_port,beta_old_dual)>0`. The first quantity is the
frozen neutral gap, the second bounds both complete charged-port sectors,
and the third bounds every occupied old or dual reference sector.
`source.py` admits all three exact rational bounds. This gap is distinct
from the selected rest frequency E*.

The actual charged source spectral bin is [L,U]. Select the energy E* in
that bin with greatest per-spin native annihilator residue, using the lower
energy to break a tie. The positive bin residue and finite dimension prove
that this selected residue is nonzero. This is a source-only finite
spectral service, not a fitted mass or an assertion that the bin is rank
two. Native charge/reflection, the actual trivial proper-cubic ground
character and the full local reference supply orthonormal states
`e_(q,s)`, q=+1,-1 and s=up,down, at excitation E*. Each is an actual odd
native state with its matching even gauge port. The two charges have
conjugate G1 matrices. A fixed pseudoreal basis change can identify these
matrices, but the conjugate bases are convenient for the field below.
Here q is the actual source Y=+1 for addition and Y=-1 for removal;
the compensating port carries -q and the boundary operator Q assigns q
to that port sector. The native e^c basis has Y=1 and integer
Peter-Weyl charge 6Y=6. No change to this normalization is made.

Set m=E*. This m is the rest-frequency parameter of the specified free
construction. It is an actual selected native excitation, constrained by
the source bin. The length scale and unit propagation speed of the macro
lattice are declared conventions. No GeV conversion, electromagnetic
identification or observed species is assigned.
The e^c-based vectorlike charge pair does not supply the electron's
Standard Model weak chirality. No second weak representation is inferred
from a four-component free Dirac readout.

## The complete doublet Fock code

Fix q and use two independent complete constituents. In the selected
three-state space of each, write 0=g and 1,2=e_(q,up/down). Define

```
E|vac>    = |g,g>,
E|up>     = |e_up,g>,
E|down>   = |e_down,g>,
E|up down>=(|e_up,e_down>-|e_down,e_up>)/sqrt2.
```

All four vectors are orthonormal. The two-occupied vector is a spin
singlet: for any two-component matrix S, its action is det(S). Native
proper rotations have det(S)=1. The one-occupied prefix carries the
actual G1 matrix S, and the vacuum is invariant. Thus the code intertwines
the full exterior action `1 direct_sum S direct_sum det(S)`. The same
identity applies to the conjugate charge representation.

The actual sum of constituent energies is 0,m,m,2m on these vectors,
above2Eg. Actual fermion parity is +,-,-,+. With the two-mode canonical
annihilators c_s, define `a_s=E c_s E†` and set them to zero on the full
orthogonal complement. They are odd, gauge invariant under each
constituent's native Gauss law, and obey

```
{a_s,a_t†}=delta_st Pcode,     {a_s,a_t}=0.
```

These are literal matrix identities, checked on all nine selected product
states. The whole native complement remains reducing. On the original
three-state code the up anticommutator is diag(1,1,0); it cannot be
identified with two-mode CAR. Two occupied constituents are necessary for
the two-particle state under this selected one-excitation-per-constituent
encoding. No additional native two-particle eigenvalue is assumed.

Use an independent pair for q=+1 and another pair for q=-1. The resulting
four-constituent code has all sixteen occupations of four canonical modes,
inside an 81-dimensional selected product space and the full native
carrier. Its number and boundary charge are

```
N=sum_(q,s) a_(q,s)† a_(q,s),
Q=N_plus-N_minus.
```

The graded tensor product of macrocell codes gives canonical anticommutation
between cells. The actual native grading supplies these signs. This is a
full-occupation code; it is not a one-particle truncation. Spatial support
of even encoded operators is the support of their constituent cells in
the declared graded source algebra. A bosonic finite-range compilation of
the grading strings is not asserted.

## One positive Hamiltonian, exact ground and full complement

Take a periodic cubic macro lattice of side N and spacing a. Let L_a be
its nonnegative nearest-neighbor Laplacian. On its momentum menu,

```
ell_a(p)=4/a² sum_i sin²(a p_i/2),
Omega_a(p)=sqrt(m²+ell_a(p)),
K_a=sqrt(m² I+L_a)-m I >=0.
```

Declare the complete source Hamiltonian

```
H_a=sum_(all constituents)(h+H_local_ref-Eg)
       +sum_(q,s,x,y) a_(x,q,s)† (K_a)_(x,y) a_(y,q,s).
```

Each coefficient is a specified finite matrix function of the declared
macro Laplacian. Each source operator is a matrix unit between actual
native spectral states and the above occupation code. All native Gauss
constraints are preserved exactly; boundary Q is conserved. The source
includes the complete native local references. The outer affine cone's
Coulomb and Peierls Hamiltonian is not part of H_a, and no interacting
affine clock or spectral result is transferred to it.

The kinetic term is positive on the entire source: write K_a=B†B and
sum it as `sum_z (sum_y B_zy a_y)†(sum_y B_zy a_y)`, for each q,s.
This argument needs no CAR outside the code. The native term has
the unique product ground G. Every a annihilates G, so G is the actual
unique physical ground of H_a. The full gap is at least the constituent
gap, uniformly in the number of macrocells. Inactive native sectors are
retained and cannot undercut the ground.

Every local code is reducing for the native term, and all a preserve the
global product code. On that code, exactly,

```
H_a=dGamma(Omega_a tensor I4).
```

The comparison is therefore an exact reducing free Fock sector of the
declared full Hamiltonian. Negative physical source energies and a changed
vacuum are unnecessary. The nonlinear square-root couplers are a disclosed
kinetic modification, not a consequence of the old uniform native source.

## Spatial nonlocality of the source law

The square-root coefficient has the integral representation

```
K_a=(1/(2sqrt(pi))) integral_0^infinity
       exp(-m²t) [I-exp(-t L_a)] t^(-3/2) dt.
```

For distinct vertices of a connected finite lattice, the heat-kernel
matrix entry is strictly positive for every t>0. Hence every off-diagonal
entry of K_a is strictly negative. The source kinetic law has full spatial
range. No finite-range microscopic light cone is claimed. The continuum
field's causal anticommutator is established separately below, using both
frequency signs and the complete occupation algebra.

## Same-clock Dirac filters

Use Dirac matrices `beta=diag(I2,-I2)` and
`alpha_i=[[0,sigma_i],[sigma_i,0]]`. Put

```
omega(p)=sqrt(m²+|p|²),
r(p)=sqrt((omega+m)/(2omega)),
U(p)=r [ I2 ; sigma.p/(omega+m) ],
V(p)=r [ -sigma.p/(omega+m) ; I2 ].
```

Direct Pauli multiplication proves

```
U†U=V†V=I2,   U†V=0,   UU†+VV†=I4,
h_D(p)U=omega U,   h_D(p)V=-omega V,
h_D(p)=alpha.p+beta m,   h_D(p)²=(m²+|p|²)I4.
```

Let a_s(p) annihilate a positive-charge code excitation and b_s(p)
annihilate a negative-charge one. Define the momentum field

```
Psi_a(p,t)=U(p) a(p) exp(-i Omega_a(p)t)
              +V(p) b†(-p) exp(+i Omega_a(p)t).
```

Every exponential is the actual H_a clock. The opposite frequency of
b† follows from its creation action; it is not a negative-energy source
state. Both terms carry boundary charge -1 and odd source parity.
Treating opposite-charge annihilators as one same-charge four-component
one-particle vector would fail this charge identity.

The equal-time field has exact canonical anticommutation on the complete
Fock code, including every occupation. Its vacuum two-point kernels are
`UU† exp(-i Omega_a t)` and `VV† exp(+i Omega_a t)`. The actual ground G
is the vacuum for a,b, and hence the Dirac quasifree vacuum for this
readout. At finite a the equation has the explicit defect

```
i partial_t Psi_a = (Omega_a/omega) h_D Psi_a.
```

For |p|<=P and aP within the first Brillouin region,

```
0<=omega-Omega_a<=a² P⁴/(24m)=delta_a.
```

The estimate follows from
`0<=p²-ell_a(p)<=a² sum_i p_i⁴/12` and division by
`omega+Omega_a>=2m`. Thus the finite field converges on a fixed momentum
and time window without a Hamiltonian or state substitution.

## Induced rotations, boosts and field locality

The frame identities define an actual free spin-one-half representation,
in addition to the inherited finite proper-cubic action. For a proper
spatial rotation with spin matrix S,

```
U(Rp)=diag(S,S) U(p) S†,
V(Rp)=diag(S,S) V(p) S†.
```

Particle one-particle states carry S and antiparticle states its conjugate;
this is consistent with their creation term in Psi. On the native code
the finite cubic subgroup agrees with the proved native G1 action.
The full native complement is not assigned an unproved continuous spin
symmetry. Continuous rotations and boosts here act on the limiting
readout sector and its test functions.

For positive on-shell p define `u(p)=sqrt(2omega) U(p)` and
`v(p)=sqrt(2omega) V(-p)`. They obey `ubar u=2m I` and
`vbar v=-2m I`. If S(Lambda) is the Dirac boost matrix, its positive
frequency intertwiner is

```
D_+(Lambda,p)=ubar(Lambda p) S(Lambda) u(p)/(2m).
```

The transformed columns solve the same positive-energy Dirac equation,
so `S(Lambda)u(p)=u(Lambda p)D_+`. Preservation of the Dirac inner
product proves unitarity; connectedness fixes determinant one. Uniqueness
of this factor proves the cocycle law. The negative-frequency formula
uses `-vbar(Lambda p) S(Lambda) v(p)/(2m)` and the conjugate action on
antiparticle states. With invariant measure `d³p/(2omega)`, these are
unitary induced one-particle actions. In the source's plain d³p normalization,
the action includes `sqrt((Lambda^-1 p)^0/p^0)` multiplying the Wigner
matrix and the amplitude at Lambda^-1 p. This Jacobian factor is required
for unitarity. Their exterior lifts preserve the
vacuum and supply the covariant Dirac field transformation. This is not
an exact finite-lattice boost symmetry. A bounded rapidity rho sends a
momentum ball P inside the ball `exp(rho)*sqrt(m²+P²)`, so a finite-window
comparison can use that enlarged menu.

The executed noncollinear control has `p/m=(3/4,0,0)` and a z boost with
cosh(chi)=5/4. It gives the exact Wigner matrix
`[[9,1],[-1,9]]/sqrt(82)`, whose unitarity and intertwiner are checked.

The limiting field anticommutator is independently fixed by both spin
sums. With p0=omega and

```
C_m(x)=integral d³p/[(2pi)³ 2omega] [exp(-ip.x)-exp(+ip.x)],
```

one obtains

```
{Psi(x), Psibar(y)}=(i gamma.partial+m) C_m(x-y).
```

C_m has zero equal-time value and initial time derivative `-i delta`.
The Klein-Gordon finite-propagation energy estimate places its support
inside the light cone: approximate its initial delta velocity by smooth
data supported in a ball of radius epsilon and integrate the positive
energy density `(|dot C|²+|gradient C|²+m²|C|²)/2` outside the expanding
ball of radius epsilon+|t|. The outward-boundary flux is nonpositive
because `2|dot C partial_r C|<=|dot C|²+|partial_r C|²`. The exterior
energy starts at zero and stays zero. Taking the distributional delta
limit gives the stated support; derivatives preserve support. Thus the Dirac
anticommutator vanishes at spacelike separation as a distribution. This
establishes locality of the limiting free CAR field. A finite momentum
projection is not exactly local. This result supplies no Maxwell Gauss-law
continuum or locality of an interacting electromagnetic dressing.

## Complete bounded charged-field instruments

For a finite four-component test f define `d=Psi_a(f,t)/C`, with
`C>=||f||_2`. On the code this is a canonical annihilator scaled by at
most one, since Psi combines mutually anticommuting a and b† modes.
On the full native carrier it remains a contraction: the local
anticommutator equals a positive coefficient times the local code
projector, so `{d,d†}<=I` and `d²=0`.

Retain a complete two-mode fermionic apparatus doublet of boundary charge
+1, raw lepton number -1, with the native spinor rotation matrix and zero
apparatus Hamiltonian. Both modes start empty. Select one mode r for the
stated test and retain its spin partner as well. Use the even invariant gate

```
W=exp[theta(r† d-d† r)].
```

It preserves total boundary charge and parity and every native Gauss law.
Its anti-Hermitian generator is a bounded finite matrix on the entire
source and reference. Read both occupations of the complete doublet,
including its zero-weight partner outcomes. For r tests use r distinct
retained reference doublets and keep all 4^r records. Measurement of an earlier reference can be deferred
because later gates do not act on it. Its complete Kraus maps are the
matrix elements of the full unitary gate sequence between the initial
empty reference and each final occupation word. Their effects sum to the
identity. No postselected source amplitude is substituted for a full law.

The charged reference is a measurement apparatus, not an erased charge
sink. These are native-boundary-charge relative records. Efficient local
implementation of spatially smeared matrix services is not asserted.

The finite source has joint macro translations and proper native/macro
cubic rotations. Retain a regular reference for
`G=(Z_N)^3 semidirect 2O`, of dimension `48N³`. On the system consisting
of source and detector doublets use its exact representation U_g. Define
`J psi=|G|^-1/2 sum_g U_g psi tensor |g>`, and replace a gate W and
each complete effect by the controlled family `U_g W U_g†` and its
effect analogue. The identities `Wtilde J=J W` and
`Htilde J=J H` prove equality of every joint record and norm error.
The regular reference is retained with zero Hamiltonian. This gives a
jointly covariant apparatus for the specified finite spatial group.
The full native trial and its retained spin references are synthesized
from the joint empty invariant state. The preparation unitaries preserve
joint cube rotations, and identical filter trees preserve macro
translations. The ideal main ground and empty detector doublets are
invariant, so their initial J state uses a product uniform regular
reference. More generally the same controlled lift applies to the complete
preparation, including every partner and flag. No asymmetric state is
silently supplied after tracing a reference.
The actual readouts are relative to that reference. Tracing it gives a
twirled marginal, not an additional detector protocol. No native internal
translation symmetry broken by the origin spectral selection is imported.

## Finite momentum, volume and complete-record errors

The frame matrix `[U V]` is unitary. Hence on a momentum window |p|<=P
and |t|<=T,

```
||Psi_a(f,t)-Psi_0(f,t)|| <= T delta_a ||f||_2.
```

No bosonic moment truncation is needed: smeared CAR fields are bounded.
For finite volume take momentum mesh h=2pi/L and embed each discrete mode
as its normalized constant function on the centered momentum cube. Let
P_h be the orthogonal cell-average projection. Successive one-dimensional
mean estimates give

```
||f-P_h f||_2 <= sqrt(3) h ||gradient f||_2.
```

The explicit frame derivatives satisfy `||partial_i[U V]||<=2/m`, and
`|partial_i omega|<=1`. Comparing the center of each cube to its points
therefore adds at most
`(3h/2)(2/m+T)||f||_2`. With support inside the retained Nyquist cells,
the total finite-mode field bound is

```
epsilon_f <= sqrt(3) h ||gradient f||_2
           +(3h/2)(2/m+T)||f||_2
           +T a²(P+sqrt(3)h/2)^4 ||f||_2/(24m).
```

For a noncompact Fourier test add its omitted L2 tail. Sobolev data give
`||1_(|p|>P) f||_2<=P^-s || |p|^s f||_2`. This handles localized smooth
position-space tests without asserting a hard cutoff is local.

For the charged swap generator, the norm of the difference is at most
the norm of the difference of d: squaring it yields the two positive
blocks `n_r delta_d delta_d†` and `(1-n_r)delta_d† delta_d`.
Unitary Duhamel comparison then gives gate error at most
`|theta| epsilon_f/C`. Telescoping a sequence and retaining every outcome
gives complete total-variation error

```
TV <= epsilon_preparation + sum_j |theta_j| epsilon_(f_j)/C_j.
```

All omitted modes are vacuum spectators in the common exterior embedding.
The native complement is covered by the full preparation bound and by
complete bounded instrument extensions. The same actual H_a, state,
references and observation times occur throughout this comparison.
Finite apparatus coefficient errors and complete preparation failure add
the terms specified below. For two observer settings related by a fixed
continuum rotation or boost, apply the bound to each transformed test;
their record discrepancy is at most the sum of the two bounds. The
boosted tests use the enlarged momentum window and their actual Sobolev
norms, including the Jacobian factor. No exact finite-box boost is used.

## An exact finite record example

Choose a periodic source with `L=8pi/(3m)`, side N=32 and
`a=pi/(12m)`. Its mode `p=(0,0,3m/4)` is exact. For
`f=(component1+component3)/sqrt2`, the actual field mode is

```
d=(2 a_up(p)+b_up†(-p))/sqrt5.
```

Use two empty charge+1 reference doublets, swap angle pi/4 on each selected
mode,
and actual elapsed source time `t=pi/(2 Omega_a(p))`. The opposite-spin
source modes, the two unused apparatus spin partners and every other
momentum remain retained vacuum spectators
for this experiment. The literal reducing source-plus-reference circuit
has sixteen states, inside the complete 256-state four-source/four-reference
mode carrier. The finite symmetry lift preserves its reference-relative
record law exactly. All gates preserve charge and parity, and exact
matrix multiplication gives

| First reference | Second reference | Ideal actual-ground probability |
| --- | --- | --- |
| 0 | 0 | `313/500+12sqrt(2)/125` |
| 0 | 1 | `137/500-12sqrt(2)/125` |
| 1 | 0 | `1/10` |
| 1 | 1 | `0` |

The four displayed probabilities sum to one. The other twelve doublet
occupation records are retained with probability zero. They are exact finite source records,
not a continuum substitution. A physical approximate-ground preparation
enlarges each interval by its complete trace-distance error, including
the zero ideal record. The complete preparation failure record is kept
separately. The scalar mass m cancels from these dimensionless coefficients
because length and time settings are explicitly tied to the same source
E*. This does not calibrate m to measured data.

## Preparation and implementation boundary

The native algebraic trial and finite energy-agreement filter apply to
each unchanged constituent h. Their ground is also the exact constituent
factor of the unique H_a ground. The retained filter-tree proof controls
the main constituents jointly with all partner sources, spin references
and filter outcomes. For B=4N³ main constituents, initial excited odds r0
and geometric contraction c from the native preparation,

```
epsilon_preparation² <= 2 B r0 c^d,
Q_preparation >= product_(j=0)^(d-1)
                   (1+r0 c^j)^[-B 2^(d-j)].
```

Depth zero means deterministic unfiltered trial preparation, with success
one and the displayed nonzero state error. Every later filter and its
failure contributes to the full law. Native reference vacua, the original
preparation references and all unused source copies are retained. Preparing
the specified state is a declared apparatus intervention on the final
source, whose autonomous observation clock is H_a. The proof does not
identify a constituent filter time with autonomous H_a evolution.

This intervention has an explicit finite, driven implementation under a
stronger apparatus premise. Retain the entire source and preparation
carrier, and call its raw Hamiltonian H_ext. Main copies have H_a, partner
copies have their complete native constituent Hamiltonian, and spin,
filter and spatial-reference registers have zero Hamiltonian. Let U be
the coherent finite unitary that synthesizes every trial and performs
the full conditional filter tree, including every flag. Measurement can
be deferred, so U is a unitary on this complete carrier. For duration
tau>0 choose the Hermitian logarithm

```
L=(i/tau) Log(exp(-i H_ext tau) U),
V_control=L-H_ext,
exp[-i tau(H_ext+V_control)]=exp(-i H_ext tau) U.
```

The logarithm uses eigenphases in (-pi,pi]; it is scalar on each repeated
eigenphase eigenspace. Then L is Hermitian and `||L||<=pi/tau`. Every
conserved Gauss, charge and parity operator commuting with both H_ext
and U also commutes with L. This is a specified programmable finite-matrix
Hamiltonian service. It supplies no autonomous local controller and no
compilation from the passive C/E law. Its coefficient evaluation and
workspace beyond the displayed registers are separate apparatus resources.
A realized L with certified norm error eta changes the pulse by at most
tau eta; these errors add to the complete record bound. An exact service
is a mathematical premise, not a claim of infinite-precision binary hardware.

Use preparation duration one in source time units. Its final common free
factor preserves the joint conditional approximation to the main product
ground, because the ideal main state is stationary and every partner is
retained. The total time and reference phases are not reset. The observation
origin labels the end of this physical preparation pulse.

The same construction implements each bounded detector gate W. For the
finite example write `T=pi/(2Omega_a(p))` and give each detector pulse width
tau_d=T/4. Use the first pulse, a bare source wait T-tau_d, and the second
pulse. Their full propagator is

```
exp(-i H_ext tau_d) W2 exp(-i H_ext T) W1.
```

The last free factor commutes with every pointer effect and therefore
leaves the complete record law unchanged. Thus the table uses the actual
separation T between pulse starts. Physical elapsed time, including the
preparation, is `1+5T/4`. This implements the two-time table with an
ordinary autonomous clock between, during and after interventions.
The zero-separation control is a single pulse targeting W2 W1 and retains
the same two pointers; it does not require two nonzero-width pulses to
start simultaneously.

For B=4N^3 main blocks, d filter levels and R retained attempts, a full
kinematic raw Hamiltonian norm bound is

```
B_ext=R [B 2^d (B_native+120)+(N/2) U B],
B_native=8*57+24+360*2^-34.
```

Here 120 is the exact largest old/dual reference cost. The kinetic term
has norm at most `||K_a|| sum_j a_j†a_j <= (4/a)B`, on every source sector;
at `a=8pi/(3Nm)`, pi>3 and m<=U give the displayed bound. H_ext is the raw
Hamiltonian, without subtracting its ground constant. Consequently the
preparation control norm is at most `B_ext+22/7`; each detector control
norm is at most `B_ext+10U`, since pi/tau_d=8Omega_a<=10m.

Control norm is not energy work. Switching a constant V on and off gives
net work equal to the change in the bare H_ext expectation, with absolute
value at most 2B_ext for either a pulse or the entire schedule. Its sum
of absolute switching contributions is bounded by 2||V|| per pulse,
so the three-pulse bound is
`2[(B_ext+22/7)+2(B_ext+10U)]`. All pointer readouts commute with H_ext;
their complete nonselective readout does not change its mean energy.
These bounds concern ideal prescribed drive and source energy. They do
not provide a thermodynamic work reservoir or charge for its efficiency.

For any finite accuracy target, finite depth gives the required conditional
error and a strictly positive explicit success bound. For unconditional
preparation accuracy, R independent complete attempts leave failure at most
`(1-Q_preparation)^R`; select the first success using a retained record and
keep all failed attempts. Finite R can make this probability arbitrarily
small. The controller and retained-source cost of this amplification is
part of its declared resource count. A limiting conditioned-state theorem
alone is not an unconditional success theorem.

For N=32 and d=3 the conditional retained-purification error is less than
2.336e-10, and one complete attempt succeeds with probability at least
0.97038569. Six independent attempts, all retained, reduce complete
preparation failure below 7e-10. The complete law has every original
filter label, the first-success index, and all sixteen detector labels.
The displayed failure category is a coarse graining of retained labels.
Selecting a successful source copy uses a controlled readout with the
same H_a on every copy. Its total record error relative to the exact
ground test is at most the conditional error plus the all-attempts-failed
probability. This permits unconditional finite-record refinement, with
finite but growing explicit storage and the declared programmable control.

One complete constituent has kinematic dimension
`D_const=2^240 5^8 119^24 (2^31+1)`. For one attempt the total displayed
log dimension, in nats, is

```
B 2^d [ln(D_const)+16ln2] + B(2^d-1)lnM + ln(48N^3)+4ln2.
```

The 16ln2 term retains the native trial's eight four-state spin references;
each filter node has its full M-state pointer. The final terms retain the
regular spatial reference and two full apparatus doublets. For R attempts
multiply the two preparation terms by R and add ln(R+1) for the materialized
first-success index. This is the full kinematic capacity, hence an upper
bound on the constrained physical dimension, and is not state entropy.
N, d and R are declared geometry and accuracy choices, not outputs of a
memory-budget scheduling law. Additional controller synthesis workspace
is excluded and labelled separately.

Spectral matrix units, source-selected E*, square-root kinetic matrices,
normalized finite test modes and their bounded gates are exact finite
services. A numerical execution on the full native carrier and efficient
local gate synthesis are not claimed. The local carrier, lattice size,
preparation tree, pointers and retained references are finite for every
member. The free comparison has no surviving interaction and no Standard
Model or GeV identification.
