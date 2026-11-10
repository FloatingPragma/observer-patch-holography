# Massive-vector records from the same interacting native blocks

An explicit change of bounded readout produces a controlled free Proca
field comparison from the three-component scalar comparison of the
frozen native quadrature source. The finite source Hamiltonian, selected
actual ground, retained reference and clock remain the same. The new
readout has a time component built from canonical momentum, an exact
polarization sum, a positive-energy spin-one induced action and a
microcausal limiting commutator. The unfiltered scalar-component labels
are not reinterpreted as a vector field.

The result concerns complete finite apparatus records and their free
massive-vector comparison. It supplies no measured Z or other particle
identity, no interacting relativistic ultraviolet limit, and no exact
boost symmetry of a finite native lattice.

## One source, state, reference and clock

The immutable parent is
`model/interacting_block_quadrature/results.json`, with digest
`64f90635637575bdaa9a24f1840ea3abb6362d6de8bc28e97de7fb275e8218e3`.
The code validates every inherited source pin and re-executes the actual
native collective-matrix and translation/cubic reference admissions.
One constituent is the complete eight-cell Hamiltonian with native
hopping `2^-34`, actual ground g and first neutral triplet with gap Delta.
Its actual currents supply the exactly dressed matrix units
`P_g O_a P_t/(60 sqrt(r))`. All other physical and kinematic native
sectors remain retained.

For each declared spacing a, periodic macro side L and integer capacity
N, use the parent's full physical Hamiltonian

```
H_(a,N)=a^-1 sum_(x,i)(h_(x,i)-E0)
          -t_a sum_<xy>,alpha X_xalpha^(N) X_yalpha^(N),
t_a=Delta(1-a^2)/(12a).
```

The state is the same selected symmetric actual global ground admitted
by the parent, with its existing regular-reference isometry J. The
source/reference clock is `exp(-it H_(a,N)) tensor I`; reference energy
is zero. The replica symmetry, full complementary-sector comparison
and actual-ground energy bound are the parent's theorems on this exact
family. A resource choice below selects N before the experiment. It
does not substitute a free Gaussian state for the actual native state.

For formulas use source-derived units

```
tau=Delta t,       xi=sqrt(6) x,
Phi=sqrt(Delta) Q, Pi=P/sqrt(Delta).
```

This is a change of labels and canonical field normalization, not a new
clock or empirical calibration. The normalized site fields obey the
same canonical comparison algebra. Their continuum normalization is
`(sqrt(6)a)^(-3/2)` at spacing `sqrt(6)a` in xi. The limiting frequency
is `Omega_0(q)=sqrt(1+|q|^2)`, with rest frequency one in tau. Restoring
the original clock gives the actual source gap Delta. No independently
selected mass or observed speed enters.

With hbar=1 and limiting speed `c=Delta/sqrt(6)` in those original
length/time units, the effective free Proca mass parameter is
`Delta/c^2`. This is a source-unit convention, not a GeV value. The
positive gap is essential; no massless limit of the singular filters
below is taken.

The free comparison frequency of the actual finite lattice is

```
Omega_a(q)^2=1+(1-a^2)/(3a^2) sum_i[1-cos(sqrt(6)a q_i)].
```

The finite physical carrier consists of the parent's replicated native
blocks. Infinite oscillator and continuum spaces are mathematical
comparisons, not extra source hardware.

## The explicit time and spatial readouts

First retain a general positive rest frequency m and write
`Omega_0(p)=sqrt(m^2+|p|^2)`. For three canonical comparison fields Phi_a
and momenta Pi_a, define Fourier readouts

```
A^0(p)=i p_a Pi_a(p)/(m Omega_a(p)),
A^i(p)=M_ia(p) Phi_a(p),
M_ia(p)=delta_ia+p_i p_a/[m(Omega_0(p)+m)].
```

Repeated a is summed. The denominator in the time component is the
actual lattice comparison frequency. The spatial filter uses the target
frequency. Both are explicit coefficients on the finite momentum menu.
The resource calculations use m=1 in the source-derived units above.

The reality relations `Phi(p)^dag=Phi(-p)` and `Pi(p)^dag=Pi(-p)` imply
`A^mu(p)^dag=A^mu(-p)`. In contrast, replacing Pi in the time component
by Phi and omitting i would make its odd-momentum Fourier coefficient
anti-Hermitian as a real field. A literal seven-state vacuum/one-quantum
control at opposite momenta checks this sign and rejects that replacement.

At the limiting frequency, the position-space formulas are

```
A^0=div Pi/(m Omega_0),
A=Phi-grad(div Phi)/[m(Omega_0+m)].
```

They are nonlocal equal-time filters of the original scalar coordinates.
Their reconstructed limiting field nevertheless has a local commutator,
as established below. The time component is constrained; it adds no
fourth oscillator or indefinite-metric state.

For a positive-frequency mode, `Pi=-i Omega_a Phi`. Consequently the
four readouts have, already at finite spacing, the target polarization
columns

```
epsilon_a^0(p)=p_a/m,
epsilon_a^i(p)=delta_ia+p_i p_a/[m(Omega_0+m)].
```

With `eta=diag(1,-1,-1,-1)` and `p_0^mu=(Omega_0,p)`, the executed
symbolic identities are

```
p_0^T eta epsilon=0,
epsilon^T eta epsilon=-I3,
epsilon epsilon^T=-eta+p_0 p_0^T/m^2,
epsilon^T epsilon=I3+2 p p^T/m^2.
```

Thus there are precisely three positive polarization degrees of freedom.
The current unfiltered spatial triplet has equal residue eigenvalues.
The filtered spatial residue is `I+p p^T/m^2`, with longitudinal
eigenvalue `1+|p|^2/m^2`. These readouts are distinct away from rest.

The Lorenz constraint follows from `dot Phi=Pi` and
`dot Pi=-Omega_0^2 Phi`:

```
partial_tau A^0+div A=0.
```

At finite lattice frequency its exact defect instead is
`i(Omega_0-Omega_a)p.Phi/m`. No exact finite-lattice Lorenz equation
is claimed. The same dispersion error that controls the covariance
controls this defect on the stated finite-moment domain.

## Positive-energy vector action and its finite spatial ancestor

On the continuum comparison one-particle space use measure
`d^3p/[(2pi)^3 2Omega_0(p)]` and three components. For a proper
orthochronous Lorentz transformation Lambda define

```
R(Lambda,p)=-epsilon(Lambda p)^T eta Lambda epsilon(p).
```

The polarization columns form an orthonormal basis of the Minkowski
orthogonal complement of p. Hence

```
Lambda epsilon(p)=epsilon(Lambda p)R(Lambda,p),
R^T R=I3,    det R=1,
R(Lambda2 Lambda1,p)=R(Lambda2,Lambda1 p)R(Lambda1,p).
```

Orthogonality follows by inserting the polarization projector; the
longitudinal p term vanishes. The determinant is positive at identity
and cannot change on the connected Lorentz group. The cocycle follows
by applying the first intertwiner twice. These identities define the
unitary induced action

```
(U(Lambda)psi)(p)
 =R(Lambda,Lambda^-1 p)psi(Lambda^-1 p).
```

The measure is invariant, the matrix is orthogonal, and the cocycle
proves the group law. Second quantization leaves the same Gaussian
vacuum invariant and does not mix creation with annihilation. Its time
and space translation generators remain Omega_0 and p. In the convention
used here the field action is

```
U(Lambda)^dag A^mu(x)U(Lambda)
 =Lambda^mu_nu A^nu(Lambda^-1 x).
```

At rest the little group acts by its ordinary three-dimensional vector
representation. This establishes spin one for this positive-energy
limiting field. It is a different induced Lorentz action from the one
that treats the three original fields as unrelated scalars. The source
spectrum alone does not choose between those observer maps.

For a spatial rotation R, the induced matrix is exactly R and
`epsilon(Rp)R=diag(1,R)epsilon(p)`. The executable checks all twenty-four
actual proper native cubic rotations. On the finite source their
ancestor rotates macro positions and native vector-current components
together. Both actions preserve H_(a,N). The internal native X-point
labels remain distinct from the macro momentum p.

The existing reference also supplies this action on the closed relational
experiment. For its actual finite group G,

```
J psi=|G|^-1/2 sum_g U_g psi tensor |g>,
R_h|g>=|g h^-1>,
(I tensor R_h)J=J U_h.
```

This is a right action on the retained reference. It differs from the
joint left action `(U_h tensor L_h)J=J` expressing common-frame covariance.
Combining the macro-position rotation with this right reference action
implements the induced native spatial/vector rotation without adding a
new reference. The lifted operator transforms as the lift of
`U_h C U_h^dag`; the same identity applies to every complete Kraus
operator. The code verifies all 576 cubic right-action blocks and uses
the parent's actual cover of order 3072 for the full source.

No exact rotation invariance is assumed for an arbitrarily selected
finite ground within a degenerate ground space. Its admitted vector
error epsilon_initial from the invariant free vacuum instead bounds
its change under these rotations by `2 epsilon_initial`, after removing
any common ground-character phase. That phase cancels from the native
current action and all records. This controls the corresponding complete
finite records; the limiting state is exactly invariant.

The boost action is not declared an exact symmetry of finite native
hardware. The code executes a noncollinear rational boost whose Wigner
matrix is nontrivial, and checks its cocycle with a second boost. The
action above belongs to the comparison field. For rapidity bounded by
rho, a momentum ball of radius P maps inside radius
`exp(rho)*sqrt(m^2+P^2)`; its transformed smooth test functions therefore
remain within a finite enlarged readout window. Applying the bounds
below to both windows controls the covariance comparison. No physical
boost gate, reset of the reference or alternate source clock is inserted.

## Limiting covariance, equations and microcausality

The same Gaussian vacuum gives

```
W^{mu nu}(x)
 =integral d^3p/[(2pi)^3 2Omega_0(p)]
    [-eta^{mu nu}+p^mu p^nu/m^2] exp(-ip.x).
```

Positivity follows from its factorization `epsilon epsilon^T`.
The measure and polarization tensor give Lorentz covariance under the
explicit induced action. Each component satisfies the Klein-Gordon
equation and the Lorenz constraint. Therefore
`partial_nu F^{nu mu}+m^2 A^mu=0` with the usual antisymmetric derivative
F. The Gaussian higher moments are those of the same three-oscillator
Fock vacuum, not new assumptions about an interacting native field.

Let C_KG be the scalar commutator distribution with its canonical initial
data. The vector commutator is exactly

```
C^{mu nu}=(-eta^{mu nu}-partial^mu partial^nu/m^2) C_KG.
```

The scalar Klein-Gordon Cauchy problem has speed-one finite propagation:
the energy density `(Pi^2+|grad Phi|^2+m^2 Phi^2)/2` and flux
`-Pi grad Phi` have nonnegative outward null-cone energy. Integrating
their conservation law outside the expanding support cone proves that
smooth compact initial data vanish there. Approximating the canonical
delta initial data by compact mollifiers proves that C_KG is supported
in the light cone. Distributional derivatives cannot enlarge support,
so the vector commutator vanishes at spacelike separation. In particular
the nonlocal equal-time filter has not imported a nonlocal target
commutator. Direct canonical checks also give

```
[A^0(x),A^i(y)]=-i partial_i delta(x-y)/m^2,
[A^i(x),dot A^j(y)]=i(delta_ij-partial_i partial_j/m^2)delta(x-y).
```

A fixed momentum cutoff has spatial tails and is not exactly microcausal.
For P>=1 in normalized units, suppose Fourier tests obey
`||hhat(p)||<=C_h(1+|p|)^(-s)` and similarly for k, with s>5/2.
The omitted covariance tail is bounded by

```
3 C_h C_k P^(5-2s)/[4 pi^2(2s-5)].
```

Indeed the polarization norm is at most `1+2|p|^2`, `Omega_0>=1`, and
the radial tail is bounded by three times the integral of
`p^(4-2s)`. Twice this covariance tail controls the commutator tail.
Smooth compact position tests admit arbitrary such s with their stated
Fourier constants. Combining this bound with the finite-source and
lattice errors supplies approximate spacelike readout commutativity;
exact locality belongs to the full limiting distribution.

## Bounded filters and exact finite instruments

On `|p|<=P`,

```
||M||<=sqrt(1+P^2/m^2),
||p/(m Omega_a)||<=P/m^2,
||epsilon epsilon^T||<=1+2P^2/m^2.
```

The lower frequency bound is `Omega_a>=m`. Real four-vector smearings
therefore give real canonical smearings f,g for
`X(h)=f.Phi+g.Pi`, and the bounded test is `W_A(h)=exp(iX(h))`.
In normalized units, if the four-vector test has L2 norm at most H,
the separate coefficient bounds are
`||f||<=sqrt(1+P^2)H` and `||g||<=PH`. They refer to the actual discrete
Parseval norm on a finite menu; the continuum comparison uses its own
corresponding L2 norm. Both bounds are required for an admitted test.

For each finite source member use exactly the parent's contraction

```
C_N=P_N W_A(h) P_N,
```

extended by zero on the entire native complement. This is a finite
matrix specification. For example a one-mode entry is the finite formula

```
<n|W(eta)|m>
 =exp(-|eta|^2/2)sqrt(n!m!)
    sum_(k=0)^min(n,m)
      eta^(n-k)(-conj eta)^(m-k)/[(n-k)!(m-k)!k!].
```

The multimode entry is a finite product on the finite admitted occupation
menu. The formula is implemented and its adjoint identity is checked.
The unbounded comparison proves the contraction norm; it is not a
continuous physical carrier hidden in the native matrix.

The same finite Halmos unitary and controlled Hadamard construction
retain four records `(s,flag)`, with s=+1,-1:

```
K_(s,0)=[I+s exp(i theta)C_N]/2,
K_(s,1)=s exp(i theta)sqrt(I-C_N^dag C_N)/2.
```

The phase on the second line can equivalently be absorbed into its
record basis. The squared effects sum exactly to identity. Both flags,
both signs and every conditional source output remain retained. A
literal exact contraction control has strictly positive leakage in
both tested phases. Applying the parent's same regular-reference lift
to each Kraus operator preserves the whole joint law and every error.
The apparatus dimension remains `3072*4^r` for r tests. No new quantum
field, boost, or orientation register is added.

The resource receipt counts all 2r retained control/flag qubits and
the shared reference separately from the full native kinematic source.
It also gives a finite upper count `6rV` of canonical filter coefficients
and the admitted matrix dimension `binomial(N+3,3)^V`. Exact filter gates
and the parent's selected-ground preparation are declared finite
services. Their efficient circuit synthesis and additional implementation
workspace are not asserted or included in the logical source/apparatus
capacity. Real matrix coefficients specify finite gates, not additional
continuous quantum registers.

The exact logical capacity of the counted full kinematic source and
retained apparatus is, in nats,

```
log dim H=NV(240 log2+8 log5+24 log119)+log3072+2r log2.
```

One block contains eight thirty-mode CAR carriers, eight five-state
scalars and twenty-four 119-state links. Tensor-product dimensions
multiply, so this logarithmic formula is exact for that counted carrier.
Restricting to the Gauss-physical subspace can only reduce the dimension.
The finite-dimensional entropy bound `S(rho)<=log dim H` follows from
nonnegative relative entropy against its uniform density. This is the
capacity-law use made here. It does not identify the pure prepared
state's entropy with that ceiling or derive N from a scheduler law.

For the two explicit experiments, the outward capacity intervals are

| a, L | Full source plus reference and readout capacity, nats |
| --- | --- |
| 1/4, 64 | `[2.697181624591e100,2.697181624592e100]` |
| 1/8, 128 | `[2.685906163883e108,2.685906163884e108]` |

The receipt also keeps much sharper exact rational interval endpoints
and the source/reference/pointer terms separately. Both totals are
arithmetically below `10^122` nats. This is a numeric comparison only;
it supplies no universe-capacity or observed-continuum identification.
N is the declared number of constituent replicas selected by sufficient
record-error inequalities, not a number inferred from a memory budget.
Unspecified extra gate-synthesis and preparation workspace remain outside
the counted capacity.

## Actual-ground and sequence errors after filtering

In normalized source coordinates the bare oscillator has `D=1/a` and
`b=(sqrt(D)Phi+i Pi/sqrt(D))/sqrt(2)`. Thus the displacement amplitude
obeys

```
||eta||^2=||f||^2/(2D)+D||g||^2/2
 <=[a(1+P^2)H^2+P^2 H^2/a]/2.
```

The executable conservatively replaces `sqrt(1+P^2)` by `1+P^2/2`
before forming a rational bound. The growth of the momentum part as
a tends to zero is retained. The detector strength is not silently
sent to zero to keep an old displacement allowance.

The parent's complete-sequence moment theorem applies with

```
R=a^-1 sum_j ||eta_j||,
F(R)=1+8R+6R^2+4R^3+R^4,
M4_sequence<=F(R)^2 M4_base.
```

It controls the coherent superpositions of histories and their retained
controls, rather than assuming a stationary Gaussian moment after each
kick. For m evolution segments, r tests and total normalized duration T,
the full comparison error is

```
epsilon_initial
 +T[epsilon_HP+2C sqrt(M4_sequence)/(K+3)]
 +(2m+4r)sqrt(M4_sequence)/(N+3)^2.
```

Here K is the mathematical intermediate total-number cutoff, N is the
actual per-site replica capacity, and C and epsilon_HP are precisely
the parent's full-complement bounds. The controlled Halmos factor is
the parent's factor one on the leaked branch; the four-record tail
factor is retained in `4r`. An uncontrolled Halmos factor sqrt(2) is
not mistakenly substituted or dropped.

`resources.py` chooses K and N using exact rational moment bounds for
the actual a,L. It increases K until the Gaussian variational cutoff
energy is at most half the allowed ground tolerance, then chooses N
to control both the full native Hamiltonian difference and all tails.
The initial actual-ground vector error is included. With requested
dynamic accuracy epsilon, its ground allowance is epsilon/5 and the
full complete-record total variation is at most `19 epsilon/20`.
This works for every positive chosen accuracy and every finite admitted
test window. Source matrices at those capacities are not numerically
enumerated; their finite algebraic program and error bounds are supplied.

## Spacing, volume and complete record transfer

In normalized unit-speed coordinates the inherited dispersion bound is

```
|Omega_a-Omega_0|<=a^2[|q|^2/2+sum_i q_i^4/4]
                <=a^2(P^2/2+P^4/4).
```

Because the A0 denominator uses Omega_a, the finite free covariance is
exactly the target polarization tensor multiplied by the actual lattice
scalar kernel. At time difference at most T, a normalized pair of tests
therefore has covariance error coefficient

```
delta_cov <=(1+2P^2)(1+T)|Omega_a-Omega_0|/2.
```

This bound includes the time and mixed components; no separate
uncontrolled derivative limit is taken.

For a box of normalized side ell, let h=2pi/ell. With smooth Fourier
tests supported in `[-P,P]^3`, put

```
F(q)=hhat(q)^dag Pi(q) khat(q)
          exp(-i Omega_a(q)t)/(2Omega_a(q)).
```

The parent's midpoint-cell estimate gives volume error at most
`(2P+2h)^3(3h/2) sup_i,q|partial_i F|/(2pi)^3`.
Writing `L_P=1+2P^2`, the uniform derivative bound is

```
|partial_i F|
 <=L_P/2 (||partial_i h|| ||k||+||h|| ||partial_i k||)
   +[sqrt(2L_P)+L_P(1+T)P/2]||h|| ||k||.
```

It uses `||partial_i Pi||<=2sqrt(2L_P)`, `Omega_a>=1` and
`|partial_i Omega_a|<=|q_i|`. These bounds are uniform in a. Adding
spacing, volume and, when used, Fourier-tail estimates controls each
finite collection of target covariances.

For completeness, covariance convergence controls the actual bounded
record law, not merely unbounded second moments. Before the final
control rotations, r ideal tests have d=2^r orthogonal control histories.
Their reduced matrix is `d^-1` times the Gram matrix of Gaussian Weyl
histories. A pair of histories contains a product of at most 2r Weyl
operators with total test norm at most `2S`, where
`S=sum_j||h_j||`. If every covariance error is bounded by
`delta_cov||h_j||||h_k||`, the Gaussian exponential identity bounds
every Gram-entry error by `2 delta_cov S^2`. Its exponential has
nonpositive real part; integrating along the segment between the two
exponents gives this bound without an exponential growth factor.

The Frobenius norm of the reduced-control difference is at most that
Gram-entry bound, and its trace norm is at most sqrt(d) times its
Frobenius norm. Hence the complete ideal sign-record variation is at
most

```
2^(r/2) delta_cov S^2.
```

All ideal internal flags are zero. Adding the full actual-source bound
above includes every possible native flag and conditional history.
No identification of different infinite Fock vectors is required for
this finite record comparison. Taking a to zero, the physical box to
infinity, and K,N from these bounds gives the stated controlled Proca
record limit. Fourier tails can then be removed on their stated test
domain. The strict limit is free; finite-N hard-core interactions have
not been promoted to an interacting relativistic continuum.

## Explicit complete finite experiments

The receipt instantiates a=1/4,L=64 and a=1/8,L=128, both of source
physical side sixteen. In normalized momentum coordinates choose

```
hhat(q)=(1-16|q|^2)_+^2,
```

first in temporal component zero at tau=0, then in spatial component
one at tau=1/2. The total allowed duration is one. The actual finite
grid contains exactly nineteen supported modes: integer triples with
squared norm zero, one or two. Their discrete test norms are below
1/8. The zero mode and both signs of every nonzero momentum are retained.

The ideal temporal variance is the finite sum of
`|hhat|^2 |q|^2/(2Omega_a ell^3)`; the spatial variance uses
`|hhat|^2(1+q_1^2)/(2Omega_a ell^3)`.
Their cross covariance and commutator vanish exactly by reflection in
q_1 at every time. Consequently their ideal Hadamard sign probabilities
factor, with marginal `(1+s exp(-variance/2))/2`.
Outward interval arithmetic encloses every sine, square root and
exponential in these nineteen-mode expressions.

The actual-source planner is instantiated with epsilon=10^-9. Its
complete variation allowance is `19/(20*10^9)`, including the selected
actual ground, evolution, both filtered instruments and all flags.
The receipt encloses all sixteen two-test outcomes, bounds the total
probability of any flag-one record by that allowance, and proves a
strictly positive actual minus-sign probability for both A0 and A1.
These intervals enclose one normalized complete law; their endpoints
are not individually renormalized. The mathematical Gaussian sum is
executed exactly to intervals; the enormous full native matrix is
controlled analytically rather than claimed to have been enumerated.

## Identity and limitations

This is a new neutral, even, native-Gauss-invariant massive-vector
readout with source-only rest frequency. Its limiting one-particle
representation is spin one under the explicitly constructed vector
action. The same source also retains its distinct unfiltered scalar
readouts. The observable map, reference-relative orientation and induced
symmetry action therefore matter to the particle description.

Neither this field's name nor its three polarizations establish an
electroweak current, a Higgs mechanism, an observed mass standard or a
Z assignment. The affine-Gauss charged source is a different Hamiltonian
and receives no automatic coupling to this field. Every finite member
and apparatus is finite; the continuum field is the controlled record
comparison with the stated momentum, time, volume and capacity errors.
