# One filtered observable algebra on the physical Wilson source

## Source and reused results

For every integer n>=0 use exactly `particle_common_contract.contract.member(n)`:
the open Wilson source has spacing `a_n=2^-n`, side `2^(2n+1)`, per-link
native word degree `n+1`, all60CAR modes per cell, electric coefficient105
and magnetic coefficient1/7. Let `P_n` be its complete Gauss projector,
including the declared zero-charge exterior endpoints. On the finite
physical carrier `H_n=ran P_n`, use the selected positive-polynomial density
`rho_n`, the source Hamiltonian `H_n^src`, and

```
alpha_n,t(A)=exp(it H_n^src) A exp(-it H_n^src),
[rho_n,H_n^src]=0,
Tr rho_n(H_n^src-E_phys,n) <= epsilon_n=2^(-n-4).
```

The superscript distinguishes the source Hamiltonian from its Hilbert
space. All operators in this construction preserve the physical carrier.
The arbitrary joint source/reference algebra is excluded: the retained
canonical purification with idle reference Hamiltonian0 need not be
stationary under the source clock. No change of clock or scalar energy
calibration is made. Subtracting `E_phys,n` in a bound changes neither
`alpha_n` nor a recorded energy difference.

The common-moment/GNS argument is the argument in
[`continuum_field_net`](../continuum_field_net/README.md#one-spacing-subsequence-for-the-complete-field-hierarchy),
applied to these actual physical finite sources. That package's CAR
relations, translation invariance, noncollapse and separate Weyl/contact
source are not transferred. The finite spectral tail and filtered norm
bounds are reused from
[`positive_energy.py`](../particle_joint_reconstruction/positive_energy.py)
and its [proof](../particle_joint_reconstruction/POSITIVE_ENERGY.md).
`Lean/QFT/GroundRecordTail.lean` proves the finite downward-weight lemma.
`Lean/QFT/ColimitGNSRepresentation.lean` records the general positive
functional/GNS interface. `ColimitSelectedStateGNS.lean` requires a coherent
finite colimit family; no such inter-level state maps are assumed here.
The sequence algebra below supplies its own positive functional and C*
domain. This proof is a mathematical extension, not a claim that a Lean
module formalizes this entire primary sequence construction.

## Countable physical probes and uniform time regularity

Use centered coordinates `x_n(c)=a_n(c-side_n/2)`. For each rational
`x in Q^3`, choose `c_n=floor(side_n/2+x/a_n)` componentwise. If this cell
is outside the interior `1<=c_j<=side_n-2`, define the probe to be zero.
Every fixed x is inside for all sufficiently large n. The coordinate
choice puts bounded x a diverging physical distance from the box walls;
it supplies no spatial translation symmetry of the state.

Two countable menus of concrete probes are available:

* For each coordinate axis and integer q, use the outgoing-link right
  endpoint spectral projection at `Q6=6Y=q`. Its exact finite spectrum is
  the native word-charge sumset from
  [`particle_physical_boundary_state`](../particle_physical_boundary_state/PROOF.md).
  An absent q gives the zero operator. A central right action commutes
  both endpoint gauge actions and all other Gauss transformations. The
  projection is even, gauge invariant and supported on the link's owner
  cell, with norm at most1.
* For each coordinate plane, use `ReTr_M15(U_p)/15`, formed on the parent
  complete plaquette and then compressed by the product word cutoff. It
  is an even gauge-invariant contraction. Its owner support consists of
  the anchor and the two forward neighboring cells. Restriction to `P_n`
  preserves its norm bound. It is an actual magnetic observable of the
  selected Wilson source.

Neither menu is claimed to generate the complete physical observable
algebra. The theorem also admits any explicitly supplied countable family
`A_n,k` of physical operators with uniform finite norm caps `C_k`. Taking
real and imaginary parts, adjoints or finite products is permitted.

For each probe choose a fixed positive rational T and the triangular filter

```
f_T(s)=(1-|s|/T)_+/T,
B_n,k=integral f_T(s) alpha_n,s(A_n,k) ds.
```

These are finite-dimensional Bochner integrals in the original source time.
`||f_T||_1=1` and `||f_T'||_1=2/T`. Therefore

```
||B_n,k|| <= C_k,
||alpha_n,h(B_n,k)-B_n,k|| <= min(2C_k, 2C_k|h|/T).
```

The real filter commutes with taking adjoints. It preserves the physical
carrier, and gauge invariance for gauge-invariant input probes. A filtered
projection is a positive contraction; it is generally not a projection.
The Bochner observable construction is not a compiled apparatus gate.
For a finite word with factor caps C_j and filter widths T_j, telescoping
gives the common-time Lipschitz cap

```
||word|| <= product_j C_j,
Lip(word) <= sum_i (2C_i/T_i) product_(j!=i) C_j.       (1)
```

Zero factors cause no division. Fixed translates and adjoints preserve
these estimates. Finite linear combinations inherit the corresponding
absolute-coefficient sums.

## One C* algebra and one state for every finite word

Let `M=product_n B(H_n)` be the bounded operator sequence C* algebra with
norm `||a||=sup_n||a_n||`, coordinatewise products and adjoints. Each
filtered probe defines a section B_k in M. Let B be the unital C*
subalgebra generated by all `alpha_q(B_k)` for rational q, all k, and
their adjoints. This is a separable C* algebra: polynomials in a countable
generator set with complex rational coefficients are norm dense.

Coordinatewise `alpha_t` is an isometric *-automorphism of M. Equation(1)
implies that any real translate of a generator is the norm limit of
rational translates. Thus `alpha_t(B)=B`. The group is norm continuous
on every generator, hence on finite polynomials by telescoping, and on
all of B by norm approximation and isometry. Operator norm continuity
here means point-norm continuity `||alpha_t(b)-b||->0` for each b; it
does not assert continuity in the norm on the space of automorphisms.

Set

```
I=B intersect {a: ||a_n|| ->0},       A=B/I.
```

I is a closed two-sided star ideal, invariant under every alpha_t. It is
proper because the unit section has norm1 at every n. A is a unital
C* algebra with a point-norm continuous induced time group. A norm-null
relation is consequently an operator relation before selecting a state.
The quotient uses the full sequence, so it can retain distinctions erased
by the state chosen on a subsequence. No faithful representation follows.

Each `omega_n(b)=Tr rho_n b_n` is a normalized positive functional on B
with norm1. Extract a diagonal subsequence on a countable norm-dense set
of rational-coefficient star polynomials. Every such expectation lies in
a compact disk. The common norm bound extends convergence to every b in B.
The resulting functional omega is linear, positive and normalized:
positivity follows by passing `omega_n(b* b)>=0`. For a in I,
`|omega_n(a)|<=||a_n||->0`; omega therefore descends to A. Exact finite
stationarity gives `omega(alpha_t b)=omega(b)` for every real t.

This single subsequence makes every finite multi-time star-word moment
converge, for arbitrary real times. The rational-time word limits extend
using(1); alternatively every fixed real-time word is an element of B.
The convergence is uniform on compact time sets for each finite word,
by the common Lipschitz bound and a finite time mesh. General elements
follow by uniform norm approximation. No effective subsequence-selection
algorithm, correlation convergence rate, unique limit, inter-level
embedding or exact finite state compatibility is inferred.

## Canonical GNS time evolution and positive energy

Construct the GNS Hilbert space from the positive state on A. Its null
space `N_omega={a:omega(a* a)=0}` is a **left** ideal and need not be
two-sided. It is used for the pre-Hilbert quotient, not for a C* algebra
quotient. Bounded left multiplication gives the star representation pi;
`ker pi` is a two-sided ideal. These ideals differ from the norm-null
sequence ideal I.

Stationarity makes

```
U(t)[a]=[alpha_t(a)],   U(t)Omega=Omega,
U(t) pi(a) U(t)*=pi(alpha_t(a)),
F_a(t)=omega(a* alpha_t(a))=<pi(a)Omega,U(t)pi(a)Omega>.
```

The point-norm continuity on A makes U strongly continuous on its dense
cyclic domain, then on the Hilbert completion. This is the canonical
implementer fixing the cyclic vector. An arbitrary implementing group
can carry a scalar phase and is not the group specified here. Write
`U(t)=exp(itL)` with this sign convention.

For any section a in B, `||a_n||<=||a||`. Let mu_n,a be the positive
energy-transfer measure of `omega_n(a* alpha_t(a))`. Its total mass is
at most `||a||^2`. The reused finite theorem gives, for every delta>0,

```
mu_n,a((-infinity,-delta]) <= ||a||^2 epsilon_n/delta.  (2)
```

This uses the ground floor of the same physical carrier. A bare charged
operator mapping out of that carrier cannot use(2) without the appropriate
compensating source and its own ground floor. Arbitrary divergent field
renormalizations are also outside the uniform norm domain.

The single extracted state gives pointwise convergence of all F_n,a(t)
to the continuous F_a(t), including t=0. The finite-positive-measure
continuity argument in the reused positive-energy proof applies: positive
limiting mass allows normalization for this proof, and zero limiting mass
implies convergence to the zero measure directly. The measures converge
weakly. Testing against compactly supported nonnegative functions in
negative energy and using(2) gives `supp mu_a subset [0,infinity)`.
Hence every negative spectral projection of L vanishes on every
`pi(a)Omega`. These vectors are dense, so **L>=0** on the entire GNS
space. No fitted energy shift, gap, nonzero spectral residue or interacting
field is inserted in this argument.

## Filtered spatial commutation on the same quotient

The all-state product-cutoff Wilson estimate from
[`particle_primary_locality`](../particle_primary_locality/PROOF.md) gives
for separated even physical observables of unsmeared owner supports X,Y

```
||[alpha_n,t(A_X),A_Y]||
 <= 2 C_X C_Y |X| exp[-((log2) a_n D-2W|t|)/a_n],
W=5232/7,  D=dist_graph(X,Y).                          (3)
```

The onsite electric interaction picture proves the ambient bound before
restriction to the physical carrier. No globally projected Hamiltonian
`P_n H P_n` is falsely decomposed into local terms. This is the actual
source clock. Filtering both probes and integrating(3) replaces |t| by
`|t|+T_X+T_Y`. This enlargement is necessary because the observables
have temporal support throughout their filter intervals.

For the centered concrete menus with assigned target separation d in the
Manhattan metric, flooring and the one-cell plaquette offsets give actual
physical owner separation at least `d-8a_n`. Suppose, more generally,
`|X_n|<=K a_n^-p` with integer p>=0 and the same separation bound is proved.
Define

```
m=d/2-2W(|t|+T_X+T_Y)>0.
```

Using `log2>=1/2`, for `a_n<=m/8` the exponent's margin is at least m/2.
Since `exp(x)>=x^(p+1)/(p+1)!`,

```
||[alpha_n,t(B_X),B_Y]||
 <= min(2C_XC_Y, 2C_XC_Y K (p+1)! a_n/(m/2)^(p+1)).  (4)
```

Thus these commutator sections belong to I. Their quotient operators
commute exactly whenever the stated conservative enlarged-cone condition
holds. The same relation holds in the GNS representation and extends to
generated local polynomial algebras by products and norm closure. A
general supplied probe menu must separately prove its support/separation
bounds; algebraic positive energy alone supplies no spatial locality.
This is an even observable algebra, not a construction of charged
fermionic fields. The bound supplies no Lorentz covariance or observed
light speed. It also supplies no spatial translation automorphism or
translation-invariant state.

## Complete finite detector records on the common algebra

For each filtered charge projection B, use `E=B`, with `0<=B<=I`. For a filtered normalized
Hermitian plaquette B, use the effect `E=(I+B)/2`. Continuous functional
calculus in the sequence algebra gives the complete two-outcome Kraus pair

```
K_0=sqrt(E),  K_1=sqrt(I-E),  K_0* K_0+K_1* K_1=I.
```

Coordinate evaluation is a unital star homomorphism, so these sequences
are exactly the finite matrix square roots of their coordinate effects.
Each pair preserves the physical carrier. The commuting square roots
give a finite source/pointer unitary with block matrix
`[[K_0,-K_1],[K_1,K_0]]`. It is a declared finite mathematical apparatus,
not a claim that a compiled efficient source implementation is executed.
A coherent linear-combination instrument with successful Kraus operator B
has effect `B*B`; it is a different instrument from this square-root
effect readout. No equality of those pointer laws is assumed.

For any fixed finite number of such instruments at specified real times,
every history Kraus product K_r belongs to B. Its original, unnormalized
probability is `omega_n(K_r* K_r)` and converges on the same subsequence.
Repeated completeness gives `sum_r K_r* K_r=I`, so all outcomes, including
the complementary outcomes, sum to one both before and after the limit.
For coherent pointer recording the density entries are
`(sigma_pointer,n)_(r,s)=omega_n(K_s* K_r)`. Every entry converges on the
same subsequence, hence the fixed finite-dimensional pointer density
converges in trace norm. Classical recording retains the complete diagonal
law. Source disturbance is included through the ordered Kraus products;
no product of independent one-time probabilities is substituted.

The finite source experiment can retain its full conjugate reference.
The displayed pointer marginal uses its mathematical trace and the same
source density; it does not claim convergence of the joint state with a
reference whose dimension grows with n. Fixed finite outcome histories
also differ from a growing complete charge PVM. Its charge-tightness
question is separate. Instrument count, filter widths and time labels
are fixed in this consequence; no uniform growing-history error rate is
derived.

## What the construction does and does not supply

The selected physical sequence has a common subsequential bounded
observable C* algebra, stationary state, continuous time group and
nonnegative canonical GNS energy, with scoped filtered spatial
commutation. Each finite construction retains its full physical carrier,
source clock and state. The existence theorem requires no exact maps
between different finite Hilbert spaces.

The unit survives, but the chosen probes can become scalar or zero in the
state representation. A fixed filter can remove spectral weight escaping
to arbitrarily high energy. The finite charge PVM is complete at each
level, while its growing outcome set need not yield a complete infinite
fixed-charge PVM: charge probability can escape to labels growing with n.
No exchange of this infinite sum with the limit is made. A retained-charge
tightness bound would justify that additional claim.

No nonzero finite-energy response, surviving interaction, faithful field
algebra, chiral species, particle mass, boost representation, scattering
limit, geometric stress tensor or Einstein dynamics follows merely from
this compactness construction. Unbounded renormalized fields require
their own norm/domain and response estimates. The comparison source in
`continuum_field_net` has a nonzero CAR identity that proves its fields
survive; that identity has not been established for this physical menu.

## Executed independent controls

`exact_two_level_control` uses the separate source
`H=diag(0,pi)`, ground `rho=|0><0|`, plus projection and filter width1.
Direct integration gives filter amplitude `4/pi^2`. At rational times0
and1/2 its two translated matrices give six finite star words. Their
complete Gram matrix factors as `V*V`; common time shifts preserve it.
The canonical two-dimensional GNS has generator diag(0,pi), spectral
weight `4/pi^4` at positive frequency pi and weight1/4 at zero. A matrix
unit supplies an explicit null left ideal element whose right product
has nonzero GNS norm. This checks the distinction between the GNS null
space and a two-sided quotient ideal.

`complete_pointer_control` executes two noncommuting two-outcome effect
instruments at times0 and1/2 on a separate two-state ground source. The
square-root Kraus operators and their two-pointer dilation are exact.
All four ordered histories are retained, their Kraus completeness is
checked, and their pointer matrix has an explicit Gram factor. The
source is traced only to display that pointer marginal.

Adversarial tests retain the escaping-energy and escaping-charge cases,
the zero-word bounds, nonphysical reference domain, changed-clock receipt
rejection and the temporal halo in(4). These controls check the theorem
mechanisms; their small matrices are not the selected spatial matrices.
