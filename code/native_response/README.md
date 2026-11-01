# One N4 process: central ports, geometry and a non-native timing response

The [N4 candidate](../../model/CANDIDATES.md#n4-a-reciprocal-coherent-exchange)
retains native port geometry but has a certified memory signal beyond native
heat dynamics. A minimal twelve-block address carrier preserves its complete
quantum channel and supports a full abstract A1–A3 construction. The resulting
native/N4 pair has the same full quantum reference, collective response and
normalized port geometry, with different clock-independent kinetic invariants.
This proves finite kinetic nonselection within the stated abstract source
interpretation. Recording an occupied port has a separate, unavoidable quantum
information cost, proved below. The address register is an explicit extension;
the native/N4 interaction and private factors are unchanged.

The live [observer-local incident-tick calculation](../../model/local_tick/README.md)
reuses this packet's certified edge bound for a different, explicitly
defined stopping rule: every attempted seam incident to one port counts as a
tick, with intervening remote events retained. It derives the corresponding
native null and a positive N4 lower bound. The global-attempt calculation
below retains its own protocol; the extension supplies no physical duration
calibration. Runnable N1 observation, record and feedback protocols live in
[model/observer/](../../model/observer/README.md), with their distinct source
and apparatus assumptions.

## Fixed inputs and experiment

The state space is the binary total-two shell (dimension 66), with one
six-dimensional private factor per port. The response matrices are exactly
the formula in [SOURCE_RECORD_GLUING.tex](../../paper/tex_fragments/SOURCE_RECORD_GLUING.tex),
in the chart of [verify_response.py](../source_selection_model/verify_response.py).
Here the first neighboring pair is **(0,2)**; (0,1) belongs to a different
distance class. The seam set is derived from this chart, not imported from
another labeling. The existing `response.json` independently binds all twelve
matrices and sixty proper carrier actions.

The candidate inputs are its binary quantum representation, response coupling
θ, equal quadrature weights, exponential exchange, the full normalized-trace
state τ, and independent uniform choices among thirty seams. The measurement
input is access to the unchanged commuting occupations N_p. These are not
consequences of A1–A3 or physical calibration. Neither mass, inertia, a bulk
boundary condition nor a new field is introduced by this calculation.

The two rational couplings reuse the preceding feedback test; neither is
fitted to this result. The discriminator is strict positivity of the
full centered-subspace memory defect after error enclosure, against exact
zero for the native and one-term controls. The compute budget uses local
72-dimensional matrices and three representative edges with protected
18-dimensional blocks, rather than constructing the full network matrix.

In independent trials, measure all initial binary loads, retain the outcome,
apply n unmodified N4 events, then measure all final loads. Use n=1 and n=2;
there is **no intermediate measurement, private reset or postselection**.
The first unconditioned measurement preserves τ. The Lüders joint law gives
`E[f(X_0)f(X_n)] = τ[f Φⁿ(f)]`, where Φ discards the event records. Thus the
following moments are operational moments of one declared instrument, not
three separately chosen geometry, response and timing fixtures. Branch label
j retains its transformed-occupation meaning; it is not relabeled a bare
current. Event count, or a supplied Poisson attempt rate, is model time.

## Geometry theorem, including all couplings

Use normalized Hilbert–Schmidt inner products throughout. On the odd
two-endpoint load block, `Z=diag(1,−1)`, private dimension is 36, and

```text
g=D_b−D_a, h=D_a+D_b−J/6, J=ΣD_p,
G=g⊗I+I⊗g, H=h⊗I+I⊗h, B=H−iG,
T=[[0,B],[-B†,0]], W=exp(θT),
R=Tr_load(W† Z W)/2, Q=Φ_e(Z)=W(I_2⊗R)W†,
β=Tr_36(R²)/36, α=1−β.
```

C and E coincide on this input: their odd-block action is the conditional
expectation `X ↦ I_2⊗Tr_load(X)/2`. The even block is fixed. Locality and
occupation conservation force each seam contribution to the bare Dirichlet
matrix to be an incidence outer product. The odd-shell probability is
`20/66=10/33`; `N_a=(N_a+N_b+Z)/2`. Proper source covariance makes α equal
on all thirty seams. Therefore

```text
K_pq=τ[N_p(I−Φ)N_q]=α L_ico,pq/396,
Cov_τ(N)=(5/33)(I−11ᵀ/12),
P(I−Φ)P on centered linear loads = α L_ico/60.
```

Here P is the orthogonal projection onto linear occupations. The last line
is a compression, not yet an invariant-subspace claim. Direct expansion of
the endpoint joint law gives `E[ΔX ΔXᵀ]=2K`; hence `60K/Tr K=L_ico`.
Four times its smallest-positive-eigenvalue projector is the source Gram
`v_p·v_q/||v_p||²`, of rank three. Selecting this band is the existing finite
port readout, not a derivation of physical spacetime dimension.

Each endpoint's protected 3+3 splitting leaves four private blocks of
dimension nine. On each, R is a traceless Hermitian contraction: contraction
follows from normalized partial trace of the involution `W†ZW`, and its trace
vanishes because W preserves that block. Nine eigenvalues in [−1,1] summing
to zero have sum of squares at most eight. Maximize the convex sum on this
polytope: at a vertex eight entries are endpoints and the last must be zero.
Consequently **`1/9 ≤ α ≤ 1` for every real θ**. No singular coupling loses
the finite port geometry.

There is also a useful quantitative bound. The block exponential gives
`R=[f(BB†)−f(B†B)]/2`, `f(x)=cos(2θ√x)`. For x≥0, `|f'(x)|≤2θ²`.
In eigenbases of Hermitian A,C,
`||f(A)−f(C)||_F²=Σ_ij |f(a_i)−f(c_j)|² |<u_i,v_j>|²`, which bounds this
by `4θ⁴||A−C||_F²`. The source identity
`Tr_6(i[g,h])²=49+19√5` and `Tr_6(i[g,h])=0` thus give

```text
α ≥ max{1/9, 1−(4/3)(49+19√5)θ⁴}.
```

Its quartic coefficient is also the exact leading coefficient of `1−α`;
the asymptotic expansion is not used as a finite-angle certificate.

## Exact memory reduction and the decisive counterexample

For centered v with `L_ico v=λv`, define `f_v=Σv_pN_p` and
`c_n=τ[f_v Φⁿ(f_v)]/τ[f_v²]`. Self-adjointness of Φ and orthogonality of P
give, without an approximation,

```text
c_1=1−r, r=αλ/60,
c_2−c_1²=||(I−P)(I−Φ)f_v||²_τ/||f_v||²_τ = b_v.
```

Under Poissonization at rate one, the initial decay is r and the initial
curvature is `r²+b_v`. The ratio `b_v/r²` is unchanged by a common clock
rescaling. A closed native heat mode has b_v=0. Positive b_v excludes that
law on the same observables; it does not rule out a larger closed system
with its private degrees of freedom retained.

The full-network norm reduces to local overlaps. Define the edge residual
`S_e=Q_e−βZ_e`, extended by zero outside odd occupation, and let ∂ have
oriented edge rows. Each diagonal load block of S_e has zero private trace.
Its normalized Gram matrix B on the full reference is:

* `B_ee=(10/33)Tr_72(S_e²)/72=(10/33)β(1−β)`.
* Disjoint edges have zero entry.
* For edges meeting at p, reduce each diagonal load block by normalized
  partial trace over the other endpoint. Write the resulting six-dimensional
  matrices as `M_e^1,M_e^0` according to p's occupation. Then
  `B_ef=[9 Tr_6(M_e^1 M_f^1)+Tr_6(M_e^0 M_f^0)]/(66·6)`.

The weights count full-shell configurations: when p is occupied there are
nine places for the other particle; when p is empty the other two endpoints
must both be occupied. Distinct-edge off-diagonal load blocks have zero
product trace. Same-edge coherences remain in the complete local norm. Thus
the leakage quadratic form is exactly `∂ᵀB∂/3600`, and dividing by the bare
covariance 5/33 gives b_v for unit coefficient vectors. This is a locality
reduction of the original network, not an independent mean-field model.

At θ=1/8 and 1/16, the interval calculation certifies strict diagonal
dominance of B with positive diagonal. Since the carrier is connected,
**every nonconstant bare carrier vector has strictly positive leakage**.
For its rank-three band, approximate `(α,b_v,b_v/r²)` are respectively
`(0.975986704560, 4.178223344e−5, 0.0206705310)` and
`(0.998236412623, 3.298779751e−6, 0.0015600348)`. Positivity is certified;
these printed decimals are numerical summaries. Both one-term exchange
controls have normal B, hence R=0, α=1 and zero leakage exactly.

## A corrected single-port readout

The same result admits a local scalar readout without reconstructing a carrier
eigenmode. Keep the tracial preparation and IID uniform attempt schedule, fix
a port p, and put `f=N_p−1/6`, `c_n=τ[f Φⁿ(f)]/τ[f²]`, and `A=I−Φ`.
The carrier identities `L_pp=5`, `(L²)_pp=30`, the covariance above and
`τ[f²]=5/36` give

```text
r = 1−c_1 = α/11,
||P A f||²_τ / ||f||²_τ = α²/110,
b_local = ||(I−P) A f||²_τ / ||f||²_τ
        = c_2 − c_1² − (1−c_1)²/10.
```

This is a corollary of the existing compression and leakage proof. The
uncorrected `c_2−c_1²` contains `α²/1210` from mixing carrier bands;
native repair already gives `1/1210`, despite closed mean dynamics. The
corrected statistic is zero for native repair and the one-term controls.
The certified edge bound `B≥mI`, together with `||∂e_p||²=5`, gives
`b_local≥m/100`. Its certified lower bounds are `186661676/10¹⁴` at θ=1/16
and `3778430107/10¹⁴` at θ=1/8. Numerical evaluations are respectively
`4.14290705794e−6` and `6.05048265887e−5`.

Only N_p needs measurement at the two endpoints of each independently prepared
trial. If E_a is its outcome projector and P_x the full load-pattern
projectors, then `Σ_{x:x_p=a} P_x τ P_x = E_a τ E_a = τ E_a`. Thus local
endpoint measurements give the same marginal joint law as the all-load
protocol. No other port's record, intermediate measurement, private reset or
postselection is needed. Preparation and measurement access remain inputs.

The stopping rule is still one or two **global IID attempts**. Stopping at
observer-local pulses, delayed messages or state-dependent events defines
another instrument and requires its own null calculation. With a supplied
constant Poisson rate ν>0, `c(t)=τ[f exp(−νAt)f]/τ[f²]` instead gives

```text
c′(0)=−νr,       c″(0)=ν²((11/10)r²+b_local),
c″(0)/c′(0)² − 11/10 = b_local/r².
```

The ratio is invariant under a common constant change of time units, not
arbitrary retiming. It does not calibrate physical duration. C and E still
give the same first two moments, so this local statistic cannot select their
completion or uniquely identify N4 among larger-state laws.

A measured test must propagate preparation, readout, scheduling and sampling
errors through the complete corrected statistic. Its quadratic terms make a
finite-sample plug-in estimate biased; a positive point estimate alone is not
evidence of leakage. Trial allocation, shared uncertainty, confidence rule and
adequate power remain prospective obligations. The regressions evaluate
unprojected incident-edge moments independently of the leakage formula and
retain the exact native false-positive control; they are finite-model checks.

## A central-port carrier that retains the quantum repair

The whole working algebra `B(H)` has scalar centre and cannot be the
twelve-central-port carrier in A1. Here H is the native load space, or for
N4 the binary load space tensored with its twelve six-level factors. Its
dimension is d. An address carrier is

```text
A = direct_sum_(p=0..11) B(H),       e_p = identity in block p,
j(X) = (X,...,X),                  decode(sigma) = sum_p sigma_p.
```

The e_p are precisely its primitive central projections; the native
occupation observables are the distinct, generally noncentral j(N_p).
The same symbol p indexes an instruction address and a native port. It
does not identify those two operators. The block readback pi_p(X)=X_p is
unital. The diagonal inclusion j is a unital faithful *-homomorphism, and
decode is its trace adjoint. The full quantum sector, including its private
coherences, is accessible within every block.

This extension is minimal among finite algebras with exactly twelve
primitive central projections and a unital inclusion of the full B(H).
Write such an algebra as a sum of M_(d_p). Each component of the inclusion
is a unital representation of M_d, so d_p=m_p d for a positive integer
m_p. Its smallest faithful Hilbert representation therefore has dimension
at least 12d. Equality gives the displayed diagonal inclusion up to
independent block unitaries. This is an algebraic minimum for that inclusion
contract; it is not a physical resource lower bound for every encoding.

### Ports select real repair operations

For each of the thirty unoriented seams e let Phi_e be the existing native
or N4 channel. Supply the existing uniform independent attempt schedule by
two operations on the same carrier:

```text
R(sigma)_p = decode(sigma)/12,
T(sigma)_p = (1/5) sum_(q adjacent to p) Phi_{pq}(sigma_p),
Tick = T R,                      Phi = (1/30) sum_e Phi_e.
```

R forgets and refreshes the classical address, without measuring or resetting
the native state. T uses that address to choose one of its five incident
seams. Both maps are CPTP; unitality of the seam channels also makes them
unital. Each seam has two endpoints, hence for every cq input, including
address/native correlations,

```text
decode Tick(sigma) = Phi(decode sigma),
decode Tick^n(sigma) = Phi^n(decode sigma).
```

The identity is at channel level and survives adjoining an arbitrary
reference system. It preserves complete native histories whenever the
inserted instruments are lifted blockwise. There is no private reset.
Replacing R by a fixed requested address gives a different, incident-star
experiment. On native occupation preparations the resulting port responses
distinguish addresses, so deleting the central address removes an actual
intervention. Omitting refresh between attempts generally changes the
channel, since the output address is correlated with the quantum state.
The regressions retain that failure control.

For a requested binary query p, a finite instrument measures N_p and
writes occupied into central port p and empty into its antipode. It retains
the corresponding quantum Lüders branch. The requested p is part of the
experiment's control contract, so the output port determines the bit.
Proper relabeling transports the request and antipode together. Reading
the output centre is nondemolition. Encoding both an unrestricted request
and its bit without that context requires at least 24 distinct classical
labels; the twelve-port carrier does not silently provide that memory.
Retaining outcomes across independent trials requires the stated external
record service, just as the endpoint protocol does.

### A full abstract source interpretation

Take specifically `H = C^66 tensor (C^6)^tensor12`, for binary N4 and its
theta=0 native control. The displayed address carrier replaces each carrier
in the [scalar-seam source](../source_selection_model/DERIVATION.md), keeping
its local-chart interpretation and midpoint refinement. The private factors
are retained at theta=0 and used by the response probes even when repair
does not act on them. This constructs an admissible abstract source with
supplied operations; it does not select a physical repair law. The general
address adapter works for other H, but this full-source proof uses these
factors. In particular, the M78 native repair is not unital and does not
inherit this binary tracial argument.

For N4 define the collective response on the same H by

```text
Delta D(v) = I_load tensor sum_(a=0..11) D(v)^(a).
```

Different factors commute, giving
`[Delta D(v),Delta D(w)]=Delta[D(v),D(w)]`. If Delta D(v)=0, differences
between tensor-product eigenvalues force D(v) to be scalar, and its sum
then forces that scalar to vanish. Thus Delta is faithful, and its image
is the same twelve-dimensional compact algebra. A positive trace form
follows. Full local effects detect the relative block phase on any one
private factor, so no response direction becomes invisible. The source
six-level generators, their ordered brackets and proper implementers are
the existing recorded operations, lifted collectively. The identity
`exp(t Delta X)=I_load tensor exp(tX)^tensor12` lifts the executed paths.
Lift each executed Cayley factor U to `U^tensor12`; applying the Cayley
transform to Delta X instead gives a different unitary.

On a proper carrier recharting a, permute the address, load labels and
private-factor locations, and conjugate every private factor by U_a.
This chart automorphism makes the readout and repair diagrams commute.
Write the actual internal chart implementer as
`Uhat_a=P_a^(load) P_a^(private slots) U_a^tensor12`. Its permutation
factor commutes with every Delta D(v), while `U_a^tensor12` belongs to the
connected collective response group. Thus Uhat has precisely the A2 form
`g_a c_a`, and the same-history response functor composes. The nontrivial
action on response labels comes from g_a. The central port permutation
is part of the algebra chart; an inner automorphism of A cannot move its
centre. The discrete location permutations do not add continuous probes.
The original slow spatial triplet and lower response triplet retain their
different marked order-five characters; no identification of those triplets
is used. This preserves the exact negative controls in #1087 and #1088.

The admitted reversible port controls are the collective response words,
their inverses, and the specified chart changes, with their control sequence
fixed independently of measured records. Nondemolition central readout is
available, but its outcome cannot select, omit or repeat a reversible word
or chart path. Feedback uses the stipulated addressed repair and reset
templates. This is an explicit grammar restriction, not a consequence of
irreversibility: unrestricted record-conditioned reversible control would
give the 144 independent directions `e_p Delta D_q` and violate the declared
twelve-dimensional response. A1 requires feedback availability, not this
larger closure under arbitrary record-controlled programs.

Native/N4 Tick, requested addressing, quantum readout and trace reset are
irreversible channels. In particular
the coherent W in a Kraus factorization is not separately exposed as an
arbitrary reversible probe. Each declared irreversible primitive loses
information on the full cq state space: refresh and requested addressing
forget the original centre, readout overwrites it, and requested-block
reset discards that block's internal state. Subsequent channels cannot
invert that loss for all states. This argument concerns compositions of
deterministic CPTP endomorphisms of the full cq algebra. It does not apply
separately to instrument branches: reading p and then restoring address p
is the identity on cq inputs, despite a fixed requested-address map having
a kernel. Admitted feedback is specified by its resulting full-input map;
normalized postselection on a restricted support is a different object.

The addressed T and individual seam maps also supply no reversible
direction, even if an average has no linear kernel. On one fixed finite
algebra, a product of CPTP endomorphisms can have a CPTP inverse only if
each factor does: a left inverse implies bijectivity, after which the
CP inverses follow successively. Such reversible maps are automorphisms.
Each native seam expectation, or its fixed N4 conjugate, kills a nonzero
operator. An incident-star average cannot be an automorphism: extremality
would force every seam summand to equal that automorphism. Thus these
fixed-control compositions also cannot create a reversible response.
The separate control-type restriction above rules out adaptive reversible
programs. Allowing independent internal
unitaries, controlled reversible probes or W as a primitive changes the
grammar and requires a fresh response-completeness proof.

A common bounded mismatch is
`M(sigma)=sum_p ||sigma_p−Tr(sigma_p) I_d/d||_1`, between zero and two.
Unital repair contracts it; address refresh contracts it by convexity.
Requested-block trace reset eliminates that block's term, so twelve such
resets repair it. It measures internal preparation mismatch; nonuniform
address weights are not silently counted as the same mismatch.

The remaining source clauses carry explicit maps:

| Clause | Construction |
| --- | --- |
| Accessible algebras, ports and observer operations | The displayed A, its twelve block interfaces, address-dependent repairs, full matrix-effect readback, requested-block reset and in-place continuation. |
| Seams and local charts | Distinct-carrier seams are scalars with normalized trace restriction `sum Tr(X_p)/(12d)`; triple restrictions are identities. Local charts use the proper automorphisms above. |
| Support and refinement | Use the existing oriented midpoint sphere, open-star nerve and identity bridge. Coarsen local states by positive descendant mixtures; lift each instrument to the tagged descendants and retain the same tag through a word. Linearity gives every unnormalized instrument square. Conditional tag weights update by Bayes' rule. |
| Complete state and constraint grammar | All `12d²−1` independent Hermitian state coordinates are observable through matrix effects. Every admitted history is a composition of their linear instrument maps and conditioning. Chart covariance and scalar seam normalization are the complete declared ontic-state constraints. |
| A3 | Specify the faithful reference with blocks `I_d/(12d)`, the full carrier cover and positive weights `1/N`. Relative entropy is nonnegative with equality exactly at that local reference, giving the unique compatible tracial family. Mixtures under refinement preserve it. |

The midpoint tower's degree-one and mesh proofs are independent of d.
If one independent coordinate is omitted, a small positive and negative
traceless Hermitian perturbation of the faithful reference agrees on every
retained coordinate and differs on an admitted effect. This supplies the
omitted-coordinate mutation for the complete cover, including coherences.
The tagged refinement proof applies to arbitrary finite-dimensional CP
instruments, including the address tick; it does not broadcast a selected
outcome to sibling carriers. Preparations, calibrated effects and the
schedule are declared operations of this source. The construction proves
their joint abstract compatibility. It does not derive their laboratory
implementation, physical duration, remote transport, coefficient selection,
or a physical continuum. Its scalar seams preserve the original remote-read
non-entailment. Occupation geometry and endogenous internal response coexist
on distinct typed observables of one carrier.

### A fixed-full-reference counterpair

Choose theta=0 and theta=1/8, with the same binary completion, C or E.
Both constructions satisfy the source clauses above. They have identical
accessible algebras, collective response, chart functor, scalar seams,
support/refinement, mismatch function, ontic-state constraints, full quantum
reference tau, cover and exact weights. Their preparation/effect menu and
fresh-address schedule are the same. Only the fixed irreversible repair
channels differ; their history/transition maps are therefore not identical.

Both repairs preserve tau. A3's ontic optimizer is exactly the same full
quantum state, not merely the same classical configuration marginal. It
does not choose between the two fixed update maps. Calibrated preparations
and observed outcomes supply explicit intervention constraints; applying
an unconstrained A3 reset after every intervention would be a different
process. Theta is fixed within each source, not an extra tunable reversible
control exposed to its port grammar.

For both sources, `60K/Tr K=L_ico` and `r=alpha/11>0`. The native source has
`b_local=0`; the existing theta=1/8 certificate proves
`b_local >= 3778430107/10^14 > 0`. Consequently the intrinsic channel invariant

```text
c''(0)/c'(0)^2 − 11/10 = b_local/r²
```

vanishes in one source and is positive in the other. A positive constant
Poisson attempt rate cancels from the ratio. Matching the initial decay by
a clock rescaling cannot match the curvature. This establishes finite
kinetic nonselection even with the full quantum reference and complete
response fixed, strengthening the configuration-only distinction in #1090.
It leaves that issue's canonical inertia theorem and its hypotheses intact.

The invariant is defined by the supplied channel and effects. Executing its
two-time estimator requires the expressly separate retained-record service;
address refresh does not preserve the first endpoint tag. This pair is an
abstract A1–A3 witness under the scalar-seam interpretation, not two physical
realizations, an all-sector independence theorem or an empirical exclusion
of OPH. Extra source constraints could select or exclude either law.

### Occupied-port records have a different information cost

An alternative request is to report a uniformly sampled occupied port.
On the binary total-two shell its POVM effects are F_p=N_p/2. They sum to
identity, but they are not orthogonal central projections. A concrete
instrument turns this into a central output by

```text
J_p(rho) = N_p rho N_p / 2,
J(rho) = direct_sum_p J_p(rho).
```

The compressed output block has dimension 11 times the private dimension,
and there are twelve such blocks. This output is a measurement result,
not a unital embedding of the full native algebra. On the address carrier,
use the same operation without compressing its blocks:
`J^A(sigma)_p=N_p decode(sigma) N_p/2`. Its output belongs to A, with each
block supported on N_p H. This declared readout forgets the old address,
so it adds no reversible direction to the source grammar. The inclusion of
its supported blocks lets Tick continue on the same carrier. Forgetting the outcome
gives a Schur channel Lambda. For two occupied pairs S,T and any private
operator X,

```text
Lambda(|S><T| tensor X) = (|S intersect T|/2) |S><T| tensor X.
```

Its fixed algebra is the sum of the 66 private matrix blocks. Repetition
multiplies one-common-port coherences by 2^(-k) and erases disjoint-pair
coherences at the first read. This distinguishes occupied-port sampling
from the address refresh R, which leaves the decoded quantum state intact.

The loss is unavoidable for every instrument with these effects, regardless
of ancillary output size. If S and T are disjoint, every Kraus operator in
outcome p annihilates at least one of |S> and |T>, since
`sum K_(p,k)^* K_(p,k)=N_p/2`. Therefore |S><T| is erased in every recorded
branch. The orthogonal states `( |S>+|T> )/sqrt(2)` and
`( |S>−|T> )/sqrt(2)` have identical outputs. Any recovery has worst-case
trace-distance error at least 1/2 on this pair, by the triangle inequality.
The forgotten Lüders output attains 1/2 on this pair. No optimality over
all input states is asserted.

More generally, suppose a channel encodes every native quantum state in a
cq carrier and has a CPTP left inverse. Read the output centre and then
decode. This is an instrument summing to the identity channel. The identity
channel has rank-one Choi matrix, so each branch is a nonnegative scalar
multiple of it. Each central probability is consequently input-independent.
Exact recovery of all quantum information cannot coexist with an informative
central occupied-port label. The addressed carrier meets this boundary by
using ports as supplied controls rather than passive copies of occupation.

### Endpoint statistics survive an explicit central measurement

For fresh occupied-port sampling at the two endpoints, tracial preparation
and n unmodified tracial-stationary attempts give the recorded joint law

```text
q_pq(n) = tau[N_p Phi^n(N_q)]/4,
c_n = (144 q_pp(n)−1)/5,
b_local = c_2−c_1²−(1−c_1)²/10.
```

This follows from `J_p(tau)=tau N_p/2`, followed by the final effect N_q/2.
Preservation of tau is needed for the centered-correlation conversion.
Proper carrier covariance makes every diagonal q_pp equal. Pooling all
equal-tag events gives `P_equal(n)=12 q_pp(n)` and
`c_n=(12 P_equal(n)−1)/5`. This uses all ports without a rare fixed-pair
filter. Without that covariance it gives the mean port correlation;
substituting the mean into the nonlinear local statistic needs a different
null. At n=0, `P_equal=1/2`, because each endpoint draws a fresh occupied port.

Thus the existing strict N4 lower bound on b_local is a statement about
central-record probabilities in a specified instrument. The address carrier
realizes the same Phi between these endpoint measurements. Even at zero
attempts, fresh sampling gives q_pp(0)=1/24, whereas one tag has probability
1/12. Treating the two samples as one persistent tag gives a false
normalization. An intermediate occupied-port read inserts Lambda and
changes the experiment; the original certificate does not describe that
modified channel. Uniform address refresh inserts no such quantum read.

## Reproduction and scope

Run `python -m pytest -q -W error code/native_response` from the repository
root with the pinned requirements. CI runs this on Linux in the
existing finite quantum-information job. The tests bind the source chart,
compare independent Kraus and source formulas, recover the existing Gram,
check full-shell overlap counts, retain the simpler controls, and replay
the finite-coupling interval certificate.

The attachment tests separately reconstruct the native 66-state Kraus
families, transport all sixty proper actions and test correlated cq inputs.
Occupied-port tests retain private coherences, lost disjoint superpositions,
a nonstationary false-positive control and actual N4 odd-sector instruments.
Small collective-response tests check exact source rank and the marked
triplet obstruction, factorwise paths and the full chart centralizer.
They do not allocate `66*6^12` matrices: the arbitrary-factor admission and
refinement arguments are analytic. Replay the inherited finite source data
with `PYTHONPATH=code python -m source_selection_model.verify`; that verifier
also distinguishes finite checks from full model membership.

The floating evaluators are diagnostics on `|θ|≤1`; they do not certify
zeros, tiny differences or arbitrary-angle forward accuracy. The separate
certificate starts from the algebraic source entries and exact rational
couplings. All-angle nondegeneracy follows from the proof above.
Address/measurement division rejects erased components. Channel and query
matrix products conservatively reject subnormal operands, whose products
can lose representable information inside BLAS; direct address moves retain
their supported subnormal values. These guards delimit numerical evaluation,
not the exact channel theorem or a general forward-error certificate.
The finite common fixed algebra, absolute quantum gap, physical clock,
source-law selection and continuum limit are not inferred. In particular
C and E share these first two occupation moments, so this test alone does
not distinguish those two completions. No M01–M12 parent or #1076 is closed
by the result.

### Finite certificate and error budget

`PYTHONPATH=code python -m native_response.certificate` prints the two
recomputed rational enclosures (set `PYTHONPATH` using the local shell).
There is no user-supplied success packet to trust. The certificate separately
reconstructs the algebraic source using directed interval operations,
including √5, in an isolated 35-digit context. It evaluates degree-24
Taylor polynomials on each protected 18-dimensional load/private block.
The interval row sums prove `||θT||∞<2`; skew-Hermiticity also gives
`||θT||₂≤||θT||∞`. The entire exponential tail is bounded by
`e²2²⁵/25! < 8·2²⁵/25! < 2⁻⁵⁰`. Each real and imaginary matrix entry is
enlarged by that error before the nonlinear trace computations. Arithmetic,
source and truncation uncertainty therefore all reach the final enclosure.
The existing broad rational endpoint export checks outward containment.

The 120 unordered pairs of incident carrier edges form two orbits of sixty
under the proper carrier action. Representatives are `{(0,2),(0,4)}` and
`{(0,2),(0,8)}`. Symmetry of the Gram and endpoint orientation signs leave two
absolute off-diagonal entries, with four of each in every row; all other
off-diagonal entries vanish. The exact permutation enumeration is tested
alongside source covariance and the independently computed thirty-edge Gram.
If m denotes the diagonal minus those eight absolute entries, the certificate
gives the strict bounds

| Exact coupling | Certified lower bound on m |
| --- | --- |
| 1/16 | 186661676 / 10¹² |
| 1/8 | 3778430107 / 10¹² |

Thus `B≥mI` and the normalized leakage matrix is at least
`(11m/6000)L_ico`. This proves positivity for the whole centered subspace,
including arbitrary superpositions, without diagonalizing an enormous
network operator. The enclosure does not assert those rational lower bounds
are optimal, nor that nonclosure holds at every nonzero coupling.
