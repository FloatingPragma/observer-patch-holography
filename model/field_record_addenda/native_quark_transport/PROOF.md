# All three quark-current families move under the native Hamiltonian

This calculation uses the unchanged eight-cell thirty-CAR/five-scalar/PW119
source `H=H0+lambda T`, `lambda=2^-34`, all onsite coefficients one, and
the complete physical preparation in `interacting_ground_preparation`.
The compensated pump and static thirty-mode reference `4 N_ref` are the
actual ones in `native_quark_response`. There is no additional hopping
term and no identification of an internal charge label with a measured
quark mass.

## A simpler exact one-body compression

The actual quark Fock character has exactly one Casimir-54 module: its
dimension is 140, its transformed occupation is twelve, and its
multiplicity is one. The source's quark spin SU(2) and color SU(3)
commute with its ten Spin(5) generators. They therefore act on the
one-dimensional multiplicity space; both actions are trivial. In
particular the entire module, not just one selected eigenvector, has
quark spin zero. The local physical ground doublet's spin belongs to
the lepton sector.

Let P denote that complete local ground doublet. Gauge invariance,
quark-spin invariance and the exact ground occupations `N_r=d_r`
for `r=u^c,d^c,Q` imply the operator identity

`P c_(s,a)^dagger c_(t,b) P = delta_st delta_ab P/2`.       (1)

This holds on every pure vector of the doublet, not just its mixed
average. Off-diagonal internal components transform nontrivially under
the gauge group, and quark-spin-vector components vanish on a singlet.
The scalar coefficient is fixed by summing the exact occupation. No
native pair has zero total hypercharge, so `P c c P=0` for every pair.
Distinct irreducible native gauge families cannot refill the same
neutral ground through an off-family one-body compression.

## Exact first-hop effect, independent of the spin background

Use the recorded atomic localization `{P0,I-P0}`, where P0 has the true
local doublet at all eight sites and every link blank. On success use
the full compensated pump `B_r=sum r_(s,a)^dagger c_(0,s,a)/sqrt(2d_r)`.
Its effect is exactly P0/2. All pump and localization failures remain.
The static reference stays at the origin.

For each of the six actual periodic incident routes, the arrival effect
P_e requires only that link excited and every nonarrival cell in its
true ground. It leaves all target matter, scalar and reference states.
The six effects are orthogonal; their full complement is retained.
The actual kinetic coefficient has magnitude 1/4 and its spin matrix
is the corresponding Pauli matrix. The literal PW119 blank-link Gram
matrix is `I/d_r` for all forward and adjoint columns. With (1), the
origin double compression and the arrival one-body compression give

`P0 B_r^dagger T P_e T B_r P0=P0/128`.                    (2)

Indeed the three one-body factors supply 1/8, the spin sum is 2,
the two internal sums give d_r^2, and pump, hopping and link factors
are `1/(2d_r)`, `1/16`, `1/d_r`. Their product is 1/128.
The code also checks all three complete two-spin matrix effects with
the quark-spin coefficient b=0. Consequently (2) holds for arbitrary
entangled local ground-spin inputs and an untouched environment.
There is no need to compare this first-hop effect with one particular
cube spin ground. CAR grading supplies a common amplitude phase only:
each unexcited local ground has fixed odd parity and the compressed
origin one-body operator is even. Distinct periodic flux routes stay
orthogonal, rather than being combined as classical paths.

Reference and flux charges make every arrival physical. Source changes
in hypercharge are `+2/3,-1/3,-1/6` and in raw B3 are `+1,+1,-1`
for `u^c,d^c,Q` removal; the static reference carries the opposites.
The origin flux cancels that reference, the destination flux cancels
the conjugate removal field, and all other Gauss constraints remain.
Raw lepton number remains unchanged. The reference is not relabelled
as the moving object.

## The actual prepared law and full source clock

Let q be physical preparation success and a be localization success
conditional on it. The original preparation supplies its actual q lower
bound and `a>=1-(r+epsilon_p)^2`, with r the invariant graph bound and
epsilon_p its retained-environment purification error. These are only
used to bound localization success. Equation (2) already holds for
every state after that recorded localization, so no ground-state
vector error is subtracted from the small transport amplitude.

For `psi=B_r P0 z`, `||psi||^2=a/2` and
`||P_e T psi||^2=a/128`. H0 preserves every link, hence
`P_e exp(-itH0) B_r P0=0` exactly. The inherited full-carrier bounds
`||T||<=360` and `||H0-8E0||<=576` imply by Duhamel

`||P_e exp(-itH)psi+i lambda t P_e T psi||`
` <=lambda t^2 360(1152+360lambda)||psi||/2`.

Writing `e_t=360t(1152+360lambda)/2`, the route probability
conditional on physical preparation, before localization/pump
normalization, is bounded by

`a(lambda t)^2 max(0,1/sqrt(128)-e_t/sqrt(2))^2`
` <=p_e<=a(lambda t)^2(1/sqrt(128)+e_t/sqrt(2))^2`.        (3)

The certificate uses the exact source time `t=2^-22`; its lower margin
is positive for each of the six routes and each of the three families.
Multiplication by q gives original-attempt probabilities. The ten
coarse records are `1-q`, `q(1-a)`, `qa/2`, the six `q p_e`, and
`q(a/2-sum p_e)`. They sum exactly to one. Every preparation pointer,
source partner, source/reference state and detector complement is
retained. Idle time or zero hopping gives exactly zero arrival.

The localization and pump are one finite complete controlled service
with three labels (localization failure, pump failure, success),
followed by autonomous waiting and a seven-outcome arrival service.
Their finite Stinespring unitaries extend on all unused pointer
inputs. A duration-one control `i Log U-H_total` implements each
service while explicitly cancelling the background during that slot;
the logarithm is chosen with norm at most pi. Failures run matching
duration no-op controls. This is a declared matrix-control resource,
not an efficient local circuit or a secretly paused source. The
two-copy raw source/reference norm is at most `1080+720lambda`;
thus each control norm is at most that number plus pi, and gross
switching work for these two pulses is at most four times that bound.
After the inherited physical preparation duration, all branches take
`2+t` more source time. Additional recorded capacity is `ln(3*7)`
nats beyond the inherited sources, reference and preparation pointers;
control-synthesis workspace is not included.

The literal native full Hamiltonian moves these compensated quark-type
excitations. The energy-bin experiments in `native_quark_response`
are different preparations/settings; their positive marginal weights
do not prove a positive energy-resolved arrival in this experiment.
No momentum band, dispersion, relativistic speed, asymptotic colored
state, physical quark identity or GeV mass is inferred from (3).
