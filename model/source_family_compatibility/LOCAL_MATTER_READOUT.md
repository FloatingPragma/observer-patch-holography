# Local interacting matter with complete observer records

The full Weyl/contact action supports odd CAR operators, internal spin
and a nonzero interacting local response on one centered ground and clock.
A uniform global repair schedule approximates its local observer histories
without sending U/c to infinity. The comparison is for a specified local
record algebra. It does not equate the global repair archive or construct
a relativistic continuum.

## Source and accessible records

Use the two-flavour full-Fock source and its declared collision-hold repair
extension, exactly as in `DEFECT_BUDGET.md`. The centered Hamiltonian is

    K=T+U D_double-(U/2)N,      c>0, U>0,

with the original nearest-neighbor Weyl links `i c sigma_j/2` and bare
`full_domain` repair in that original basis. Each full-Fock bond is even,
supported on its seam and has norm c. Each bare repair Kraus operator is
even and supported on its seam. These are all the action/repair hypotheses
used by the local-history comparison. In particular, no alignment sector
or covariance of the repair under hidden spin is required.

A second direct instance uses the signed Hubbard links and `full_domain`
repair in that chart. The onsite frame in `model/field_symmetry` transports
the two actions, including the pi-flux signs. Unitary equivalence of the
whole protocol additionally conjugates every repair Kraus map and CAR probe.
The bare repairs in the two charts are not identified. Both direct instances
satisfy the local theorem independently, because onsite conjugation preserves
the Hamiltonian bond norms and both bare instruments have even seam support.
Neither the signs nor U are selected by this comparison. Centering fixes
one clock for the state and all probes. The aligned-hole construction has
its separate same-chart invariance requirement.

The experiment uses at most k local Stinespring interventions, with every
outcome retained, on s sites in a cube of radius r. Ordinary source time
runs from zero to T. A final accessible field output on these sites is
allowed, as are arbitrary observer memory and untouched references.
Record-adapted controls have the same policy in both comparisons and keep
this spatial support. No global controlled Hamiltonian measurement is
treated as a local intervention.

There are m repair attempts at declared times. Each seam is uniform on
the E source edges, conditional on the available history. The law and coin
records can be retained, as can the raw records of seams meeting the
experiment's sites. Distant raw values are not accessible to this protocol:
they are traced out, or replaced by one declared outside-record symbol.
They are retained in the physical source archive and no outcome is
postselected. Controls and later schedule choices cannot depend on those
discarded raw values. All source initial states and references are admitted
in the finite comparison.

The ideal experiment replaces each repair by identity and writes the same
seam/law/coin labels, raw zero on accessible seams, and the same outside
symbol. Physical service errors for either implementation are additional.

## A local complete-history bound

The full-CAR locality theorem in `model/full_fock_field/README.md` applies
to this centered action with identical constants. Removing the onsite
terms by an interaction picture leaves even two-site terms of norm c.
For a norm-one observable supported on the s sites, and R>r,

    B_R(T)=min(2,6c s(2R+1)^2 T 2^(-(R-r))(exp(48cT)-1))     (1)

bounds the difference between full and radius-R dynamics. It holds for
odd as well as even local CAR fields, and uniformly after adjoining
arbitrary memory/reference systems. The source proof instantiates the
fermionic commutator-chain theorem of
[Nachtergaele, Sims and Young, Section 3](https://arxiv.org/html/1705.08553v3#S3).

For any ideal suffix with at most k interventions, its dual output
observable has a radius-R approximation with norm error at most

    epsilon_R=min(2,(2k+1)B_R(T)).                           (2)

To prove this, defer the local instruments to their Stinespring unitaries
and keep their memory. In the interaction picture each intervention is
one evolved local unitary, whose full-versus-cut difference is at most
B_R(T). Telescoping their ordered product costs k B_R(T). The dual
observable has that product on both sides, contributing twice this
amount. The final local field observable contributes B_R(T). Each
operator norm is at most one. Taking the supremum over output observables
and reference dimensions therefore supplies the complete output bound;
it is not an inference from finitely many selected scalar expectations.

Every repair Kraus operator is even and belongs to its seam's CAR algebra.
The exact signs in the full-domain source are essential here. A seam
disjoint from the radius-R cube therefore fixes the cut dual observable,
including its local memory and reference, after its raw record is traced.
For a far seam, repair minus identity has dual norm at most 2epsilon_R.
For a seam meeting the cube the elementary bound is two. There are at
most 6(2R+1)^3 meeting seams. Put

    q_R=min(1,6(2R+1)^3/E).

Telescope the actual prefixes and ideal suffixes at each repair. Averaging
the conditional uniform seam at that step gives, in half-trace distance,

    D(output_actual,output_identity)
          <= min(1,m(q_R+epsilon_R)).                     (3)

All preceding actual repair dynamics is inside a contractive prefix;
no defect expectation or invariant occupation subspace is required.
Every local fine outcome is retained. For the special case of only a
final local readout, k=0. Intervening local readouts and feedback use (2).
The error is independent of initial entanglement and of the onsite U.

For fixed c,U,T,r,s,k,m and target eta>0, first choose R so that
m epsilon_R<=eta/2. Then choose the outer source box with
E>=12m(2R+1)^3/eta. Equation (3) is at most eta. Both choices have an
explicit rational finite search: bound exp(48cT) by 2^ceil(96cT), and
use the exact open-cube edge count E=3(L-1)L^2. The source box must also
contain the radius-R cube. Counts can be evaluated without allocating
its Fock space. This proves a local observer comparison; it does not
claim execution of the resulting large native tape.

For growing time, support or repair budgets, the same inequalities are
the sufficient conditions. A fixed positive rate per seam is a different
schedule and does not satisfy the volume-dilution argument by assertion.
No physical duration per repair attempt or capacity trajectory is derived.

## One centered ground and surviving interaction

Use the centered thermodynamic state construction in
`model/centered_source_ground/README.md`: normalized full finite-box
ground projectors, tracial exterior, translation averaging and one weak-star
subsequence. Its local ground inequality, stationary positive GNS clock,
gauge symmetry, hidden SU(2), particle-hole symmetry and one-site covariance
I/2 hold on the same state. There is no protected-spin compression.
The odd local addition/removal preparations each have original weight 1/2,
and the internal spin observables are even bilinears in those same fields.
Their finite-energy response and physical interpretation require their
own estimates; a local odd field is not a selected Standard Model particle.

All fixed complete local experiments converge along that state subsequence
by the norm-local dynamics construction and finite-dimensional local
instrument approximation. Assigning (3) a vanishing error gives the same
limits for the specified repaired experiments. This is an iterated or
chosen diagonal construction, without an effective modulus for selecting
the ground-state subsequence. It retains the prepared-source and coherent
control premises.

The same limiting field is quantitatively interacting at fixed lattice
spacing. For a local spin component a=c_xsigma, put

    A(t)={tau_t(a),a*},
    F(t)=omega(A(t)*A(t))-|omega(A(t))|^2.

The centered source has no onsite kinetic term and
`[K,a]=-(hopping annihilator)-U(n_opposite-1/2)a`. Since
omega(n_opposite)=1/2 and n_opposite is a projection,

    F(t)/t^2 -> U^2/4.

The source bounds `||[K,a]||<=3c+U/2` and
`||[K,n_opposite]||<=6c` give
`||[K,[K,a]]||<=9c^2+9cU+U^2/4<=M`, where
M=9c^2+12cU+U^2. The same Taylor/centering argument as
`model/retarded_spatial_locality/README.md` yields

    sqrt(F(t)) >= U|t|/2-2Mt^2,
    F(t) >= U^2 t^2/16,   0<|t|<=U/(8M).                  (4)

F is the sum of four connected CAR four-point functions, with the
ordinary fermionic Wick signs. At least one consequently has magnitude
at least U^2 t^2/64. Gauge symmetry removes odd one-point functions;
the disconnected contractions cancel algebraically, without assuming
a Gaussian state. These bounded local multi-time records are among the
experiments controlled by (3). Their implementation uses the declared
coherent local probe services and their full failure/error budgets.

Equation (4) excludes a Gaussian full time-history law on these same
fields. It establishes neither interacting scattering nor a sharp pole.
The lattice spacing is fixed: bounded-U ballistic refinement has its
separate free-limit theorem. A common relativistic chiral gauge action,
physical spin/statistics, a source-selected scaling and laboratory particle
identification are not consequences of this local construction.

## Exact locality and archive controls

The producer verifies complete dual identity on every local CAR basis
element for a disjoint repair seam. Its crossing-seam control retains the
actual CAR sign across an intermediate site; omitting that sign destroys
the odd-field identity. Every padded coin/raw slot remains.

The original-source control reconstructs the actual centered Weyl action
on three-site chains along all three spatial axes. It checks the bond
spectrum and the local charged second-moment operator identity on every
Fock input. A bare repair on the disjoint physical edge fixes all sixteen
operators in the observed site's CAR algebra. No spin-frame conjugation
of that instrument occurs in this control. The corresponding centered
non-Gaussian derivative is checked as a full operator identity as well.

A second exact control entangles the local up/down alternatives with a
remote ordinary/collision configuration. Discarding the remote raw record
leaves the complete local density unchanged, while the retained global
archive has original collision probability 1/2. Thus the local theorem
cannot be promoted to equality of the complete global archive. A source
observer that reads those values needs the full-archive defect budget or
another justified comparison.
