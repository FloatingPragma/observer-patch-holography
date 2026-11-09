# Local calibration of the cubic Hubbard action

The source has actual signed cubic hopping
h_(x,x+e_j)=(c/2)(-1)^(sum_(k>j)x_k), two flavours and onsite contact U.
The physical units of c,U and time are not selected here. The construction
identifies their ratio from the original probability of the unchanged local
repair collision record, with a volume-independent error. It consumes a
known finite calibration preparation, rather than the interacting ground.

## Exact bulk coefficient from a local source neighbourhood

Take adjacent interior vertices x,y, with every nearest neighbour of both
present. Start with one up fermion at x and one down fermion at y. In this
invariant two-particle sector, the centered clock differs from
H=T+U D only by the scalar -U. The positive defect convention
T+(U/2)sum_z(N_z-1)^2 differs by U(V-2)/2. Both scalar changes cancel from
the experiment. No volume-dependent scalar is included in a norm bound.

Let Q be the projector on a doublon at x or y. It is exactly the event
raw8 for either full-domain repair law on edge(x,y). All law, coin and raw
records remain in the complete instrument; only the statistic is grouped.
On the full source its effect is D_x+D_y-D_xD_y. On the two-particle input
it reduces to D_x+D_y. Put p(t)=<psi,e^(itH)Qe^(-itH)psi>.
The coherent clock acts between the supplied preparation and this single
readout. Intervening repairs are not silently removed from a history.

The source is real in this chart, so p is even. Bipartite parity gives
QT²psi=0. Since a hop from a doublon leaves the doublon subspace,

    QHpsi=QTpsi,
    QH²psi=U QTpsi,
    QH³psi=QT³psi+U²QTpsi.

Consequently, with A=||QTpsi||² and B=<QTpsi,QT³psi>,

    p(t)=A t²-(U²A/12+B/3)t⁴+R6(t).                     (1)

These coefficients are determined by a one-step graph collar of(x,y).
The one-up/one-down CAR sector is unitarily represented by two coordinate
slots; in that representation T=h tensor I+I tensor h. The unitary includes
the exact sign for converting spin-first creation order to the source's
site-first occupation order. Hence A=c²/2. Write a=h/c. Then

    B/c⁴ = sum_(z=x,y) a_xy [ (a³)_xy+3a_zz² a_xy ],      (2)

where a_zz² denotes (a²)_zz and the two terms use their corresponding
oriented endpoints. Both endpoints have degree six, so (a²)_zz=3/2.
There are exactly fifteen length-three walks from x to y. Three use only
the edge direction. For each of the four transverse directions there
are two retracing walks and one plaquette detour. The actual pi-flux signs
give positive relative sign to the eleven retracing walks and negative
sign to the four detours. Thus (a³)_xy=7 sign(x,y)/8 and

    B/c⁴=25/8,
    p(t)=c²t²/2-c²(U²+25c²)t⁴/24+R6(t).                 (3)

Every such walk is in the union of the endpoint neighbourhoods. Thus the
calculation applies to every open cubic source containing that collar,
independently of its exterior or volume. The code exhausts all fifteen
walks in each of the eight coordinate parities and three directions. The
source sign depends only on those parities. This is a finite local
exhaustion proof, not an extrapolation from a sequence of boxes. Removing
the pi-flux signs gives33 in place of25 and fails the source check.

## Uniform finite-time inference

For degree at most six, ||h||<=3c by the row-sum bound. Second quantization
on the two-particle sector gives ||T||<=6c and ||UD||<=U. If specified
bounds give c<=c_max,U<=U_max, put M=6c_max+U_max. The sixth derivative of
p has modulus at most ||ad_H^6(Q)||<=64M^6. Taylor's theorem therefore gives

    |R6(t)|<=4M^6|t|^6/45.                               (4)

Define C_hat=(8p(t)-p(2t)/2)/(3t²) and
Q_hat=(8p(t)-2p(2t))/t⁴. Then

    |C_hat-c²|<=32M^6t⁴/27,
    |Q_hat-c²(U²+25c²)|<=544M^6t²/45.                    (5)

Probability error epsilon at each time adds17epsilon/(6t²) and
10epsilon/t⁴ respectively. Positive rational enclosures C and Q yield

    (U/c)² in [Q_-/C_+²-25,Q_+/C_-²-25] intersect[0,infinity).

With t->0 and epsilon=o(t⁴), these intervals contract uniformly in volume.
U>0 is the source-domain assumption. Common rescaling(c,U,t)->(sc,sU,t/s)
leaves the entire experiment unchanged, so absolute units require an
independent clock calibration. No zero-time derivative is measured.

The finite witness evolves the actual invariant one-up/one-down sector by
a sparse rational Taylor polynomial of order n. Its vector error is at
most((6c+U)|t|)^(n+1)/(n+1)!, since the exact evolution is unitary. The
collision probability error is at most2r+r². All coefficients are actual
source coefficients. A separate occupation-mask CAR calculation checks
the spin-first representation and its signs.

## A finite empty collar suffices for arbitrary distant matter

Global vacuum preparation is unnecessary for this local statistic. Let
Lambda_R be a source cube of radius R containing(x,y) in its radius-one
interior. Require its initial marginal to be within half trace distance
delta of the known two-particle product with all other collar sites empty.
The exterior and any exterior reference are otherwise arbitrary. The
exact pure marginal case factorizes from that exterior automatically;
the approximate case here only asserts a bound on local probabilities.

Let H_R omit bonds crossing the collar. The actual CAR locality bound in
`local_matter.boundary_upper`, for the even norm-one collision effect with
support size two, gives

    |p_full(t)-p_R(t)|<=delta+B_R(t).                     (6)

This follows by first replacing the evolved effect with its cut-clock
effect in operator norm, then comparing the collar marginals. No trace of
the distant physical archive is equated with a full-record channel bound.
The remaining fine records and their quantum outputs are retained by the
actual complete repair instrument. Inequality(6) concerns the collision
statistic alone.

The rational planner uses the inherited bound

    B_R(t)<=min(2,12c(2R+1)²|t| 2^(ceil(96c|t|)-(R-1))).  (7)

For every positive error target, increasing R terminates. Choose
delta+B_R(2t)=o(t⁴), along with sampling and coherent-service errors of
the same order, and apply(5). This is a local preparation accuracy
requirement with a computable finite collar. A supplied local reset can
consume2|Lambda_R| blank mode registers and export the overwritten
registers; no such reset tape or global-vacuum producer is executed here.
The preparation premise cannot be dropped: a completely filled collar
has collision probability one at t=0 and does not obey(3).

The collision effect depends only on local doublon occupation and is
invariant under the actual onsite Weyl/signed frame. The known collar
product maps to another known occupation product with local phases.
Thus(6) also applies to the original Weyl action and its bare collision
effect. It does not claim that all fine repair maps are frame covariant.

Independent identically prepared trials can estimate the complete
collision indicator with the two-time Hoeffding budget of
`source_action_identification.sampling_budget`. Failures remain in the
denominator and complete implementation errors are added to epsilon.
If the selected scheduler chooses uniformly from E source edges, use the
indicator of the joint event "seam=(x,y), raw8". Its unconditional mean
is p(t)/E. Multiplying its sample mean by the known E estimates p(t)
without conditioning on the seam. Applying the same Hoeffding budget at
precision epsilon/E gives a sufficient number of global attempts; its
quadratic E cost is retained. All other seam and raw outcomes stay in the
original sample denominator. Fairness without the uniform seam law does
not supply this normalization. A known nonuniform seam weight requires
its corresponding factor instead. This defines the clock as a global
attempt count and does not identify that count with coherent action time.
The construction tests the configured c,U on growing sources. It supplies
no capacity-to-action law, thermodynamic ground preparation, physical time
unit, particle species or mass assignment.
