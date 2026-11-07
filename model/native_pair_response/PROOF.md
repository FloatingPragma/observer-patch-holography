# Spatial pair responses in one native occupation sector

## One source, reference and counted clock

Use the fair connected binary federation with V ports and M seams. The
selector draws each seam independently with probability 1/M, repair has
probability one, and the endpoint coin is fair. Fix one conserved total
2<=q<=V-2 and the working reference tau_q=I/binomial(V,q). The source may
use its fixed C/E mixture during an unread body; recorded probes use E.
The uniform reference is a declared preparation, or a separately bounded
preparation target. Neither the density nor a laboratory time unit is
selected here. The default biased H66 source is outside this theorem.

The owning [correlation-sector proof](../native_correlation_sectors/PROOF.md)
gives P_q=I-G_q/(2M), with G_q=sum_e(I-S_e), and the primitive degree-two
inclusion I_(2,q). On this subspace its squared uniform norm factor is

    c_q=binomial(V-4,q-2)/binomial(V,q),       v_q=4c_q.

The [spatial source proof](../native_spatial_scaling/PROOF.md) binds the
combinatorial geometry, not a physical metric. All finite program examples
here retain the admitted levels zero through two. Statements for arbitrary
finite graphs or refinement levels are mathematical kernel statements with
the stated geometry and limit hypotheses, not unrestricted host execution.

## A spatial matrix of complete reflected probes

For each setting i choose two disjoint original seams (a,b),(c,d), and set

    h_i=e_ac-e_ad-e_bc+e_bd,
    phi_i=(n_a-n_b)(n_c-n_d)=I_(2,q)h_i.

Each h_i is primitive, has squared norm four and coefficient L1 norm four.
Probe settings may overlap each other; no orthogonality is assumed. Read
the two seams of setting j, perform T uniform attempts, then read the two
seams of setting i. Retain every raw flag, coin, current and FIFO export.
The score is the outgoing contrast product at the past read and the
incoming contrast product at the future read:

    s_past=(d1-2j1)(d2-2j2),        s_future=d3*d4.

The complete native two-read maps obey F_j K=s_past K and
K F_i=s_future K, where F_i=diag(phi_i). Each pair is unital and trace
preserving on the complete binary carrier. Its 625 formal slots include
589 zero maps. These are unchanged owning-source identities on every
matrix unit, including all spectator and reference extensions. Unitality
at tau_q and the unread population action therefore give

    C_ij(T)=E[s_future s_past]=c_q h_i^*P_2^T h_j,
    B_ij(T)=C_ij(T)/v_q=h_i^*P_2^T h_j/4.                 (1)

This is one fixed q, one reference and one global-attempt clock. Evaluating
the right side on two-subset coefficients is the exact intertwiner; it
does not change the working preparation from q to two. Equation (1) extends
the owning autocorrelation to a spatial matrix and distinct probe supports.
The matrix is real symmetric and positive semidefinite because P_2 is a
positive operator on the actual federation. B(0) is its probe Gram matrix,
which need not be the identity. These are positive relaxation responses,
not particle masses or Lorentzian propagators.

For finite setting weights a_i,b_j, linear combinations have response
a^*C(T)b. They can be estimated using separate experiments with the retained
setting identifiers, or randomized settings with their declared probability
weights. Scores are bounded by ||a||_1||b||_1 after that normalization.
This procedure does not simultaneously measure noncommuting settings or
silently discard their measurement disturbance. A stationary covariance
formula is not a claim that all microscopic fields have Gaussian laws.

The admitted program contains four fixed E commands and one repeated
uniform body command. Its attempt count is T+4<=2^32; its one-slot FIFO
exports every oldest record. Exact sparse response evaluation sums the
specified complete kernel without enumerating exponentially many histories.
The expanded evaluator can return RESOURCE_STOP, and its stop is retained
as a different returned result. No source admission or evaluator cap changes.

## Local windows and preparation errors

Let S contain both endpoints of all four observed seam commands, and let
B be its radius-r graph ball. Compare the actual output with the control
that keeps the same global selector and clock, uses the original repair
inside B, and acts as identity on B for other choices. The initial control
state is the exact restriction of tau_q; it retains every local number
sector. It is not a smaller sharp-q source.

The owning [local-response theorem](../native_local_response/PROOF.md)
applies with these roots. Probe commands have both endpoints in S and do
not enlarge the conservative backward support. Only the T uniform body
attempts contribute to its escape bound. Thus, with maximum degree Delta,

    p_escape <= min(1, |S| binomial(T,r+1)(Delta/M)^(r+1)),
    |C_ij-C_ij,B| <= 2 p_escape,
    |B_ij-B_ij,B| <= 2 p_escape/v_q.                       (2)

The underlying bound is on the complete selected local classical-quantum
output, including probe backaction and arbitrary untouched references.
Initial compared cone/reference states must have the same reduction, or
their joint half-trace error must be added. A bound on an unreferenced
canonical marginal does not supply that joint premise. The observation
traces remote raw values; retaining those values, using an adaptive global
archive controller, or observing an extra seam requires its dependencies
as additional roots. The observer is allowed to read the source and
scheduler specifications. Equation (2) does not erase that access.

For a complete scalar preparation-score allowance eta, the normalized
response error adds eta/v_q. Full half-trace error epsilon supplies
eta<=2epsilon, while a direct score-preparation theorem can be sharper.
For a finite smeared menu multiply these allowances by the score's L1
weight bound. In particular fixed small q has v_q of order V^-2; neither
the normalized response nor its preparation/sample cost is volume free.

## Exact exclusion versus a collision-retaining independent pair

The native two-subset generator has off-diagonal entries -1 for a legal
single endpoint hop. At {a,b} its diagonal is
deg(a)+deg(b)-2*1_(a~b). Let L be the one-port graph Laplacian. Compress
L tensor I+I tensor L to the symmetric distinct-pair space. Its diagonal
is deg(a)+deg(b), with the same legal off-diagonal entries. Consequently

    G_ind,compressed=G_2+2A,
    R=I-G_ind,compressed/(2M)=P_2-A/M,                     (3)

where A is the adjacent-pair diagonal projector.

There is a complete probabilistic comparison at precisely this clock:
select one seam uniformly and one of the two labelled walkers uniformly,
then swap that walker's endpoint on the selected seam. This is an auxiliary
independent-generator process, not a changed native rule. It updates at
most one walker per attempt, and is not P_1 tensor P_1. If the walkers
coincide, retain an absorbing collision flag. Before collision its unordered
distinct-pair kernel is exactly R. An adjacent pair has collision probability
1/M per attempt, comprising the two placements of weight 1/(2M). All other
pairs have zero immediate collision probability. R is nonnegative and
substochastic; the absorbing flag restores total probability one. It is
never normalized on surviving pairs.

P_2 and R have nonnegative entries, and the exact identity is

    P_2^T-R^T = sum_(k=0)^(T-1) P_2^(T-1-k)(A/M)R^k.      (4)

Thus the difference is entrywise nonnegative, and its column sum at pair
x is exactly the collision probability p_x(T)=1-1^*R^T e_x. This is not
an assertion that (4) is positive semidefinite as an operator.
For h_i,h_j of the probe menu, put mu_j(x)=|h_j(x)|/4. Equations (1),(4)
give the finite normalized comparison

    |h_i^*(P_2^T-R^T)h_j|/4 <= E_(mu_j) p_x(T).            (5)

This compares coefficient responses in the two-subset transfer dictionary.
It does not replace a q-occupation quantum instrument by independent
particles. The reference factor c_q in the native observable remains the
same throughout. Native raw records are not identified with the auxiliary
walker's random labels or collision flag.

The collision hazard is at most 1/M at each surviving attempt, giving
p_x(T)<=min(1,T/M). If the two initial ports have distance d, collision
requires at least d single-walker moves. Their conditional probability per
attempt is at most Delta/M. Union over d marked attempt positions gives

    p_x(T) <= min(1,T/M,binomial(T,d)(Delta/M)^d).          (6)

Use the minimum d on the support of mu_j for a uniform bound. Exact sparse
propagation retains the collision probability and can be smaller. The
bound can be large at a growing diffusive time and is not an interacting
continuum estimate. A four-port path supplies a strict finite control:
for the contrast on edges (0,1),(2,3), the one-attempt difference between
the normalized native and killed-pair autocorrelations is 1/(4M)=1/12.
Removing the exclusion correction therefore changes this observed finite
response. It cannot be justified by matching the one-body transfer.

## Smooth primitive fields on the same occupation sector

The following continuum statement concerns bounded macroscopic smearings,
not a fixed microscopic seam contrast. For real mean-zero port profiles
f,g put A_f=sum_i f_i(n_i-q/V), chi=q(V-q)/(V(V-1)), and

    Y_f=A_f/sqrt(V chi),       a=(V-2q)/(V-2).

The pair coefficient b_ij=f_i g_j+f_j g_i has lowering sum
(D_2 b)_i=-2f_i g_i. Define

    u_i=[-2f_i g_i+<f,g>/(V-1)]/(V-2),
    h_ij=b_ij-u_i-u_j.                                    (7)

Since D_2 U_1 u=(V-2)u+(sum_i u_i)1, D_2h=0 exactly. Including this
primitive vector in the same q sector yields the observable identity

    W_(f,g)=I_(2,q)h/(V chi)
      =Y_f Y_g-<f,g>/V
       -a sum_i f_i g_i(n_i-q/V)/(V chi).                 (8)

There is no selected two-particle working state in (7),(8). The subtraction
removes the exact constant and one-body components on the actual q slice.
The squared L2 error from dropping the last term is exactly

    a^2 ||P_0(fg)||_2^2/(V^2 chi)
      <= a^2 ||f||_infinity^2 ||g||_infinity^2/(V chi).     (9)

This is a finite projection estimate. At half occupation a=0; at other
occupations the one-body subtraction is essential, even at finite size.

These are accessible occupation functions, not proposed arbitrary quantum
effects. A final native spanning-tree E detector with the public conserved
q decodes the body-terminal mask. Reversing all retained original body
currents reconstructs each earlier mask, hence every profile and (8) at
each body cut. The [owning detector proof](../native_sector_readout/PROOF.md)
is valid at every global number and retains its post-detector disturbance.
The present adapter attaches it to a pure-E independent uniform body in the
same admitted program, checks original `_outcomes` word probabilities and
the actual FIFO/archive, and gives the decoder no initial mask. Initial
masks in test controls are preparations used to verify the independent
record reconstruction. Complete zero words and host resource stops remain.

This readback uses the entire body tape and global final tree. It is a
different observation map from the local four-E menu in (2); applying that
local boundary bound to the global tree would be invalid. The mathematical
family of all finite body words is also distinct from a completed expanded
host run at arbitrarily large capacity or time.

## Moment control and the inherited second-chaos limit

Use precisely the spatial subsequence and bounded mean-zero test profiles
of the [native density limit](../observer_time_fields/REFINEMENT.md#actual-native-density-observations-at-several-continuum-times):
V_l=240*4^l, q_l=V_l/3, global count floor(1260*16^l*tau), and f_l,g_l
the repeated cell-average profiles. The owning theorem proves joint weak
convergence of actual native occupation snapshots Y_l to the stationary
Gaussian field with covariance

    K_(tau,u)(f,g)=<f,exp(-|tau-u| A)g>.                   (10)

Uniform moments needed for products can be proved directly from the same
canonical state. If sum_i f_i=0 and |f_i|<=C, then for a uniform q-subset X,

    E exp(t sum_(i in X)f_i)
      =e_q(exp(tf_1),...,exp(tf_V))/binomial(V,q)
      <=(V^-1 sum_i exp(tf_i))^q
      <=exp(q t^2 C^2/2).                                 (11)

For completeness, fix the sum of positive variables z_i. Holding all but
z_a,z_b fixed, e_q(z) is z_a z_b e_(q-2)(rest) plus terms depending only
on z_a+z_b. Replacing the pair by its average increases e_q. Repeated
averaging converges to equal variables, proving the first inequality.
The chord bound on [-C,C] and mean zero give mean exp(tf_i)<=cosh(tC),
and cosh x<=exp(x^2/2) proves the second.

Thus Y_f is subGaussian with variance allowance
sigma^2=qC^2/(V chi), at most 3C^2/2 for q=V/3. Chernoff and integration
give, for every positive integer m,

    E|Y_f|^(2m) <= 2(2 sigma^2)^m m!.                      (12)

Every native snapshot has this same canonical marginal. Hölder therefore
bounds all required joint products uniformly, regardless of temporal
dependence. In particular eighth moments bound fourth moments of Y_fY_g.
The correction in (8) is another centered bounded-profile linear statistic,
with prefactor a/sqrt(V chi); (11),(12) give its vanishing L2 norm and
uniform higher moments. Consequently (8),(9), joint Gaussian convergence
and the continuous mapping theorem imply joint convergence

    W_l(tau;f,g) -> :Y(tau;f)Y(tau;g): .                  (13)

The uniform fourth moments make pairwise products uniformly integrable,
so this conclusion includes covariance convergence, not only weak laws:

    Cov(W(tau;f,g),W(u;h,k))
      =K_(tau,u)(f,h)K_(tau,u)(g,k)
       +K_(tau,u)(f,k)K_(tau,u)(g,h).                     (14)

This is the symmetric second Gaussian chaos of the inherited diffusive
field, with transfer generator A tensor I+I tensor A. It is not a newly
selected interacting continuum field: finite exclusion corrections coexist
with this Gaussian collective limit. No claim is made for unbounded L2
tests, unsmeared pointwise products, fixed microscopic contrast profiles,
or a quantitative spatial convergence rate. A different non-Gaussian
continuum sector needs its own observables, scaling and limit proof.

The finite matrix, collision comparison, readback and collective limit use
the stated same-q source contracts. They supply no Lorentzian action,
physical preparation or clock, chiral particle sectors, measured spectrum,
or completed source-to-field identification.
