# Connected native contact response

The fair connected source supplies a controlled dressed local resolvent.
Its independent-pair comparison threshold also tends to zero with capacity.
These two results distinguish a measurable finite contact response from a
proof that the isolated twelve-port pole persists as a bound particle.

## Source, reference and observation

Use the original fair federation with V=240*4^l ports, M=630*4^l seams,
independent uniform seam selection and repair probability one. Each cell
contains the original thirty icosahedral seams. Each of its three exterior
seams uses a different port and reaches a different neighbouring cell.
The maximum port degree is six. The pair action and attempt transfer are

    G=sum_e(I-S_e),       P=I-G/(2M),       0<=G<=24I.       (1)

The last bound follows from the pair-configuration graph's maximum degree
twelve. These are the actual exclusion maps derived from the fair C/E
instrument. No spatial hopping Hamiltonian, coherent real-time operation
or physical clock is inserted. The working reference for this result is
the canonical q=2 population state I/binom(V,2).

Take the two actual cell seams (0,1),(2,4), translated to a chosen cell.
Their contrast product h has four coefficients +/-1, squared norm four
and zero lowering sum at every port. Its zero extension remains in the
global primitive degree-two sector W2. The complete four-E experiment
uses outgoing contrasts on the past pair and incoming contrasts on the
future pair, separated by n actual unread attempts. Its raw mean is

    C(n)=v_V <u,P^n u>,   u=h/2,   v_V=8/[V(V-1)].         (2)

The source-derived [pair response](../native_pair_response/PROOF.md) proves
this identity including all zero labels, holds, destructive back-action,
and recorded currents. The past instrument is not replaced by a chosen
prepared vector. Its scored canonical output is diag(h)/binom(V,2).
`records.py` independently evaluates this output through the original
fixed-E kernels using all local occupation classes and their exact exterior
multiplicities, then checks the actual level-zero unread kernel and future
records. A fixed E repair has 25 slots; a pair has 625, of which 589 are zero.
The actual uniform attempt has 24M+1 slots. None is dropped or conditioned on.

The declared native programs use a one-slot FIFO with oldest-record export.
Admission and expanded execution are separate: the explicit branch-cap
control returns RESOURCE_STOP. Sparse scored factorization does not claim
that the expanded interpreter completes every admitted experiment. The
analytical family in n also does not override the attempt or level caps.
Source specifications are available to the observer; no assertion of
restricted specification access is made.

## The connected comparison threshold

Let Q2 project unordered pair coefficients onto W2. Compressing two
independent one-body Laplacians to distinct pairs and then W2 gives

    H0=Q2(G+2A_adj)Q2 | W2,                              (3)

where A_adj marks pairs whose ports share a source seam. This is the
same signed-contact comparison used by the
[finite contact proof](../native_contact_response/PROOF.md). On one
isolated cell its lowest observed pole is E>5, below the cell's centered
tensor threshold 2(5-sqrt(5)). That threshold cannot be transferred to the
connected family.

Here is an upper bound on the actual compressed threshold in (3), not an
inference from a vanishing lower bound. Give all twelve ports in a cell its
auxiliary affine triangle-centroid coordinate on the flat regular
icosahedron of circumradius one. These are comparison coordinates, not a
measured spatial dictionary. Shared midpoint identifiers define their
combinatorial refinement. Signed coordinate symmetries and cyclic axis
permutations preserve the centroid multiset and its cell adjacency at every
formal level. The centroid mean is zero and its coordinate Gram and energy
matrices are scalar. The face-plane inradius satisfies
r^2=(5+2sqrt(5))/15>3/5, while every centroid has radius at most one.
For any unit coordinate f, therefore,

    sum f_i=0,  |f_i|<=1,  V/5<A=sum f_i^2<=V/3,
    E_f=sum_edges(f_i-f_j)^2=E_l/3<4,
    E_l=20-4sqrt(5)+(-40+16sqrt(5))/(3*2^l).               (4)

The exact E_l formula and all-level combinatorial attachment are proved in
[native spatial scaling](../native_spatial_scaling/PROOF.md). Distinct
floating port choices do not alter (4), since the coordinate is constant
inside each cell. Successful routing must still satisfy the owning source
postconditions. This argument does not prove that a floating constructor
returns at every mathematical level. Literal controls evaluate levels 0–2.

Put b_ij=f_i f_j and

    c_i=(-f_i^2+A/[2(V-1)])/(V-2),
    w_ij=b_ij-c_i-c_j.                                   (5)

Since sum_j b_ij=-f_i^2 and sum_i c_i=-A/[2(V-1)], the lowering sum of
w is zero. This is exactly Q2 b, using D2 D2*=(V-2)I+J. With B=sum f_i^4<=A,

    ||w||^2=A^2/2-B/2-[B-A^2/(2(V-1))]/(V-2)
            >=A^2(1/2-5/[2(V-2)])>=A^2/3>=V^2/75.        (6)

The second inequality uses V>=240. Each transposition commutes with Q2,
so <w,Gw><=<b,Gb>. Directly exchanging a port along an edge gives

    <b,Gb>=sum_edges(f_i-f_j)^2(A-f_i^2-f_j^2)
           <=A E_f<=4V/3.                               (7)

Also |c_i|<2/(V-2), hence |w_ij|<=1+4/(V-2)<21/20.
Since M=21V/8, the contact energy obeys

    2sum_edges w_ij^2 <=9261V/1600<6V.

Combining (3),(6),(7) proves the variational bound

    0<min spectrum(H0) <=550/V=55/(24*4^l).              (8)

Strict positivity holds at each finite connected level because G has only
the constant zero mode and W2 excludes it. The isolated E>5 is above this
connected comparison bottom at every level, including level zero. The
bound does not exclude a bound band at fixed momentum, a dressed resonance,
other localized modes or a different source limit. Those require their own
definitions and persistence proof. It invalidates only the proposed use of
the isolated-cell threshold as a global binding certificate.

`threshold.py` evaluates (5) without storing all pair coefficients. It
checks exact Q(sqrt(5)) norms and Rayleigh numerators on the actual routed
graphs. The projected G energy subtracts
sum_edges(f_i^2-f_j^2)^2/(V-2) from (7), as follows from the same inclusion
intertwiner. Neither a numerical eigenvalue fit nor a quadratic trial outside
W2 is used.

## Controlled exterior dressing

Let Pi retain the 66 pair configurations wholly inside the chosen cell and
let Q=I-Pi. Write the exact global action in blocks

    G = [ A  B* ; B  D ],
    A=G_cell+N_boundary,     B*B=N_boundary<=2I.          (9)

Here N_boundary counts the occupied ports among the three crossing ports.
An outward move sends {a,b} to {outside(a),b}. Distinct crossing endpoints
and distinct neighbouring cells give distinct target pairs, even when both
ports are crossing ports. Thus the images do not interfere and (9) follows
exactly. All source exits are retained. D=QGQ obeys 0<=D<=24I.

For s>0 define the Schur self-energy and its finite exterior-path expansion

    Sigma(s)=B*(D+sI)^(-1)B,
    Sigma_J(s)=sum_(k=0)^J B*(24I-D)^k B/(24+s)^(k+1),
    F(s)=A+sI-Sigma(s),   F_J(s)=A+sI-Sigma_J(s).          (10)

The standard block inverse gives Pi(G+sI)^(-1)Pi=F(s)^(-1).
Its norm is at most 1/s, so F(s)>=sI. The expansion has positive operator
terms because 0<=24I-D<=24I. With r=24/(24+s),

    0<=Sigma-Sigma_J<=2r^(J+1)/s I,
    F_J>=F>=sI,
    0<=F^(-1)-F_J^(-1),
    ||F^(-1)-F_J^(-1)||<=2r^(J+1)/s^3.                 (11)

Inverse ordering proves the one-sided interval; the resolvent identity and
both inverse bounds prove its width. No primitive projection is inserted
inside (10): an exterior excursion can return a core component outside
cell W2. The zero-extended fixed four-E probe h is global W2, and the full resolvent
preserves that global subspace. These are compatible statements.

The matrices B*(24I-D)^k B depend only on exterior configurations reachable
in k single moves from the exit images. `dressing.py` enumerates these paths
with exact integer coefficients, then inverts only the 66-dimensional F_J.
The finite-depth calculation includes all diagonal holding degrees and
all signs; deleting an unvisited path is justified by polynomial locality,
not by a reflecting boundary or survival postselection. The bound in (11)
is independent of capacity. More depth supplies a certified approximation
for every s>0, although the implemented exact evaluator has its own finite
degree budget.

At s=24 and J=2 the normalized scalar response satisfies

    s<u,F_J^(-1)u> <= s<u,(G+sI)^(-1)u>
       <=s<u,F_J^(-1)u>+1/2304.                         (12)

For cell0 at each literal levels 0–2 the lower endpoint is between
0.660081931825 and 0.660081931826. The upper endpoint is below 0.660515959604.
The isolated-cell value is exactly33205121978357/49638960845760,
between 0.668932657182 and 0.668932657184. Their separation is greater than
1/120, even after the exterior tail allowance. This is a source-attached
finite observable discrepancy, not an inferred pole location. The finite
controls use 425 reachable exterior configurations at level 0 and 414 at
levels 1–2. Agreement of their first moments is not asserted as agreement
of the full connected response across capacities or all cells.

At positive spectral parameter near E, D-E may have spectrum arbitrarily
close to zero. Equations(10),(11) do not supply its inverse or a separating
gap there. Off-spectrum dressing therefore supplies no positive-energy
pole continuation, localization length or particle residue theorem.

## Complete native Abel readout and finite error budget

For the same native attempts and source reference, set q=2M/(2M+s).
The elementary geometric series of (1) yields

    s<u,(G+sI)^(-1)u>
       =(1-q)sum_(n>=0) q^n C(n)/v_V.                   (13)

The discount is an analysis parameter determined by s and M, not an added
scheduler or physical time unit. Finite truncation at N has error at most
q^(N+1), since |<u,P^n u>|<=1. Each term is a separately specified complete
native four-E program. A randomized choice of lag would need its own admitted
selector; (13) alone does not create one. A finite example at level 0 has
s=24,N=512 and tail less than1/10000. Its largest lag program uses 516
attempts, before any preparation prefix. One may compare this finite
statistic to (12) with the sum of the certified errors.
The exact rational tail evaluator separately caps its cutoff at 4096 and
its power-size allowance at 12000 bits before exponentiation. These are
analysis-evaluator budgets, not changes to native program admission.

The canonical reference in (2) is explicit. The existing
[scored preparation theorem](../native_correlation_preparation/PROOF.md)
can instead bound each normalized four-E score after actual unread warmup,
including arbitrary initial reference correlations in its stated domain.
If every normalized score differs by at most epsilon_p, its finite Abel
combination differs by at most epsilon_p because its nonnegative weights
sum to at most one. The preparation archive is retained and ignored by
this fixed score; no product approximation for that archive is implied.
All requested warmup and lag attempts must fit the native program caps.

For independent complete repetitions the raw score is in [-1,1]. Normalized
accuracy epsilon requires raw accuracy v_V epsilon. The inherited exact
Hoeffding budget is therefore at least 2log(2/alpha)/(v_V^2 epsilon^2), with
an explicit integer logarithm allowance in the producer. Applying failure
probability alpha/(N+1) per lag gives a union bound for the finite sum.
Repeated preparation and independent trials are premises, not free source
resources. Timeout or abort trials cannot be silently removed. The source
q=2 signal shrinks as V^(-2), so normalization does not eliminate its
sampling cost. Source state error with an arbitrary correlated reference
requires a joint bound or the stated scored-preparation theorem.

These constructions yield a controlled connected response and a precise
failure of the isolated binding comparison. They neither construct an
interacting relativistic particle sector nor identify an observed mass.
