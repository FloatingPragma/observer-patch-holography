# Gaussian fields as limits of native occupation records

Written source theorem; no new simulation or physical-time identification.
This extends the one-unit resolvent interpretation by deriving the Gaussian
higher-moment law from long native record blocks. It does not impose Wick's
rule as a new field dynamics.

## 1. Finite source and observations

Let G be a finite connected undirected graph with M vertices and positive seam
weights w_e. Put W=sum_e w_e and L=B diag(w_e) B^T. One integer load evolves by
the unchanged native endpoint event: choose e with probability w_e/W and a fair
placement coin. The exact transition on its location is

    P=I-L/(2W).

The uniform law pi_i=1/M is stationary. Connectivity makes eigenvalue 1 simple.
For each edge the fair update is (I+Swap_e)/2, an orthogonal projection. Their
probability mixture is symmetric positive semidefinite; hence on the zero-sum
subspace its largest eigenvalue r satisfies 0<=r<1. We do not assume that the
global multi-load source is symmetric.

The observer holds a disclosed initial location (or a successful current-only
census with its separately charged error) and every seam/current record.
Continuity reconstructs the current location X_t. Define the vector record

    Z_t=e_(X_t)-(1/M)1,
    Y_n=sqrt(M/(4W n)) sum_(t=1)^n Z_t.

All coordinates are functions of the retained tape. The reference uniform
preparation, weights, counting window, normalization and linear readout are
declared apparatus. This does not read a hidden working state. The scalar
fluctuation Y_n is not identified with the original moving load as a particle.

Let Pi=I-11^T/M and L^+ be the inverse of L on Pi, zero on the constant mode.

## 2. Exact finite covariance and an explicit error

Stationarity and reversibility give, for k>=0,

    E[Z_t Z_(t+k)^T]=(1/M) Pi P^k Pi.

Summing over all ordered pairs of times proves

    C_n=E[Y_n Y_n^T]
       =(1/(4W))[Pi+2 sum_(k=1)^(n-1)(1-k/n) Pi P^k Pi].       (1)

The convergent infinite sum is

    C_infty=(1/(4W)) Pi(I+P)(I-P)^(-1)Pi
           =L^+ - Pi/(4W).                                    (2)

Both equations describe the actual records. On a P-eigenvector of eigenvalue
0<=z<1, their difference is exactly

    C_infty-C_n = z(1-z^n)/(2W n(1-z)^2).                     (3)

The identity follows by summing the geometric series and its derivative. It
includes z=0. In particular the matrix difference is positive semidefinite and

    ||C_infty-C_n|| <= r/[2W n(1-r)^2].                       (4)

This is a finite-window population bound, not a confidence interval for an
estimate from finitely many blocks. Sampling uncertainty is additional.

For the twelve-port unweighted graph, W=30 and the existing exact native
spectrum gives r=(55+sqrt(5))/60<23/24. Consequently the conservative rational
bound ||C_infty-C_n||<=46/(5n) holds. No experimental value fixes that bound.

## 3. Gaussian limit and higher moments

For any real zero-sum test vector g, let S_n=sum_t g(X_t). For complex theta
near zero, put D_theta=diag(exp(i theta g_i)). The exact characteristic function
of S_n is

    pi^T D_theta (P D_theta)^(n-1) 1.                         (5)

The simple isolated eigenvalue 1 of P continues to an analytic eigenvalue
lambda(theta) of P D_theta near zero. This follows in finite dimensions from
the simple root of its characteristic polynomial and the analytic implicit
function theorem. The spectral projection is analytic by the resolvent formula;
all other eigenvalues remain uniformly inside a circle of radius less than 1.
Thus their contribution to (5) tends to zero exponentially, uniformly on a
sufficiently small complex neighborhood. At zero the leading prefactor is 1.

Differentiating the eigenvalue equation with normalization pi^T v(theta)=1
gives

    lambda'(0)=i pi^T g=0,
    v'(0)=i (I-P)^(-1) P g,
    lambda''(0)=-<g,(I+P)(I-P)^(-1)g>_pi.                    (6)

Insert theta=u/sqrt(n). Taylor expansion of log lambda gives convergence of
(5) to exp(-u^2 sigma_g^2/2), where sigma_g^2 is the nonnegative expression in
(6). The continuity theorem for characteristic functions therefore gives the
one-dimensional centered normal limit. Applying the same argument to every
linear combination of the finite vector gives

    Y_n -> Normal(0, C_infty) in distribution.                (7)

The argument is uniform on compact complex u sets: n log lambda(u/sqrt(n))
has a uniformly vanishing cubic remainder there, and the prefactor tends to 1.
Uniform convergence of these analytic moment-generating functions implies
convergence of every fixed mixed moment by differentiation at zero. Hence the
limiting even moments obey Gaussian pairings and the odd moments vanish.
They follow from the native record law, rather than a separately declared
Gaussian field measure.

This theorem fixes the graph first and then takes n to infinity. It supplies
no uniform central-limit rate for a simultaneously growing graph. Equation (4)
controls covariance uniformly once W and r are specified. A joint refinement
of higher moments needs a diagonal sequence or an additional quantitative
bound; finite resource claims must charge that choice explicitly.

## 4. Connection to the existing slice theorem

On the declared four-axis product routing of E0/E1, L is the source operator
L_int+a^(-2)L_external. The reviewed slice-refinement theorem derives the
massive Green covariance L^-1, native internal gaps and common limiting speed.
Equation (2) differs from L^+ only by Pi/(4W), a local contact term plus the
explicit finite-volume constant-mode projection. For spatially separated
supports, the identity part has no cross term. For mean-zero observables the
constant-mode part vanishes as well.

At fixed finite graph the Gaussian limit therefore has the same cross-slice
Green correlations between disjoint supports as the native resolvent, provided
the zero-mode projection and its boundary conditions are kept explicit. It is
not legitimate to discard the finite-volume constant subtraction for arbitrary
test functions: a reflection-positivity claim must use an admissible zero-mode
treatment or prove the infinite-volume limit. The scalar source's fourth
coordinate remains a supplied Euclidean slice, not native event time.

In an infinite-volume limit where the mean-mode projection vanishes locally,
the Green covariances and their derived Gaussian moments can feed the usual
Euclidean field reconstruction. The physical time/intervention correspondence,
actual reflection-positive domain, same-source matter/light comparison and
unit attachment remain separate source statements. The previous two-load
reflection-positive product experiment is a different observable family; its
failure is neither erased nor transferred to this collective-fluctuation limit.

## 5. What this result changes

The old E0/E1 one-unit covariance did not determine higher correlation
functions without an added Gaussian rule. For the explicitly defined native
long-block occupation fluctuation, the finite-source theorem above now selects
Gaussian higher moments. Its source is native repair, not a scalar Hamiltonian
installed and evolved. Its normalization and limiting regime are visible.

There is no claim here that all native records, all preparations, individual
loads, or interacting particle observables are Gaussian. Declared topology,
the extended scheduler and the collective decoder make a new candidate when
executed beyond registered A-cont. No full OpenWave criterion or accepted cell
is asserted by this written theorem.

## 6. Native coin completion of the contact term

An explicit enlarged observation protocol removes the contact term. For each
vertex i and each block step t, execute an auxiliary native event on a disclosed
zero-load pair and retain its returned coin c_(t,i). Both endpoints remain zero
and the returned current is zero. Independent fair source coins give independent
xi_(t,i)=2c_(t,i)-1, also independent of the occupation trajectory. These events,
their preparation, scheduler and records are part of the enlarged candidate.
No Gaussian draw is supplied. The auxiliary pair may be reused serially; no
reset is needed because every auxiliary endpoint remains zero.

Define

    V_n=(1/sqrt(4W n)) Pi sum_t xi_t,   F_n=Y_n+V_n.

Independence gives exactly Cov(V_n)=Pi/(4W) for every n. Its multivariate
central limit follows directly from the product characteristic function
prod_i cos(u_i/sqrt(n))^n, with Pi and the displayed normalization applied.
Its moments converge by the same compact analytic argument. Independence
and section 3 imply jointly

    F_n -> Normal(0,L^+), with every fixed mixed moment converging.

Moreover L^+-Cov(F_n) is exactly the matrix in equation (3), so the explicit
finite-window bound (4) is unchanged. The contact completion coefficient
1/(4W) is an explicit decoder normalization obtained from the native transfer
identity. It is not a measured coupling or a physical energy normalization.
This is a stronger protocol than the occupation-only one in sections 1–5;
its auxiliary records and cost must never be omitted from an execution claim.

## 7. Rooted observables give the finite Dirichlet Green law

Choose a distinguished vertex o and let D have rows e_i^T-e_o^T for i!=o.
Let L_D be the principal minor of L obtained by deleting row and column o.
Connectivity implies L_D is positive definite. Then

    D L^+ D^T=L_D^(-1).                                      (8)

To prove (8), for any vector b on the remaining vertices, extend it by
b_o=-sum_i b_i. The solution u=L^+ b satisfies Lu=b. Subtract u_o from every
component to obtain a solution with value zero at o. Its other components
are D L^+ D^T b and solve L_D v=b, proving the identity.
Thus the actual rooted fluctuation records D F_n have the full Gaussian
Dirichlet Green law in their block limit. No zero-mode removal is implicit.
For fixed graph their finite covariance error is at most
||D||² r/[2W n(1-r)²]; ||D||²=M for this choice of root.

For a finite lattice box, join each interior site to its usual interior
neighbors and join each missing boundary-neighbor slot to the single root
with that slot's seam weight. Parallel addresses implement integer weights.
The principal minor is precisely the Dirichlet lattice Laplacian. This
declares a graph and an observation reference, not a Gaussian update law:
the only evolving source event remains integer endpoint repair, including
the auxiliary coin events of section 6. A source unit may visit the root;
the observer's fluctuation difference is what has the zero boundary value.

## 8. Finite scalar reflection positivity on a symmetric box

Take a reflection-symmetric nearest-neighbor box and a bond reflection
through a plane between two slices. The Gaussian limit of the rooted
records has density proportional to exp(-phi^T L_D phi/2). This density is
derived from sections 6–7, not used to generate the records. Each edge
crossing the plane pairs a variable x_i on one side with its reflection
y_i and contributes exp(w_i x_i y_i), after putting its diagonal quadratic
terms into the two equal side factors. Expanding

    exp(sum_i w_i x_i y_i)
      =sum_(multi-index k>=0) prod_i(w_i^k_i/k_i!) x^k y^k

shows that the reflected expectation of conjugate f on one side times f
on the other is a sum of nonnegative coefficients times squared absolute
side integrals. Absolute values are dominated by replacing each cross
product with its absolute value; on every orthant this is the reflected
positive Gaussian form. This supplies an integrable Gaussian majorant
for polynomial f and justifies termwise integration. Each coefficient is
nonnegative because every crossing weight is positive. Edges within a
slice, including internal port edges, only change the identical side
weights. A nonnegative mass term also changes only these weights.

This is a finite mathematical reconstruction statement for the specified
native-record limit. It does not identify the slice coordinate with native
event time or realize a physical measuring instrument for the reconstructed
unitary. Maxwell vector reflection needs its own gauge-invariant positivity
proof; scalar positivity alone does not prove that assertion.
