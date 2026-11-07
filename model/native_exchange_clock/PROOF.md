# Native attempt time and the local exchange action

The fair native source fixes a local exchange operator and a discrete record
transfer on the same occupation space. Their exponential responses have a
controlled comparison in a clock computed from native attempt counts. This
comparison holds in every occupation sector. Convergence of the generators
themselves has a narrower bounded-correlation-degree domain.

## Source, state and observed response

Use the unchanged binary federation graph with V ports, M seams and maximum
degree Delta<=6. Every attempt selects a seam independently and uniformly,
uses a fresh fair endpoint coin, and repairs with probability one. The C/E
mixture is fixed. All raw records, zero outcomes and exports belong to the
native instrument. Its unread population transfer on the q-occupation sector
is

    G_q=sum_e(I-S_e),        P_q=I-G_q/(2M),
    tau_q=I/binomial(V,q).

Here S_e swaps the two endpoint occupations. These identities and the
uniform bound (10/21)I<=P_q<=I are supplied by the owning
[federation construction](../federated_observer/SCALING.md). The source graph,
selector and global occupation are declared specifications. Reading these
specifications does not require a physical calibration. The default weighted
H66 member is a different source and is rejected by this attachment.

On L2(tau_q), P_q and G_q are self-adjoint and the constant function is a
ground vector. The native working state tau_q is distinct from that ground
vector and from the transfer reference P_q^L/Tr(P_q^L). No equality of these
states is used.

The [correlation-sector construction](../native_correlation_sectors/PROOF.md)
supplies primitive vectors h_i from two disjoint native seams and
phi_i=I_(2,q)h_i=(n_a-n_b)(n_c-n_d). With 2<=q<=V-2 their common variance is

    v_q=4 binomial(V-4,q-2)/binomial(V,q)>0.

The complete four-E experiment scores the outgoing contrasts of the first
pair and the incoming contrasts of the last pair. Its stationary expectation
after n intervening attempts is exactly

    C_ij(n)/v_q = h_i^* P_2^n h_j/4.                         (1)

The two probe pairs may differ; each pair has disjoint endpoints. The fine
Kraus identities F_i K=s_out K and K F_j=s_in K prove (1) by unitality and
diagonal preservation, including the two probes' backaction. Their complete
flags are retained. A numerical sparse evaluation of (1) is not an expanded
execution of every history. The native interpreter's finite program and
RESOURCE_STOP boundaries continue to apply.

## Uniform endpoint comparison without a Poisson scheduler

Set a=10/21. For p in [a,1], put delta=1-p and
r=-log(p)-delta. The convergent positive logarithm series gives

    0<=r<=delta^2/(2a).

For every integer n>=1,

    0<=exp(-n delta)-p^n
      =exp(-n delta)(1-exp(-nr))
      <=n delta^2 exp(-n delta)/(2a)
      <=2/(a e^2 n)<3/(5n).                                (2)

The maximum of x^2 exp(-x) on x>=0 is 4/e^2. The rational last bound uses
e^2>sum_(j=0)^4 2^j/j!=7. Spectral calculus therefore proves, on every q,

    ||P_q^n-exp[-nG_q/(2M)]||<=min(1,3/(5n)).               (3)

At n=0 the difference is zero. Define native rescaled time t_n=n/(2M).
The factor 2M is fixed by the fair endpoint rule and the global selector;
it is a convention in counted attempts, with no seconds or physical energy
unit assigned. Equation (3) compares the actual deterministic n-attempt
endpoint transfer with the local exchange response at t_n. The latter is
computed from the same G_q and the same stationary reference.

For a specified t>=0, let n=floor(2Mt), kappa=2Mt-n. If n>=1,

    ||P_q^n-exp(-tG_q)||<=min(1,(3/5+kappa/2)/n).            (4)

Indeed the rounding difference is bounded by
kappa delta exp(-n delta)<=kappa/(en)<kappa/(2n). For n=0,
the bound one is valid; at t=0 the error is zero. These statements are
uniform in q, the sector dimension and the number of modes. At a fixed
positive t the bound tends to zero as M grows. They compare two operators
on each finite graph. Convergence of the graphs, observables and states
to one continuum is a separate condition.

For arbitrary normalized L2 vectors f,g, (3) bounds their matrix element
error by the same number. It therefore bounds every normalized entry in
(1). For a family u_i=h_i/2, a linear combination with coefficient vector c
has quadratic-form error at most the bound times ||sum_i c_i u_i||^2.
Equivalently, the matrix error is bounded by the bound times the norm of
the probes' Gram matrix. No orthogonality of different probes is assumed.

## Local action on a bounded correlation degree

The Boolean-slice decomposition identifies the primitive degree-j block in
every sector j<=q<=V-j with the same reducing block of G_j. On j-subsets,
G_j is the graph Laplacian of the exclusion configuration graph. Each
configuration has at most j Delta available exchanges. Its quadratic form
is a sum over crossing edges of |f(X)-f(Y)|^2, so

    0<=G_j<=2j Delta I.

Thus the orthogonal sum of degrees zero through r has spectral ceiling
Lambda=2r Delta, independent of q and V. Let M>r Delta and define the exact
rescaled transfer logarithm H_M=-2M log(P_q). The positive series gives

    0<=H_M-G_q
       <=G_q^2/[4M(1-r Delta/M)]
       <=epsilon_r I,
    epsilon_r=r^2 Delta^2/(M-r Delta).                      (5)

All inequalities in (5) are restricted to those reducing degree blocks.
Consequently

    ||P_q^n-exp(-t_n G_q)||<=min(1,t_n epsilon_r),
    ||P_q^floor(2Mt)-exp(-tG_q)||
       <=min(1,t_n epsilon_r+(t-t_n)Lambda).                (6)

The certificate uses the smaller of (3) and (6). The same generator
estimate bounds the abstract unitary continuations by |t|epsilon_r;
this operator identity does not execute either unitary on the native
working carrier. In particular G_q and H_M are distinct from the working
reference's modular generator -log(tau_q), which is scalar on fixed q.

Equation (5) supplies a local action approximation from the actual discrete
transfer on bounded correlation degree. Its coefficient is fixed by native
counts. The exchange operator contains the density interaction
-2 n_a n_b at fixed hopping, as proved in the owning
[effective-field construction](../native_effective_fields/README.md).
This does not identify a chiral gauge action or a relativistic Hamiltonian.

## Why generator and response bounds have different domains

The full-sector generator difference need not tend to zero. Under uniform
q-occupation, the mean number of crossing seams is

    Tr(G_q)/binomial(V,q)=2M q(V-q)/[V(V-1)].

The largest eigenvalue is at least this mean. Since
-2M log(1-x/(2M))-x>=x^2/(4M), positivity and spectral calculus give

    ||H_M-G_q||>=M [q(V-q)/(V(V-1))]^2.                    (7)

For even V and q=V/2 this is at least M/16. The full-sector difference in
generators therefore grows with the number of seams even though (3)
gives a vanishing full-sector endpoint response error at n proportional
to M. Rapidly decaying modes account for this distinction. Equation (7)
does not obstruct the full-sector semigroup comparison.

## Exact finite response enclosures and preparation error

Sparse applications of the literal local G give the rational moments
mu_k=h^*G^k h/||h||^2. For even m, the degree-(m-1) and degree-m Taylor
polynomials bracket exp(-x) for x>=0. Therefore

    sum_(k=0)^(m-1) (-t)^k mu_k/k!
      <= h^*exp(-tG)h/||h||^2
      <=sum_(k=0)^m (-t)^k mu_k/k!.                         (8)

The enclosure width is at most (t Lambda)^m/m!. The same odd/even Taylor
argument applied to (1-x/(2M))^n brackets the actual native response by

    sum_k (-1)^k binomial(n,k) mu_k/(2M)^k,                 (9)

truncated at the same orders; terms beyond n vanish. The derivatives have
the required alternating signs on 0<=x/(2M)<=1. Clipping these bounds to
[0,1] uses the positive spectral measure. All retained intervals and moments
are exact fractions; no floating eigenvalue or fitted pole is a premise.

If preparation supplies an unnormalized complete-score error eta relative
to (1), the final normalized response error is at most eta/v_q plus (3)
or (6). In particular a half-trace preparation error epsilon gives
eta<=2epsilon for the score in [-1,1]. A direct score-preparation theorem
may give a smaller eta without preparing the full quantum state. The
division by v_q is explicit; fixed small q does not give a volume-independent
statistical or preparation cost. Arbitrary-reference preparation estimates
keep their own reference interface. The endpoint operator bound alone does
not upgrade a scalar response comparison to a quantum-instrument comparison.

## Complete-history and physical boundaries

The proof of (3) is spectral, and introduces no random attempt count.
Although the exponential can also be represented by a Poisson sum, that
representation has a different full record-length law. At deterministic
n=1 the length equals one with certainty; a Poisson mean-one length has
probability exp(-1) of length one. Their complete length laws have total
variation 1-exp(-1)>1/2. Thus an endpoint approximation cannot be asserted
as equivalence of complete native histories.

The complete C/E quantum process remains its specified dissipative
instrument. A field-clock identification must preserve the accessible
observations it claims to represent, including any additional coherent
services. The result here fixes a finite local exchange response and its
relation to counted native attempts. It supplies no physical clock,
source-selected unitary implementation, common Lorentz cone, continuum
state/action limit, chiral matter theory or particle mass identification.
