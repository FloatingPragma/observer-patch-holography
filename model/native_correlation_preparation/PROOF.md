# Preparing native correlation records

Independent fair repair attempts prepare the expectation of a specified
four-read correlation experiment from any quantum state in a fixed binary
occupation sector. The error is controlled in a space of polynomials of
degree at most four. Its dimension is polynomial in the number of ports,
including sectors whose full Hilbert dimension grows exponentially.
This is an observable preparation theorem; it does not assert preparation
of the full density matrix or identify a physical vacuum.

## Source, reference and complete observation

Take a connected graph with V ports and M seams. Each warmup and lag attempt
selects a seam independently and uniformly, uses a fresh fair placement coin,
and repairs with probability one. A fixed C/E mixture is allowed. These are
the unchanged [fair native grammar](../native_program/interpreter.py) and
[federation endpoint law](../federated_observer/source.py). The biased
twelve-port law, adaptive seam selectors and a changed repair probability
are outside this theorem.

Fix 2 <= q <= V-2. Let tau_q be the uniform q-sector density and let P be
the population transfer. The
[native correlation decomposition](../native_correlation_sectors/PROOF.md)
gives P=I-G/(2M), with G the exclusion Laplacian, and exhausts its orthogonal
correlation degrees. The [native mixing theorem](../observer/QUANTUM_MIXING.md)
imports the exclusion consequence of the Caputo--Liggett--Richthammer gap
theorem and proves the nonconstant contraction factor

\[
 r=1-\gamma,\qquad \gamma=\lambda_1(L)/(2M)>0.                 \tag{1}
\]

No numerical eigenvalue is needed. The existing
[graph certificate](../native_equilibration/construction.py) supplies
gamma >= 2/(M V D) when D bounds the graph diameter. This follows from the
range/path lower bound lambda_1 >= 4/(V D). The primitive twelve-port source
also has the exact rational certificate gamma >= 11/240.

Choose two disjoint existing seams e=(a,b), f=(c,d). Define
phi=(n_a-n_b)(n_c-n_d), F=diag(phi), and

\[
 v_q=\mathbb E_{\tau_q}\phi^2
     ={4(q)_2(V-q)_2\over(V)_4}>0.                            \tag{2}
\]

After n warmup attempts, read E on e and f, perform t lag attempts, and
read E on e and f again. Retain all raw records, including the zero-weight
addresses. Score the past pair by its outgoing contrasts
(d1-2j1)(d2-2j2), and the future pair by its incoming contrasts d3*d4.
The product S lies in [-1,1]. Warmup records are ignored without conditioning
when forming this observation. They can remain in the complete physical
archive; the comparison does not assert a product law for that archive.
The observer uses the declared native working-attempt clock throughout.

## Exact low-degree effect, including backaction

For a seam g, write H_g u=(u+u composed with endpoint swap_g)/2.
The [squarefree exclusion dual](../observer_time_fields/REFINEMENT.md)
identifies this with the unread native adjoint on every diagonal polynomial.
Both C and E have this same diagonal adjoint. Every raw Kraus map is a
scaled partial injection, so its adjoint sends a diagonal effect to a
diagonal effect, including arbitrary off-diagonal input matrix units.

The complete two-E identities F K_h=s_out(h)K_h and
K_h F=s_in(h)K_h give the exact initial scored effect

\[
 A_t=H_e H_f[\phi\,(P^t\phi)].                               \tag{3}
\]

The two past readouts are included in (3); they are not nondisturbing probes.
For example, on the four-port path at q=2 and t=0, A_0=phi^2 takes both zero
and one, while its stationary mean is 2/3. Replacing A_t by its stationary
mean on arbitrary input would omit preparation.

The dual P preserves squarefree degree, multiplication uses n_i^2=n_i,
and H_e,H_f average site permutations. Thus A_t has degree at most four.
Writing C_q(t)=E_tau A_t, the complete reflected-read theorem gives
C_q(t)=E_tau[phi P^t phi]. Since |phi|<=1 and P,H_e,H_f contract the
uniform-sector L2 norm,

\[
 \|A_t-C_q(t)\|_2\le\|A_t\|_2
       \le\|\phi P^t\phi\|_2\le\sqrt{v_q}.                  \tag{4}
\]

This argument also covers different disjoint-seam contrasts phi_i and phi_j
at the past and future cuts: replace (3) by
H_(past1) H_(past2)[phi_i P^t phi_j]. Each contrast has the same v_q and
the bound (4) is unchanged. It applies entry by entry to a finite declared
cross-response menu; arbitrary smearing coefficients require their explicit
sum-of-absolute-coefficients error allowance.

## A polynomial-size pointwise mixing bound

Set s=min(4,q,V-q). The space of polynomials of degree at most four restricted
to the q-slice has dimension binomial(V,s), by the existing exhaustive
harmonic decomposition. Its mean-zero subspace W therefore has dimension
N=binomial(V,s)-1. Although the native graph need not be vertex-transitive,
W is invariant under every site permutation.

Choose an orthonormal basis u_1,...,u_N of W in L2(tau_q). The diagonal of
its reproducing kernel, sum_i |u_i(X)|^2, is constant on the slice: site
permutations preserve W, its orthogonal projector and the uniform measure,
and act transitively on q-subsets. Its uniform average is N, so its value
at every X is N. Cauchy--Schwarz therefore proves
|g(X)| <= sqrt(N)||g||_2 for every g in W. P preserves W and contracts it
by r. Applying this fact to (4) gives

\[
 \|P^n A_t-C_q(t)\|_\infty
       \le\sqrt{N v_q}\,r^n.                               \tag{5}
\]

The graph gap in (1) is the only graph-dependent mixing input. No assertion
of a graph-independent mixing rate or a degree-four gap larger than the
usual exclusion gap enters (5). It holds for every finite lag t; the
implementation evaluates only its explicitly bounded sparse controls.

## Arbitrary quantum input and untouched reference

Let rho_QR be any normalized positive density on the fixed q source and an
arbitrary untouched reference R. Let Z_R be the signed sum of the terminal
R blocks weighted by S, after summing every preparation and experiment
record. Pulling the effect back through the n warmup attempts gives

\[
 Z_R-C_q(t)\rho_R
   =\sum_X [(P^n A_t)(X)-C_q(t)]\,\rho_{R,XX}.                \tag{6}
\]

Every rho_(R,XX) is positive and their traces sum to one. The trace norm of
(6) is at most the supremum in (5), with no assumption that Q and R are
independent. Source coherences do not affect this signed observation because
its exact effect is diagonal. Consequently

\[
 {\|Z_R-C_q(t)\rho_R\|_1\over v_q}
       \le\sqrt{N/v_q}\,(1-\gamma)^n.                        \tag{7}
\]

The scalar expectation follows by tracing R. More generally, any one
bounded reference observable has the same error. Equation (7) does not
bound the full pointer/reference probability law, a conditional law given
warmup records, or an adaptively chosen future experiment. At q=0,1,V-1,V
the contrast variance vanishes; no normalized statement is made.

## Other complete low-degree observations

The kernel argument and (6) apply to any real diagonal effect A of degree
at most four. If a separate exact instrument attachment identifies A as its
scored effect, a proved bound ||A-E_tau A||_2^2 <= u and a positive declared
normalization z give normalized error sqrt(N u)/z times r^n. This statement
does not infer an instrument from a supplied polynomial. The function
`polynomial_budget` consumes these explicit norm and normalization premises.
For a bounded sparse polynomial A=sum_a c_a n_a, `polynomial_moments`
derives them exactly from

\[
 E_\tau n_a={(q)_{|a|}\over(V)_{|a|}},\qquad
 E_\tau A^2=\sum_{a,b}c_a c_b
                     {(q)_{|a\cup b|}\over(V)_{|a\cup b|}}.
\]

The zero falling factorial terms are retained. These formulas keep the
fixed-number correlations; independent occupation is not assumed.

For a retrospectively read native trajectory, two degree-two snapshots
F_i(X_0), F_j(X_t) have effect A=F_i P^t F_j. This is also degree at most
four. A final native sector readout and reverse-current reconstruction
must establish that observation map; it is not the four-E apparatus in
(3), and no past seam averages are inserted. If |F_i|<=a and
||F_j||_2<=b, then ||A-E_tau A||_2<=ab. The sum of the absolute polynomial
coefficients is one valid bound for a. The helpers `snapshot_effect` and
`polynomial_moments` supply this algebra without changing the source or
introducing a readout service.

## Integer budgets, resource limits and sampling

Let gamma_0 be a proved positive rational lower bound in (1),
b=ceil(1/gamma_0), and c=ceil(sqrt(N/v_q)). Bernoulli's inequality gives
(1-gamma_0)^b <= 1/(1+b gamma_0) <= 1/2. With
B=ceil(log2 c)+p, n=bB native warmup attempts make the right side of (7)
at most c/2^B <= 2^-p. The executable computes all these integers using
rational comparisons and integer square roots. It never substitutes an
uncertified floating gap or logarithm.

For fixed q=2, the normalized prefactor grows as O(V^2). At a fixed density
strictly between zero and one it also grows as O(V^2), since s=4 and v_q is
bounded away from zero. The gap and therefore the warmup duration can
depend strongly on capacity. The separately proved global-state comparison
would use binomial(V,q)/v_q as prefactor. The scalar improvement does not
assert rapid full-state mixing.

The source program includes all n warmup repairs, t lag repairs and four
E reads. Its complete cost is n+t+4, with FIFO overflow set to export.
A single program requires n+t+4 <= 2^32 and the unchanged description and
record-storage admission bounds. Each admitted program retains every
selector, law, coin, raw zero and committed record. A mathematical sufficient
budget beyond this ceiling is reported as outside single-program admission.
It is not executed by extending the interpreter's limit.

The exact H66 control uses the proved rational gap 11/240 and all sixty
population counts imported from the actual complete native instrument. It
propagates the four-E adjoint through the certified warmup count with integer
numerators, retaining all 66 input columns. Thus its reported directed
80-bit rational error envelope covers every initial source density and
untouched reference by (6). This is an evaluated marginal calculation;
it is not a completed expansion of the full record tree.

The expanded-tree evaluator has separate branch and step limits. Its
RESOURCE_STOP is retained and never called a completed correlation run.
If an actual execution is stopped at prefixes with total probability p_stop,
give these outcomes score zero and compare with the same native kernel's
uninterrupted completion. Complete positive branch extension bounds the
additional unnormalized reference-score trace norm error by p_stop. The
normalized bound is therefore (7) plus p_stop/v_q. This term can be large;
success conditioning does not remove it. The proof concerns exact native
maps; a changed numerical implementation needs its separate error bound.
For a separately implemented generic score bounded by S_max, the stop
allowance is S_max p_stop divided by its declared normalization.

Independent experimental repetitions with score S in [-1,1] have sample-mean
variance at most 1/R for R runs. Thus Chebyshev's inequality bounds a
normalized sampling error exceeding eta by 1/(R eta^2 v_q^2). This sampling
term is separate from (7). Independence of repetitions and their fixed
preparation protocol are explicit hypotheses, and no practical sample-count
or physical mass conclusion follows from the expectation theorem.
