# Native transfer correlation sectors and reflected records

The fair connected native source has an exact decomposition by correlation
order. A four-read experiment isolates one degree-two component, including
its measurement backaction. Its normalized record correlation has the same
transfer poles at every allowed global occupation on a fixed graph. This is
a decomposition of the population/history transfer, not a particle assignment
or a symmetry of each recorded quantum branch.

## Owning source and reused results

Let the declared native graph have V binary ports and M seams, with one
independent uniform seam selection, a fresh fair coin and repair probability
one per attempt. The C/E mixture is fixed. Native endpoint enumeration gives

\[
 P_q=I-G_q/(2M),\qquad G=\sum_{ab\in E}(I-S_{ab}),                 \tag{1}
\]

on the q-occupation configuration space, where S exchanges endpoint bits.
All hold, raw zero-address and terminal records remain in the instrument;
(1) sums them without conditioning. The graph, total occupation and selector
are declared inputs. The default biased H66 law is a different source.

The graph is supplied by the unchanged fair federation admission. Evaluated
geometry in this package is levels zero, one and two. Algebraic statements
below apply to every finite graph with this literal endpoint law. No claimed
all-level geometry, infinite-volume limit or host execution is inferred from
those three evaluations.

Existing results are reused rather than reproved as additional claims:
[federated-observer scaling](../federated_observer/SCALING.md) supplies density
closure, the full-sector bound 10I/21<=P<=I on its federation, and the uniform
logarithm series. [Native equilibration](../native_equilibration/README.md)
supplies the connected-sector stationary preparation and its finite error;
its exclusion-gap input is explicitly external. The
[squarefree dual](../observer_time_fields/REFINEMENT.md#exact-exclusion-dual-for-the-entire-four-bit-law)
supplies the monomial identity under endpoint exchange. The
[local effective action](../native_effective_fields/README.md) identifies
its occupation interaction. The result here is the orthogonal separation
of those monomials, its exact occupation-independent normalized response,
and the complete native record experiment that reads it.

The Boolean-slice harmonic decomposition is standard; see
[Filmus, Electronic Journal of Combinatorics 23(1), P1.23 (2016)](https://arxiv.org/abs/1406.0142).
The finite proof below is self-contained and attaches that algebra to the
native transfer and its recorded probe.

## Orthogonal correlation degrees

Let H_r have the orthonormal basis of r-subsets of the V ports. Define

\[
 U_r|A\rangle=\sum_{i\notin A}|A\cup\{i\}\rangle,
 \qquad D_{r+1}=U_r^*,\qquad W_r=\ker D_r.                       \tag{2}
\]

Take D_0=0. On each basis vector, additions and deletions at different sites
cancel in the commutator, and the same-site counts give

\[
 D_{r+1}U_r-U_{r-1}D_r=(V-2r)I.                               \tag{3}
\]

Thus ||U_r v||^2=||D_r v||^2+(V-2r)||v||^2, so U_r is injective
for r<V/2 and D_r is surjective for 1<=r<=V/2. Consequently

\[
 \dim W_r={V\choose r}-{V\choose r-1},\quad 0\le r\le V/2,      \tag{4}
\]

with the lower binomial zero for r=0. For h in W_r, induction in (3) gives
D U^k h=k(V-2r-k+1)U^(k-1)h. Define the inclusion
I_(r,q)=U^(q-r)/(q-r)!, which sends a coefficient vector h to the
squarefree field sum_(A subset X,|A|=r)h(A). It follows that

\[
 I_{r,q}^*I_{r,q}|_{W_r}={V-2r\choose q-r}I,
 \quad r\le q\le V-r.                                        \tag{5}
\]

Images of different W_r are orthogonal: move the adjoint lowering powers
across their inner product and use D_s g=0 for the higher primitive degree
s. Their dimensions telescope to binomial(V,q), giving the exhaustive
orthogonal decomposition

\[
 H_q=\mathop{\bigoplus}_{r=0}^{\min(q,V-q)} I_{r,q}W_r.           \tag{6}
\]

Every site permutation commutes with adding or deleting a site. Hence each
S_ab, G and P intertwines U,D and I. Since they are self-adjoint, W_r and
each image in (6) reduce the unread population transfer. The normalized
inclusion from (5) is unitary onto its image and identifies that block with
P_r restricted to W_r. This is an exact spectral nesting theorem, including
multiplicities. Different degrees can share eigenvalues; primitive degree
two does not imply a additional numerical eigenvalue or particle.

The existing positivity bound permits the same decomposition for
H_q=-log P_q, the faithful references P_q^L/Tr(P_q^L), and analytic
functional calculus. These references live on the transfer Hilbert space.
They are not the prepared working density I/binomial(V,q). The local
operator I-P is the generator only after the specified Poissonization of
attempts; its logarithm and unitary continuation do not become the native
working clock. U,D are algebraic intertwiners, not admitted particle-changing
commands. No SU(2) gauge or physical spin identification is made.

## A primitive degree-two observable

Choose two disjoint existing seams (a,b) and (c,d), with their ordered native
endpoints. Put

\[
 h=e_{ac}-e_{ad}-e_{bc}+e_{bd},\qquad
 \phi(X)=(n_a-n_b)(n_c-n_d)=I_{2,q}h(X).                         \tag{7}
\]

Each one-site lowering sum cancels, so D_2h=0 and ||h||^2=4. Under the uniform
q-sector law, phi is orthogonal to constants and every one-body density
field. Its variance is

\[
 v_q=4\,{\binom{V-4}{q-2}\over\binom Vq}
     ={4(q)_2(V-q)_2\over(V)_4}.                               \tag{8}
\]

This is positive precisely for 2<=q<=V-2. Equations (5) and (6) give, for any
nonnegative integer t,

\[
 C_q(t)=\mathbb E[\phi(X_0)\phi(X_t)]
       ={\binom{V-4}{q-2}\over\binom Vq}\,h^*P_2^t h,
 \qquad {C_q(t)\over v_q}={h^*P_2^t h\over4}.                  \tag{9}
\]

Thus the entire normalized curve is independent of global q on this fixed
graph and common attempt clock. It is not independent of graph capacity or
of the clock convention. At q=0,1,V-1,V the observable vanishes and the
normalization is undefined. The fair federation transfer's positive
spectral theorem writes the right side of (9) as a positive probability
mixture of lambda^t, with 10/21<=lambda<1. Strict inequality at one follows
from connectedness and orthogonality to the unique constant ground. Equal
poles in degree one and two are allowed. The coefficients are weights of
this finite scalar history statistic, not universal quantum effects.

## Complete four-E readout with backaction

The incoming E record is d=n_a-n_b, and the signed current obeys
n_a'=n_a-j,n_b'=n_b+j. Its outgoing contrast is therefore d-2j.
Read the two disjoint seams, retain their complete raw flags and coins, then
run t uniform native attempts, then read the same two seams again. Score the
first pair by s_out=(d1-2j1)(d2-2j2), and the last pair by s_in=d3*d4.
Disjointness ensures that each product belongs to one time cut: the outgoing
cut after the first pair, or the incoming cut before the last pair.

For the two-read Kraus operators K_h and F=diag(phi), the complete native
address maps give

\[
 F K_h=s_{out}(h)K_h,\quad K_hF=s_{in}(h)K_h,
 \quad\sum_h K_h^*K_h=\sum_h K_hK_h^*=I.                        \tag{10}
\]

Every fixed-seam E command has all 25 interpreter slots, including inactive
C and hold slots. The pair has 625 slots. Zero and abort-address maps have
zero weight and are retained; no successful-record renormalization occurs.
The exact local proof checks all 16 binary configurations on the four
ports. Its partial-injection identities extend by the identity on every
exterior and untouched reference, and restrict to every global q sector.

Starting with tau_q=I/binomial(V,q), the first scored instrument produces
sum_h s_out K_h tau_q K_h^*=F tau_q by unitality. The final scored effect is
sum_h s_in K_h^*K_h=F. The intermediate native channel preserves diagonal
operators and acts there as P_q. Therefore its complete observed score is

\[
 \sum_{\text{all records}}p(h)\,s_{out}(h_{past})s_{in}(h_{future})
 =\operatorname{Tr}[F\Phi^t(F\tau_q)]=C_q(t).                  \tag{11}
\]

This proves (9) for an actual backacting instrument. The two initial E reads
are included in the apparatus, not silently removed from a stationary
classical path. Using incoming contrasts on the past side is different:
on a four-port path at q=2,t=0 it gives zero, whereas the reflected score
gives 2/3. The exact enumeration retains this countercontrol.

If a source preparation is within half trace distance epsilon of tau_q
(with its stated reference interface), any common exact complete instrument
preserves that bound. Since the unnormalized score lies in [-1,1], its
expectation changes by at most 2epsilon. Dividing by the exact v_q amplifies
this to 2epsilon/v_q. Fixed q has v_q of order V^-2, so the normalized formula
is not a capacity-independent preparation or sample-cost guarantee. A
sample mean over independent repeats is an estimator with separate sampling
error; (11) is the population expectation. No numerical coefficient is
fitted to measured particle data. No survivor conditioning, arbitrary
adaptive selector, continuum state, physical vacuum, chirality or Lorentzian
dynamics is supplied.

The executable constructs and admits the full four-E program with exported
records and a nonadaptive uniform body. Its finite program domain requires
4+t<=2^32 attempts. The mathematical kernel identity holds for every finite
t; concatenated prefixes have their explicit record interface. The generic
expanded-tree evaluator can return RESOURCE_STOP before completing even the
short program, and the certificate retains an exact identity-channel stop
control. Neither a sparse expectation evaluation nor an admitted program
is reported as a completed expanded-tree run.
