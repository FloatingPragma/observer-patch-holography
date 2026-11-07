# Complete record response of the weighted native source

The source is the default program in `model/native_program/interpreter.py`:
66 binary occupation configurations with total load two, twelve ports with
weights `(25,25,8,8,8,8,8,8,8,8,8,8)`, uniform selection among the thirty native
seams, equal tagged C/E law weights, and the endpoint placement coin fixed by
those weights. One unit of time is one attempted repair, including zero
current. The one-record FIFO exports its oldest entry. No source pulse,
spatial coordinate chart, Hamiltonian, Higgs field or target spectrum is added.

The spatial labels in this calculation are the actual twelve ports and thirty
seams. They are not asserted to be points in three-dimensional physical space.
The calculation concerns population relaxation in this source, with its
actual recorded measurements and their backaction.

## Complete process and its reference

Write a configuration as `x`, the occupation of port `a` as `n_a(x)`, and the
port weights as `theta_a`. For seam `e=(a,b)`, coin one places the sole endpoint
load at `a` and has probability `theta_a/(theta_a+theta_b)`. Coin zero places it
at `b`. Empty and doubly occupied endpoints hold. The invariant probability is

```text
pi(x) = product_a theta_a^n_a(x) / 7505.
```

Its values are `625/7505`, `200/7505` and `64/7505`, with multiplicities 1, 20
and 45. Every seam is a conditional heat-bath update for this same measure.
Thus its transition satisfies detailed balance; their uniform mixture does
too. Connectivity gives a simple stationary eigenvalue.

For each complete raw flag `h`, let `V_h` be its positive partial-injection
matrix and `w_h` its squared coefficient. Its quantum operation is

```text
Phi_h(rho) = w_h V_h rho V_h^dagger.
```

The coefficient includes seam selection, the tagged law and the biased coin.
There are 720 raw C/E slots, including 420 zero slots, of which 240 are
reserved abort addresses. The interpreter also retains one global-hold slot of weight zero,
so this package retains 721 outputs. Empty slots have score zero. Complete
histories use the products of these operations, without dropping any flag or
renormalizing any branch. The full source/reference output is

```text
sum_h |h><h| tensor (Phi_h tensor id_R)(rho_SR).
```

The same definition applies to any finite word, with the ordered product of
its branch maps. Every word is a partial injection. Its effect is diagonal in
the occupation basis. Therefore every complete raw-record probability depends
only on the initial source diagonal. With an untouched reference, its joint
record/reference output depends on the initial diagonal source blocks
`rho_ii` on the reference. The final retained source/reference state can
contain off-diagonal matrix entries within a history fiber. Population
observability below is not tomography of the full quantum state.

Let `T_h[j,i]=w_h` when branch `h` maps input `i` to output `j`, and zero
otherwise. Then `P=sum_h T_h` is column stochastic. The code extracts every
`T_h` from the literal interpreter using classical input-reference markers,
and independently checks its sum against the endpoint transition. The matrix
formulae below are a contraction of the complete instrument, not a different
unrecorded process.

## Native measurements and exact backaction

The E flag reports the pre-repair contrast `d_e(x)=n_a(x)-n_b(x)`. Define an
unconditional score on every raw flag by

```text
Y_e(h) = 60 * 1[law=E and seam=e] * d_e(h).
```

The factor 60 is the reciprocal of the actual tagged-seam probability. All
other records, including C, zero currents and reserved addresses, have their
assigned score, usually zero. This statistic involves no conditioned success
sample. It uses only exported native records; it never inspects a hidden
occupation register.

Define `B_e=sum_h Y_e(h) T_h`, and let `D` be the matrix with rows `d_e`.
The native score is an unbiased measurement of the incoming contrast:

```text
1^T B_e = d_e^T.
```

It does not leave that incoming contrast intact. Put

```text
b_e = (theta_a-theta_b)/(theta_a+theta_b),
a_e(x) = 1[n_a(x)+n_b(x)=1],
g_e(x) = b_e a_e(x).
```

Conditional on the spectator occupations and endpoint total, the heat-bath
update independently resamples the outgoing placement from its equilibrium
law. In equilibrium the conditional incoming mean contrast is `b_e` on
endpoint total one and zero otherwise. Hence

```text
B_e pi = diag(pi) g_e.
```

The code checks this identity exactly for all thirty seams and every one of
66 output configurations. For each of the 22 equal-weight seams, `b_e=0`, so
`B_e pi=0`. Its recorded contrast has zero correlation, under a fixed continuation, with every later
function of future native records or final source state after the repair,
under the stationary source law. This means zero columns of the ordered
future/past contrast response. A zero row is not implied: an equal-weight
contrast can respond to an earlier biased-seam observation.

For example, seam zero is `(0,1)` and seam one is `(0,5)`. The same-seam
next-attempt covariance at seam zero is zero, while the covariance of seam
zero after seam one is exactly `5440/49533`. The same-time variance of seam
zero is positive. A passive expression `d_0^T P diag(pi) d_0` is also positive,
but it omits the measurement-producing repair and is not this record
correlation. Replacing the actual `B_e` with a passive diagonal insertion is a
failing control.

## Multi-time response and its complete finite resolvent

Let `m=D pi` and let `G` have rows `g_e`. At lag zero only one tagged seam can
occur, so

```text
C_ef(0) = 60 * delta_ef * E_pi[d_e^2] - m_e m_f.
```

For integer lag `k>=1`, direct composition of the same complete process gives

```text
E[Y_e(t+k)Y_f(t)] = 1^T B_e P^(k-1) B_f pi,
C(k) = D P^(k-1) A,
A = diag(pi) G^T - pi m^T.
```

There are `k-1` unread repairs between the two measured repairs. For increasing
attempts `t_1<...<t_r`, all higher ordered moments are

```text
1^T B_(e_r) P^(t_r-t_(r-1)-1) ... B_(e_2)
    P^(t_2-t_1-1) B_(e_1) pi.
```

These equations retain readout backaction at every measured event. The
package checks three-time instances against direct independent endpoint-path
enumeration. Indicator scores in place of the displayed contrast scores give
any specified complete raw-record event probability by the same composition.

Since `1^T A=0`, the stationary pole cancels. For `|z|>1` the exact matrix
response is

```text
R(z) = sum_(k>=1) z^(-k) C(k) = D (z I-P)^(-1) A.
```

The rational expression supplies its continuation away from its poles. The
finite spectrum belongs to the discrete repair-attempt process. It is
separate from the modular flow of the diagonal reference `pi`, and from the
analytic continuation of a reconstructed transfer Hamiltonian. Neither a
Lorentzian dispersion relation nor a particle rest mass follows from a pole
of this response.

## Every population relaxation mode is observable

The occupation contrast measurements supply an exact finite observation
matrix

```text
O = [1^T; D; DP; ...; DP^12].
```

Its rank is 66. The centered record insertions supply

```text
K = [A, PA, ... , P^10 A],
```

whose rank is 65. Exact elimination modulo the prime 65521 supplies nonzero
minors of these sizes. All rational denominators are invertible modulo this
prime, so these are lower bounds on rational, real and complex ranks. The
ambient dimension bounds the first rank above by 66; `1^T P^j A=0` bounds the
second above by 65. These matching bounds prove the stated ranks. Intermediate
modular ranks are lower bounds; no claim about their equality to every
intermediate rational rank is required.

Thus the record response is observable on the whole population space and
reachable on its centered subspace. Equivalently, every nonstationary
population eigenspace appears in some measured response component; none is
removed by restricting to the native E record contrasts. This is an exact
algebraic identification result. It does not assert that inversion of noisy
correlations is well conditioned, or that a finite experiment resolves every
pole.

The integer matrix `1980 P` has a squarefree characteristic polynomial. Its
irreducible factors have degrees 1, 13, 15, 15 and 22; the linear factor is
`x-1980`. The certificate retains every coefficient. An independent check
uses power traces and Newton identities modulo 65521. Detailed balance makes
`P` similar to a real symmetric matrix. Its minimum diagonal is `791/990`,
so stochastic Gershgorin bounds put every eigenvalue at least `296/495`.
Connectedness gives every nonstationary eigenvalue strictly below one. There
are therefore 65 distinct positive real nonstationary eigenvalues, each an
actual pole of the complete centered record response.

For a right and left eigenvector pair `v_lambda,w_lambda` normalized by
`w_lambda^T v_lambda=1`, the pole residue is
`(D v_lambda)(w_lambda^T A)`. Both factors are nonzero by the two rank
certificates, so each residue has rank one. These ordered stochastic-response
residues are not asserted to be positive quantum spectral weights.

They are relaxation factors in native attempts. In particular, they are not
65 particles. The decay of a positive reversible population channel supplies
no undamped relativistic propagation or spatial momentum dictionary.

## The weighted occupation hierarchy

Every raw record carries a decoded current `j`, and the actual source update
obeys the pathwise continuity law

```text
n_k(output)-n_k(input) = j * (1[k=b]-1[k=a]).
```

For a selected seam the conditional mean current is

```text
E[j|x,e] = (d_e(x)-b_e a_e(x))/2
          = (theta_b n_a - theta_a n_b
             + (theta_a-theta_b)n_a n_b)/(theta_a+theta_b).
```

Averaging seam selection yields the exact occupation evolution. Unlike the
fair source's closed linear mean law, the fixed biased source involves pair
correlations on its eight unequal-weight seams. The span of the twelve
occupation functions has dimension 12; adjoining their one-attempt images
raises it to 19. An explicit witness is the drift of port zero on pairs
`(0,1),(5,6),(0,5),(1,6)` with signs `+,+,-,-`: its rectangle contrast is
`17/990`. Every affine single-occupation function has zero contrast on this
rectangle. Therefore a twelve-component linear mean propagator cannot
represent the full biased source for all initial population distributions.
The 66-state response above supplies an exact closed description without
inserting a guessed kinetic law.

## Preparation, records and execution scope

The fixed tagged mixture is admitted by the finite preparation theorem in
`model/local_repair_family/mixing.py`. After 1600 uniform attempted repairs,
its reduced source/reference replacement error in half trace distance is at
most `542159/10^12`, uniformly for every initial source/reference state. The
comparison traces the preparation records. Those records are retained by the
source and can remain correlated with it; the result does not compare the
entire preparation archive to an independently prepared product state.

A fixed continuation and any complete future record event contract this
error, so its probability differs from the stationary prediction by at most
that bound. An r-fold contrast moment differs by at most
`2*60^r*(542159/10^12)`, since its score is bounded in absolute value by `60^r`.
Conditioning the continuation on selected preparation records is outside this
stationary prediction. The program keeps the fixed independent uniform seam
selector throughout.

Matrix contractions compute finite response laws exactly. They do not expand
or sample the entire `721^n` history tree. The interpreter's expanded-tree
resource stop is retained as such; it is not reported as a completed long
run. A finite response estimate has its own sampling error, and the rank
certificate alone supplies no stable pole-estimation bound.

These results produce a same-source native measurement dictionary, its
complete finite response, an exact continuity relation and a necessary pair
correlation hierarchy. They do not supply a continuum field, canonical
position and momentum, an interacting particle sector or a mapping to
measured particle masses. No physical premise is marked discharged.


## Reproduction

From the repository root, run:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -W error -m pytest -q model/native_record_response/test_response.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.native_record_response.response --check
```

The receipt binds its code, proof, tests, actual native interpreter, native
geometry and preparation sources. It stores the complete polynomial
coefficients and rank-minor row indices, with no floating spectral fit.
