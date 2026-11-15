# Joint time and field laws on matched source screens

The positive incident-clock transfer in [the field construction](README.md)
and the [stopped-clock screen comparison](../native_action/INCIDENT_REFINEMENT.md)
apply to one common finite-density preparation and one common field protocol.
The resulting comparison concerns a complete four-bit record distribution,
including record-controlled native feedback. It does not infer a joint law
from a matching gap or covariance.

The preparation is an independently supplied Bernoulli product state at
density `rho=1/3`. This extends the positive-transfer dictionary beyond the
earlier fixed-total-two executions. On the two full screens it gives exactly
matching restrictions to the common ball and nonzero local fluctuations.
No physical duration per repair, fitted Hamiltonian, or independent-particle
approximation is used.

## One executable observer domain

This is a diagonal recorded branch of the
[canonical primitive observer](../primitive_observer/README.md). C and E
have the same diagonal native transition law. Their coherent continuations
are not identified by this experiment.

| Interface | Instantiation |
| --- | --- |
| State | Binary occupations on the actual full source screen. The certified common-ball calculation retains the exact exclusion action on field polynomials. |
| Preparation | Independent Bernoulli occupations of density `rho`, with occupations at root ports zero and one initially disclosed as `a,b`. Both disclosures are retained, including every label. |
| Native operation | Uniformly select a source seam and execute the existing fair-coin native repair. The unread diagonal map on a seam `e` is `Q_e=(I+SWAP_e)/2`. |
| Observer region and clock | The original twelve-port carrier is the clock region. Every attempted seam meeting it is a tick, including attempts with zero current. |
| Field access | Two counters begin at the disclosed occupations. They update using the authentic signed currents on incident seams. Subsequent field bits are processing of these records, with no new occupation measurement at each tick. |
| Read protocol | Retain `c=n_0` after the first scheduled incident hit and `d=n_1` after the `k`th scheduled incident hit, with `k>=2`. The success record is `(a,b,c,d)`. |
| Feedback | If enabled and `c=1`, apply one additional native repair on the internal root seam `(0,1)` immediately after the first read. The policy uses the committed read bit. |
| Clock accounting | Feedback is itself one additional incident tick. The final read is at total clock value `k+c` with feedback and `k` without it. The first read precedes the intervention. |
| Continuation | Remote evolution, field-conditioned state and the result of feedback all continue into the next interval. There is no reset or fresh preparation at a read. |
| Finite approximation | The common-ball calculation caps each scheduled waiting interval at `h` attempts and retains every timeout prefix. Positive integer arithmetic supplies enclosing masses; discarded numerical mass is charged. |
| Output | The four-bit success law and all timeout prefixes are a declared coarsening of the richer current/controller history. A complete law here means all outcomes of this declared output, not every source mark or terminal quantum state. |

The source identification is the one checked by the stopped-clock
comparison:

| Quantity | Value |
| --- | ---: |
| Full screen levels | 3 and 5 |
| Full screen ports | 15,360 and 245,760 |
| Full screen seams `M` | 40,320 and 645,120 |
| Root clock region | 12 ports |
| Incident clock seams `m` | 33 |
| Support of incident marks `S=B_1(root)` | 15 ports |
| Common induced ball `Q=B_9(root)` | 222 ports, 561 internal seams |
| Outgoing seams of the common ball | 45 |
| Common-ball clock fraction `p` | `33/561=1/17` |

Retaining incident marks requires the fifteen-port support `S`. The
positive-transfer dictionary below has the narrower condition that its
field projectors and intervention commute with the remote transfer. The
two root-port fields and the internal repair `(0,1)` satisfy both
requirements. Arbitrary interventions elsewhere in `S` need not do so.

## The same positive transfer for a finite-density preparation

On the complete binary configuration space of any one of these finite
graphs, use column distributions and set

```text
A = (1/M) sum_(remote e) Q_e,
B = (1/M) sum_(incident e) Q_e,
D = I-A,                  p = m/M,
K = B D^(-1),
T = D^(-1/2) B D^(-1/2).
```

Each `Q_e` is a nonnegative symmetric orthogonal projection. Consequently
`0<=A< I`, `p I<=D<=I`, `0<=T<=I`, and `T` is symmetric stochastic and
entrywise nonnegative. The completed native clock has transfer `K`, with
the remote repairs before the terminating incident repair. The exact
similarity is `K=D^(1/2) T D^(-1/2)`.

Let `pi_rho(x)=rho^N(x)(1-rho)^(V-N(x))`, where `N(x)` is the total
occupation, and let `Omega_rho(x)=sqrt(pi_rho(x))`. Every native repair
preserves `N`, including each marked branch. Therefore its matrix commutes
with the diagonal matrix `W=diag(Omega_rho)`. Field projectors commute with
`W` as well. For every complete branch product `C`,

```text
1^T C pi_rho = Omega_rho^T C Omega_rho.
```

This follows by writing `pi_rho=W Omega_rho`, moving `W` through the
number-preserving product, and using `1^T W=Omega_rho^T`. The identities
`A Omega_rho=(1-p)Omega_rho`, `B Omega_rho=p Omega_rho` and
`D^s Omega_rho=p^s Omega_rho` also hold.

Remote seams have no endpoint in the root. Root occupation projectors and
the internal repair `(0,1)` thus commute with `D` and its real powers.
Insert the similarity for each `K` in a field cylinder. Adjacent factors
cancel through these projectors and repairs; the remaining boundary
factors cancel against `Omega_rho`. This proves equality of every complete
local field word under the native `K` process and the effective `T`
process, with the same projectors and internal intervention dictionary.
The proof is branchwise and so permits the stated record-dependent choice
to apply or omit the repair. Initial conditioning is represented
by its retained projectors; it does not require a new stationary state
after a read.

If more incident marks are retained, the existing marked construction uses
`T_m=D^(-1/2) B_m D^(-1/2)`, with `sum_m T_m=T`. Its nonnegative branch
matrices have the same complete cylinder probabilities after the same
cancellation. This identifies a marked effective instrument, not the
native process's unobserved attempt times or hidden trajectories.

The transfer `T` is determined by the graph and clock. The density changes
its boundary preparation, not the transfer itself. The argument also
follows by mixing the existing uniform fixed-number-sector identities
with their Bernoulli weights, without restricting to the two-particle
sector.

Strict positivity holds on the entire binary configuration space. For
each triangular face of the root carrier,

```text
sum_(e in triangle) Q_e >= (3/2) I.
```

Fix the spectator bits. On a triangle with zero or three occupied ports,
all three maps are identity. On its one- or two-occupation sector the
projection sum has eigenvalues `3` and `3/2`. Summing the carrier's twenty
triangles, with every one of its thirty edges counted twice, gives
`sum_(root-carrier e) Q_e >= 15 I`. All these seams are incident, so

```text
B >= (15/M) I,        T >= (15/M) I,
H = -log T.
```

The producer checks this premise against its actual root seam menu. It
enumerates the twenty source triangles and verifies the double cover of
all thirty incident carrier seams. On all eight binary triangle states it
also calls the original native endpoint/coin rule, checks identity-plus-swap,
and verifies the exact symmetric count-matrix identity
`(C-3I)(C-6I)=0` for `C=2 sum_e Q_e`. Thus the strict lower floor is not
inferred from a floating spectrum or assumed from the carrier's name.

The certified lower floors are `1/2688`, `1/43008` and `5/187` for levels
three, five and the common ball respectively. Each finite `H` is therefore
well defined, and the existing logarithm series gives, for a floor `alpha`,

```text
0 <= H - sum_(j=1)^r (I-T)^j/j
  <= (1-alpha)^(r+1)/((r+1)alpha) I.
```

The floors shrink with the full screen's seam count. This does not supply
a uniform global operator-norm estimate, a limiting infinite-screen
Hamiltonian, or a coherent native implementation of `exp(-i s H)`.
The comparison established below is of specified finite field laws.

## Exact exclusion dual for the entire four-bit law

[The backward operator](dual.py) works with squarefree monomials
`n_X=product_(i in X) n_i`. A native seam satisfies exactly

```text
Q_e n_X = (n_X + n_(swap_e X))/2.
```

It preserves the cardinality of `X`. The degree-two action is the
two-label exclusion lift: when both endpoints belong to `X`, swapping
does not create independent walks or a repeated port. This is the exact
action on binary configurations.

The observable dual of the completed distribution kernel is

```text
L = K^T = (I-A)^(-1) B.
```

For a cap of `h` attempts in each scheduled interval, it is

```text
L_h = sum_(j=0)^(h-1) A^j B,
r = 1-p,                  q = 1-r^h.
```

The incident action comes first when propagating backwards. Reversing the
order is not a harmless implementation choice, even though a stationary
commuting field dictionary can conceal that reversal in some cylinder
probabilities. Every column and row sum of `L_h` is `q` on each monomial
degree.

For `k` scheduled intervals define the following polynomials, evaluated
on the initial configuration:

```text
v    = L_h^(k-1) n_1,
F    = L_h n_0,
L0   = L_h v,
J0   = L_h (n_0 v),
J1   = L_h (n_0 Q_(0,1) v)     with feedback,
J1   = J0                      without feedback.
```

Here multiplication by `n_0` occurs after propagating through the
intervention backwards. Since `n_0^2=n_0`, the products `n_0 v` and
`n_0 Q_(0,1)v` have degree at most two. Their coefficient masses before
the first interval are `q^(k-1)`; those of `L0,J0,J1` after it are `q^k`.
The coefficient mass of `F` is `q`.

Condition on the initial retained bits `(a,b)`. The expectation of a
squarefree monomial is the product of `a` or `b` for its disclosed ports
and one factor `rho` for each other distinct port. Write this conditional
expectation as `E_ab`. The four successful conditional masses are exactly

```text
Pr(c=1,d=1, success | a,b) = E_ab[J1],
Pr(c=1,d=0, success | a,b) = q^(k-1) E_ab[F] - E_ab[J1],
Pr(c=0,d=1, success | a,b) = E_ab[L0] - E_ab[J0],
Pr(c=0,d=0, success | a,b) = q^k - q^(k-1) E_ab[F]
                                      - E_ab[L0] + E_ab[J0].
```

Multiplying by `pi_ab=rho^(a+b)(1-rho)^(2-a-b)` gives all sixteen joint
success probabilities. Feedback appears only on the `c=1` branches;
replacing both `J0` and `J1` by a single averaged controlled expression
would implement a different observer.

Thus the entire four-bit joint law, including its fourth connected
cumulant, is obtained from exact degree-one and degree-two source actions.
Conditioning accounts for the two initial bits before taking expectations.
This is not a first/second-moment closure of a general multi-time process.
The common-ball calculation uses 222 singleton and 24,531 pair monomials
instead of enumerating its `2^222` configurations.

## Timeouts are output branches

Let `j` count successful scheduled intervals before a timeout, with
`0<=j<k`. Before the first hit the only retained labels are `(a,b)`:

```text
Pr(a,b,timeout after 0) = pi_ab (1-q).
```

For `j>=1`, the observed `c` remains in the timeout record:

```text
Pr(a,b,c=1,timeout after j)
    = pi_ab q^(j-1)(1-q) E_ab[F],
Pr(a,b,c=0,timeout after j)
    = pi_ab q^(j-1)(1-q) (q-E_ab[F]).
```

Subsequent success and timeout probabilities depend only on the scheduler,
not the state or first bit. The native intervention is trace preserving,
so these formulas apply with or without feedback. Individual prefix
probabilities retain their first-read dependence; replacing them merely
by their aggregate `q^j(1-q)` would lose output information.

There are `16+4+8(k-1)` records, hence 44 at `k=4`. Their total mass is one,
the success mass is `q^k`, and the total timeout mass is exactly
`tau=1-q^k`. Later timeout branches retain the additional native feedback
tick when `c=1`. No successful record distribution is renormalized.

## Directed arithmetic and complete error accounting

[The field-law implementation](refinement.py) computes coefficient lower
bounds on a grid `G=2^bits`. Each sparse subset action has integer entries
and row/column sum `2*number_of_subset_seams` before division by `2M`.
Division is always downward. Positive maps preserve componentwise
domination, including through successive capped intervals. The feedback
average of the two degree-one coefficients is also rounded down.

All supplied coefficient masses are at most one. For the retained ball
and `bits=52`, `2M G=1122*2^52<2^63` bounds every integer multiplication
and positive sparse accumulation. Partial completed sums remain at most
`G`. The implementation checks the corresponding bound for its declared
graph before using signed 64-bit arithmetic. No floating normalization,
eigensolver tolerance or unmeasured matrix residual enters this proof.

If a nonnegative exact coefficient vector has known mass `s` and its
computed lower vector is `u/G`, then

```text
delta = s - sum_i u_i/G,
E_ab[u/G] <= E_ab[exact] <= E_ab[u/G] + delta,
```

because each conditioned monomial lies in `[0,1]`. Here `s` is `q` or
`q^k`, including any error inherited by a previous lower vector. Exact
rational interval addition and subtraction then enclose every displayed
success and timeout formula. Negative terms use the opposite endpoint;
known probability bounds may be intersected with the interval.

Let `[ell_z,u_z]` be the intervals for every record of the complete capped
law `P_h`. The construction proves `ell_z<=P_h(z)<=u_z` and

```text
delta_cap = 1 - sum_(all z) ell_z
          = ||P_h - ell||_1.
```

This is a one-sided subprobability error, not a claim that half this L1
quantity bounds every event against an unequal-mass vector. Its exact
timeout total is stored separately from its arithmetic interval widths.

For the completed sixteen-record law `P`, every successful capped lower
mass remains a lower bound. Set

```text
mu = 1 - sum_(successful z) ell_z.
```

Then `mu` includes both unresolved waiting and arithmetic losses. For any
subset of the four retained bits, summing `ell_z` over records where those
bits are all one gives a lower indicator moment `l`; its completed moment
lies in `[l,min(1,l+mu)]`. This charges missing mass once and does not
condition on success. In the mathematical finite-source model, sending `h`
to infinity removes the waiting tail; for each fixed finite `h` and `k`,
an unbounded integer grid makes the directed arithmetic error vanish.
The retained evaluator supports `2<=k<=8`, `h<=2048` and `bits<=52`,
subject to its signed-64-bit bound. It certifies the displayed finite
budgets. Arbitrarily small numerical error would require an arbitrary-precision
integer implementation or another certified evaluator.

Receipt endpoints may be exported by exact outward rounding to a fixed
dyadic grid. The complete deficit must then be recomputed from the
exported lower masses, or explicitly enlarged by the sum of their downward
rounding losses. Likewise the completed missing mass uses the sixteen
exported success lower masses when it certifies that exported vector.
Rounding one aggregate deficit cannot in general cover all
separate row-rounding errors.

## Composing with the actual-screen locality theorem

The supplied Bernoulli preparations have identical restrictions to `Q`.
Both initial field disclosures, every counter update and the internal
record-controlled repair are permitted local operations for the existing
complete-output theorem. The extra repair is counted in the observer's
clock and stays in the local support; it is not an additional random
scheduled hit in the stopped-word estimate.

Let `epsilon_Q(k)` be the existing one-screen first-exit bound. It is
computed directly from actual integer walk counts in the common ball,
using `m=33` and summation through degree 48 with its rational remainder.
For the completed field law on either full screen `G`, postprocessing
contractivity gives

```text
TV(P_G, P_Q) <= epsilon_Q(k),
TV(P_level3, P_level5) <= 2 epsilon_Q(k).
```

These are normalized classical laws. They are a consequence of the
stronger complete quantum-output comparison, not a new diagonal-only
locality assumption. In particular the spatial allowance is added after
the field event is extracted; it is not multiplied by a density or guessed
amplitude factor.

The retained outward pair allowances from that theorem are:

| Scheduled incident hits | Complete pair allowance |
| ---: | ---: |
| 2 | `2.731e-10` |
| 4 | `6.656e-9` |
| 8 | `4.296e-7` |

Thus a common-ball completed indicator interval `[l,u]` becomes
`[max(0,l-epsilon_Q), min(1,u+epsilon_Q)]` for either full screen.
The positive-transfer dictionary identifies these same completed field
laws with their own source-derived finite effective transfers and boundary
vectors. It does not assert that the three global transfer matrices are
equal, or close in operator norm.

The fourth joint cumulant is evaluated from its fifteen set partitions:

```text
kappa(X0,X1,X2,X3)
  = sum_(partitions pi) (-1)^(|pi|-1)(|pi|-1)!
                         product_(C in pi) E[product_(i in C) Xi].
```

All moment intervals lie in `[0,1]`; multiplication is monotone there and
negative partition coefficients reverse the endpoints. The calculation
uses all sixteen success labels and the completed missing-mass budget.
Spatial error may be added to the individual moments before this interval
calculation. Alternatively, normalized laws at TV distance `epsilon`
have fourth-cumulant difference at most `75 epsilon`: each block moment
changes by at most `epsilon`, each product changes by at most
`|pi| epsilon`, and the coefficient sum is `1+14+36+24=75`.

## Executions, resources and remaining obligations

The full-source execution reuses the actual-screen
[current-history runtime](../native_action/clock_histories.py), with both
initial field disclosures made before evolution. It retains the same
first-bit-controlled native repair and reconstructs the second field from
its own incident currents. The probability certificate is deterministic
and does not interpret a small set of seeded histories as a statistical
confidence interval.

The retained executions use seeds zero through seven with feedback both
disabled and enabled, for sixteen paired full-screen executions. Seeds five
and six exercise a feedback repair with signed current `+1` on `(0,1)`,
moving occupation from root port zero to root port one. Every execution
also resumes a serialized intermediate checkpoint and compares the full
states, records and random-generator continuation. These examples check
actual field access and continued work; the analytic law and error bounds
do not depend on their observed frequencies. A separate one-slot run
retains an externally forced timeout.

The serialized observer checkpoint wraps the existing full-source
checkpoint with both initial disclosures, both current counters, and the
first-read memory for each screen. Loading checks these original values
against the source prefix, including exact types and record ancestry.
This source replay is a validation operation, not another observation
available to the controller. After restoration, the observer advances the
saved counters and first-read memory using only the subsequent record
suffix. It neither resets the counters from terminal occupations nor
reprocesses the old prefix as part of this operational continuation.
Independent full-history reconstruction and comparisons with terminal
occupations are audit checks of the result, not hidden-state inputs to
the feedback policy.

The expected number of scheduled attempts for the completed target is
`k M/33` on each full screen and `17k` on the common ball. With feedback,
the expected additional native work is `rho`: the first field has marginal
`rho` before any controlled repair. Thus the expected complete native
work is `k M/33+rho` or `17k+rho`. These are operation counts, not physical
durations. A cap of `h` per common-ball interval permits at most `kh+1`
native attempts including feedback; the computational cost of propagating
24,753 dual coefficients is a separate resource.

The common-ball cap approximates completed-clock laws. It is not a common
cap on the two full screens' global attempts: identity-slot erasure changes
those counts. An externally capped full-source simulation retains its
own timeout and must charge its separately certified cap allowance when
compared with a completed target.

For four scheduled hits, cap 512 and a `2^52` coefficient grid, the
reproducible producer gives the following outward summaries. All quantities
are unconditioned on success.

| Quantity | Baseline | First-bit-controlled repair |
| --- | ---: | ---: |
| Complete capped record count | 44 | 44 |
| Complete capped lower-mass L1 deficit | `<1.209e-9` | `<1.212e-9` |
| Total interval-cap timeout | `<1.324e-13` | `<1.324e-13` |
| Fourth connected cumulant on either full screen | `[0.0008638252,0.0008638837]` | `[-0.0049753855,-0.0049753193]` |

The cumulant intervals include the one-screen complete-output allowance
`epsilon_Q(4)<3.328e-9`, applied at the indicator moments. In particular,
their opposite signs survive the complete numerical and spatial budgets.
The final field mean increases under the native feedback by an amount in

```text
0.0836900722 <= E[n_1(final)]_feedback - E[n_1(final)]_baseline
            <= 0.0836900880
```

on either full screen. This is a difference of the two independently
enclosed protocols; it does not assume feedback preserves the stationary
preparation. The [retained receipt](refinement_results.json) contains all
44 row intervals per protocol, lower-mass deficits recomputed after dyadic
export, completed missing-mass and timeout allowances, joint-four and
connected-four intervals, and the actual full-source execution ledgers.

The independent [small-source oracle](test_refinement.py) propagates all
binary configurations through the original native endpoint/coin rule.
It checks every success and timeout record against the dual intervals at
multiple caps and densities, including coarse grids. Its complete terminal
states distinguish reversed waiting/incident order even where a stationary
local dictionary hides that mistake. Separate controls distinguish removed
remote motion and omitted record-controlled repair. The fourth-cumulant
oracle uses an explicit partition formula independent of the producer's
partition generator.

Reproduction uses the same producer for generation and checked replay:

```text
python -m model.observer_time_fields.refinement --check model/observer_time_fields/refinement_results.json
python -m pytest model/observer_time_fields/test_dual.py model/observer_time_fields/test_refinement.py
```

The proof supplies a finite source-relative clock and effective-field
comparison for two actual screens with a compatible, non-diluting
preparation. It leaves an all-level refinement limit, physical preparation
and detector error, a native coherent realization of the effective
Hamiltonian, and dictionaries not commuting with the remote transfer as
explicit further obligations. A movement-only clock retains its negative
eigenvalue obstruction; changing the clock is not a free application of
the positive-transfer theorem.

## Canonical equilibrium without product preparation

The [canonical-sector construction](canonical_sector.py) removes independent
Bernoulli occupations as a mandatory preparation for this clock/field
connection. It reuses the fixed-number native equilibrium and the same
exclusion dual; its additional result is a complete finite-density field law
with **exact conservation correlations**, attached separately to each actual
screen. Nothing is factorized after evolution or conditioned on successful
completion.

For a connected binary source with `V` ports and conserved occupation `N`,
every native seam map is the symmetric stochastic projection
`(I+SWAP)/2`. Adjacent swaps on a connected graph connect every pair of
binary configurations of total `N`. Every configuration also has a positive
self-loop. Hence the diagonal process is irreducible and aperiodic in that
sector, with unique stationary law

```text
pi_(V,N)(x) = 1/binomial(V,N),     sum_v x_v = N.
```

This is a finite-chain proof, not an assumption of initially independent
ports. The producer calls the original two-port rule to check identity plus
swap and verifies connectivity of all 15,360 and 245,760 source ports.
The empty and full sectors are the trivial single-state cases. This argument
concerns the recorded diagonal process; it does not assert a new general
coherent-state uniqueness theorem.

Repair preserves `N`, so it cannot select the sector or density. The displayed
non-diluting examples use `N=5120` and `81920`, respectively, giving density
`1/3` on each screen. The equilibrium law is derived conditional on that
conserved boundary. Choosing to prepare it is an operational service.
An exact sampler uses successive uniform integer choices: with `r` positions
remaining and `n` occupations left, the next bit is one when a uniform
integer from `0,...,r-1` is below `n`. Every final `N`-subset has probability
`1/binomial(V,N)` under the stated random-draw service. No finite relaxation
time or exact equilibrium after an arbitrary short native run is claimed.

The restriction to the retained ball `Q` is computed, not replaced by a
fixed number inside that ball:

```text
pi_Q(x) = binomial(V-|Q|, N-sum(x)) / binomial(V,N).
```

It is uniform conditional on its local occupation total, but that total
fluctuates. Its weights depend only on local total number, so every retained
native seam preserves it. The positive-transfer argument above applies
with `Omega_Q=sqrt(pi_Q)` and with the full-screen uniform sector boundary.
Number preservation supplies the required boundary commutation. Independent
product structure is unnecessary.

Let `(z)_r=z(z-1)...(z-r+1)` denote a falling factorial. For initial disclosed
bits `(a,b)` with nonzero probability, the exact preparation identities are

```text
Pr(a,b) = (N)_(a+b) (V-N)_(2-a-b) / (V)_2,
E[n_X | a,b] = (factors a,b present in X)
              * (N-a-b)_r / (V-2)_r,
```

where `r` counts the distinct undisclosed ports in `X`. Replacing the product
expectation in the existing dual by this identity evaluates the same complete
44-record protocol. Known coefficient mass bounds each expectation
because every monomial lies in `[0,1]`. All directed rounding, timeout,
cumulant and outward export allowances therefore remain valid. Impossible
initial disclosures retain explicit zero-mass branches.

The two finite screens have **different** restrictions, despite equal
density. Distinct initially disclosed ports have covariance
`-rho(1-rho)/(V-1)`: it is exactly `-2/138231` at level three and
`-2/2211831` at level five. The producer evaluates each restriction
separately. The complete-output locality theorem compares each full screen
to its own exactly restricted common-ball preparation, with the same
one-screen allowance below `3.328e-9`. It does not assume that the two
restrictions coincide.

The [retained results](canonical_sector_results.json) certify:

| Full screen | Baseline connected fourth cumulant | Feedback connected fourth cumulant | Feedback increase of final mean |
| --- | --- | --- | --- |
| Level 3 | `[0.0008630090,0.0008630675]` | `[-0.0049747133,-0.0049746471]` | `[0.0836955211,0.0836955370]` |
| Level 5 | `[0.0008637742,0.0008638327]` | `[-0.0049753435,-0.0049752773]` | `[0.0836904127,0.0836904286]` |

These intervals include numerical, stopping and spatial allowances. The
receipt also bounds total variation between the two normalized completed
sixteen-word field laws. For each word, its capped interval enlarged by
the timeout mass encloses its eventual probability. Interval differences
bound their complete L1 distance; adding the two separate spatial allowances
then bounds the full-screen distance. Finite-number corrections are retained
in this comparison rather than absorbed into an assumed identical state.
The resulting TV interval is `[0.00003106953,0.00003108527]` without feedback
and `[0.00003102585,0.00003104159]` with it; both resolve a nonzero finite-size
effect after the locality allowance.

Eight paired seeded executions initialize the actual full screens with their
exact totals and run the unchanged native repair engine, including remote
repairs and record-controlled feedback. The initial screen preparations are
separate; agreement of their local initial configurations is not assumed.
Both current counters, first-read memory, full states, record histories and
random-generator continuation survive serialized restart. Moving and
zero-current feedback occur, and a forced one-slot timeout remains an output.
These seeded pseudorandom executions test the implementation; their small
seed ensemble is not statistical evidence for a uniform distribution over
all configurations. The analytic law uses the specified uniform-draw model.

Independent small-source tests enumerate every fixed-number configuration
and every native branch. They check all records on both full small graphs
and restrictions of larger sectors, including empty/full and coarse-grid
cases. Controls reject replacement by a product law, replacement of the full
sector by a fixed number in the retained ball, and a disconnected graph's
claim to a unique global-sector equilibrium. This removes a preparation
independence premise, while leaving the conserved density, preparation
service, coherent continuation and physical identification explicit.

```text
python -m model.observer_time_fields.canonical_sector --check model/observer_time_fields/canonical_sector_results.json
python -m pytest -q -W error model/observer_time_fields/test_canonical_sector.py
```

## Native temporal correlations and the spatial Green field

[The temporal implementation](spatial_dynamics.py) connects the actual
[source Green construction](../physical_space_time/OCCUPATION_FIELD.md#controlled-green-elimination-on-the-actual-source-screens)
to correlations at separated native attempts. It uses the same binary
repair, uniform IID seam selection and stationary canonical number sector.
For collective fields the observer admits the initial occupation disclosure
and complete seam/current tape; these reconstruct the occupations through
actual repairs. A single primitive port tape does not supply all these reads.
The sector can have any `0<q<V`; the continuum field uses `q=V/3`, while
the exact elementary oracle executes all 3,960 branches of `V=12,q=2`.

Let `L_l` be the full port Laplacian, `M_l=630*4^l` its seam count,
`s=4^l`, `P0` remove the conserved constant, and
`chi=q(V-q)/(V(V-1))`. For `eta=x-q/V`, native repair gives exactly

```text
E[eta_(n+1) | eta_n] = K_l eta_n,     K_l=I-L_l/(2M_l),
Cov(eta_n,eta_0) = chi K_l^n P0.
```

The first identity follows on every seam from the mean of the two native
endpoint outcomes. The second uses the exact canonical covariance
`E eta_0 eta_0*=chi P0`; it is not an independent-particle assumption.
These are two-point laws of bounded native occupations, not a claim that
their instantaneous joint distribution is Gaussian. Initial occupation
disclosures used for readback retain every label and are not conditioned
away in these expectations.

Set `N_l=2M_l s`. For any rational discount rate `z>0`, the geometric series
is an exact processing of correlations at integer native attempts:

```text
r=N_l/(N_l+z),
(1/N_l) sum_(n>=0) r^(n+1) K_l^n = (s L_l+zI)^-1.       (T)
```

The omitted tail after `h` terms has norm at most `r^h/z`. No extra Poisson
scheduler or fractional native event is needed. The constant mode is
`1/z` on both sides and disappears for centered fields. Its mean-free
`z->0` limit is `(1/s)L_l^+`, exactly the Green operator in the spatial
fluctuation construction. Thus that static covariance and (T) come from
one repair law and readout dictionary. The observation-count Brownian
limit concerns cumulative fields on long histories; its fluctuation-time
parameter is different from the native attempt separation in (T).

### A uniform temporal bound for eliminating carrier interiors

Let `U` repeat a cell value over its twelve ports with factor `1/sqrt(12)`.
Decompose the actual full Laplacian between cell constants and their
orthogonal complement:

```text
L_l = [[A,C],[C*,D]],       S_l=A-CD^-1C*.
```

The internal native gap is greater than `lambda0=2`. The cross-carrier
seams form a port matching, so their Laplacian has norm at most two.
The existing Schur comparison gives
`S_l >= lambda0/(lambda0+2) A`. Consequently

```text
CD^-2C* <= (A-S_l)/lambda0 <= 2S_l/lambda0^2.
```

The Schur formula at positive frequency is

```text
U*(sL_l+zI)^-1U = (sS_l+zI+E_z)^-1,
E_z=sC[D^-1-(D+z/s)^-1]C*,
0 <= E_z <= [2z/(lambda0^2 s)] sS_l.
```

For `R=(sS_l+zI)^-1`, positivity and the resolvent identity yield
`0<=R-(R^-1+E_z)^-1<=R E_z R`. Since
`max_(a>=0) a/(a+z)^2=1/(4z)`, this proves at every source level

```text
0 <= (sS_l+zI)^-1-U*(sL_l+zI)^-1U,
||(sS_l+zI)^-1-U*(sL_l+zI)^-1U|| <= 1/(8s).             (R)
```

The constant mode cancels in the difference. The result holds for every
admitted port matching, independently of how geometric directions label
its ports, and survives restriction to the fixed ancestry observations.
The executable verifies exact full 240-port residuals for `z=1/2,1,4`;
independent full-matrix solves and changed valid port routings test (R).

Cell means do not form a closed Markov process at finite level. Their
two-attempt discrepancy from resetting to cell means is
`CC*/(2M_l)^2`; the base cell-zero/minus-cell-one contrast gives exactly
`13/19051200`. Weak convergence of compressed static Green operators
alone does not establish a limiting resolvent identity. The
[source compactness theorem](../physical_space_time/OCCUPATION_FIELD.md#a-common-positive-spatial-generator-from-the-native-green-field)
supplies the stronger convergence used below; it does not make the finite
cell process Markovian.

### Two continuations of the same native density spectrum

On the scaled low band `a=s lambda<=Lambda`, the positive attempt transfer
has `H_density=-log(I-L_l/(2M_l))`. Its rescaled eigenvalues obey

```text
0 <= N_l H_density-a <= Lambda^2/[2(N_l-Lambda)].        (L)
```

This is a bound on the declared band, not on the full growing spectrum.
For the root-carrier incident clock, write `m=33`, let `L_i` contain its
incident seams, and let `L_r=L_l-L_i`. The density transfer has positive
similarity representative
`T=(2mI+L_r)^-1/2(2mI-L_i)(2mI+L_r)^-1/2`.
Generalized min-max gives ordered eigenvalue bounds

```text
lambda/(2m+lambda) <= eigenvalue_j(I-T) <= lambda/(2m),
|2ms eigenvalue_j(-log T)-s lambda_j|
    <= Lambda^2/[2(2ms-Lambda)].                        (I)
```

Only the low band with `Lambda<2ms` is used. This is an ordered spectral
comparison, not a scalar conversion of every native field or marked clock
history. Complete incident laws use the earlier noncommuting transfer
dictionary and retain their waiting and readout correlations.

There is also a positive wave construction. Reuse the
[coherent native-arc continuation](../observer/COHERENT_CONTINUATION.md):
coherently prepare the uniform seam/coin alternatives with isometry `A0`,
reverse each native current arc with involution `J`, and reflect about the
prepared subspace. For `W=J(2A0 A0*-I)`, its exact compression is

```text
A0* W^n A0 = T_n(K_configuration),
K_configuration G = G(I-L_l/(2M_l)).
```

Here `T_n` is a Chebyshev polynomial and `G[x,p]=x_p-q/V`. These identities
hold on each fixed binary number sector: native branch reversal is a
bijection on every seam, and each density column obeys the preceding
endpoint identity. They do not replace exclusion by independent particles.
Therefore the density-mode compressed amplitudes are
`cos(n arccos(1-lambda/(2M_l)))`. Put `Q_l=M_l s` and define the coherent
step parameter `t=n/sqrt(Q_l)`. For `a=s lambda<=Lambda<4Q_l`,

```text
Omega_l(a)=sqrt(Q_l) arccos(1-a/(2Q_l)),
0 <= Omega_l(a)/sqrt(a)-1 <= Lambda/[6(4Q_l-Lambda)] = e_l,
||cos(n theta)-cos(t sqrt(sL_l))||_band <= t sqrt(Lambda) e_l.   (W)
```

The zero mode is obtained by continuity. To prove the bound, set
`x=sqrt(a/(4Q_l))`. Then `arccos(1-2x^2)=2 asin(x)` and
`asin(x)-x<=x^3/[6(1-x^2)]`, by integrating
`(1-u^2)^-1/2-1<=u^2/[2(1-x^2)]`. Thus (W) gives an explicit
square-root spatial dispersion for an admitted continuation of the same
native density transfer. The finite oracle executes all 3,960 native arcs
through eight coherent steps and retains both Ramsey outcomes.

These are coherent amplitudes and an apparatus clock, not the passive
occupation process. Coherent access to uncommitted records, their reversible
re-encoding, the reflection, preparation of the stated amplitude and
controlled Ramsey readout are required services. Opposite reflection gives
the same first coherent native event but a different second Ramsey result.
A phase on a complete committed current label preserves each CP instrument
branch while changing the coherent transfer's third spectral moment by
`-1/900`. The committed native instrument therefore does not force this
continuation. No native pulse compiler for the growing walk is asserted.
The amplitude proportional to a centered occupation is a declared coherent
preparation; it is not the stationary classical mixture used for (T).
The common source operator and field intertwiner connect these protocols,
without identifying their preparation or measurement contracts.

Its relation to the density log is explicitly nonlinear:
`theta=arccos(exp(-H_density))`. The
[active three-class modular clock](../transfer_modular_probe/README.md)
uses a different observable subspace of the same elementary source; its
two nontrivial eigenvalues `(25+/-sqrt(3))/30` do not equal these density
bands. Its compiled gate does not by itself implement (W), supply the
growing field generator, or establish a cross-sector clock calibration.

Finally, a shared scalar clock conversion cannot turn the density log
into the wave dispersion. On the elementary carrier its three Laplacian
bands are `5-sqrt(5),6,5+sqrt(5)`, with log rates
`h(lambda)=-log(1-lambda/60)`. An arbitrary mass and speed would require
`h(lambda)^2` to be affine in `lambda` after a common positive scale.
This function is strictly convex. Rational logarithm enclosures give a
three-band affine-fit defect in `[0.910524428210,0.910524428211]` after
normalizing by the first squared frequency; any affine fit has maximum
normalized squared-frequency error at least `0.329430832885`.
This rejects that specific identification while retaining the constructive
coherent wave route (W). The time-average Green covariance `a^-1` is also
different from a free scalar vacuum covariance `(2sqrt(a))^-1`.
Neither spectral construction alone supplies a physical quantum field
state, interactions, Lorentz locality, conformal invariance or measured
spacetime. Those require the corresponding common observable dictionary
and controlled regime, beyond the finite and bandwise results here.

```text
python -m model.observer_time_fields.spatial_dynamics --check model/observer_time_fields/spatial_dynamics_results.json
python -m pytest -q -W error model/observer_time_fields/test_spatial_dynamics.py
```

## One continuum operator for fields, native decay and coherent waves

[The continuum dynamics implementation](continuum_dynamics.py) uses the
actual source compactness result, the same ancestry readouts and the
native resolvent identity (R). No continuum Laplacian is substituted for
the finite repair operator. Write `H=L2_0(mu)`, where `mu` is the normalized
flat-polyhedron area transported to the carrier sphere. Let `J_l` identify
cell coefficients with normalized step functions, `P_l=J_l J_l*`, and

```text
B_l=4^l S_l,   G_l=J_l B_l^+ J_l*,   D_l=4^l L_l.
```

The coarse-patch Poincare estimate gives compact Green images and upgrades
the existing common subsequence to `||G_l-G||->0`. Smooth-test variational
bounds make `G` injective on `H`; it is positive, compact and bounded by
`C=4096/3`. Its inverse `A=G^-1`, with domain `range(G)`, is therefore a
densely defined positive self-adjoint operator with compact resolvent.
This defines an actual continuum spatial action
`E(f)=||sqrt(A)f||^2` on `D(sqrt(A))`, independently of a guessed metric.
The [full-port recovery argument](../physical_space_time/OCCUPATION_FIELD.md#locality-and-spectral-dimension-of-the-source-continuum)
shows that adjoining the constant zero mode gives a conservative strongly
local Dirichlet form; its heat-trace spectral dimension is two.
The selected subsequence is common to every readout below; neither a
whole-sequence limit nor its numerical value is inferred from finite runs.

### The full native dynamics have the same limit

For every `z>0`, functional calculus of the bounded Green operators gives

```text
J_l(B_l+z)^-1 J_l* = G_l(I+zG_l)^-1 -> (A+z)^-1
```

in operator norm on the mean-zero space. The identity acts on the
level-`l` cell subspace; the left side is zero on its orthogonal complement.
Bound (R) shows that replacing the left side by
`J_l U_l*(D_l+z)^-1 U_l J_l*` changes it by at most `1/(8*4^l)`.
Consequently the **full native port dynamics**, including their carrier
interiors, have the same limiting resolvent. Finite elimination memory
has not been discarded as an exact identity.

The functional-calculus step requires control of high spectral mass.
For any fixed `f in H`, the positive full-native spectral measure of
`U_l J_l* f` has total mass `||P_l f||^2->||f||^2`. Its Stieltjes transform
converges by the preceding resolvent result to that of `A` at `f`.
The limiting measure has the same total mass because `A` is densely
defined. More explicitly, for any `R,z>0`,

```text
mass_l([R,infinity)) <= (R+z)/R
    [ ||P_l f||^2-z <f,J_l U_l*(D_l+z)^-1 U_l J_l* f> ].
```

Take the level limit, then `z` large, then `R` larger than `z`. This makes
the bound arbitrarily small. Stieltjes-transform uniqueness and this
tightness give convergence for bounded continuous spectral functions.
Polarization gives weak vector convergence of compressed propagators;
the inequality
`||J_l U_l* F(D_l) U_l J_l* f||^2
 <= <U_l J_l* f,|F|^2(D_l) U_l J_l* f>`
gives the matching norm upper bound and hence strong convergence.
On a compact spectral interval the heat and cosine functions are
uniformly continuous in their parameters; the controlled tails extend
this to compact parameter intervals. Thus, for each fixed input `f`,

```text
J_l U_l* exp(-tau D_l) U_l J_l* f -> exp(-tau A) f,
J_l U_l* cos(t sqrt(D_l)) U_l J_l* f -> cos(t sqrt(A)) f,
```

uniformly for `0<=tau,t<=T`. This is strong convergence on fixed inputs,
not operator-norm convergence of unbounded-frequency wave propagators.

Both functions are limits of the **actual integer-step protocols**. Put
`N_l=1260*16^l`, `Q_l=630*16^l`, choose native attempt count
`floor(tau N_l)` and coherent step count `floor(t sqrt(Q_l))`, and retain
all native outcomes or both Ramsey outcomes as appropriate. On a band
`0<=a<=Lambda<N_l`, the respective errors against the preceding finite
propagators are at most

```text
T Lambda^2/[2(N_l-Lambda)] + Lambda/N_l,
T sqrt(Lambda) Lambda/[6(4Q_l-Lambda)] + sqrt(Lambda)/(25*4^l).
```

The last term includes integer coherent-step rounding, using
`sqrt(630)>25`. For a unit input with `D_l` energy at most `E`, each
remaining spectral tail contributes at most `2 sqrt(E/Lambda)`. For
arbitrary fixed continuum inputs, the preceding tightness argument
provides the qualitative tail control. These are finite-step error
budgets, separate from the unquantified spatial subsequence limit.
Native occupation correlations decay with `exp(-tau A)`; the admitted
coherent amplitudes propagate with `cos(t sqrt(A))`. Their parameter
normalizations differ and are not calibrated laboratory seconds.

### Actual native density observations at several continuum times

The same spatial subsequence also gives the joint laws of **actual native
occupation snapshots**, with attempt count and source size tending to
infinity together. Put `V_l=240*4^l`, `q_l=V_l/3`,
`chi_l=q_l(V_l-q_l)/[V_l(V_l-1)]` and `a_l=1260*16^l=2M_l 4^l`.
Start in the declared uniform fixed-number equilibrium and retain the
uniform IID seam schedule and fair repair. For a bounded real mean-zero
`f in L2(mu)`, let `f_l` repeat its cell averages on the twelve ports of
each cell. Define

```text
Y_l(tau;f) = sum_p f_l(p)[n_floor(a_l tau)(p)-1/3]/sqrt(V_l chi_l).
```

For `tau>=0`, every finite collection of these variables converges jointly
to a centered stationary Gaussian family with

```text
Cov(Y(tau;f),Y(u;g)) = <f, exp(-|tau-u| A) g>.                (OU)
```

The argument establishes Gaussianity of the joint process; convergence of
its two-point function alone would not suffice. Each unread native repair
has exactly the law of an independent fair identity-or-endpoint-swap.
Consequently its occupation path can be represented as `n_j=Pi_j n_0`,
where `Pi_j` is a product of independent seam permutations, independent of
the uniform initial `q_l`-subset. This is an equality of the existing
occupation path laws, not an added carrier process or preparation. Repair
marks themselves are not identified with the auxiliary swap coins.

Here is the needed bound on the randomness remaining in those permutations.
On any finite source with `V` ports and `M` seams, let
`K=I-L/(2M)` be the native one-particle kernel. It is a positive contraction,
being an average of `(I+swap)/2`. For deterministic real profiles `f,g`,
with `||g||_infinity<=C`, independent stirring gives, uniformly in `n`,

```text
E[V^-1 <f,Pi_n g>] = V^-1 <f,K^n g>,
Var(V^-1 <f,Pi_n g>) <= 4 C^2 ||f||_ell2^2/V^2.             (OU1)
```

To prove (OU1), use the Doob martingale
`M_j=V^-1 <K^(n-j)f,Pi_j g>` and set `h=K^(n-j)f`.
Its conditional increment variance is at most

```text
(2 C^2)/(M V^2) sum_edges (h_a-h_b)^2
    = (4 C^2)/V^2 <h,(I-K)h>.
```

Indeed an actual swap changes the pairing by
`(h_a-h_b)[(Pi_(j-1)g)_b-(Pi_(j-1)g)_a]/V`, and occurs with
probability `1/(2M)` per seam. Positivity of `K` gives
`<K^k f,(I-K)K^k f> <= ||K^k f||^2-||K^(k+1)f||^2`.
Summation telescopes, proving the bound independently of the number of
attempts. No temporal mixing estimate is used here.

Condition on all stirring permutations at finitely many snapshot times.
A Cramer-Wold combination is a canonical subset sample with coefficients
`b_l=sum_i c_i Pi_(n_i)^T f_i,l`, whose sum is zero and whose supremum is
bounded independently of `l`. Its conditional variance after the displayed
normalization is exactly `V_l^-1 ||b_l||^2`. For each pair of times, the
corresponding overlap depends on the intervening permutation product.
(OU1) makes that overlap converge in probability to
`V_l^-1 <f_i,l,K_l^|n_i-n_j| f_j,l>`. The full-native heat convergence
above identifies its limit as the entry in (OU), including the integer
rounding of both attempt counts. Thus the entire conditional covariance
matrix converges to the deterministic matrix in (OU).

For a positive limiting Cramer-Wold variance, the finite-population CLT
applies: bounded coefficients and `q_l=V_l/3` give
`max_p b_l(p)^2/[v_l min(q_l,V_l-q_l)]->0`, where
`v_l=||b_l||^2/(V_l-1)`. This is
[Hájek's criterion, stated as Li–Ding Theorem 1](https://arxiv.org/pdf/1610.04821).
For random conditional populations, every subsequence has a further
subsequence on which its covariance converges almost surely. Apply the
deterministic CLT there; bounded conditional characteristic functions
then converge in expectation. This proves the unconditional conclusion
for the original sequence. If the limiting variance is zero, its uniformly
bounded conditional variance converges to zero in expectation, and
Chebyshev supplies the degenerate Gaussian limit. Cramer-Wold proves the
joint claim. The exact identity
`E|Y_l(tau;f)-Y_l(tau;g)|^2=||P_l(f-g)||_L2(mu)^2`
extends the result by bounded-test approximation to all mean-zero `L2`
tests. The path theorem below supplies continuous distribution-valued
versions and path convergence. No spatial convergence rate is asserted.
The process is stationary Ornstein–Uhlenbeck in scaled global-attempt units.

The executable checks the native identity-or-swap dictionary on every
total-four carrier branch, all 3,600 two-attempt stirring branches and
their conditional canonical variances by independently summing all 495
initial configurations. The conditional variances fluctuate at
finite size. Reusing one randomly chosen seam at the second attempt changes
the mean overlap from `K^2` to `K`, although the uniform one-time law is
preserved. A native snapshot's fourth moment also differs from the Gaussian
value. These controls keep the IID schedule and the actual conditional CLT
in the proof.

Since the gap bounds make `A^-1/2` bounded on the mean-zero space, the
declared spatial readout `Phi(tau;f)=Y(tau;A^-1/2 f)` has covariance
`<f,G exp(-|tau-u|A)g>`. This connects the spatial Green field to actual
native density dynamics through a declared spectral filter. Finite readouts
can use the cell coefficients of `G_l^1/2 P_l f`: norm convergence of
`G_l` gives convergence of these tests to `A^-1/2 f`, and the exact `L2`
covariance identity above controls their approximation. Thus no uniform
supremum bound on these filtered tests is needed. The unfiltered
equal-time density covariance is the identity; neither field has the
cosine time covariance of the separately synthesized wave ensemble.

There is also a direct Euclidean positivity statement. For positive
`s,t`, the reflected covariance is
`<f,exp(-(s+t)A)g>=<exp(-sA)f,exp(-tA)g>`.
Stationarity supplies the two-sided extension used in this reflection.
The Gaussian family has conditionally independent past and future given
its time-zero field: subtracting their conditional means leaves zero
cross covariance. For each bounded positive-time cylinder functional `F`,
this gives `E[conjugate(F after reflection) F]=E|E[F|Y(0)]|^2>=0`.
On the one-particle reconstructed Hilbert space, a positive time shift has
semigroup `exp(-tau A)`, with generator **A**, not `sqrt(A)`.
This connects the actual density law to a reflection-positive Gaussian
reconstruction. It supplies no coherent preparation in that reconstructed
state, physical relativistic vacuum, Lorentz symmetry or laboratory clock.

The Green-filtered normalization `phi=Phi/sqrt(2)` makes this reconstruction
explicit. In a mode of `A` with eigenvalue `lambda>0`, its covariance is
`exp(-lambda*|tau-u|)/(2lambda)`. Its temporal Fourier transform is
`1/(omega^2+lambda^2)`. Thus the finite-mode Gaussian precision, and the
corresponding continuum Cameron–Martin quadratic form, are

```text
S_E(phi) = (1/2) integral [||partial_tau phi||^2 + ||A phi||^2] d tau.
```

The expression specifies the form on its finite-energy domain; it is not
a density relative to an infinite-dimensional Lebesgue measure. Gaussian
Wick decomposition gives the symmetric Fock reconstruction, vacuum energy
zero and Hamiltonian `dGamma(A)`. The same statement follows mode by mode
from the oscillator covariance, and extends by its orthogonal chaos
decomposition. Its spatial precision is `A^2`, not `A`: the source form
is equivalent to a second-order elliptic energy, so this is the associated
spatial precision of spectral order four. No smooth bi-Laplacian or
`H2` operator-domain identification is required. A single
rescaling of the clock cannot turn all unequal frequencies `lambda` into
`sqrt(lambda)`. This is a free field reconstructed from the declared
classical native readout, with diffusive dispersion; an observed
relativistic quantum field or a native preparation of its vacuum does not
follow.

### Continuous paths of the actual native density and Green field

Along the same spatial subsequence, linearly interpolate actual snapshots
at the attempt times `k/a_l`. For every finite `T`, jointly with the
finite Green-filtered snapshot readout above, they converge in

```text
Y_l -> Y       in C([0,T],H^-3(S2)),
Phi_l -> Phi   in C([0,T],H^-2(S2)).                         (OP)
```

These are simultaneous source-refinement/attempt-time limits, with
`a_l=1260*16^l`. The stationary preparation, IID seams, fair repair and
full snapshot reconstruction are the same as in (OU). They are not the
count-first long-block construction of the separate cosine wave.

Here is a discrete proof of tightness that keeps the stopping-time and
reversal hypotheses explicit. Write `s=4^l`, `P=I-L/(2M)`, `a=2Ms`,
`F_f(x)=sum_p f_p(x_p-1/3)/sqrt(V chi)` and
`E_l(f)=s<f,Lf>/V` for a real mean-zero full-port profile.
The canonical variance is `||f||^2/V`. Set
`g=2(I+P)^-1 f`, so `E_l(g)<=4 E_l(f)`, and fix the finite reversal
horizon `N=ceil(aT)`, `T_l=N/a`. Define

```text
m_j = F_g(X_j)-F_Pg(X_(j-1)),
r_j = F_g(X_(j-1))-F_Pg(X_j),
M_k = sum_(j=1)^k m_j,       Mhat_r = sum_(j=N-r+1)^N r_j.
2[F_f(X_k)-F_f(X_0)] = M_k-[Mhat_N-Mhat_(N-k)].              (OP1)
```

The first process is a martingale; stationarity and reversibility make
the second a martingale in the reversed filtration. The last identity is
pathwise: subtract the increments and telescope `(I+P)g=2f`. Using `f`
in place of `g` leaves a nonzero discrete endpoint term. No independence
of the two martingales is used.

A swap on seam `{b,c}` changes `F_g` by
`(g_b-g_c)(x_c-x_b)/sqrt(V chi)`. Thus, **at every current state**,

```text
E[m_j^2 | X_(j-1)] <= E_l(g)/(a chi) <=4 E_l(f)/(a chi),
max_j |m_j|^2 <=4 E_l(g)/(s chi) <=16 E_l(f)/(s chi).        (OP2)
```

For the first inequality bound variance by the second moment of the raw
increment and sum the seam swaps, each of probability `1/(2M)`. For the
second, a raw jump squared is at most `E_l(g)/(s chi)`; centering increases
its amplitude by at most two. The same bounds hold in reverse. Hence a
martingale increment over a bounded stopping time and a further interval
of length `delta` has second moment at most
`4 E_l(f)(delta+1/a)/chi`. For uniformly bounded-energy test profiles,
this tends to zero as `delta->0` and `l->infinity`, and all jumps vanish.
Doob's inequality supplies stochastic boundedness. The
[stopping-time tightness criterion and continuous-limit criterion in
Whitt, Theorem 3.5 and Lemma 3.10](https://arxiv.org/pdf/0712.1929)
therefore give continuous subsequential limits for both martingales.
Their joint tightness, reversal on `T_l->T`, and (OP1) give tightness of
every finite collection of interpolated field modes. Deterministic-time
increment variances alone would not justify the stopping-time step.

There is a uniform supremum estimate. Stationarity gives exactly
`E m_j^2=<g,(I-P^2)g>/V=4<f,(I-P)(I+P)^-1 f>/V`
`<=4 E_l(f)/a`. If
`sigma^2=E M_N^2<=4T_l E_l(f)`, Doob bounds the `L2` norm of the forward
maximum by `2sigma` and of the reversed terminal-minus-suffix maximum by
`3sigma`. Equation (OP1), followed by the initial variance bound, gives

```text
E sup_(t<=T) |Y_l(t;f)|^2 <=2||f||^2/V+50T_l E_l(f).       (OP3)
```

To control infinitely many modes, use the real spherical harmonics
`Y_jm` normalized for spherical probability area, with weights
`w_j=1+j(j+1)`. Their addition identities are
`sum_m Y_jm(x)^2=2j+1` and
`sum_m |grad Y_jm(x)|^2=(2j+1)j(j+1)`. Let `a0` be the base polyhedron's
side and `r` its inradius. The radial map has derivative norm at most
`1/r`, and any two points in adjacent fine triangles can be connected
within their union by length at most `2a0*2^-l`. Jensen's inequality for
the cell averages, summed over `30*4^l` cross seams, therefore gives

```text
sum_m Var Y_l(0;Y_jm) <=2j+1,
sum_m E_l((Y_jm)_l) <=[a0^2/(2r^2)](2j+1)j(j+1)
                    <=2(2j+1)j(j+1).                     (OP4)
```

Here `a0^2<2`, `r^2>9/16` prove the final conservative constant; internal
seam differences vanish. Subtracting each `mu` mean only decreases its
variance and leaves energy unchanged. Conservation removes degree zero.
The exact telescoper
`(2j+1)/w_j^2<=1/j^2-1/(j+1)^2`, together with (OP3), implies

```text
E sup_(t<=T) ||Pi_>J Y_l(t)||_-3^2
    <= (2+100T_l)/(J+1)^2.                                (OP5)
```

For the actual Green-filtered test `h_l=G_l^1/2 P_l h`, the existing source
comparison `L_cell/24<=S_l<=L_cell/12` gives
`||h_l||^2<=C||h||^2` and, for its cell-constant full-port lift,
`E_l(h_l)<=2<h_l,B_l h_l><=2||h||^2`, where `C=4096/3`.
Consequently (OP2) applies without a supremum assumption on the filtered
tests, and (OP3) gives

```text
E sup_(t<=T) ||Pi_>J Phi_l(t)||_-2^2
    <= (2C+100T_l)/(J+1)^2.                                (OP6)
```

Finite-mode tightness and these uniform supremum tails prove compactness
in the stated continuous-path spaces: choose a sequence of harmonic
cutoffs with summably small tail probabilities and compact coordinate
path sets; uniform tails make their intersection relatively compact.
The joint laws (OU) and their filtered versions identify every limit,
proving (OP). The piecewise-constant snapshot paths have the same limit
in the Skorohod `J1` topology. Indeed every raw snapshot jump has squared
`H^-3` norm at most `4/(V chi)`, by the addition identity and Jensen;
the filtered `H^-2` bound is `2/(s chi)`, by its aggregate energy bound.
Both vanish, so the interpolation error tends to zero deterministically.

The executable checks (OP1), zero conditional means, the state-uniform
bracket and stationary variance on all 29,700 native carrier branches.
It also retains the nonzero uncorrected-profile defect and a fixed-initial-
state control for which the stationary reverse conditional mean fails.
These path results do not identify global-attempt time with proper time,
change the supplied spatial topology, or turn the diffusive Gaussian
reconstruction into a Lorentzian vacuum.

Changing only the native count scale does not produce a propagating field.
At the ballistic attempt scale `b_l=2M_l*2^l=1260*8^l`, define
`Y_l^bal(t)=Y_l(t*2^-l)` and likewise for the filtered field, using the
same linear interpolation. Before adding its initial variance, (OP1)
and Doob give `E sup_(k<=N)|F_f(X_k)-F_f(X_0)|^2<=25(N/a_l)E_l(f)`.
Set `N=ceil(b_l T)` and sum the same harmonic energy bounds as in
(OP5)–(OP6). Then

```text
E sup_(t<=T)||Y_l^bal(t)-Y_l(0)||_-3^2
 <=50(T*2^-l+1/a_l),
E sup_(t<=T)||Phi_l^bal(t)-Phi_l(0)||_-2^2
 <=50(T*2^-l+1/a_l).                                      (OP7)
```

Both native paths freeze to their initial random field in this scaling.
The separate coherent-walk step parameter is not a rescaling of these
passive attempts. This diagnoses the proposed change of count units for
this density sector; it does not exclude a different microscopic causal
theory whose effective long-scale density is diffusive.

### Finite propagation on the source comparison space

The same full-port convergence also supplies a support bound. Restore
the constant zero mode and write `Abar=A direct-sum 0` on **all**
`L2(mu)`. A support statement for arbitrary local inputs must use this
extension: subtracting a mean is itself nonlocal. Let `d` be intrinsic
distance on the unit-circumradius flat comparison polyhedron. A bounded
real `kappa`-Lipschitz function `psi` is sampled at cell centroids and
repeated on the twelve ports. It has zero difference on internal seams;
on a cross seam its difference is at most `kappa*2^-l`. Cross seams form
a port matching. Expanding the Hermitian part edge by edge therefore gives

```text
Re <u,e^psi_l D_l e^-psi_l u>
 >= -4^l [cosh(kappa*2^-l)-1] ||u||^2.                  (FS1)
```

For an edge with weight difference `delta`, its real contribution is its
unweighted squared difference minus
`2(cosh(delta)-1) Re(conj(u_a)u_b)`. Bound the latter by
`(cosh(delta)-1)(|u_a|^2+|u_b|^2)` and sum. The matching, not the full
native degree, gives the coefficient in (FS1). The differential norm
inequality for the finite matrix semigroup yields

```text
||e^psi_l exp(-tau D_l)e^-psi_l||
 <= exp(tau*4^l [cosh(kappa*2^-l)-1]).
```

Cell compression commutes with these cell-repeated multipliers. The
centroid weights converge uniformly to `psi`, and the full-port heat
convergence above extends exactly to the constant mode. Taking the same
subsequence gives
`||e^psi exp(-tau Abar)e^-psi||<=exp(tau*kappa^2/2)`.
For inputs supported in sets `E,F` at distance `r>0`, use
`psi(x)=kappa*min(d(x,E),r)` and optimize at `kappa=r/tau`. Thus

```text
||1_F exp(-tau Abar)1_E|| <= exp(-r^2/(2*tau)), tau>0.    (FS2)
<g,cos(t sqrt(Abar))f> = 0 when |t|<sqrt(2)*r.            (FS3)
```

The last implication is the nonnegative self-adjoint
[Davies–Gaffney/finite-propagation equivalence of Coulhon and Sikora,
Theorem 3.4](https://arxiv.org/abs/math/0609429), applied with distance
`sqrt(2)*d`. Its normalization is `exp(-distance^2/(4*tau))` for unit
wave speed. Consequently (FS3) bounds speed by `1/sqrt(2)` in the declared
comparison distance per coherent wave parameter. Integrating the cosine
also gives this support bound for the sine-over-generator propagator of
a localized initial velocity.

The exact source control verifies (FS1) using rational exponential
weights. Allowing weights to vary inside a carrier defeats its
cross-seam-only penalty. At level zero a cross-seam matrix element of the
**finite** continuous cosine starts with `t^2/2` (with `4^l*t^2/2` at
level `l`), so exact finite-level
propagation is not asserted. Nor does (FS3) make the stationary Gaussian
field's long-range correlations vanish: it concerns propagation of local
changes. This is a finite-speed continuum wave law on the supplied
two-dimensional comparison topology. It does not identify a measured
metric, light speed, Lorentz symmetry, or laboratory clock.

### Heat response and the boundary of a finite signal cone

The passive heat dynamics have a different support law from (FS3).
The full-space extension `Abar=A direct-sum 0` has
`ker(Abar)=span{1}`, since the mean-zero generator has the proved positive
gap. Its conservative positivity-preserving heat semigroup is irreducible.
Indeed a closed subspace invariant under all bounded multiplications is
`L2(B)` for a measurable set `B`. If it is also heat invariant,
self-adjointness makes its orthogonal complement invariant; multiplication
by `1_B` therefore commutes with the heat semigroup. Conservation gives
`exp(-tau*Abar)1_B=1_B`, so `1_B` belongs to the constant kernel.
Consequently `B` has zero or full measure. The
[positivity-improvement equivalence, Keller–Lenz–Vogt–Wojciechowski,
Theorem A.2](https://arxiv.org/pdf/1101.0373), then proves

```text
<g,exp(-tau*Abar)f> > 0
for every tau>0 and nonzero nonnegative f,g in L2(mu).        (HR1)
```

This is an almost-everywhere positivity statement, with no assumed smooth
heat kernel or point-detector readout. It cannot be applied directly to
the mean-zero projection: subtracting a mean destroys positivity and
local support.

The support result also applies to an intervention mean.
For nonzero nonnegative bounded Lipschitz profiles `f,g` with supports at distance
`d>0`, use the [existing activation controller](../predictions/README.md#a-controlled-continuum-mean-with-an-executed-event-clock)
with `h_l=cell_average(f)/sqrt(V_l chi_l)` and retain the density readout
against `g`. Subtracting either profile's `mu` mean changes neither the
activation differences nor its response. At a fixed number `q` of
**selection slots**, the already verified exact response is
`2chi_l(I-Q_l^q)h_l`, where
`Q_l=(I+K_l)/2=I-L_l/(4M_l)`. Put
`q=floor(4M_l*4^l*tau)=floor(2520*16^l*tau)`.
This formula holds for every profile: on each seam the activation slope
is antisymmetric under its endpoint swap, so the uniform preparation's
incoming execution derivative cancels. The activity derivative is
`-h_l^T(I-K_l)n`. Uniform reversibility, the exact linear density transfer,
and `sum_(j=0)^(q-1) Q_l^j(I-K_l)=2(I-Q_l^q)` give the stated response.
In the normalized density readout this is exactly

```text
d E[eta_l(g)]/d epsilon|_0
    = (2/V_l) g_l^T(I-Q_l^q) f_l
    --> 2<g,(I-exp(-tau*Abar))f>.                           (HR0)
```

Cell averaging preserves positivity and converges strongly to the test
profiles (uniformly for these Lipschitz tests); their cell supports are
disjoint once the mesh diameter is below `d`. The same full-port heat
convergence applies with the doubled slot normalization, including the
integer rounding. Thus `R_tau=2(I-exp(-tau*Abar))` is the limiting linear
response map for this test class. This is a derivative at zero bias;
the coordinate-specific fixed-amplitude error budget is not asserted for
these arbitrary localized profiles.
No state-dependent executed-clock stopping correction is substituted for
this fixed-slot calculation. The existing executed-clock derivative has
the same limit, but is a separate finite formula.

Equations (HR1) and (FS2), and the disjoint supports, now give

```text
<g,R_tau f> = -2<g,exp(-tau*Abar)f> < 0,
0 < |<g,R_tau f>| <=2 exp(-d^2/(2tau)) ||f|| ||g||,
for every tau>0.                                           (HR2)
```

Thus this particular response has no exact finite support cone, even
though its distant short-time tail is exponentially bounded. At each
finite source the controller changes activation probabilities only on
seams meeting the profile support; outside it the differences vanish.
Conservation explains the negative remote density response. Its
nonvanishing is a derivative of an intervention mean, not an inference
of signalling from equilibrium correlations.

If the same support dictionary, control operation, readout and parameter
were identified with physical local interventions and elapsed time,
(HR2) would contradict an exact finite-cone response: any `c<infinity`
would require zero response when `0<tau<d/c`. A positive scalar clock
conversion leaves the contradiction intact. The result does not supply
that physical identification. Uniform seam selection, activation control,
initial/current disclosure, field reconstruction and the aggregate slot
clock are declared services; their physical locality and accessible
duration have not been derived. It therefore establishes no observed
faster-than-light signal, excludes no alternative observer dictionary or
microscopic causal completion, and does not invalidate the separate
comparison-wave support theorem (FS3). A hydrodynamic diffusion limit
need not preserve the causal cone of another microscopic theory.

The exact source control uses
`exp(-tau*L)=exp(-D*tau) sum_k tau^k(DI-L)^k/k!`, with
`D` the maximum native degree. Every summand is nonnegative, and a
shortest native path gives a strictly positive coefficient. In contrast,
the fixed finite selection-slot kernel is exactly zero before that many
steps. Deleting the cross seams leaves disconnected carriers and zero
cross-component heat, while projecting away constants gives a negative
off-diagonal entry even at zero time. These controls preserve the
connectedness, constant-mode and finite-versus-continuum distinctions.

### A continuous wave readout from one native record history

The existing native functional CLT supplies independent Gaussian limits
for two disjoint consecutive long observation blocks of the **same**
stationary history. No independent carrier replicas are required. At
fixed level their normalized cell fields `xi_l,zeta_l` have covariance

```text
C_l=B_l^+ - c_l P0,       c_l=1/(2520*16^l).
```

Define the declared spectral processing of these two recorded fields by

```text
X_l(t)=cos(t sqrt(B_l)) xi_l + sin(t sqrt(B_l)) zeta_l.
```

The sine coefficient is the velocity divided by `sqrt(B_l)`, so the two
initial Gaussian fields have the same covariance. The source bounds
`0<=B_l<=4^l I/2` and `B_l C_l=P0-c_l B_l` imply
`0<=B_l C_l<=I`. Both covariance and conserved constant mode are checked
exactly on the actual twenty-cell source. Independence of the two CLT
increments is essential: reusing one field changes the time-dependent
variance and fails the stationary covariance identity.

The joint law is centered Gaussian with

```text
Cov(X_l(t),X_l(u)) = C_l cos((t-u) sqrt(B_l)).
```

Because `g cos(t/sqrt(g))` extends continuously by zero at `g=0`, norm
convergence of `G_l` and `c_l->0` gives every finite space-time covariance
along the same subsequence. Its limit is

```text
Cov(X(t),X(u)) = G cos((t-u) sqrt(A)).                    (CW)
```

The spherical harmonic topology used for the native spatial field makes
this a continuous distribution-valued process. With harmonics normalized
for spherical probability area, `dmu/domega<2`, and with `Pi_>J` omitting
degrees through `J`, the previously proved harmonic summation gives

```text
E sup_(0<=t<=T) ||Pi_>J X_l(t)||_-2^2
    <= 4(C+T^2)/(J+1)^2,
E ||X_l(t)-X_l(u)||_-2^2 <= 2 |t-u|^2,
E ||X_l(t)-X_l(u)||_-2^4 <= 12 |t-u|^4.
```

The first bound follows from the initial tail plus the integral of the
velocity tail; its covariance is bounded by the identity. The second
uses the same velocity bound and Cauchy–Schwarz; the third is the Gaussian
fourth-moment inequality. Finite harmonic paths are tight by the increment
bound, and the uniform supremum tails control the remaining modes.
Hence `X_l` converges in `C([0,T],H^-2(S2))` to the nonzero Gaussian wave
readout (CW). The bounded covariance extends its Gaussian test variables
to all `L2(mu)` smearings. For `f in D(A)`, these variables satisfy
`d^2 X(t;f)/dt^2 + X(t;Af)=0` in mean square, distributionally in time.
This does not assert pathwise pairing of an `H^-2` distribution with every
operator-domain test; that would require additional spatial regularity.
Its displacement covariance is
`A^-1` and its velocity covariance is the identity, which specifies a
classical Gaussian wave ensemble. It is not the quantum vacuum covariance
`(2sqrt(A))^-1` and asserts no canonical commutation relations.

Observation count tends to infinity first at each finite source. The
source subsequence then tends to the continuum. These limits do not give
a simultaneous finite-count/refinement rate. The construction processes
native records through a declared spectral readout; passive occupation
waiting retains its different heat-type correlation law. Coherent native
amplitudes converge independently to the **same** cosine operator, with
their admitted preparation, reflection and readout services.

### Modular phase matching and its order of limits

The [native return reference](README.md#a-native-return-reference-for-the-coherent-wave-generator)
has positive phase generator `Theta_l=acos(K_l)`. Its finite-cap version
`Theta_l,h=-log(P_h(K_l)/P_h(1))` satisfies the existing uniform logarithm
bound. The sufficient joint budget `h_l=64*(4^l)^4` gives

```text
||sqrt(Q_l)(Theta_l-Theta_l,h_l)|| <= 26/4^l.
```

For bounded `t`, its reconstructed modular cosine therefore has the same
limit `cos(t sqrt(A))`; the additional error is at most
`min(2,26|t|/4^l)`. This compares spectral phase magnitudes and cosine
amplitudes. The arc walk has both phase signs and does not implement every
vertex modular unitary or modular intervention. The
[full-reference limit below](#the-full-reference-state-and-its-quantum-fluctuation-boundary)
shows why the rescaled phase match does not supply the desired quantum
state: these full-configuration references approach normalized trace, and
the density-wave sector has vanishing reference weight.

The cap order matters in a stronger way than absence of a uniform error
certificate. For a **fixed** `h`, positivity of the polynomial coefficients
gives `P_h(x)/P_h(1)>=x^h`, hence, on any fixed band `a<=Lambda<2Q_l`,

```text
sqrt(Q_l) Theta_l,h(1-a/(2Q_l))
    <= 26 h Lambda 4^l/(1260*16^l-Lambda) -> 0.
```

A fixed capped reference freezes its low-band modular phase under spatial
refinement; the uncapped wave phase tends to `sqrt(a)`. The executable
checks this control against the actual positive length polynomial.
The increasing sufficient budgets are analytic resource bounds, not
executed samplers at those large lengths.

This construction connects a spatial covariance, a continuum action,
native temporal correlations, continuous wave observations and matching
modular phase magnitudes through one source operator. The source sphere
is a two-dimensional readout domain. Neither its three ambient coordinate
functions nor the internal rank-three band add a physical radial variable.
The comparison wave support cone is proved above. A Lorentzian metric,
an observed `1+3` identification, interacting CFT/QFT and atomic-clock
calibration require additional proved observable relations. They are not
conclusions of this continuum limit.

```text
python -m model.observer_time_fields.continuum_dynamics --check model/observer_time_fields/continuum_dynamics_results.json
python -m pytest -q -W error model/observer_time_fields/test_continuum_dynamics.py
```

### The full reference state and its quantum-fluctuation boundary

The native return construction normalizes over **all** fixed-number
configurations, not just its density-mode subspace. This distinction has
a quantitative consequence on the same refining source. Write
`V=240*4^l`, `q=V/3`, `d=binom(V,q)`, `M=21V/8`, and let
`tau(X)=Tr(X)/d` on the full configuration algebra. For each native seam,
let `S_e` swap its two binary occupations. The actual transfer is

```text
K=(1/M) sum_e (I+S_e)/2,       kappa=tau(K)=1-chi,
chi=q(V-q)/[V(V-1)].
```

This uses the unchanged uniform seam law and fair native repair.
Normalized permutation traces count fixed binary subsets. One swap,
two swaps sharing an endpoint, and disjoint swaps have respective traces
`1-2chi`, `1-3chi`, and `1-4chi+4psi`, where
`psi=q(q-1)(V-q)(V-q-1)/[V(V-1)(V-2)(V-3)]`.
There are `M` equal ordered edge pairs, `45V/2` distinct ordered pairs
sharing a vertex, and `M(M-1)-45V/2` disjoint pairs: exactly `3V/4`
ports have degree five and `V/4` have degree six. Consequently

```text
v_l := tau((K-kappa I)^2)
 = 2(365V^2-1631V+1560)/[11907(V-2)(V-1)^2] <= 1/(16V).    (S1)
```

The executable verifies the census on the actual source levels and
independently sums all native endpoint/coin branches in the twelve-port
total-two and total-four sectors. The inequality is all-level algebra:
subtracting (S1) from `1/(16V)` has positive denominator and numerator
`227V^3+4564V^2+9615V-23814`, positive for `V>=3`.

Set `f(x)=exp(-acos x)` and `sigma=f(K)/Tr f(K)`. For the present density,
`2/9<=chi<=1/4`. Convexity of `f` gives the relative-to-`kappa` bound
`|f(x)-f(kappa)|<=|x-kappa|/chi` for `0<=x<=1`; its largest secant
slope is bounded by the secant to one. Also `f>=exp(-pi/2)>1/5`.
Cauchy–Schwarz applied to the spectral trace therefore gives

```text
(1/2)||sigma-I/d||_1 <= (5/chi) sqrt(v_l)
                           <= 3/(8*2^l) = epsilon_l.      (S2)
```

The estimate holds on the full algebra, uniformly over bounded readouts.
It concerns the reconstructed reference; it does not claim that coherent
`sigma` has been prepared on the native working source. For the admitted
growing cap `h_l=64*(4^l)^4`, its additional half-trace error is at most
`delta_l=1/[2*(4^l)^2]`, so the capped reference obeys (S2) with
`epsilon_l+delta_l`.

There is a corresponding statement about **fixed dimensionless modular
parameter** `u`. Let `Theta=acos(K)`, `theta0=acos(kappa)` and
`||X||_(2,tau)^2=tau(X*X)`. The concavity of `acos` gives
`|acos(x)-theta0|<=8|x-kappa|`. For `||B||<=1`,

```text
||exp(-iuTheta) B exp(iuTheta)-B||_(2,tau)
 <= 16|u| sqrt(v_l) <= 4|u|/(15*2^l).                    (S3)
```

The cap contributes at most `2|u|/(4^l)^2`. Relative reference densities
are bounded by five times normalized trace, so the same vanishing holds
in the reference-state `L2` norm. This is not operator-norm convergence.
Nor does (S3) imply triviality for the wave scaling
`u=t sqrt(630)*4^l`, which diverges. Rescaled low-mode phases can remain
nontrivial while fixed-parameter bounded-state observations become tracial.

More generally, any compatible product-preserving bounded algebra readouts
whose moments have a limit inherit a tracial state: cyclicity of `tau`
and (S2) make every bounded commutator expectation vanish. On the faithful
support of such a tracial limiting state its modular automorphism is the
identity. This conditional statement does not supply the missing algebra
embeddings or preclude a differently scaled dynamical limit.

The exact native density-mode subspace has dimension `V-1` and is invariant
under `K`. Its orthogonal projection `P_density` satisfies, for the ideal
reference and every polynomial cap,

```text
Tr(sigma P_density) <= 5(V-1)/binom(V,V/3) -> 0.           (S4)
```

Here the scalar density relative to normalized trace is at most five:
for the cap polynomial, `P_h(1)/P_h(0)<=exp(pi/2)<5`.
The bound `binom(V,V/3)>=3^(V/3)>=2^(V/3)` makes the decay explicit without
allocating the configuration space. Renormalizing this subspace would
therefore be a new sector-conditioning operation with a vanishing success
probability under the full reference. The nonzero density-wave generator
alone does not identify its state as the unconditioned source vacuum.

This also identifies a precise boundary for a quantum fluctuation claim.
Recorded occupations and their deterministic processing commute exactly.
Allow instead coherent local observables `A_j,B_j` of norm at most `B`,
at most `aV` terms per family, with at most `b` overlapping terms of the
other family for each term. Assume disjoint supports commute and every
operator preserves the binary fixed-number sector. For centered fluctuations
`F=V^(-1/2)sum A_j`, `H=V^(-1/2)sum B_j`,

```text
||[F,H]|| <= 2ab B^2,
|Tr(sigma_h [F,H])| <= 4ab B^2 (epsilon_l+delta_l) -> 0.   (S5)
```

Thus a quantum central limit with moment convergence on this specified
fluctuation domain cannot have a nonzero scalar canonical commutator.
No quantum CLT for the refining source is assumed or proved here. The
[tracial quantum CLT of Kuperberg](https://arxiv.org/abs/math-ph/0202035)
is relevant context, but its IID-copy premise is different from this
source refinement; (S5) follows directly from trace distance and locality
of the tested operator family.

The finite references are not exactly tracial. With `A=n_0`,
`B=i[A,S_(0,1)]` and the cap-one reference `(I+K)/Tr(I+K)`, the native
oracle gives `-i Tr(sigma_1[A,B])=chi/[M(2-chi)]`, equal to `2/435`
in the total-four carrier. A separately supplied coherent pure state on
two swapped configurations gives one instead. These controls exclude an
overbroad no-go: different nontracial preparations, different observable
scalings, unbounded overlap or a justified sector-conditioned construction
need their own analysis. The present full reference does not obtain the
quantum vacuum covariance `(2sqrt(A))^-1` by its existing unconditioned
bounded-algebra limit.

### Fractional wave covariance and its quantum boundary

The declared filter `A_l^(1/4)` changes the synthesized classical Gaussian
wave to `W_{1/4,l}`, with stationary covariance

```text
C_l(t,u)=A_l^(-1/2) cos((t-u)sqrt(A_l)).                    (W4)
```

It equals twice the real part of a canonical free-field vacuum two-point
function. Its equal-time covariance is the same `G_l^(1/2)` as the
[half-cylinder comparison](../physical_space_time/OCCUPATION_FIELD.md).
These identities supply neither quantum commutators nor a physical bulk.
Velocity and gradient-energy variances both equal `sqrt(lambda)` in each
mode, so their finite-mode total energy is `Tr sqrt(A_l)`. The executable
checks this against independent Schur elimination and matrix square roots.

The path statement uses **fixed spherical harmonics**, not eigenvalue-index
weights. For their centered cell averages `h_jm,l`, (OP4), Schur
minimization and Cauchy–Schwarz give

```text
sum_m <h_jm,l,sqrt(A_l) h_jm,l>
 <= sqrt(sum_m ||h_jm,l||^2 * sum_m <h_jm,l,A_l h_jm,l>)
 <= sqrt(2) (2j+1) sqrt(j(j+1)) <= 12 j^2,  j>=1.
```

The stationary increment covariance is
`2 A_l^(-1/2)[I-cos((t-u)sqrt(A_l))] <= |t-u|^2 sqrt(A_l)`;
both independent Gaussian quadratures contribute. With
`w_j=1+j(j+1)`, for every `s>3/2` the geometric `H^-s` second-moment
increment is bounded by `C_s |t-u|^2`, where
`C_s=12[1+1/(2s-3)]`. Its fourth moment is at most
`3 C_s^2 |t-u|^4`. The initial covariance is bounded by `sqrt(C) I`,
`C=4096/3`; the velocity tail after degree `J>=1` is at most
`12 J^(3-2s)/(2s-3)`. Integrating the velocity bounds the expected squared
supremum tail by twice the initial tail plus `2T^2` times this vanishing tail.
Finite harmonic tightness and these tails prove convergence along the same
spatial subsequence in `C([0,T],H^-s(S2))`. The executable retains rational
budgets for `s>=7/4`; its separate eigenvalue-index diagnostics are not
labeled geometric Sobolev norms.

The sine kernel `D_perp(t)=A^(-1/2)sin(t sqrt(A))` is the candidate
Pauli–Jordan kernel of a quantum completion. For separated **mean-zero**
supported tests its pairing vanishes outside the existing comparison cone.
For unrestricted local tests the correct kernel is
`D_bar(t)=D_perp(t)+t P_const`. Omitting the constant mode leaves a nonlocal
`-t P_const` term, so the projected kernel alone does not define a local net
for all smearings. Restoring this massless zero mode gives a free particle,
which has no normalizable stationary ground state. A mass, restriction to
zero-mean observables or a nonstationary zero-mode state is an additional
choice. The code bounds finite-level dipole leakage without asserting an
exact finite-level cone or a native commutator from classical covariance.

### Source-coded quantum fluctuations and an ordered detector law

The [source compiler](../observer/source_scalar/README.md) already realizes
proper-code Pauli controls, blank registers, coherent code gluing and complete
export outcomes. In `N` independently prepared `|0>` codes define
`Q_N=sum_j X_j/sqrt(2N)` and `P_N=sum_j Y_j/sqrt(2N)`. Then

```text
[Q_N,P_N]=i sum_j Z_j/N,
E Q_N^2=E P_N^2=1/2,       E Q_N P_N=i/2,
E Q_N^4=3/4-1/(2N).
```

The finite commutator is not `iI` in operator norm. For
`v=(s,t)`, set `W_N(v)=exp(i(s Q_N+t P_N))`. Tensor independence gives
its **ordered** product expectation as the `N`th power of the one-code
matrix element. Taylor expansion, or the IID quantum CLT of
[Jakšić, Pautrat and Pillet, §§3–4](https://jaksic.xyz/papers_pdf/qcltiid.pdf),
therefore gives, for every fixed finite list,

```text
lim_N E product_j W_N(v_j)
 = exp(-|sum_j v_j|^2/4 - i sum_(j<k) sigma(v_j,v_k)/2),
sigma((s,t),(s',t'))=s t'-t s'.                            (QF)
```

This is the CCR vacuum state, including its noncommuting phase, on the
chosen fluctuation algebra. The finite source execution uses four codes and
one pointer, retains every work/pointer outcome, and reads both quadratures
of the ordered product by controlled source pulses. The phase-space rotation
angles use a fixed `2^-44` grid; telescoping bounds their operator error by
the sum of angle errors. Reversing the order
conjugates the phase; the tracial product reference removes it. These
controls distinguish (QF) from a classical covariance fit. Growing `N` in
this IID preparation is separate from source spatial refinement.

For a native density eigenvalue `lambda` at level `l`, the coherent walk
fixes `theta_l=acos(1-lambda/(2Q_l))`, `Q_l=630*16^l`, in scaled-Laplacian
units. The declared mode encoding applies `exp(i theta_l Z/2)` to each
code per incident trigger. Its positive oscillator generator has frequency
`Omega_l=sqrt(Q_l) theta_l`, approaching `sqrt(lambda)` on bounded bands
by (W). With displacement normalized by `1/sqrt(Omega_l)`, the encoded
ordered covariance is exactly `exp(-i Omega_l(t-u))/(2 Omega_l)` even at
finite `N`; higher products approach the free-field state by (QF).
For a fixed spectral band its phase error obeys the same clock budget as
the wave, since `|exp(ix)-exp(iy)|<=|x-y|`.
The existing [active incident-clock theorem](../relative_observer_geometry/README.md#active-coherent-waves-in-a-native-incident-clock)
gives `t=(2d/sqrt(630)) tau`, `d=5` or `6`, under its protected-register
and atomic-service contract. This is a common native frequency and event
clock dictionary, with pulse-service costs still charged separately.

A source kick `exp(i epsilon phi(u))` followed by the same mode readout
has derivative `i E[phi(t),phi(u)]=sin(Omega(t-u))/Omega` at zero kick.
Thus the antisymmetric two-point function also determines an intervention
response, not merely an equilibrium correlation. The elementary execution
retains both pointer outcomes for positive and negative kicks; the central
difference error is bounded by the exact sine remainder. Summing mode
responses gives the sine kernel above under the declared mode dictionary.

The executable checks source pulses at all three elementary native bands.
It does not execute a growing spectral encoder or a growing detector net.
Independent polarized preparation, mode labels, angle synthesis, controlled
readout and protected storage remain explicit services; passive repair has
not selected them. Finite-mode approximation followed by the common spatial
limit defines the mean-zero free-field completion for admissible smearings,
with sine-kernel response on that domain. Field smearings can lie in
`L2_0(mu)`; momentum smearings require `D(A^(1/4))` for finite variance.
The CCR normalization sets the encoded action unit; it does not calibrate
Planck's constant in laboratory units. Its locality has the zero-mode
qualification above. Interactions, conformal symmetry, a physical spatial
bulk and atomic-clock identification do not follow from this encoding.

### Faithful state-clock compatibility

The finite inner-flow result in
[`PrivateInner.lean`](../../Lean/Dynamics/PrivateInner.lean) and the
[finite Gibbs/KMS construction](README.md#the-same-three-class-reference-as-the-compiled-modular-probe)
apply to the encoded mode algebra once its state and generator are specified.
For a code with `H=Omega(I-Z)/2`, `Omega>0`, set
`alpha_t(A)=exp(itH) A exp(-itH)`. A faithful state `rho` has
`sigma_s(A)=rho^(is) A rho^(-is)`. Equality
`sigma_s=alpha_(-b s)` for every real `s` is equivalent to
`log(rho)+bH` commuting with the full matrix algebra. Its center is scalar,
so normalization proves the necessary and sufficient relation

```text
rho=exp(-bH)/Tr exp(-bH)=(I+rZ)/2,
g=log((1+r)/(1-r))=b Omega,    r=tanh(b Omega/2).          (FT)
```

Here `b>0` and `0<r<1`. The polarization axis as well as the eigenvalues
must match the clock. The finite executable tests the whole Pauli basis,
an equally polarized but off-axis state, and equal `r` at unequal native
frequencies. The latter gives unequal inferred `b` and excludes a common
modular time conversion. Matching only integer steps is insufficient:
gaps differing by `2pi` agree at integer parameters and disagree at a half
step. The tracial endpoint `r=0` has trivial modular flow; the pure
reference `r=1` is not faithful on the full code algebra. Neither endpoint
supplies this finite-temperature clock relation.

For `N` independent copies of this faithful state, the canonical
normalization becomes `Q_N=sum X/sqrt(2Nr)`,
`P_N=sum Y/sqrt(2Nr)`. Their commutator expectation is `i`,
their variance is `1/(2r)`, and
`E Q_N^4=(3/4-1/(2N))/r^2`. The same ordered product calculation (QF)
therefore yields the thermal CCR state

```text
lim_N E product_j W_N(v_j)
 = exp(-|sum_j v_j|^2/(4r)-i sum_(j<k) sigma(v_j,v_k)/2).
```

For `phi_N(t)=(Q_N cos(Omega t)+P_N sin(Omega t))/sqrt(Omega)`,
even the finite-copy two-point function is exactly

```text
C(t-u)=[cos(Omega(t-u))-i r sin(Omega(t-u))]/(2r Omega),
C(t-i b)=C(-t),       n_b(Omega)=(1-r)/(2r).
```

The positive and negative frequency weights have ratio
`exp(-b Omega)`, proving the two-point KMS boundary identity. Ordered
Gaussian products in the copy limit give the full thermal free-field
state. The antisymmetric part and linear kick response remain
`sin(Omega(t-u))/Omega`; the symmetric noise changes with `r`. Thus
temperature is inferable from noise relative to the same response and
frequency, after fixing the encoded quadrature normalization.

The execution uses two system codes, two retained purifying references
and a pointer. A source `Y` rotation and source CNOT prepare each
`sqrt((1+r)/2)|00>+sqrt((1-r)/2)|11>`; restriction to the system codes
gives the faithful product state. All 32 ordered Ramsey outcomes and all
16 outcomes for each signed kick are retained. No success branch is
discarded. The state on the enlarged pure system/reference algebra is
not faithful. Frequencies and polarizations are evaluated at 80 digits and
pinned to the existing `2^-44` dyadic synthesis grid; their exact receipt
rationals do not depend on platform floating-point library rounding.
Rotation angles use the same grid;
the receipt separately checks preparation error, full Pauli modular
action and finite-kick remainder. Its rational polarization fixes an
exact logarithmic `b`; numerical diagnostics compare this with the
requested common value. The ideal relation (FT) is exact, while the
finite source preparation has its stated synthesis error.

For the mean-zero continuum operator, the existing
[native spectral bound (L4)](../physical_space_time/OCCUPATION_FIELD.md)
gives `lambda_n(A)>=n/49920`. With `Omega=sqrt(A)`,

```text
Tr exp(-b Omega) <= integral_0^infty exp(-b sqrt(x/49920)) dx
                 = 99840/b^2,
log Z_b = sum_(m>=1) Tr exp(-mb Omega)/m <= 149760/b^2.    (FT2)
```

The last inequality uses `sum m^-3<=3/2`. Consequently the mean-zero
bosonic Fock space carries the faithful normal state
`rho_b=Z_b^-1 exp(-b dGamma(Omega))`, with modular flow
`sigma_s=alpha_(-b s)` on its full encoded algebra. Its mode covariance
is `coth(b Omega/2)/(2 Omega)`; the lower spectral gap bounds the
thermal factor, so the field and momentum smearing domains above still
apply. Finite-mode copy limits followed by the existing spatial
subsequence define this state and clock together. This is a conditional
continuum construction, not execution of a growing source compiler.

Using the active incident parameter `t=kappa tau`,
`kappa=2d/sqrt(630)`, gives `s=-kappa tau/b`; inverse temperature in
`tau` units is `b/kappa`. At refinement level `l` the continuum clock
uses `Omega_l=sqrt(Q_l) theta_l`, so its preparation requires
`r_l=tanh(b Omega_l/2)` and `g_l/theta_l=b sqrt(Q_l)`.
Holding polarization fixed while changing this clock normalization
would change temperature. A common `b` is a measurable preparation
compatibility condition. Carrier architecture and passive repair have
not selected it. The encoded full-field algebra is not a spatial patch
algebra: restriction to a patch need not preserve this modular generator.
The excluded constant mode and its locality qualification remain as
above. Geometric boosts, an Unruh temperature, interactions and
laboratory seconds do not follow from (FT).
The screen construction uses `sqrt(A)`. The following encoding joins the
radial product field to an active clock under the same source services.

### One product field, state and active clock

Fix the [radial product dictionary](../physical_space_time/OCCUPATION_FIELD.md#independent-native-records-give-a-radial-product-field):
an interval `(0,L)` with Dirichlet ends, a radial grid of `q` interior
points and spacing `h=L/(q+1)`, and the actual full native angular generator
`D_l=4^l L_l`. Its radial frequencies satisfy
`c_k,h^2=4h^-2 sin^2(k pi/(2(q+1)))`. A native angular mode of scaled
eigenvalue `lambda` has the already established coherent phase
`theta_l=acos(1-lambda/(2Q_l))`, `Q_l=630*16^l`, and frequency
`a_l=sqrt(Q_l) theta_l`. Use one protected code per mode and copy.
Each zero-inclusive incident event triggers two existing source pulses,
in chronological order:

```text
Z(theta_l/2), X(c_k,h/(2 sqrt(Q_l))).
```

These pulses implement `exp(i c X/(2sqrt(Q))) exp(i a Z/(2sqrt(Q)))`.
The Pauli axes anticommute, so

```text
(a Z+c X)^2=(a^2+c^2) I,
omega_l,k=sqrt(a_l^2+c_k,h^2),
H_l,k=[omega_l,k I-a_l Z-c_k,h X]/2.                      (PC)
```

Thus the product dispersion is determined by composition of the native
angular clock and declared radial apparatus. No independent fitted product
frequency is applied as a pulse. The finite split word is **not** exact
evolution by (PC). Duhamel's formula and
`||[iaZ/2,icX/2]||=ac/2` bound one step by `ac/(4Q)`;
telescoping bounds `n` steps by `n ac/(4Q)`. The retained execution
measures the split-order difference and the finite state's stationarity
defect. Replacing the radial `X` axis by `Z` changes the gap to `a+c`
and fails the independent exponential check.

For the exact joint generator (PC), write
`N=(cX+aZ)/omega`, `E=(aX-cZ)/omega`, and retain `Y` as the other
transverse axis. The faithful state is
`rho=(I+rN)/2`, `r=tanh(b omega/2)`. A source `Y` rotation of the
existing retained purification prepares this axis; it also rotates every
field readout and intervention. The relation (FT) now gives the modular
flow of **this same product generator**. IID copies have transverse
fluctuations `sum E/sqrt(2Nr)` and `sum Y/sqrt(2Nr)`, so (QF) gives
the joint thermal CCR state. Ordered products at two different times
are measured by controlled source words with all 32 system, reference
and pointer outputs retained. Signed kicks retain all 16 system/reference
outputs. Independent matrix calculations check the actual finite split
word; its ideal modular action and sine response have explicit split
error budgets. They are not asserted to hold exactly before the limit.
The modular check uses the actual reduced source state at `s=-t/b`;
its budget includes the rounded polarization's gap error relative to the
requested `b` and the compiled preparation/axis error.

The finite nonconstant compiler fixture is an exact full 240-port mode:
`(0,0,0,0,1,-1,-1,1,0,0,0,0)` on one carrier and zero elsewhere has
eigenvalue `6` and vanishes on every external seam endpoint. It has zero
cell compression, so it is an internal control, **not** a nonconstant
mode surviving in the continuum readout. Under identical local refinement
its scaled eigenvalue is `6*4^l`, outside any fixed low band. The constant
angular mode is separately executed and does survive. Extension to bounded
full-port bands uses the spatial compression/resolvent theorem, not
equality of a finite full-port matrix with its cell Schur complement.
No growing nonconstant-mode encoder is executed.

The finite cell-visible encoder uses the actual level-zero 240-port
Laplacian directly. Starting from the disclosed cell contrast
`w0=(1 on cell0,-1 on cell1,0 elsewhere)`, it forms
`w=(12I-L)^1024 w0` by exact integer arithmetic and gcd reduction.
The normalized packet `v=w/||w||` has exact zero mean, Rayleigh quotient
`lambda=<v,Lv>` in `[.044442251696,.044442251697]`, residual
`||(L-lambda)v||<.000304`, and cell-visible squared weight
`||V* v||^2>.996458773594`, where `V` is the normalized cell embedding.
This deterministic polynomial fixes a packet without a numerical
eigensolver or an externally supplied mode. Its nonzero residual is kept;
it is not declared an exact eigenvector or an independently identified
continuum mode.

Two Householder reflections map logical addresses zero and one to the
native constant and this packet. Each reflection is compiled as a source
binary-tree state preparation, a zero-address phase, and the inverse
preparation. The resulting eight-qubit word is executed on all 256 input
columns; all sixteen padding addresses remain outputs. A source-prepared
Bell reference, the encoder, its source phase and inverse decoder give a
512-address interference readout. Each incident trigger applies the
phase `theta=acos(1-lambda/(2M))`, `M=630`, to the encoded packet.
The whole-register pulse word implements this restricted phase operation;
it does not replace the full native propagator by a target matrix.

For the actual native compression `K=I-L/(2M)`, the spectral theorem and
`sup_[-1,1]|T_n'|<=n^2` give

```text
||T_n(K)v-cos(n theta)v|| <= n^2 ||(L-lambda)v||/(2M).
```

The retained 25-trigger native discrepancy is about `.000149026`, below
this bound; source synthesis and phase rounding are charged separately.
Its wave parameter is `t=n/sqrt(630)`, with the same angular dispersion
correction as (PC). The native one-attempt census retains all 1,260
seam/coin/current records, including incident events with zero current.
Those records decode the same packet field through
`F(after)=F(before)+current*(v_b-v_a)` and trigger its reference readout.
The 25-trigger encoded experiment costs 1,747,068 source pulses, each
with one proper-code load and export, on eight work qubits and one
reference qubit. Supplied storage, controller access, initial preparation
and paused coherent services remain premises. This is a cell-visible
finite native packet encoder, separate from the prepared twelve-state
Gibbs marginal; it does not construct the full growing mode family,
thermal selection, radial apparatus or laboratory units.

The native fixed-port incident count has mean `2d*4^l` per density
parameter unit `tau`, for `d=5` or `6`. Hence the product word has wave
parameter `t=(2d/sqrt(630)) tau`, just as the angular clock. On a bounded
joint band, `product_clock_bound` adds the native count error, angular
dispersion, two-axis splitting and pulse synthesis error. For a code with
`n` triggers it charges two pulses per trigger. Count tolerances may tend
to zero with `4^l tolerance_l^2 -> infinity`; complete-step synthesis
errors must satisfy `sqrt(Q_l) error_l -> 0`. A fixed numerical precision
does not certify unlimited refinement. Timing repairs pause during atomic
services on the protected code; these pulse costs are not laboratory
seconds. The finite joint clock/pointer law sums all native seam/coin
branches, including incident events with zero current.

Taking the copy limit, bounded spectral approximations, radial refinement
and the common native spatial subsequence gives one encoded free field with

```text
B=-d_r^2+(A direct_sum 0),
C_b(t-u)=[coth(b sqrt(B)/2) cos((t-u)sqrt(B))
                         -i sin((t-u)sqrt(B))]/(2sqrt(B)),
D(t-u)=B^(-1/2) sin((t-u)sqrt(B)),
sigma_s=alpha_(-bs),       s=-(2d/sqrt(630))*tau/b.        (PC2)
```

The angular constant is included, and Dirichlet radial boundaries give
`B>0`, avoiding the removed-constant nonlocal correction of the screen-only
field. The full product sine kernel has support in the comparison cone
`|t-u|^2 >= |r-r'|^2+2d_screen^2`. The product counting bound
`N_B(x)<=C_L x^(3/2)`, `C_L=L(49920+L^2/9)/3`, gives
`Tr exp(-b sqrt(B))<=6C_L/b^3` and
`log Z_b<=8C_L/b^3`, using `sum m^-4<=4/3`.
Thus `Z_b^-1 exp(-b dGamma(sqrt(B)))` is a faithful normal state on
the full encoded Fock algebra. Field smearings lie in `L2`; momentum
smearings require `D(B^(1/4))`. The corresponding classical symmetric
thermal covariance is obtained from native radial block records by the
specified filter; only the ordered code construction supplies the CCR.

This joins the product operator, state, intervention response and active
clock conditionally. Radial apparatus, spectral mode access, independent
preparation, common `b`, controlled readout, protection and sufficient
analogue precision remain declared services. Passive repair has not
selected them. The comparison cone does not identify physical radial
transport, a smooth observed metric, an interacting QFT, a CFT, spatial
patch modular boosts or a laboratory time unit.
