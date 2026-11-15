# Coherent observer outputs on matched source screens

This construction joins the existing complete quantum locality theorem to
an executable coherent observation on the actual level-three and level-five
source screens. It retains the local quantum state, an untouched reference,
two incident-clock records, record-controlled feedback and every timeout.
The surrounding source has Bernoulli density `1/3`; its density does not
vanish as the screen grows. The calculation uses a finite operator dual,
not the exponentially large density matrix of the whole screen.

The prior [incident-clock comparison](../native_action/INCIDENT_REFINEMENT.md)
proved a quantum output bound but executed the large screens on diagonal
states. The [quantum observer correspondence](../quantum_observer_correspondence/README.md)
retains reference coherences on one total-two carrier. Here the preparation
and protocol are chosen so that the same complete quantum output can be
computed on the matched 222-port ball and compared with both larger sources.
This is a finite coherent refinement result. It supplies no continuous
physical clock, Lorentzian action, selected repair law or continuum field.

## One state, operation and record contract

| Component | Specification |
| --- | --- |
| Source | The full binary occupation tensor space of each actual source graph, with the existing C or E instrument on every seam. Uniform independent seam choices and fair native placement coins are supplied. |
| Matching region | The root twelve-port carrier `A`, its 33 incident seams, observation support `S=B_1(A)` of 15 ports, and the actual common induced ball `Q=B_9(A)` of 222 ports and 561 seams. Source address maps identify the native instruments. |
| Logical input | An arbitrary qubit on the pair `a=0,b=1`, encoded by `|0_L>=|10>` and `|1_L>=|01>`, possibly entangled with an arbitrary finite untouched reference. |
| Environment | Every other port is independently occupied with probability `rho=1/3`, independently of the logical input/reference. This is a product preparation, not a fixed global total-two sector. |
| Local clock | Every attempted seam meeting `A` ticks, including zero-current attempts. Remote repairs continue on the same source while the observer waits. |
| Public mark | The intersection mask of the selected incident seam with `{a,b}`: `0` for neither, `1` for `a` only, `2` for `b` only, `3` for both. Coin and current records are coarse-grained as CP alternatives. |
| Protocol | Request two scheduled incident ticks. After first mark `1`, apply the declared ideal local phase `Z_a`; otherwise apply no phase. No preparation or environment reset occurs between ticks. |
| Output | The complete record and full reduced pair/reference state. The finite-ball implementation caps each requested interval at `h` native attempts and halts on timeout. |
| Services | Logical/reference preparation, the product environment, the mask readback, timeout service and feedback phase are supplied. The phase contributes a separate gate cost and no native clock tick. |

In the full-screen comparison the target is the completed, two-tick
experiment. The finite attempt cap belongs to the calculation on `Q`.
It is not identified with the differently normalized global attempt count
of either full screen. The controller consults its first mask only; it does
not see remote attempts, hidden occupations or private Kraus labels.

There are exactly 21 capped record labels: 16 two-hit words, four
first-hit/second-timeout words, and one initial timeout. All remain in the
output, including zero quantum blocks. Successful outputs are never
renormalized. Retaining a quantum state is distinct from assigning
simultaneous classical values to incompatible local fields.

## Exact low-degree operator closure

Let `n_i=|1><1|`, `sigma_i^+=|1><0|` and
`sigma_i^-=|0><1|`. On disjoint sets `P,M,U`, consider

```text
O_(P,M,U) = (product_(i in P) sigma_i^+)
            (product_(j in M) sigma_j^-)
            (product_(k in U) n_k).
```

For the unread C and E maps on one binary seam, the following Heisenberg
rules follow directly from the native Kraus maps.

* A seam disjoint from `P union M` averages the occupation monomial with
  its endpoint-swapped monomial. Its occupation degree is preserved.
* A seam meeting exactly one signed site multiplies the monomial by
  `1/2` for C and by zero for E.
* A seam joining two opposite signs kills the monomial under either law.
  A seam joining two equal signs leaves it unchanged.

For the second rule the other endpoint may belong to `U`. The surviving
native branch leaves that occupation unchanged. These statements use C's
equal waiting multipliers `g=h=1/2`, E's `g=h=0`, their common `k=1,m=0`,
and the actual classical exchange law. They are not assertions for the
entire classified record-law family.

The signed support therefore never moves. Diagonal monomials of degree at
most two close under every restricted seam sum, remote waiting power,
incident mask and timeout. This supplies the full pair populations: for
each record, propagate the three terminal functions `n_a`, `n_b`, and
`n_a n_b`, together with the constant function. If their unnormalized
expectations are `x,y,z,w`, respectively, the pair populations are

```text
p_11=z,  p_10=x-z,  p_01=y-z,  p_00=w-x-y+z.
```

For either logical basis input, the initial expectation of a squarefree
occupation monomial is the product of its factors: the two logical sites
are pinned to `10` or `01`, and each other factor contributes `1/3`.
This evaluates a dual occupation law exactly without listing the
`2^222` configurations. At degree two there are only
`1+222+binom(222,2)=24,754` squarefree monomials.

The only initial logical off-diagonal block is
`|10><01|` tensored with the diagonal product environment. A remote seam
leaves the pair unchanged. A seam avoiding both pair sites preserves the
Bernoulli environment. On a one-endpoint seam, the surviving C coherence
is multiplied by `1/2` without moving that environment; E removes it.
Thus the complete off-diagonal output needs one scalar for each word.
The feedback phase changes its sign and leaves every population unchanged.
No moment truncation or Gaussian closure is used.

## The complete logical-input channel

Write `c_m` for the number of root-incident seams with mask `m`,
`pi_m=c_m/33`, and define

```text
chi_C=(1,1/2,1/2,0),       chi_E=(1,0,0,0),
p=33/561,  r=1-p,         q=1-r^h.
```

The mask counts are read from the actual matched source graph, not fitted
to a field law. Since the selected seam is independent of the state,
the complete record probabilities are

```text
w_timeout = 1-q,
w_(m,timeout) = q pi_m (1-q),
w_(m,n) = q^2 pi_m pi_n.
```

They sum to one for both laws and every logical/reference input. The
off-diagonal multipliers of the unnormalized output are exactly

```text
c_timeout = 1-q,
c_(m,timeout) = w_(m,timeout) (-1)^[m=1] chi_L(m),
c_(m,n) = w_(m,n) (-1)^[m=1] chi_L(m) chi_L(n).
```

Let `u_(w,z)` and `v_(w,z)` be the dual-computed population at pair
configuration `z` for logical input zero and one, respectively. For an
arbitrary logical/reference state with reference blocks `rho_ij`, the
entire output block on word `w` is

```text
sum_z |z><z| tensor [u_(w,z) rho_00 + v_(w,z) rho_11]
 + c_w |10><01| tensor rho_01
 + conjugate(c_w) |01><10| tensor rho_10.
```

This formula is an exact channel identity on every logical matrix unit.
Linearity and tensoring with the identity prove its arbitrary-reference
scope. Complete positivity comes from the native instrument composition,
including the environment preparation and partial trace; it is not inferred
from separately plausible population and coherence coefficients. The
retained Bell execution instantiates a two-dimensional reference. The
theorem does not cover arbitrary initial correlations with the environment.

## Quantum witnesses that a classical record fit misses

For Bell input
`(|10>|0_R>+|01>|1_R>)/sqrt(2)`, C and E give identical record probabilities
and identical pair/reference diagonal entries. Their coherent blocks differ.
Define `lambda_L=sum_m pi_m chi_L(m)`. Summing the singular values of the
off-diagonal differences gives the exact complete CQ half-trace distance

```text
( q^2 (lambda_C^2-lambda_E^2)
  + q(1-q) (lambda_C-lambda_E) )/2.
```

Record-controlled signs do not change this distance because the public
words remain separate blocks. Erasing the public record before comparison
can cancel those signs and is a different observation.

There is also a direct entanglement witness. On the success word `(1,0)`,
no executed seam meets `b`. In the reference partial transpose, the
principal block indexed by `|10>|1_R>` and `|01>|0_R>` has zero diagonal
and off-diagonal entries of magnitude `w_(1,0)/4` for C. Its determinant is

```text
-w_(1,0)^2/16 < 0.
```

For E that word has no coherence and is diagonal and separable. A
negative partial-transpose block proves entanglement of the retained CQ
output itself: it is a direct-sum block of that complete output. This
argument does not discard failures or divide by the selected word's mass.
The determinant alone is not a perturbation budget; quantitative transfer
of a witness to a full screen uses its expectation and the complete-output
error bound below.

The Bell input has a diagonal reduced logical state. A phase analyzer
acting only on the pair cannot read its reference coherence. Such an
operational test needs a joint pair/reference correlation readout, or a
separately declared unreferenced coherent preparation. Any complete terminal
instrument can be appended to the proven channel comparison by trace-norm
contractivity, with its own service and resource budget.

## Spatial, cutoff and arithmetic budgets

The existing [stopped-clock locality theorem](../native_action/INCIDENT_REFINEMENT.md)
applies to the same native C/E maps, initial matching restriction and
untouched reference. All mask readouts and the phase feedback are supported
inside `S`. It therefore bounds the complete reduced output of either full
screen against the completed process on `Q` by the source-derived
one-screen allowance `epsilon_space(2)`. This includes the classical
record and its quantum disturbance; the estimate is charged after the
specified output map. It is not an independent error per density entry.

Couple completed and capped ball histories until a requested interval times
out. There are two successful intervals with probability `q^2`, independent
of the logical state and the allowed feedback. Adding the common timeout
labels to the completed output space gives

```text
d(full-screen completed output, capped-ball output)
 <= epsilon_space(2) + (1-q^2).
```

Here `d` is half trace norm on the full CQ output, including arbitrary
untouched references. Comparing the two completed full screens through
the common completed ball uses twice the one-screen spatial allowance;
the cap is an additional computational comparison, not a physical feature
inserted into the completed targets.

The finite dual computation may use directed integer rounding. Entrywise
lower coefficient arrays remain lower under the nonnegative dual maps;
their missing mass bounds expectation errors because each initial
occupation monomial lies in `[0,1]`. Linear population combinations then
use interval arithmetic and the independently known probability range.
No lower coefficient array is renormalized, and a resulting interval is
not presented as an exact CP approximation channel.

For a reference-uniform numerical channel budget, suppose a chosen
population table has coefficient errors bounded by `e_(w,z,i)` for
logical input `i=0,1`, and its coherence coefficient error is at most
`f_w`. Positivity of `rho_00,rho_11` and
`||rho_01||_1 <= sqrt(tr rho_00 tr rho_11) <= 1/2` give

```text
eta_arithmetic <= (max_i sum_(w,z) e_(w,z,i) + sum_w f_w)/2.
```

Exact coherence coefficients make the second sum zero. This bound, the
spatial allowance and the cap tail are distinct terms. Any preparation,
feedback-phase or terminal-service errors must be added under their own
complete-output contracts. No host floating-point diagnostic supplies a
hardware accuracy claim or an outward roundoff certificate.

## Resources and limits

On the capped ball there are at most `2h` native attempts and one phase
gate. The exact expected scheduled attempt count is
`q(1+q)/p`, the expected number of completed incident ticks is `q+q^2`,
and the expected phase-gate count is `q pi_1`. These expectations follow
from geometric waiting and the public mask law; they do not replace the
worst-branch cap in an arbitrary error-composition theorem. Initial
logical/reference and product-state preparation costs are explicit
services. A final measurement, if chosen, contributes its own cost.

The physical pair output has dimension four; the Bell instance has joint
pair/reference dimension eight on each of the 21 records. The polynomial
dual is a method of calculating this channel, not extra physical observer
memory and not an assertion that the global source state has only that
many quantum degrees of freedom.

This result extends the coherent observation domain and attaches it to
actual matched source screens with a complete error budget. It does not
identify the incident count with seconds, make fractional powers of the
Euclidean transfer into CP dynamics, or select C over E as a physical law.
The existing [action audit](../repair_action_audit/README.md) and
[physical audit](../PHYSICAL_AUDIT.md) retain those separate obligations.
Only the stated finite level pair and preparation are compared; an
all-level limit, other environment states and more general interventions
require further estimates or a larger dual closure.

## Retained execution and reproduction

The actual incident-mask counts are `(23,5,4,1)`. Consequently the unread
exchange-coherence multipliers per completed incident tick are `5/6` for C
and `23/33` for E. The complete Bell-output distance after two completed
ticks is `101/968` (approximately `0.10434`); the producer retains the
exact capped value separately. Both laws have the same populations and
clock probabilities throughout this comparison.

The dual runs with a cap of 384 attempts per interval and a 52-bit integer
grid. Its directed arithmetic allowance is below `4.862e-10`, the complete
timeout probability below `1.552e-10`, and the one-screen spatial allowance
below `1.366e-10`. Including outward receipt export, the midpoint linear
table is within **`7.779e-10`** half trace norm of either full-screen
completed output, uniformly over the logical input and reference. The
table is an interval-based computational representation; it is not asserted
to be a separately executable CP instrument. The two completed full-screen
channels differ by at most **`2.731e-10`** under their matching input contract.

An independent five-port oracle enumerates full native configurations with
all four mask categories and genuine remote waiting. It checks all 1,344
entries of the 21 pair/reference output blocks for each law. The canonical
primitive runtime independently composes the stopped native Kraus maps,
checks full source/reference outputs, and serializes, replays and continues
its first-event checkpoint. These small full-state executions audit the
same algebra used by the 222-port operator calculation; they are not
substitutes for that larger calculation.

The tests also reject phase corruption with unchanged probabilities,
incorrect terminal populations, omitted timeout mass, damaged budgets and
malformed exact inputs. Retained runtime diagnostics must satisfy the same
acceptance gates as fresh execution before tolerant replay comparison. An
independent complex three-dimensional reference control checks the channel
allowance against the original native Kraus maps. Cache keys use immutable
integer density ratios; mutating a returned `Fraction` cannot change the
parameters of a cached calculation. No receipt is accepted solely because
its own producer digest matches.

```sh
python -m model.coherent_screen_observer --check model/coherent_screen_observer/results.json
python -m pytest -q -W error model/coherent_screen_observer
```
