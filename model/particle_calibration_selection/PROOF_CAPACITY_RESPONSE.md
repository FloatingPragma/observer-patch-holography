# Finite record capacity and the native family response

The literal reset-port2 memory with native probabilities `p=1/12, q=1/6`
has a definite response at every finite attempt budget. Its first-order
coherent coefficient has rank two after one attempt and rank three after
every integer number `n>=2` of attempts. Uniform bounds on the latter
coefficients give

```math
\lambda_2/\lambda_1 < 1171058/6633 < 177,\qquad
\lambda_3/\lambda_1 < 8616/11 < 784.
```

The inequalities hold for all such integers, without a numerical scan. If
these ordered coefficient eigenvalues are attached linearly to physical
charged-lepton masses, varying the retained attempt budget cannot produce
the measured lepton ratios. This statement concerns the specified native
controller, reset preparation and coefficient readout. It excludes no other
controller, interacting self-energy or physical particle construction.

## Source, capacity and the complete response

The owning [native family-memory source](../common_source/native_family_record_memory.py)
constructs 128 raw C/E labels, including 104 labels of zero weight. Its
native branch law has 24 positive labels and is independent of the old
memory and the working field state. All 12 old-memory inputs are retained.
An enabled C record updates the port by the source half-turn; a simultaneous
E event applies the second half-turn. Other records reset to port2. The
declared feedback service applies `exp(-ih X tensor P_i)` at the updated
address on C-enabled records and identity on every disabled record.
The service's coherent pulse, its angle and the two-component register
are explicit architectural inputs. A compiled native implementation of
this pulse is not supplied by the capacity calculation.

Starting from reset port2, the reached addresses are `{2,3,9,10}`.
The actual rank-one family projectors satisfy `P9=P10`. Hence their unread
response has an exact three-ray quotient, with rays `(P2,P3,P10)`.
Write `a=p(1-q)=5/72`, `b=pq=1/72`, and use column probability vectors.
The full transition and its C-enabled part reduce to

```math
T=\begin{pmatrix}1-a&1-b&1-a\\a&0&0\\0&b&a\end{pmatrix},\qquad
E=\begin{pmatrix}b&a&b\\a&0&0\\0&b&a\end{pmatrix}.
```

`capacity_response.source()` reconstructs the 1,536 original
old-memory/raw-label rows before checking this quotient. It preserves the
difference between disabled and zero-probability labels. The quotient is
used for an unread calculation, with complete fine records retained in the
underlying instrument. It is not an identification of the public record
algebras at addresses 9 and 10.

Let `eta_0=(1,0,0)`, `eta_k=T^k eta_0`, and

```math
w_n={1\over n}\sum_{k=0}^{n-1}E\eta_k,
\qquad Y_n=\sum_{i=1}^3 (w_n)_i P_i.
```

Consider the exact complete unread channel with `n` attempts and pulse
angle `h=t/n`. On the X-plus block, its derivative at `t=0` is
`-i[Y_n,rho]`. Differentiate each pulse at zero and use trace preservation
for every other step. The intermediate memory law is precisely `T^k`;
the result does not replace correlated memory by independent samples.
The full auxiliary-register derivative uses `X tensor Y_n`.
Every column of E sums to p, so `tr Y_n=p` for every n. The first examples
are exact:

| Attempts | P2 weight | P3 weight | P10 weight |
| ---: | ---: | ---: | ---: |
| 1 | 1/72 | 5/72 | 0 |
| 2 | 41/2592 | 695/10368 | 5/10368 |

For fixed reset input and complete raw-record retention, each of the 24
positive raw labels is possible at every step. The deterministic memory
path adds no independent label once the fine raw history is specified.
There are exactly `24^n` positive histories and `128^n` allocated histories.
If the declared continuations are halt and exact reread, the compound
confusability graph has no edges. The existing correctable-record theorem
therefore gives `M0=24^n` and direct archive capacity
`N_record=n ln24` nats. Retaining this exact archive requires at least that
capacity; the literal padded allocation uses `n ln128` nats. The 12-address
controller, source and pulse apparatus have separate resource accounts.

This is a concrete capacity-dependent response through a finite record
budget, once a policy choosing n is specified. The geometric pixel
coordinate P, global physical capacity N, FIFO slots and this archive's
correctable-record capacity are distinct objects. The source supplies no
map identifying them. The [FIFO theorem](../particle_identification/NATIVE_CAPACITY.md)
also proves zero complete-channel change when the same admitted prefix is
stored with a different `export_oldest` FIFO size and the archive is kept.
Specification reads can distinguish the storage metadata; that does not
alter the repair instrument in that theorem's fixed-program domain.

## All-horizon spectral bounds

The three source rays satisfy `tr(Pi Pj)=1/5` for distinct rays. Their
spanning Gram determinant is

```math
d=2/5-2\sqrt5/25>1/5.
```

For positive weights `(u,v,w)`, their coefficient has invariants
`tr Y=u+v+w`, `e2(Y)=4(uv+uw+vw)/5` and `det Y=d u v w`.

The transition gives the invariant bounds

```math
(\eta_k)_1\ge 1-p\quad(k\ge1),\qquad
0\le(\eta_k)_2\le a,\qquad
0\le(\eta_k)_3\le ab/(1-a).
```

The first bound follows from the reset component, and the second from
`(eta_{k+1})_2=a(eta_k)_1`. The last follows by induction from
`(eta_{k+1})_3=b(eta_k)_2+a(eta_k)_3`, starting with zero.
The instantaneous enabled weights are

```math
u_k=b+(a-b)(\eta_k)_2,\quad
v_k=a(\eta_k)_1,\quad
w_k=b(\eta_k)_2+a(\eta_k)_3.
```

Thus `u_k>=b`, `v_k>=a(1-p)` including k=0, and
`w_k>=ab(1-p)` for k>=1. The last bound uses
`(eta_1)_2=a` and `(eta_k)_2>=a(1-p)` thereafter. For every n>=2,

```math
Y_n\succeq L=bP2+a(1-p)P3+\tfrac12ab(1-p)P10.
```

All three lower weights are positive, hence Y_n has rank three. Their
values are `l=(1/72,55/864,55/124416)`. If the eigenvalues of L are positive,
its least eigenvalue is `det L/(lambda2 lambda3)>=det L/e2(L)`.
The strict bound on d therefore gives

```math
\lambda_1(Y_n)\ge\lambda_1(L)>
{l_1l_2l_3\over4(l_1l_2+l_1l_3+l_2l_3)}={11\over103392}.
```

Since `tr Y_n=p`, `lambda3/lambda1<p/(11/103392)=8616/11`.
For the middle eigenvalue, separate the rank-one term `vP3`. The min-max
principle gives `lambda2(Y_n)<=||uP2+wP10||<=u+w`. The invariant bounds imply
`u+w<=b+(a-b)a+ab/(1-a)`. Division by the preceding strict lower bound gives
`lambda2/lambda1<1171058/6633`.

The independently stored [CODATA comparison](../charged_lepton_response/calibration_reference.json)
has `mu/e=206.7682827` and `tau/e=3477.23`, with quoted standard uncertainties
`0.0000046` and `0.23`. Both exceed these source-only bounds by large margins.
No comparison value enters the transition, coefficient or inequalities.
This excludes the interpretation of finite-horizon Y_n eigenvalues as the
three measured masses with one common factor. It does not exclude a
different, derived relationship between the coefficient and measured
threshold or resonance masses.

## Convergence and the finite-step noisy response

The source's reset decomposition `T=(1-p)R_reset+pS` gives
`||T^k(eta_0-pi)||_1<=2p^k`. Since every E column sums to p,

```math
\|Y_n-Y_\infty\|\le {2p(1-p^n)\over n(1-p)},\qquad
Y_\infty=\sum_i(E\pi)_iP_i.
```

This is an explicit transient correction to the configured native
coefficient. Increasing archive capacity at fixed protocol sends it to the
fixed native stationary response. It does not select a different limiting
attenuation law. Stationary product preparation removes this particular
first-order transient for every n. The owner's
[complete memory response theorem](../particle_memory_response/PROOF.md)
supplies the distinct all-lag second-order coherent and dissipative
coefficients under its pulse-clock hypotheses.

There is a direct finite-n noisy witness on the same reset source. Prepare
the family state P2 in the X-plus block. The all-disabled coarse history
has weight `w0=(1-p)^n` and leaves P2 unchanged. The history with only its
last C event enabled, E false, has weight `w1=(1-p)^(n-1)a` and applies
`exp(-ihP3)`. Their squared state overlap is
`1-(16/25)sin²(h/2)`. For an ensemble of pure states,
`1-tr(rho²)=sum_{r,s} w_r w_s(1-|<psi_r,psi_s>|²)`.
Every term is nonnegative, so the complete output satisfies

```math
1-\operatorname{tr}(\rho_n^2)\ge
{32\over25}a(1-p)^{2n-1}\sin^2(h/2)>0
\quad(h\notin2\pi\mathbb Z).
```

All other histories, including zero histories, remain in the source. This
bound uses two coarse groups to lower-bound their complete sum; it does not
condition the output on either group. A complete unitary dilation may keep
coherence jointly with an archive, while the unread working-register
channel is mixed. Consequently its finite-n response cannot be replaced by
a unitary mass correction on every input. `complete_unread_density` verifies
the exact density and positive bound at rational half-angle tangents.

## Physical selection boundary and useful existing bridges

The finite carrier fixes the rays and their overlaps. The source law fixes
p and q. The complete memory theorem fixes the transient and noisy response
once the feedback service, initial memory and clock are given. None of these
steps reads the geometric pixel coordinate P.

The [same-metric Lean attachment](../../Lean/Screen/SeamCurrentPhysicalMetricAttachment.lean)
derives seam and vertex lengths from P only after a joint physical
identification is supplied. Its structure explicitly requires the port-dual
cell area, UV area element, response frame and propagation action to share
one metric. It supplies neither the feedback rate p(P), q(P), nor a family
attenuator. This theorem is useful for a physical spatial scale after its
attachment is proved; inserting that metric scale into the memory rates
would add a separate dynamics law.

The [source-selection construction](../interacting_extension/source_selection/PROOF.md),
section 6, also separates attempt count, record capacity and P, and proves
fixed-prefix independence from unused tape capacity. Its private counterpart
in `pro/Particles/05/oph_interacting_extension/source_selection/PROOF.md`
has the same boundary. The private
`pro/Particles/04/RESULTS.md` classifies its coherent Fourier readout and
mean-count-to-phase assignment as extra premises. Those constructions
supply no hidden P-to-controller selector.

The minimal missing map specifies how one physical capacity changes the
actual initial memory, event instruments, pulse gain and duration, and the
field observable, on the same source. A policy that changes only the
retained attempt count is excluded as a route to measured lepton ratios
through Y_n. The complete C-duty inverse can instead realize prescribed
positive coefficients by changing its probes and gain; that is a
[conditional implementation](../mcpr_native_memory_bridge/PROOF.md), with
selection of those controls as an additional dynamics premise.

## Reproduction

```sh
python3 -m pytest -q model/particle_calibration_selection/test_capacity_response.py
```

The tests reconstruct all native rows, differentiate the independent
12-address channel, check exact short-horizon coefficients, verify the
invariant inequalities and strict spectral constants, and calculate the
complete finite-step density. Experimental comparisons are separate
arguments to `comparison`, after the source bounds have been evaluated.
