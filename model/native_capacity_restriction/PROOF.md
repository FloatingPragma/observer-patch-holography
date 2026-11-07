# Canonical restriction of the native source

This is a full-number-mixture corollary of the canonical source used in
`model/observer_time_fields/canonical_sector.py`. Its two-port weights and
falling-factorial helper are reused directly. The additional interface gives
the entire retained number law, exact nesting, state-distance bounds and
integer-population evaluation without constructing the global Hilbert
space. Hypergeometric sampling and trace-distance contraction are the
standard mathematical ingredients.

[Independent native products](PRODUCT.md) have a separate exact dynamic
restriction for repair/hold programs. That result preserves the weighted
H66 law without identifying it with the fair connected federation.

## The admitted source parameters

`model/completion/program_class_v2.json` defines a class of executable
programs. `model/native_program/interpreter.py` fixes their admission and
complete finite instruments. A member requires the following inputs.

| Input | Admitted domain and role |
| --- | --- |
| Local carrier | One H78 total-two integer cell, or its H66 binary shell. |
| Independent carrier factors | A positive finite number of H78 or H66 cells; each repair addresses one literal cell. This is distinct from the connected federation. |
| Connected binary carrier | The pinned geodesic `port_pair` federation, with `V_l=240*4^l` ports and `630*4^l` seams, integer `l>=0`, and one declared conserved global number `0<=m<=V_l`. |
| Local port weights | Twelve primitive integers in `1..32`, maximum/minimum at most four. They determine endpoint heat-bath coin probabilities. |
| Local seam weights | Thirty integers in `1..16`. They determine the uniform-menu weighted scheduler. |
| Local probabilities | Exact rational C/E mixture and repair probabilities with denominator at most 256, repair probability at least `1/16`; record-only exploration probability is at least `1/16`. |
| Federation law | Fair endpoint coins and equal seam weights. The interpreter does not admit an endpoint-biased federation lift. |
| Control | A finite literal schedule of repair and hold commands. Selectors are fixed, uniform, last-current, or the admitted exploratory record-only controller. |
| Initial state | A supplied source/reference operator. The interpreter checks the finite carrier and coefficient format; positivity and normalization of an actual preparation are separate requirements. |
| Memory | FIFO capacity `K=1..8`, with complete oldest-record export or explicit halt on overflow. The archive and all zero, hold, abort and stop records are retained. |
| Coherent services | The separately counted proper-code X/Y/Z/I/H/J source words and their supplied initialization/load/export services. These are not outputs of passive repair. |

For fixed inputs, each complete history has a uniquely defined instrument
block. At a repair attempt, its scalar probability is the product of the
repair probability, the record-selected seam probability, the declared law
weight and the endpoint coin probability. Its operator is the specified
native partial injection. Composing these operators with holds and the
actual FIFO update gives the source, reference and complete-record output.
Reading the simulator specification therefore determines the configured
law. It does not supply a missing map from a physical capacity to these
inputs.

The named `record_clock_member` uses the biased H66 source with relative
weights `(25,25,8,...,8)`, a one-slot export FIFO, a retained record inference
and a dimensionless modular clock protocol. No theorem identifies this
weighted source and its complete-record response with a growing fair
federation. In particular, the fair federation is a separate admitted
branch, not a proved capacity continuation of the biased clock member.

The resource envelope and local evaluator are also distinct. The class has
finite request bounds on attempts, description size, pulses and retained
bits. The local evaluator stops on federation levels above two, more than
64 independent factors, or its expanded-history limits. A closed-form
canonical restriction at a larger level is not an executed or
resource-admitted large-history program.

The FIFO parameter `K`, factor count, federation level, conserved number
`m`, local record capacity and proposed cosmological capacity `N` are
different quantities. The source contract provides no equation identifying
them. The spherical count `V_l` is not a derived three-dimensional physical
volume. A physical capacity family requires a specified map from `N` to
one source, state, observer, clock and comparison interface.

## Exact restriction and nesting

Fix integers `V>=1`, `0<=m<=V`. On the full binary occupation tensor product,
let

\[
 \rho_{V,m}=\binom Vm^{-1}\sum_{|S|=m}|S\rangle\langle S|.
\]

The unread native fair endpoint channel is unital and conserves total
occupation. Consequently every fixed seam and every source-independent
mixture of seams preserves this canonical source state. The complete
retained-record instrument also outputs records correlated with the source.
A fixed selector with positive weight on every seam of a connected graph
has a unique canonical stationary population law within each number sector.
This result does not select the conserved sector, a mixture of sectors, or
a finite exact preparation. Record-dependent control requires its own
stationary-state argument.

For the admitted federation set `V=V_l`. The formulas below also hold on
an abstract `V`-port binary tensor product; arbitrary integer `V` is a
combinatorial domain and does not create an additional admitted geometry.

Let `W` be any `v` distinct retained ports, `0<=v<=V`. The actual occupation
partial trace gives

\[
 \rho_W=\operatorname{Tr}_{W^c}\rho_{V,m}
       =\sum_{k=0}^{v}h_{V,m;v}(k)\rho_{v,k},\qquad
 h_{V,m;v}(k)=\frac{\binom vk\binom{V-v}{m-k}}{\binom Vm}.       \tag{1}
\]

An impossible binomial coefficient is zero. Every allowed exterior pattern
contributes the same diagonal coefficient. Counting its `m-k` occupied
ports proves (1), independently of the retained-port placement. Thus the
probability of one specified retained pattern with `k` ones is

\[
 \frac{(m)_k(V-m)_{v-k}}{(V)_v},                                \tag{2}
\]

where `(a)_r=a(a-1)...(a-r+1)` and `(a)_0=1`. This is exactly the owning
two-port formula when `v=2`. The implementation extends its finite degree
and population domain, preserving exact integer arithmetic.

For nested regions of sizes `u<=v<=V`, actual partial traces compose:

\[
 h_{V,m;u}(j)=\sum_{k=0}^{v}h_{V,m;v}(k)h_{v,k;u}(j).            \tag{3}
\]

One proof conditions on the intermediate occupation `k` and counts the
same global configurations. Equivalently, trace the two disjoint exterior
sets in either order. Equation (3) retains every intermediate number
sector. Selecting only one sector and renormalizing it is a different,
conditional operation.

The numerator of (2) uses at most `v` integer factors. The complete number
law is generated from its first supported value using

\[
 \frac{h(k+1)}{h(k)}
 =\frac{(v-k)(m-k)}{(k+1)(V-m-v+k+1)}.                          \tag{4}
\]

It requires `O(v)` rational arithmetic operations, with operand size
depending on `v log V`; it does not enumerate `binom(V,m)` states. A
description by number weights is complete for this permutation-invariant
diagonal state, since each sector is uniform. Emitting every basis entry
or arbitrary global observation can cost much more.

## What can be compatible across sizes

Orthogonality of distinct number sectors gives the exact half trace
distance to the canonical sharp-`q` coarse state:

\[
 D(\rho_W,\rho_{v,q})=1-h_{V,m;v}(q).                           \tag{5}
\]

For an arbitrary normalized density `sigma_q` supported in that sector,
the correct statement is `D(rho_W,sigma_q)>=1-h(q)`. The number measurement
proves this lower bound. Equality in (5) uses the matching canonical
conditional state, not an arbitrary density within the sector.

For `V>1`, the retained number has variance

\[
 \operatorname{Var}N_W
 =v\frac mV\left(1-\frac mV\right)\frac{V-v}{V-1}.             \tag{6}
\]

For `1<=v<V` and `0<m<V` this is strictly positive. Thus a nontrivial
uniform **global** fixed-number family cannot be exactly compatible with
one sharp number on each smaller carrier under these tensor partial
traces. This statement does not concern independent products whose local
factor numbers are fixed, and it does not exclude approximate convergence
on fixed local observable menus.
For `V=1`, the allowed number and every retained number are deterministic,
so the variance is zero.

As an admitted static control, take successive federation sizes `V=4v`
and the declared density `m/V=q/v=1/3`. For `v=240,960,3840`, the exact
distances in (5) lie respectively in
`(0.937004677738,0.937004677739)`,
`(0.968471222490,0.968471222491)` and
`(0.984231719274,0.984231719275)`.
These are full coarse-state comparisons; they do not contradict the
owning finite-window convergence results.

Equation (1) itself supplies the compatible static coarse state. There
is no need to choose one of its sectors. On the full-Fock mathematical
extension, the binomial sector mixture with weights
`binom(V,m)p^m(1-p)^(V-m)` is `diag(1-p,p)^tensor V` and is exactly
restriction-compatible for any declared `p`. The interpreter's declared
global fixed-number member is not silently replaced by this mixture.
The fair native law selects neither `p` nor a function `m(V)`.

## Finite-size and complete-readout bounds

Let `p=m/V`, and compare (1) with the Bernoulli product state
`beta_p^tensor v`. Sampling `v` labels with replacement from a population
with `m` occupied labels produces that product law. Couple it to sampling
without replacement until the first repeated label. The collision union
bound gives

\[
 D(\rho_W,\beta_p^{\otimes v})
 \le 1-\frac{(V)_v}{V^v}
 \le\min\{1,\ v(v-1)/(2V)\}.                                 \tag{7}
\]

Both densities are diagonal and uniform within each number sector, so
their full state distance equals the total variation distance of their
number laws. No independent factorization of the actual exterior and
retained region is asserted. The case `v=0` has distance zero.

The exact distance to local vacuum and its union bound are

\[
 D(\rho_W,|0\rangle\langle0|^{\otimes v})
 =1-\frac{(V-m)_v}{(V)_v}
 \le\min\{1,vm/V\}.                                          \tag{8}
\]

For a fixed window and fixed global occupation, the canonical restriction
therefore converges to vacuum. A homogeneous nonvacuum limit needs a
different occupation trajectory. A state conditioned on a moving occupied
core is a different preparation/readout; (8) makes no claim about it.

The analogous finite conditional-product measure with positive weights
of ratio at most `kappa` obeys
`P(n_i=1)<=kappa*m/(V-m+kappa*m)`. To see this, use
`m e_m=sum_{|S|=m-1} w_S sum_{j outside S} w_j` on the `V-1` weights
excluding `i`; its inner sum is at least `(V-m)w_min`.
For `0<m<V`, the odds at `i` are therefore at most `kappa*m/(V-m)`;
the cases `m=0,V` are deterministic and satisfy the bound directly.
The local biased twelve-port source is such a finite example. A biased
growing federation is outside the admitted source class, so this bound
does not supply one.

Any common complete quantum instrument acting on the retained region
contracts the distances in (7) and (8). It may retain output quantum
states, all records, zero outcomes, failures and further adaptive local
operations. The same statement holds after tensoring both inputs with an
identical independent reference. For correlated reference extensions a
joint input-distance bound is required; a source marginal alone does not
determine the joint output. If repairs cross the retained boundary, the
instrument acts on a larger dependency region. That region must be
included, or an independent locality error must be charged. Conditional
successful-output normalization does not inherit (7) without its
conditioning bound.

At the declared level `l=200`, `V=240*4^200` is a 408-bit integer.
With the supplied density `1/3`, all 13 number outcomes of a twelve-port
restriction are computed exactly and sum to one. Equation (7) bounds the
full state error by `66/V`, approximately `1.065e-121`. This demonstrates
large-population local evaluation. It supplies no relation between `V`
and `10^122` nats, no global history execution, no exact finite
preparation, and no physical clock, field or particle-mass assignment.

## Verification

`restriction.py` reuses the owning canonical helpers and compares the
general formula directly with the interpreter's actual sparse occupation
partial trace. The bounded replay enumerates all fixed-number states up
to eight ports, every retained size, every coarse canonical number and
every nested size. It checks normalization, moments, the exact distances
and both error bounds. These small abstract populations test the tensor
identity; they are not declared geodesic geometries. The source builder
checks the pinned level-zero port count separately.

The tests include every retained subset at a smaller size, an explicit
arbitrary-sector-state strict inequality, malformed domains, independent
reference extension, retained source/reference correlations, fair native
stationarity with every endpoint record slot retained, the two-port
owning interface, and the 408-bit integer evaluation. The source receipt
binds the implementations, proofs, tests and adopted source contract.
