# Finite local simulations and capacity dependence

A finite simulation can approximate a complete local experiment without
simulating its entire ambient universe. The existing
[local field theorem](../full_fock_field/LOCAL_SIMULATION.md) establishes this
for the specified full-Fock Weyl/contact action at fixed spacing. The result
below adds imperfect preparation, differences in the local action and clock,
and imperfect instruments. It gives sufficient conditions for comparing
different capacities. It does not derive those conditions from OPH's
capacity law or identify a particle mass.

## The capacity coordinates are different

The [native record construction](../observer/CAPACITY.md) uses
`N_record=ln M0`, where `M0` is the largest correctable set of reachable
public records for a specified continuation family. The proposed global
closure adds `N_global=ln D=ln M0(U_N)` on the same physical carrier.
Neither is a number of lattice sites, a Fock occupation number, the number
of repeated trials, a scalar occupation cap, or a gauge-flux cutoff.

The [capacity source family](../../code/capacity_readback/complete_packet_capacity_lift.py)
has exact all-rung capacity formulas. Its universal membership in the
complete source class and the physical selector remain open. Its status is
`NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT`. Those formulas supply
no map from global capacity to the field couplings, spacing, clock,
reference state or regulator cutoffs used here.

The [finite Higgs construction](../interacting_extension/higgs_join/PROOF.md)
likewise has specified couplings, graph, scalar cap and reference. An `N`
parameter is absent from its response calculation. This establishes what
that calculation takes as input; it is not a theorem that physical masses
are independent of global capacity.

## An error bound for the whole local experiment

Consider the existing Weyl/contact field on any ambient box containing the
required buffers. Its target Hamiltonian is

```
H(c,U) = c H_hop + U H_contact,     c,U >= 0.
```

Fix a chronological list of target durations `t_i`, an observer region, and
instruments on that region and retained memory. All Kraus operators are
even in the CAR algebra. Instruments may depend on earlier retained
records. All their outcomes, including failures, remain in the output.

The [existing planner](../full_fock_field/local_simulation.py) constructs
nested buffers backwards. At stage `i` its local operator error is
`delta_i`, uniform at every untouched-reference matrix level. The initial
simulation input is a state on the largest buffer together with the actual
observer memory and any retained reference. Assume its trace distance
from the restriction of the target joint state is at most `epsilon_0`.
Agreement of field marginals alone does not establish this premise when
the reference is correlated with the field.

In each finite buffer, let the simulator use coefficients `c',U'`, duration
`s_i`, and a split schedule with unitary operator error at most `eta_i`.
Let its complete instrument differ from the target instrument in
half-diamond norm by at most `zeta_i`. This must hold uniformly over all
incoming states, references and earlier records. A successful-branch
average is insufficient.

Write `V_i` and `E_i` for the buffer's site and spatial-edge counts. The
actual [source term catalogue](../full_fock_field/field.py) has two hopping
terms per spatial edge, each with norm `c/2`, and one contact projector per
site. Therefore

```
a_i = E_i |t_i c - s_i c'| + V_i |t_i U - s_i U'|
    >= ||t_i H(c,U) - s_i H(c',U')||.
```

For Hermitian `A,B`, integrate the derivative of
`exp(-i(1-r)A) exp(-irB)` over `0<=r<=1`. Its norm is bounded by
`||A-B||`, giving `||exp(-iA)-exp(-iB)||<=||A-B||` without assuming
commutation. The half-diamond distance between the associated unitary
channels is at most this operator difference. Thus `a_i` bounds both
coupling and clock errors. In particular, rescaling both couplings by
`kappa>0` and every duration by `1/kappa` gives `a_i=0` exactly.

Replace the target evolution by its finite buffer, then its generator by
the simulator's, then its exact evolution by the split schedule, and
finally its instrument. Telescope the compositions. Completely positive
trace-preserving maps contract trace distance; the corresponding dual
maps are complete contractions. The nested buffers make each locality
replacement valid on the whole support of the later approximate effect.
Consequently the joint final state of local field, memory and complete
record satisfies

```
D(final_target, final_simulator)
    <= min(1, epsilon_0 + sum_i (delta_i/2 + eta_i + a_i + zeta_i)).
```

The total variation distance of the classical record marginal obeys the
same bound. There is no outcome-count multiplier or division by a success
probability. Feedback through earlier records is part of the complete
instrument; taking a uniform error over its classical controls preserves
the bound. The executable helper fixes the durations and coefficients at
each stage. An extension with record-controlled durations needs bounds
uniform over those branches as well.

`local_history_budget` reuses the target planner's buffers and step counts,
but recomputes the split error using the simulator's actual coefficients
and durations. Its `epsilon` is a planning budget. The total can exceed it
when the supplied preparation, instrument or action errors are large.
`certificate()` gives a rational example with radius four and total error
below `1/100`, including these extra errors. It describes a supplied
schedule, not an executed large Fock-state array.
The helper inherits the planner's finite domain: radii through 4096 and
`48c|t_i|<=64`, with time subdivision required outside the latter range.
It rejects a request when its admitted radii cannot meet the budget.

## What would make the approximation uniform in global capacity?

Suppose a source theorem identified a common local algebra and clock for
a family indexed by global capacity. Suppose also that its local couplings
converged, remained uniformly bounded, and its initial joint restrictions
and instruments converged in the norms above. For every fixed finite
experiment and error tolerance, choose the nested buffers and schedules
using uniform coefficient bounds. Their sizes depend on the local
experiment and tolerance, not on ambient volume. With those finite buffers
fixed, the preparation, Duhamel and instrument errors tend to zero.
This proves convergence of complete finite local experiments under the
stated premises. The source theorem and its convergence rates are still
required before a particular manageable simulation certifies a particular
physical capacity.

Volume independence alone cannot replace these premises. A simple
counterexample uses zero hopping and alternates the local contact
coefficient between zero and one as a formal family index grows. Give the
observer a declared two-mode compensating reference with zero Hamiltonian,
and prepare

```
|psi> = (|pair,empty_reference> + |empty,pair_reference>)/sqrt(2).
```

Both terms have total number two. The projector onto `psi` commutes with
total number and is an even CAR observable. After time one its probability
is `p_+=(1+cos U)/2`. Disconnected ambient sites do not affect it, yet the
two subsequences differ by `(1-cos 1)/2 > 0.2298`. Alternating Taylor bounds
give the rational enclosure checked by `coefficient_counterexample`.
The test includes the reference explicitly and retains both measurement
outcomes. This is a counterexample within the declared even-CAR model
class; neither alternate law is asserted to be a complete admissible OPH
universe.

## A local record theorem is not yet a mass theorem

Convergence at each fixed time does not imply a uniform long-time error,
convergence of every spectral gap, or an isolated relativistic particle
pole. For instance, `cos(t/k)` converges to one on every fixed finite time
interval while each finite `k` has a nonzero frequency. A physical mass
identification additionally needs the source-selected reference, common
energy-momentum dictionary, physical clock, the relevant pole or spectral
sector with controlled errors, and the appropriate volume and regulator
limits. Charged probes need a consistently calibrated reference as well.

Reading the simulator specification can provide exact laws and enable
analytic spectral calculations. It is not restricted by a bound on what
a finite detector sample alone can resolve. But that access does not turn
an absent capacity-to-action or state map into a derived one. No mass,
mass ratio, physical capacity value, or capacity independence of the
observed particle spectrum is asserted by this error budget.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/particle_identification/test_capacity.py
```
