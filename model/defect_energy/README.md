# Repair cost and action of the same persistent memory defect

This packet composes with `model/defect_protection`, on its **same** binary
66-state source, complete native C/E branches, and one coherent orientation
qubit. It derives the action of that specified repair/memory law and tests
which repair-cost summaries retain it. No separate Fock action, spectrum,
momentum dictionary, or particle mass is introduced.

Run `python3 -m model.defect_energy --check` and
`python3 -m pytest model/defect_energy -q -W error`.
All receipt calculations are integer or rational; the largest Markov carrier
has 132 states. Word enumeration is bounded to eight attempts. The receipt
lists direct implementation bindings; `defect_protection` owns the wider
source-execution and interface audit. Expanded compiler words are counted,
not newly executed in this packet.

## Constructive repair-to-action identity

Write the unweighted source partial injection of raw branch `h` as
`M_h=sum_x |f_h(x)><x|`, with the sum only over its admitted inputs. Actual
ordered-particle transport gives a bit `p(x,h)` on each column. Holds have
zero bit. With `P_minus=(I-X)/2`, the specified memory update satisfies

```
X^p = exp(-i*pi*p*P_minus),
A_h = pi sum_x p(x,h) |f_h(x)><f_h(x)| tensor P_minus,
K_h = exp(-i*A_h) (M_h tensor I) / sqrt(2).
```

The output-supported action is well-defined because each raw branch is a
partial injection. The amplitude `1/sqrt(2)` belongs to the fair native coin;
uniform random seam selection contributes an additional `1/sqrt(30)`.
Zero and abort slots retain zero support. This is an identity of the
**subnormalized complete branch**, not a claim that repair is unitary.
The implementation checks all 7,920 C/E source columns, including 416
nonzero parity columns. Linearity gives source superpositions and arbitrary
untouched references without an extra assumption about input states.

Along an admitted source history, the memory gate is
`exp(-i*pi*sum_t p(x_t,h_t)*P_minus)`. Thus repair supplies an exact relative
action, in units where the gate angle is dimensionless. It is an operator
and depends on the source state reached before each update. A vertex
rechart changes an open-path parity by its endpoint bits; closed holonomy
and the experiments below are unchanged. No quantum phase is inferred from
classical topology alone: it comes from the **specified coherent memory
update**, whose inclusion is explicit in the joint carrier.

## A matched repair-cost experiment

Both following paths start and end with occupied pair `(0,11)`. Every arrow
is one actual moving native attempt; the other occupied port is transported
by the same law.

* Even path: move `11->5->1`; move `0->5->11->0`; move `1->5->11`.
* Odd path: move `0->5`, `11->0`, `5->11`; move `11->5->1->5->11`.

These histories have the **same multiset of complete
`(law,seam,coin,raw)` labels**. Consequently every additive scalar cost
`sum_t c(law,seam,coin,raw)` is identical, for arbitrary `c`. This includes
all the supplied native and memory-feedback compiler counts:

| Resource on either seven-move history | C | E |
| --- | ---: | ---: |
| Complete native attempts / movements | 7 / 7 | 7 / 7 |
| Native pulses including orientation feedback | 3,561,005 | 3,561,026 |
| Orientation feedback pulses | 27,908 | 27,908 |
| Identity-phase pulses (included in native total) | 1,820 | 1,820 |
| Fresh fair coin bits | 7 | 7 |
| Retained native tape bits | 21 | 28 |

Yet the logical orientation toggles **twice** on the even path and **once**
on the odd path. The repair-derived relative action is respectively
`2*pi` and `pi`; the final memory gates are `I` and `X`. Initial memory
`|0>` therefore gives opposite computational readouts. The full 128 coin
outcomes of each fixed seven-seam word are retained, with the displayed
branch carrying probability `1/128`; it is not silently renormalized. The
receipt gives each outcome's resource counts, including load/export and
complement slots. Source preparation and storage duration are supplied
services, not included in these counts. The final memory read uses one bit
and no basis-change pulse. A separate pair of three-move loops also supplies
the simpler equal-attempt/equal-movement counterexample.

This is a constraint on the **form** of a repair-cost energy derivation.
A scalar function of the total raw-branch cost cannot reproduce the full
relative action. The actual state-dependent operator action above does.
This does not exclude a repair-cost Hamiltonian on the full process, an
energy depending on more state information, or an ordered sequence of
different generators. Coarse operation counts discard precisely the
source-conditioned action that this example retains.

## Persistence and the continuing-source memory response

Every complete branch commutes with `I66 tensor X`; hence its two charge
sectors persist under the specified memory law. This is the same conserved
observable certified by `defect_protection`. The joint defect here is an
observer memory sector, not yet a source-localized particle.

For initial `I66/66 tensor rho`, fair coins, independent uniform seams, and
the continuing source traced only at the observation boundary, exact
propagation gives a memory channel

```
Lambda_n(rho) = (1-p_n) rho + p_n X rho X,
p_1 = 26/495,
p_2 = 733/7425.
```

The formula also holds with an arbitrary reference attached to `rho`.
The classical ordered-cover lift computes it exactly: source-diagonal
blocks remain diagonal, each carrying `rho` or `X rho X`. The independent
test pushes a Bell state through direct native label tracking and checks
the full rational Choi state, including coherences.

`Lambda_1(|0><0|)` has smaller eigenvalue `26/495`, so its trace distance
from any pure state is at least `26/495`. This lower-bounds the normalized
channel distance from every memory-only unitary by the same amount. Also,
the two-attempt Z multiplier exceeds the square of the one-attempt
multiplier by `398/245025`. Therefore the reduced memory channels are not
powers of a fixed memory-only channel: continuing source correlations
matter. Neither result excludes a Hamiltonian part plus noise, a Hamiltonian
on a larger system, or an energy for the full source/memory/record process.

The uniform source is stationary and has mean native movement count
`5n/33` in `n` attempts for every initial memory state, including the two
sharp X sectors. Persistence and maintenance counts are thus composed on
one carrier. The full process retains its growing record tape; stationarity
is asserted only for the reduced uniform source, not the archive.

## Remaining repair-cost-to-mass join

The exact action fixes an event's relative gate angle, not a unique energy
per unit duration: endpoint exponentials admit logarithms differing by
`2*pi` on spectral projectors, and an event duration has to come from the
same specified clock. Reading the simulator, memory, and scheduler
specifications is permitted; all quantities above come from those inputs.
If a fuller implementation fixes durations and within-event generators,
those are further readouts to bind, rather than inaccessible assumptions.

For a Lorentz-compatible completion, the same derived rest repair energy
`E0` and transport speed `c` must satisfy `m=E0/c^2`. The independent
transport check is then `d^2 E(p)/dp^2 at p=0 = c^2/E0`, with the momentum
and motion belonging to this same defect. `InternalEnergyInertia.lean`
proves inertial coupling conditional on the chosen proper-time action;
`ProperTimeInternalAction.lean` states the refinement and proper-time
conditions selecting that action. This packet supplies the finite native
repair/action side of the proposed join. It does not yet supply spatial
motion or a momentum for this observer-memory sector, nor discharge those
Lorentz/proper-time premises. Transferring a different field's dispersion
would not check `E=mc^2` for this defect. The matched-cost experiment tells
that future construction exactly which action information its repair-cost
functional must preserve.
