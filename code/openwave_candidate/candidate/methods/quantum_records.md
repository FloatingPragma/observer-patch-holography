# Integer records and a monitored finite quantum instrument

This method joins a port/transit Hamiltonian to integer half-seam records
in a finite one-load sector. It defines an executable discrete protocol.
It does not establish a continuum theory or a complete OpenWave criterion.

For each oriented seam `e=(a,b)`, introduce two integer registers
`F[e,a], F[e,b]`. A basis state is `|z,F>`, where `z` is either a port or
a transit seam. With a supplied cutoff `K`, retain `|F[h]| <= K` and the
states reachable from `|initial_port,0>` by the half-hops below. There is
no reduction modulo an integer.

The half-hop from endpoint `x` to transit seam `e` increments `F[e,x]`
by one. Its reverse decrements it. The matrix element is `+1` at `a`
and `-1` at `b`. Hops exceeding the cutoff are omitted in both directions,
which gives a Hermitian hard-wall regulator. It does not produce a
stochastic coherent-hop rejection record.

```text
D |x,F> includes s(e,x) |e,F+unit(e,x)>;
Gamma = +1 on ports, -1 on transit seams;
H_m = D + m Gamma;
D Gamma + Gamma D = 0;
H_m^2 = D^2 + m^2 I.
```

Here `D` is the lifted operator on position and integer records. Its square
must not be replaced by the bare graph Laplacian on this larger space.
On the one-edge tree sector the implementation agrees with the bare
port/transit Hamiltonian after reordering its three states.

The conserved integer Gauss charges are

```text
G_x = 1[z=x] + sum_{e incident on x} F[e,x];
G_e = 1[z=e] - F[e,a] - F[e,b].
```

A forward half-hop removes one unit of endpoint occupancy and adds one
unit to its register, preserving `G_x`; transit occupancy and the negative
register term change together, preserving `G_e`. The reverse follows by
inversion. The mass term is diagonal. Hence all matrix elements of `H_m`
stay within the same exact integer Gauss sector.

For a native committed event on a port basis state, let `n` be its one-hot
port occupancy and use the literal integer rule

```text
d = n_a - n_b;
j = (d - (2 coin - 1) (abs(d) mod 2)) // 2;
n'_a = n_a - j; n'_b = n_b + j;
F'[e,a] = F[e,a] + j; F'[e,b] = F[e,b] - j.
```

The instrument measures position and all half-seam registers. It selects
seam, fair coin, and record law C or E uniformly. C retains `(coin,j)`;
E retains `(coin,d)`. E here is a pre-repair contrast record, distinct
from the separate seam-offset dynamics. The native event preserves every
`G_x` and `G_e`. An attempted target outside the finite basis is instead
retained as `boundary_rejected`, with unchanged state, zero applied current,
nonzero proposed current when applicable, and no emitted native record.
A transit source produces a retained `transit_noop`. Zero-current native
events also retain their coin and C/E records.

For source basis label `s`, event choice `alpha`, probability `w_alpha`,
and the declared target map `f_alpha`, the Kraus operators are

```text
K_(s,alpha) = sqrt(w_alpha) |f_alpha(s)><s|;
sum_(s,alpha) K_(s,alpha)^dagger K_(s,alpha) = I.
```

Distinct sources and event choices remain distinct outcomes even when their
target states coincide. Thus the instrument is completely positive and trace
preserving. Its additional source/record measurement is a supplied monitoring
choice, with no claim that it is the least disturbing native instrument.

One protocol step first applies `U=exp(-i H_m t)` and then applies the
measured instrument with supplied probability `p`:

```text
rho' = (1-p) U rho U^dagger + p C(U rho U^dagger).
```

The ledger retains a no-commit outcome of probability `1-p`. This is a discrete
clock in inverse hop-amplitude units. At `p=0` it gives the joint coherent
law. At `t=0,p=1`, from a port state away from the boundary, its port
population response is the native `I-L/(2M)`. No continuous-time,
fast-commit, Poisson-event, or physical-clock equivalence is asserted.

| Equation or observable | Portable implementation |
| --- | --- |
| Reachable bounded sector, half-hops | [quantum_records.py](../oph_openwave/quantum_records.py), `JointRecords.__post_init__`, `half_hops` |
| Integer Gauss charges and lifted generator | same file, `gauss`, `hamiltonian` |
| Native event, explicit rejections and C/E records | same file, `native_target`, `commit` |
| Discrete coherent/commit protocol | same file, `step` |
| Regulator comparison and scoped checks | same file, `cutoff_experiment`, `run` |
| Independent Kraus reconstruction, square and event record deltas | [test_quantum_records.py](../tests/test_quantum_records.py) |

The executed regulator experiment uses the oriented triangle
`((0,1),(1,2),(2,0))`, initial port 0, zero records, `m=0.23`, and eight
steps with `t=0.37`, `p=0.29`. Cutoffs 1 and 2 have dimensions 13 and 25.
The maximum difference in final position probabilities is `0.0208933517`.
The cumulative expected native rejection counts are `0.00573201863` and
`0.000000361110879`; traces differ from one by less than `6e-15`.
These counts concern monitored native attempts only. They do not bound
coherent hard-wall error. Two cutoffs provide a sensitivity measurement,
not a regulator-removal theorem.

Run `python -m oph_openwave.runner quantum_records` from the candidate directory.
The receipt records this experiment alongside scope, supplied choices and
local gates. Supplied choices include the graph, half-seam fields, one-load
preparation, mass, cutoff, hop clock, monitoring, transit noops and commit
probability. Many-load joint records, exchange statistics, removal of the
cutoff, physical calibration and platform acceptance remain open.
