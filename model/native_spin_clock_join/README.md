# Selection-field and exchange clocks on one prepared state

The finite ordered-spin preparation and its zero-field exchange dynamics
have a quantitative common-state comparison. Every native outcome, including
failures and untouched reference systems, is retained. The bound applies to
growing finite observation times; its field requirement is explicit. It
does not replace a finite field-selected ground by an infinite vacuum.

## Complete finite histories

On a finite bipartite graph with N spins, use the existing source operators

```
H_J = J sum_edges (S_x.S_y-1/4),
Z = sum_x epsilon_x S_x^1,       H_B = H_J-B Z,       B>=0.
```

The selection signs are those used by the
[native ground preparation](../native_ordered_ground/README.md). Since
`||Z||<=N/2`, Duhamel gives
`||exp(-itH_B)-exp(-itH_J)||<=B N |t|/2`.
For a coherently controlled evolution, take the largest absolute duration
over the controller blocks. This preserves the operator bound with an
arbitrary spectator reference. A common scalar shift may be compared by
applying the same compensating controller phase to both clocks.
Independently re-centering them is allowed only when the relative controlled
phase is explicitly compensated; a clock interferometer can detect a
branch-dependent phase.

Interleave these evolution slots with the same complete quantum instruments.
Their records, zero outcomes and all failure branches are included in the
output. For each slot use an upper bound t_k on its absolute duration, valid
on every adaptive branch. Trace-norm contractivity and a telescoping sum give

```
D(output_B,output_J) <= min(1, B N sum_k t_k/2),                 (1)
D(rho,sigma) = ||rho-sigma||_1/2.
```

The output can include the entire finite field and all records and reference
systems. A conjugation counts its forward and backward evolutions separately.
Unbounded stopping times require their own tail bound. The input state is the
same for both clocks, including when it is a prepared H_B ground.

An input half-trace error delta, a product-formula error eta and a complete
native-dilation error nu add to (1). The existing
[integrated-port construction](../native_spin_observer/README.md) gives
`nu<=29 P 2^-b` for P pulses at b fractional bits. This uses its specified
coherent load/export service and completed-gate repair schedule. P includes
the actual preparation, controlled evolutions, every readout and its failure
records. The comparison constructs no unspecified controller and executes
no additional large native tape.

## Local histories with a volume-independent cutoff

Let the interventions and final accessible field lie on s sites in a cube
of radius r. There are m interventions by time T, each with a local
Stinespring implementation and arbitrary memory and reference systems. The
initial state on the common infinite algebra is identical for both clocks.
For R>r the existing spin commutator-chain estimate is

```
B_R(T) = min(2,6 J s (2R+1)^2 T 2^(-(R-r)) (exp(48JT)-1)).
```

The same estimate holds for H_B. In the interaction picture of its onsite
field, every bond keeps its support and norm J; the time-ordered
commutator-chain estimate uses only those norms. Local onsite rotations
preserve the intervention supports and operator norms.

For each clock the complete-history cutoff has unhalved trace error
`(2m+1)B_R(T)`. The radius-R models start from the same field/reference
marginal. Applying (1) between them and the triangle inequality yields

```
D(accessible output_B,accessible output_J)
 <= min(1,(2m+1)B_R(T) + B (2R+1)^3 T/2).                    (2)
```

The finite-history timing here is ordinary evolution from 0 to T. More
general forward/backward or controlled schedules must use their actual
absolute time budget and a proved local-intervention representation.
In particular a global controlled-energy measurement on a box is not a
local intervention on a smaller fixed region. Such a measurement uses (1)
with its full box size, without silently substituting (2).

For any rational T and target epsilon>0, `local_clock_plan` chooses R with
`(2m+1)B_R(T)<=epsilon/2`, and then the positive rational field
`B<=epsilon/((2R+1)^3 T)`. Both clocks therefore agree within epsilon.
The radius calculation uses `exp(48JT)<=2^ceil(96JT)` and integer logarithmic
inequalities, without allocating an exponentially large integer. It works
for arbitrary finite T. This also supplies a field choice for a sequence
of increasingly long experiments. If a packet has success bounded below by
`c p^4`, assigning clock error `o(p^4)` controls its unnormalized readout at
that scale. Preparation and native errors must satisfy their own budgets.

This construction specifies how small B must be for the clock comparison.
It provides no effective box-size or field-size modulus for convergence of
the finite H_B ground to a particular infinite ordered-state marginal.
The ordered-state thermodynamic limit retains its stated order of limits;
an arbitrary diagonal sequence B(L) is not identified with that state.

## Finite ground energy and a complete exact control

Let Omega_B be an exact finite H_B ground and E_J the finite H_J ground
energy. The variational principle and `||B Z||<=BN/2` imply

```
0 <= <Omega_B,(H_J-E_J)Omega_B> <= BN,
Var_Omega_B(H_J) = B^2 Var_Omega_B(Z) <= (BN/2)^2.             (3)
```

The equality follows from `H_J Omega_B=E_B Omega_B+B Z Omega_B`.
For a prepared state within half-trace delta of this ground, the energy
upper bound is `BN+J |edges| delta`: H_J is a sum of negative singlet
projectors and has spectral range at most `J |edges|`. These are finite
energy and stationarity estimates. A positive-energy ground spectral measure
centered at the prepared state does not follow from them.

An exact two-spin control uses J=1 and B=3/8. If s is the singlet and
`t=(|11>-|00>)/sqrt(2)`, its selected ground is `(3s+t)/sqrt(10)` with
H_B energy -9/8. The complete H_J energy measurement has outcomes

| H_J excitation above its own ground | Probability |
| --- | --- |
| 0 | 9/10 |
| 1 | 1/10 |

Thus its H_J excitation mean is 1/10 and its variance is 9/100. The state
is not H_J-stationary. Evolve for time pi, undo this state's preparation,
and retain both return outcomes. H_B gives `(1,0)`; H_J gives
`(16/25,9/25)`, with complete-record total variation 9/25. The full state
half-trace distance is 3/5. All matrices, projector completeness, weights
and the selection-field torque are checked exactly. The finite H_B gap is
9/8 while the finite H_J gap is 1; neither number is a relativistic rest
mass. The preparation and inverse are within the existing finite native
control construction, with their actual compilation errors charged in (1).

## Energy currents of the same spin action

Assign half of each bond energy to each endpoint and retain the declared
onsite selection energy:

```
e_x = (1/2) sum_y h_xy-B epsilon_x S_x^1,
H_B = sum_x e_x,              j_xy = i[e_x,e_y].
```

These currents are Hermitian and antisymmetric. Their support has graph
range at most two. Direct Heisenberg differentiation gives the exact local
and weighted conservation identities

```
d e_x/dt + sum_y j_xy = 0,
d (sum_x phi_x e_x)/dt = sum_(x<y) (phi_y-phi_x) j_xy.          (4)
```

For two adjacent exchange bonds,
`i[h_xy,h_yz]=J^2 S_x.(S_y cross S_z)`. These are the energy currents of
the same ordered-spin action and clock used in the readout comparison.
For a bipartite staggered field, the pair currents equal their zero-field
expressions: the additional bond terms are proportional to
`epsilon_x+epsilon_y` and cancel. The energy densities and Hamiltonians
nevertheless differ. Tests check the nonzero torque obtained by evolving
the zero-field density with the selected-field Hamiltonian.

An interleaved complete instrument contributes the energy jump
`sum_r Tr(K_r rho K_r^* H)-Tr(rho H)`, retaining every r. Equation (4)
applies between instruments and does not assert conservation through their
actions or construct an apparatus energy operator. A continuum stress tensor
requires scaling control of these current matrix elements, spatial-momentum
Ward identities and relativistic covariance. Neither a geometric stress
tensor nor Einstein dynamics is supplied by (4).

Run `python -m model.native_spin_clock_join --check` and
`python -W error -m pytest -q model/native_spin_clock_join`. The receipt
contains exact arithmetic and analytic resource certificates. It does not
claim a native pulse replay, Lorentz covariance, an interacting continuum,
a sharp mass shell, or an effective infinite-ground selection modulus.
