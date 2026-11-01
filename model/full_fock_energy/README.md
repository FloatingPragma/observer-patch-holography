# Energy transport on the complete field algebra

The same [full-Fock Weyl/contact action](../full_fock_field/README.md) now has
local energies, oriented energy currents and a complete native readout energy
ledger. This closes the finite carrier mismatch between the previous
[two-particle energy calculation](../field_stress/README.md) and the
all-number reference and charged-excitation calculations. It supplies no
physical maintenance-to-rest-energy law or gravitational stress tensor.

## Exact local and regional continuity

Write the supplied finite action as `H=sum_b H_b+U sum_x n_x0 n_x1`, with
each unordered hopping bond counted once. Specify the local energy convention

```
h_x = U n_x0 n_x1 + (1/2) sum_(b incident x) H_b,
j_(x->y) = i[h_x,h_y].
```

Then `sum_x h_x=H`, each `h_x` and `j_(x->y)` is Hermitian,
`j_(x->y)=-j_(y->x)`, and

```
d h_x/ds = i[H,h_x] = -sum_y j_(x->y),
d H_R/ds = -sum_(x in R,y outside R) j_(x->y).
```

Here `s` is Hamiltonian model coordinate with supplied `hbar_model=1`.
These are operator identities, so they hold for every normalized state,
every occupation sector, coherences between sectors and untouched quantum
references. They require neither a Gibbs state nor stationarity.
The regional identity follows by cancelling the two orientations of every
internal current. Even CAR operators on disjoint mode sets commute, hence
these currents vanish beyond graph distance two. Distance-two terms can be
nonzero because adjacent bond contributions occur in both local energies;
omitting them gives an incorrect regional boundary current. The onsite
interaction changes the distance-one currents and cannot be dropped.

`energy.py` checks these statements with exact Gaussian-integer sparse
arithmetic, including every site subset, and independently verifies that
the local sum equals the executable full-Fock Hamiltonian. The reuse of
the existing integer-arithmetic checker changes its carrier, not its
arithmetic assumptions. Checked multiplication bounds exclude integer
overflow. The explicit certificate admits one through four sites, all
`2^(2V)` Fock states, binary-rational nonnegative `c,U`, local-energy
denominator at most `2^16` and bounded integer entries. It rejects
unrepresentable inputs instead of rounding them onto this domain.
The algebraic derivation applies to every finite box and finite real
coefficients of the stated self-adjoint action. The choice of half-bond
local energy is supplied; another partition can change local currents.

## Chemical potential and charged excitations

The faithful reference above uses `H` at zero chemical potential. The
half-filled particle benchmark instead uses `K=H-mu N`, with `mu=U/2`.
These generators must not be silently identified for a charge-changing
comparison. With `n_x=n_x0+n_x1`, define

```
k_x = h_x-mu n_x,          K = sum_x k_x = H-mu N,
j_N(x,y) = i([h_x,n_y]+[n_x,h_y]),
j_K(x,y) = i[k_x,k_y] = j_H(x,y)-mu j_N(x,y).
```

Since `[h_x,N]=0` and `[n_x,n_y]=0`, summing the first current over `y`
gives `-i[H,n_x]`; its two orientations cancel in a region. For the
declared half-bond convention it has only nearest-neighbor support.
Expanding the second commutator proves the grand-energy current identity
without a state or reference assumption. All subsets obey their complete
number and grand-energy continuity equations. The executable checks use
exact rational `mu` with the same bounded integer-arithmetic admission.

Although `i[K,h_x]=i[H,h_x]` for number-conserving local energies, charged
probes distinguish the two evolutions. For `[N,O]=q O`,

```
[K,O] = [H,O]-mu q O,
Delta K = Delta H-mu Delta N.
```

The second equality follows directly from `K=H-mu N` for any two normalized
states, and for a charge-resolved transition `Delta N=q`. The certificate
checks the commutator identity for every CAR creation and annihilation
operator. On the actual dimer at `U=1,c=1/8,mu=1/2`, the lowest particle
and hole gaps are both `Delta K=sqrt(17)/8-1/16`, approximately
`0.4528882032`. Their bare `H` costs are respectively `0.9528882032`
and `-0.0471117968`. An independent complete-matrix spectral check retains
both members of each degenerate doublet. Equality of grand-canonical gaps
therefore does not imply equal bare field-energy costs, nor does it select
a physical chemical potential or reservoir energy account.

## Readout changes energy

The binary occupation read at mode `a` has the two projectors `n_a` and
`1-n_a`. On the full input algebra its energy change is exactly

```
E_read^*(H)-H = -[n_a,[n_a,H]].
```

This follows by expanding `n_a H n_a+(1-n_a)H(1-n_a)-H` and using
`n_a^2=n_a`. Conservation of total particle number does not make the read
energy preserving: hopping coherence is disturbed. The onsite interaction
commutes with this read, but the Gibbs state of the interacting action
typically gains energy after it. Such backaction is recorded separately
from the closed Hamiltonian evolution.

`instrument.py` runs the literal source continuation followed by its actual
number-pointer word, preserving both unnormalized outcomes, the untouched
reference, and every first failed proper-code export. The integer identity
predicts the ideal read's change on the actual evolved success state.
The code compares it with the complete source pointer's two weighted energy
outputs. No postselection or success renormalization is used. Preparation
and read resources are separately counted.

For the ideal symmetric split, the existing complete-history bound is
`delta <= min(1,B^3 |s|^3/(3 n^2))`, where `B=c E+U V` and `n` is the
number of split steps. Since `||H||<=B`, the ideal unitary energy drift
is at most `2 B delta`; this is also valid with an arbitrary untouched
reference. The output gives that rational bound separately from numerical
residuals. No floating or hardware error is outward-certified here.

A failed proper-code export terminates with a different output space, on
which this field Hamiltonian is undefined. Its probability and complete
Kraus family remain retained. Its physical energy is **unassigned**, not
zero. The ledger therefore reports the success-sector field energy and
unresolved abort energy separately. It does not infer total conservation
for the field, source, controller, memory and apparatus together, nor a
thermodynamic cost from pulse counts. Those require an extended energy
observable on the actual services and failure spaces.

## Reproduce

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.full_fock_energy
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q -W error model/full_fock_energy
```

The executable prepares the complete sixteen-state `beta=1, mu=0` Gibbs
reference with the existing source compiler, then evolves and reads that
same preparation. All preparation aborts and the purifier remain accounted.
Independent tests assemble Jordan--Wigner matrices on all three bond axes,
compare a four-site complete algebra, test nonzero interaction and
distance-two currents, check the sign against time derivatives, and inject
source leakage that must fail the ideal-instrument comparison. Exact
number/grand-energy checks and independent charged-transition costs also
prevent confusing a particle gap with bare `H` energy. The finite
continuity identity does not establish Lorentz covariance, a physical
clock, defect protection or Einstein's geometric premises.
