# Fixed-volume electric cutoff and complete records

## One parent Hamiltonian

Fix a finite spatial carrier, spacing, boundary condition and positive
electric coefficient. Retain every CAR occupation. The uncut Wilson
coordinate source of [FIELD_ROUTE](../particle_observer_outputs/FIELD_ROUTE.md)
has `H=E+B`, where E is the positive electric Laplacian on a finite product
of compact gauge groups and B is bounded. Its compact resolvent follows
from the finite-multiplicity Casimir spectrum tending to infinity and the
bounded perturbation. The Fock factor is finite dimensional. The spatial carrier and its field action are held fixed.

The nested native word projections P_l commute with E and converge strongly
to identity. The faithful native matrix coefficients and their conjugates
separate gauge-group points; their unital star algebra is dense in continuous
functions, hence in Haar L². Products give the same result for all links.
Because P_l commutes with E, `P_l psi -> psi` in the electric form norm for
every psi in D(sqrt(E)). Thus their union is a form core for H. This covers
both total-degree and product-per-link spaces, while their finite matrices
and states remain distinct. For general form convergence background, see
[Simon, J. Functional Analysis28 (1978),377–385](https://doi.org/10.1016/0022-1236(78)90094-0).
Here commutation with E supplies the form-core approximation directly.

The Gauss Haar projector P_G commutes with E, H and every P_l. Applying it
to the approximating vectors gives a form core on its invariant subspace.
All arguments below therefore apply separately to the full carrier or the
specified Gauss carrier, with their respective ground energies. They do
not equate these two minima.

## Ground energies and state subsequences

Let H_l be the compression of H to P_l. The variational principle gives
`E0(l)>=E0`, monotonically decreasing as the subspaces grow. Approximate a
normalized ground vector in the form norm by P_l and normalize; its energy
converges to E0. Hence `E0(l)->E0`.

Use the finite polynomial cooling construction with mean energy excess
epsilon_l tending to zero relative to E0(l). Embedded into the parent
space, its density rho_l has `Tr(rho_l H)<=E0(l)+epsilon_l`. Since
`H>=E-||B||`, the electric means have a common finite bound C. For the
finite-rank electric spectral projection R_M below M,
`Tr rho_l (I-R_M)<=C/M`. Finite-dimensional compactness and this uniform
tail bound imply trace-norm relative compactness. Every convergent
subsequence has a normalized positive limit. Lower semicontinuity of
the positive form H+||B|| gives mean energy at most E0, so that limit
is supported in the full parent ground eigenspace. This selects no
unique state when that eigenspace is degenerate.

Canonical purifications retain the conjugate reference. Their mean of
`E tensor I + I tensor E^T` is `2Tr(rho_l E)`, including when rho_l and E
do not commute. That joint positive operator also has compact resolvent.
The same compactness proof therefore gives joint trace-norm subsequences,
without discarding the reference. This joint positive operator is used
only to bound spectral tails; the reference's idle Hamiltonian remains zero.

## The same clock and complete apparatus

Extend the compressed evolution by the parent electric evolution on its
complement: `Htilde_l=E+P_l B P_l`. On P_l this is H_l. Bounded strong
convergence `P_l B P_l -> B` and the Duhamel formula give strong convergence
of propagators on each fixed vector, uniformly on compact time intervals.
The integrand is dominated by2||B||||psi|| and tends to zero; the argument
uses the common electric domain. This does not imply operator-norm or
unrestricted diamond-norm convergence.

For a fixed finite complete bounded Kraus menu K_j, compress each operator
to `S_j(l)=P_l K_j P_l`. Complete it with the retained failure operator
`S_fail(l)=[P_l-sum_j S_j(l)* S_j(l)]^(1/2)` on P_l. Its radicand is
positive by completeness. Each S_j converges strongly to K_j; the failure
operator converges strongly to zero. The finite Stinespring operators,
propagators and finite compositions therefore converge strongly, with
their adjoints where required. Together with joint input trace-norm
convergence, this gives trace-norm convergence of complete histories and
their subnormalized postmeasurement states, with the reference retained.
The failure label is never removed or renormalized. This is a fixed-state,
fixed bounded-apparatus assertion; detector menus depending arbitrarily
on cutoff are outside it.
On a Gauss-restricted carrier the menu must act within that carrier.
Charged probes require their explicitly compensated apparatus and the
corresponding enlarged physical domain before this argument is applied.

## A quantitative electric-complement bound

For a projection P commuting with E, suppose `PHP>=L P`,
`QHQ>=b Q`, b>L and `||PHQ||<=c`, Q=I−P. For u=||Ppsi||, v=||Qpsi||,

```
<psi,Hpsi> >= L u²+b v²−2cuv
           >= (L−c²/(b−L))(u²+v²).
```

The second inequality follows by completing the square; it is formalized
in [ElectricCutoffFloor.lean](../../Lean/QFT/ElectricCutoffFloor.lean).
Together with a retained trial upper bound U it gives
`L−c²/(b−L) <= E0(H) <= U`. An electric spectral cutoff at M has
b at least M−||B||. A word cutoff needs its own omitted electric floor;
its degree alone is not substituted for M.

## Executed native rotor comparison

The directional source is the same full four-mode CAR/flux rotor as
`field_route.exact_rotor(K)`, with flux shifts of six. The uncut flux set is
the integers. Its bounded matter part has norm at most4: each unitary
flux phase gives a four-row one-body symbol of norm at most2 and zero
trace, whose complete Fock lift has norm at most4. Finite flux matrices
are compressions of this parent. For |flux|>K, electric energy is at least
`(K+1)²/2`, so b=`(K+1)²/2−4` and c=4 are valid.

The code builds every retained state and exact reducing component directly
from the unchanged native onsite and hopping matrices. At K=6 it agrees
entrywise with the inherited208-state source in its real phase chart.
Exact rational positive-definiteness bisection bounds every occupation
sector minimum for K=6,12,18. The two-block inequality then encloses the
uncut sector minima and their adjacent-charge thresholds, with no missing
flux states silently set to zero. At K=18 the uncut ground is enclosed
approximately by[-2.249980,−2.160425] and either adjacent-charge threshold
by[1.017679,1.197349]. The machine receipt supplies exact rational endpoints.
The finite K=6 and K=12 ground enclosures are disjoint, certifying an actual
cutoff change rather than a rounding discrepancy.

These bounds concern one infinite-flux directional rotor at fixed volume.
Its finite K=6 detector residues belong to that stated cutoff; they are not
transferred to the uncut rotor by equality of labels. Physical spatial
refinement, chiral fields, observed particles and Einstein dynamics are
independent claims. No joint volume/spacing/electric-cutoff limit is proved.
