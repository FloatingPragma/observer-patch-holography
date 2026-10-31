# Complete finite loop-flux response

This package puts an external Peierls loop phase on the same four-site
periodic Weyl/contact field used by the
[particle response](../particle_response/README.md). It computes conserved
number currents, a nonzero equilibrium loop response and its full spectral
curvature. All 256 Fock states and every spectral weight are retained. This
supplies finite current and response identities for the stated field; the
loop phase and field action are supplied observer operations.

No source pulse execution is added here. The calculation does not identify
a dynamical photon, laboratory electromagnetic flux, Lorentz spin, physical
clock or thermodynamic conductivity. The three-dimensional and physical
identification requirements of the particle programme remain separate.

## Same action and an actual loop

For four sites with two CAR modes each, let

```
T_x(phi) = exp(i w_x phi) (i c/2) a_x^* sigma_z a_(x+1),
H(phi) = sum_x (T_x + T_x^*) + U sum_x n_x0 n_x1,
K(phi) = H(phi) - (U/2) N,       sum_x w_x = 1.
```

Indices wrap modulo four. Here `c=1/8`, `U=1`, and the reference is the
unique ground of the **complete** `K(phi)` whenever its gap is numerically
resolved. It is not the faithful Gibbs reference used by other interfaces.
At zero flux this matrix is exactly the existing periodic field matrix.
`phi` is dimensionless, and all energy derivatives use the supplied
Hamiltonian model units.

The oriented number current is `j_x = i(T_x-T_x^*)`. The CAR relations give,
without choosing a state,

```
d n_x/dt = i[K,n_x] = j_(x-1) - j_x,
J = dK/dphi = sum_x w_x j_x,
D = d²K/dphi² = -sum_x w_x² (T_x+T_x^*).
```

Thus current conservation holds for every state and untouched reference.
In a stationary state all four bond expectations coincide, so their common
value equals `<J>`. This is the response conjugate to the declared phase;
converting it to amperes needs a physical charge, clock and flux convention.

Redistribute the phase by replacing `w` with `v`, both of total weight one.
Choose `chi_0=0` and `chi_(x+1)=chi_x+(w_x-v_x) phi` and let
`G=exp(i sum_x chi_x n_x)`. Direct substitution gives
`K_v=G K_w G*`. The total loop phase cannot be removed by this transformation.
The full spectrum is `2*pi` periodic. Because `G` depends on `phi`, the
individual derivative operators and diamagnetic/paramagnetic contributions
are generally **not** separately related by ordinary conjugation. Their
combined energy derivative and curvature are gauge independent.

## Closing the finite response premise

Let `|0>` be a normalized nondegenerate ground state and `E_n,|n>` the
complete eigensystem. Differentiate `K|0>=E_0|0>` in a local phase convention
with `<0|0'>=0`. Projecting onto every orthogonal eigenvector gives
`<n|0'>=<n|J|0>/(E_0-E_n)`. Differentiating the energy once and twice gives

```
E_0'  = <0|J|0>,
E_0'' = <0|D|0> - 2 sum_(n != 0) |<n|J|0>|²/(E_n-E_0).
```

The ground's isolation, finite Hermitian action and differentiable phase
law are the exact premises. No infinite-volume limit, relaxation-time
ansatz, missing-current assumption or fitted velocity enters. Degenerate
excited levels cause no problem: their equal-denominator weights sum over
the full eigenspace. A degenerate or numerically unresolved **ground** is
rejected before division; the code does not delete it or invent a gap.

The complete weights also satisfy

```
sum_n |<n|J|0>|² = <0|J²|0>,
sum_n (E_n-E_0)|<n|J|0>|² = (1/2)<0|[J,[K,J]]|0>.
```

These follow by inserting the identity on the complete carrier. The
producer checks both identities and retains all 256 gap/weight entries,
including the elastic term and numerically tiny dark-channel amplitudes.
No visibility threshold enters any sum. The JSON rounds numbers for
portable display/replay; the in-memory arrays and computations retain the
unrounded values. Replay compares complete spectral-subspace weights
`<0|J P_E J|0>` and multiplicities, because individual eigenvector residues
inside a degenerate eigenspace depend on the eigensolver's arbitrary basis.
Every original eigenvalue and weight remains in the receipt. Only this
comparison projection groups levels within `2e-10`; an intergroup separation
below `2e-8` is rejected as numerically ambiguous. These are numerical replay
rules, not a proof of exact degeneracy. The retained fixtures' distinct gaps
are separated by at least `2.7e-6`.

Energy comparisons allow two eleven-decimal serialization units; group-weight
comparisons allow the sum of the two lists' rounding budgets plus two units
for floating reconstruction. No bound on individual basis-dependent residues
is relaxed. Nonspectral fields, aggregate current/curvature, source hashes
and claim flags remain exact comparisons, with exact JSON types and finite
values required throughout. The CLI rejects duplicate keys at every depth.
Tests rotate complete degenerate
eigenspaces by random complex unitaries and reject shifted interlevel weight,
changed energies or multiplicities, missing modes, nonfinite weights and
scope promotion. The CLI uses this same verifier. Exact identities are distinguished from floating
numerical evaluations; no outward numerical error certificate is claimed.

## Results and controls

At zero flux, the ground has number four, energy `-2.04517902015`, and a
resolved first excitation gap `0.01395766521`. Its diamagnetic contribution
is `0.00545072671`, the paramagnetic sum is `0.00643576067`, and the net
curvature is **`-0.00098503397`** in model units. The negative finite-size
curvature is allowed: zero flux is not assumed to minimize this finite
ring's ground energy. It is not a thermodynamic instability conclusion.

At `phi=0.4`, the stationary loop-current response is
**`-0.00038343508`** and curvature is `-0.00090614394`. Redistributing all
flux onto the winding bond leaves these values invariant. It changes the
diamagnetic term from `0.00546916870` to `0.02187667479`, offset by the
corresponding paramagnetic change. Dropping either term therefore gives
the wrong curvature and destroys this gauge check.

For the noninteracting control with `0<phi<2*pi` and uniform phase,
`alpha=phi/4` gives the complete half-filled ground energy
`E_0=-2c(cos(alpha)+sin(alpha))`. Its curvature is
`c(cos(alpha)+sin(alpha))/8`; at `phi=0.4` this is `0.01710683722`.
The zero-flux free control and the atomic `c=0,U=1` control both have
16-fold full-Fock ground degeneracy, retained explicitly. The atomic
current and diamagnetic operators are identically zero.

Tests independently construct Jordan–Wigner matrices, differentiate the
full matrix and ground energy with five-point differences, check local
continuity, gauge redistribution, flux periodicity, the complete
time-domain current correlator, and both analytic controls. An omitted
winding bond, either omitted curvature term, or discarded degenerate
ground cannot pass those tests. Finite current response is useful evidence
for the supplied charged field, but neither proves particle localization
under native repair nor identifies observed particle species.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q -W error model/flux_response
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.flux_response --check model/flux_response/results.json
```

`response.report()` produces `results.json`, binding the implementation,
field and observation helpers, and this derivation by SHA-256. To refresh
after a reviewed change, remove only that generated receipt and use
`--write model/flux_response/results.json`. The writer refuses to overwrite
an existing file. The calculation is bounded to four sites and three
interacting phase values; it is not an empirical particle comparison.
