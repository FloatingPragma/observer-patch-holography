# Exact energy Ward cells and their coarse-graining

The supplied complete-Fock Weyl/contact action now has an exact nearest-face
energy conservation column on a three-dimensional lattice, and that column
commutes with spatial and temporal block coarse-graining. This attaches an
actual conserved source observable to the **energy column** of the finite
Ward interface used by the Einstein contracts. It does not supply the other
three momentum columns, a symmetric Lorentz tensor, a physical clock/length,
or an inhabitant of the full Einstein premise register.

The local energy convention and action are the existing ones in
[full_fock_energy](../full_fock_energy/README.md):

```
h_x = U n_x0 n_x1 + (1/2) sum_(b incident x) H_b,
H = sum_x h_x,       J_xy = i[h_x,h_y].
```

Their continuity law holds for every state, number sector and untouched
reference; no Gibbs or stationary reference is needed. The new result is
not another proof that a finite Hamiltonian conserves its total energy.
It routes the actual range-two local currents to standard lattice faces,
forms exact finite spacetime cells, and proves that successive coarse cells
retain the same source conservation identity.

## Source-bound nearest-face routing

For each nonzero oriented `J_xy`, distribute it equally over every shortest
nearest-neighbor path from x to y. Along an edge use the path orientation;
reverse orientation negates the current. Each path's incidence is
`delta_x-delta_y`, since every intermediate vertex cancels. Averaging paths
therefore leaves the divergence exactly unchanged:

```
sum_i [F_(x,x+e_i) - F_(x-e_i,x)] = sum_y J_xy,
d h_x(t)/dt + div F_x(t) = 0.
```

Outside edges of the open box have zero flux. All paths stay in the box.
At distance two, a straight pair has one shortest path and a bent pair has
two, each with weight 1/2. The prescription treats coordinate permutations
and reflections equally; it does not choose a preferred axis ordering.
The half-bond local energy and this routing remain stated conventions.
Neither is claimed to be a unique physical stress prescription.

`ward.py` computes the same source action using exact rational Majorana
monomials. Their multiplication is the Clifford relation, with Python
unbounded integers and rational coefficients; no floating zero test is
used. The `2^(2m)` canonical monomials are a faithful basis of the complete
operator algebra on m CAR modes. Thus the calculation never projects to a
particle sector. On the actual `3x3x3`, `c=1/8,U=1` action it checks all 351
site pairs, 81 nonzero currents and 27 exact nearest-face Ward equations,
routing onto 54 spatial faces. This is a symbolic identity on all `2^54`
Fock states, not an execution or allocation of that many state vectors.
The bounded symbolic evaluator admits at most 64 sites and binary-rational
source coefficients with denominator at most `2^16`; the algebraic proof
applies to every finite rectangular box and finite real self-adjoint action
coefficients. Independent small full-Fock Jordan--Wigner matrices bind all
three hopping axes and the contact to the existing source action. Tests
also exhaust all 65,536 products of the dimer's Majorana basis.

## Exact time cells and the Ward column

For `t_n<t_(n+1)`, define the **integrated** oriented face current

```
I_(x,i)(n) = integral_(t_n)^(t_(n+1)) F_(x,x+e_i)(t) dt.
```

The fundamental theorem of calculus gives the operator equation

```
h_x(t_(n+1))-h_x(t_n)
  + sum_i [I_(x,i)(n)-I_(x-e_i,i)(n)] = 0.
```

Taking any normalized state's expectation, including an arbitrary untouched
reference, preserves the equality. On the standard four-step chart use
`T_00(n,x)=<h_x(t_n)>` and `T_i0(n,x)=<I_(x-e_i,i)(n)>` for i=1,2,3.
Then the forward-difference definition of `ddiv` in
`Lean/ObserverPatchHolography/EinsteinBranch/Tensor.lean` has
`ddiv(step,T)(n,x)(0)=0`. These are cell energy and integrated face flux;
no continuum density, metric factor or pointwise relativistic stress is
being inferred. The other columns are **unspecified**, not filled with
zeros and counted as physical conservation laws. On a finite box the
statement uses its interior chart and zero external face currents.

The independent complete 64-state, three-site spectral diagnostic evaluates
the integral analytically in the energy eigenbasis, retaining degenerate
frequencies through the sinc formula. At cell width 3/4 its operator Ward
residual is below `5e-15`. Replacing the integral by the endpoint current
times the width gives a nonzero operator residual about `0.0235723`.
These floating residuals are diagnostics; the exact operator identity above
supplies the all-time conclusion.

## Coarse-graining is a commuting operation

For any finite block map `b`, including ragged rectangular blocks, sum the
energies in a block and the oriented face fluxes crossing each block pair.
Every internal face appears with both signs and cancels. In incidence form,

```
D_coarse P_faces = P_sites D_fine.
```

Consequently the same finite-cell Ward equation holds for every block.
Two successive block maps give exactly the same energy and boundary flux
as their composition. Merging successive time cells sums their integrals
and telescopes the energy endpoints. The certificate checks all 27
rectangular width choices on the 27-site box; independent tests check
non-geometric, non-monotone block labels and their compositions.
`Lean/QFT/CoarseEnergyWard.lean` proves the generic finite-sum divergence,
continuity, block composition and time telescoping statements without
admissions or additional axioms. The concrete CAR action is bound by the
exact Python certificate, rather than represented as a Lean value.

## Controlled source approximation and finite current samples

This attachment reuses the already proved
[complete-field locality and source schedule](../full_fock_field/LOCAL_SIMULATION.md).
It does not introduce a second thermodynamic-limit claim. Write
`A=U+3c`. The bond norm is at most c, so `||h_x||<=A` and
`||J_xy||<=2A^2`. A fixed nearest face receives total absolute path weight
at most seven: one direct pair, two straight pairs, and eight bent pairs
with weight 1/2. Thus `||F_face||<=14A^2=:B` uniformly in ambient volume.
Its support lies in a radius-three cube, with at most 343 sites; local
energy has support at most seven sites in radius one.

Let `L(R,r,s,t)` be `local_tail_boundary_bound` for radius R, support radius
r, support size s and time t, and let `epsilon` be the existing ideal
Strang unitary error bound on the radius-R patch, uniform for `0<=s<=t`.
Then uniformly in all normalized states and references,

```
energy endpoint error <= A min(2, L(R,1,7,t)+2 epsilon),
integrated face error <= t B min(2, L(R,3,343,t)+2 epsilon).
```

Only Hamiltonian terms meeting the face support enter its derivative.
Their summed norms are bounded by `343(U+6c)`. Therefore
`||dF/dt||<=2 B 343(U+6c)`. Using m midpoint samples adds at most
`t^2 B 343(U+6c)/(2m)` to the integrated face error, by integrating the
Lipschitz bound in each sample cell. This is an outward rational bound,
not a claim that endpoint sampling is exact. Supplied endpoint measurement
error e_h adds e_h; supplied instantaneous face measurement error e_F adds
`t e_F`. A cell residual is bounded by twice the endpoint error plus six
times the integrated-face error. Block errors add over the stated measured
boundary faces and endpoints; there is no hidden volume-uniform claim for
an arbitrarily large extensive block.

The receipt evaluates `R=18,t=1/8,n=10^7` Strang steps and `m=10^12`
midpoint samples with exact rational arithmetic. That bound describes a
`37x37x37` source patch in a larger ambient lattice, separately from the
27-site exact identity check. The analytic identity and norm estimates apply
to both under the same action law. These are resource/error
bounds, **not executed source or measurement counts**. The source factor
compiler, arbitrary local input restriction and readout access are admitted
services. State preparation, finite-shot statistics, floating/hardware
errors, source export failures, controller energy and failure-space energy
require their own explicit accounts; none is silently assigned zero here.
The function accepts additional observable readout error bounds when
available. Its inherited tail domain requires `24 c t<=32` and `3<R<=4096`.

An exact conserved energy column stable under these coarse operations closes
this particular finite attachment. It does not derive Lorentz symmetry,
the full conserved symmetric stress tensor, geometric modular kernels,
area entropy, generalized-entropy stationarity, a Newton constant, a
continuum scaling limit or 3+1 Einstein dynamics.

```sh
python -m model.coarse_energy_ward --check model/coarse_energy_ward/results.json
python -m pytest -q -W error model/coarse_energy_ward
cd Lean && lake env lean QFT/CoarseEnergyWard.lean
```
