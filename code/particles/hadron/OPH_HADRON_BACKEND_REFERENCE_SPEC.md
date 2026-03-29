# OPH Hadron Reference Backend and Export Contract v1

This document is a concrete reference backend specification for the current OPH
stable-channel hadron lane. It is designed to remove the remaining operational
underspecification called out in the RHMC/HMC handoff dossier and to make one
thing implementable right now:

- a real backend/export path that can honestly emit
  `backend_correlator_dump.production.json`
- on the frozen seeded `2+1`, QED-off geometry
- for the first promotable stable channels `pi_iso` and `N_iso`

It does **not** claim that the papers uniquely force this backend choice.
Rather, it fixes one honest, explicit, production-capable reference path so the
engine can be coded without hidden conventions.

## 1. Branch and scope

- branch: seeded `N_f = 2+1`, QED-off
- targets: `pi_iso`, `N_iso`
- export granularity: per ensemble, per cfg, per source, per channel
- exported channels:
  - `pi_iso`
  - `N_iso_direct`
  - `N_iso_exchange`
- repo-side nucleon combination rule:
  - `N_iso = N_iso_direct - N_iso_exchange`

This contract is **only** for the stable-channel branch. `rho` remains outside
scope and stays a separate finite-volume scattering problem.

## 2. Gauge action

Use the isotropic Wilson plaquette gauge action

\[
S_g(U;\beta)
=
\beta \sum_p \left(1 - \frac{1}{3}\Re\operatorname{Tr} U_p\right).
\]

Concrete conventions:

- gauge group: `SU(3)` fundamental representation
- isotropy: same coupling in spatial and temporal directions
- link storage order in backend: `mu-major`, then `site-major`
- gauge boundary conditions: periodic in all four directions

## 3. Fermion discretization

Use clover-improved Wilson fermions with a single explicit profile choice:

\[
D_f = D_W(am_f) + \frac{i}{4} c_{SW} \sigma_{\mu\nu} F_{\mu\nu},
\qquad c_{SW} = 1.0.
\]

Concrete conventions:

- light sector: degenerate `u/d` doublet at mass `am_l`
- strange sector: one flavor at mass `am_s`
- preconditioning: even-odd preconditioned operator for all solves
- spin basis: Euclidean Dirac gamma matrices with
  `gamma_4` Hermitian and `gamma_5 = gamma_1 gamma_2 gamma_3 gamma_4`
- color-spin index order in propagator objects: `color-major`, `spin-minor`

The choice `c_SW = 1.0` is deliberate: it is not the only possible production
choice, but it removes ambiguity and is straightforward to implement and audit.

## 4. RHMC/HMC factorization

Use the following exact determinant representation.

### 4.1 Light doublet

For the degenerate light doublet, use standard HMC with one pseudofermion field
for the positive Hermitian operator

\[
Q_l^\dagger Q_l,
\qquad
Q_l = \gamma_5 D_l.
\]

Pseudofermion action:

\[
S_{pf,l} = \phi_l^\dagger (Q_l^\dagger Q_l)^{-1} \phi_l.
\]

### 4.2 Strange flavor

For the strange flavor, use RHMC on

\[
(Q_s^\dagger Q_s)^{-1/2},
\qquad
Q_s = \gamma_5 D_s.
\]

with a minimax rational approximation of order `12` over the spectral interval

- `lambda_min = 1.0e-6`
- `lambda_max = 25.0`

The backend **must** persist the actual partial-fraction coefficients used in the
raw export manifest under `solvers.rhmc_strange.rational_coefficients`. Those
coefficients are part of execution provenance, not an implicit runtime choice.

No Hasenbusch splitting is used in the reference profile. That is conservative
but removes another source of ambiguity.

## 5. Molecular-dynamics integrator

Use a single-timescale second-order Omelyan integrator with

- `trajectory_length = 1.0`
- `n_steps = 16`
- `lambda = 0.1931833275037836`
- momentum refresh: fresh Gaussian momenta before every trajectory
- accept/reject: Metropolis test on the full Hamiltonian at end of trajectory

This choice is intentionally modest and explicit. The engine may later add a
more elaborate tuned profile, but the reference path does not depend on that.

## 6. Linear solvers

Use conjugate gradient on the even-odd preconditioned normal equations with
float64 arithmetic throughout.

- solver id: `cg_normal_eq_eo`
- force-solve residual tolerance: `1.0e-8`
- measurement-solve residual tolerance: `1.0e-12`
- max iterations: `20000`
- failure policy: abort the cfg measurement and record failure in manifest

If an implementation later upgrades to mixed precision or multigrid, that is a
new profile and must not silently reuse `oph_reference_backend_v1`.

## 7. Boundary conditions

- gauge links: periodic in x, y, z, t
- fermions: periodic in x, y, z
- fermions in time: anti-periodic

The backend must record these exact values in
`boundary_conditions.gauge` and `boundary_conditions.fermion`.

## 8. Source and sink construction

Use deterministic local point sources only.

- source coordinates: exactly the `src0` / `src1` coordinates already fixed by
  the receipt/payload contract
- gauge fixing: none
- source smearing: none
- sink smearing: none
- dilution: none
- stochastic noise: none

Every exported dataset must correspond to one exact source coordinate from the
receipt/payload contract. No source averaging is allowed inside the backend
export.

## 9. Stable-channel correlator definitions

### 9.1 Pion channel

Export `pi_iso(t)` as the real zero-momentum projected connected correlator

\[
C_\pi(t; s)
=
\sum_{\mathbf x}
\Re\operatorname{tr}_{c,\mathrm{spin}}
\left[\gamma_5 S_l(x;s)\gamma_5 S_l(s;x)\right].
\]

This is exactly the array the repo expects as `pi_iso` for one cfg and one
source.

### 9.2 Nucleon channel

Use the local proton interpolating field

\[
N_\alpha(x)
=
\epsilon^{abc}
\big(u^{aT}(x) C\gamma_5 d^b(x)\big)
 u^c_\alpha(x),
\]

with positive-parity projector

\[
P_+ = \frac{1 + \gamma_4}{2}.
\]

Export the projected zero-momentum direct and exchange pieces separately:

- `N_iso_direct(t)`
- `N_iso_exchange(t)`

with the exact convention that the repo performs the subtraction

\[
N_{\mathrm{iso}}(t) = N_{\mathrm{iso,direct}}(t) - N_{\mathrm{iso,exchange}}(t).
\]

The backend must **not** fold the minus sign into the exported exchange array.
The subtraction happens only once, in the repo pipeline.

## 10. Numeric export semantics

For every dataset:

- dtype: little-endian float64
- shape: `(T,)`
- units: lattice units
- content: raw projected correlator sequence in Euclidean time
- no effective masses
- no source averaging
- no cfg averaging
- no normalization beyond the operator definitions above

## 11. Raw export bundle format

The raw export bundle contains exactly two files:

- `backend_run_manifest.json`
- `correlators.h5`

The HDF5 dataset layout is canonical:

```text
/ensembles/{ensemble_id}/cfgs/{cfg_id}/sources/{src_id}/pi_iso
/ensembles/{ensemble_id}/cfgs/{cfg_id}/sources/{src_id}/N_iso_direct
/ensembles/{ensemble_id}/cfgs/{cfg_id}/sources/{src_id}/N_iso_exchange
```

The JSON manifest carries:

- backend/build/run provenance
- the fixed profile id
- exact physics / solver / integrator / BC choices
- cfg ids and `trajectory_stop`
- source ids and source coordinates
- the HDF5 dataset paths for each source/channel pair

The repo-side adapter then converts this raw bundle into the existing frozen
JSON production dump.

## 12. Required manifest provenance

The backend manifest must populate these top-level fields before any publication
attempt:

- `backend.family`
- `backend.name`
- `backend.version`
- `backend.git_commit`
- `backend.run_id`
- `backend.build_id`
- `backend.machine`
- `profile_id`
- `solvers.rhmc_strange.rational_coefficients`
- `integrator`
- `boundary_conditions`
- `sources`
- `contractions`

## 13. Current publication boundary

This reference backend spec closes the **exporter/backend-definition** gap. It
does not by itself fabricate production arrays, and it does not erase the
existing requirement for:

- real backend execution on the seeded family
- nonempty `cfg_source_corr_t` arrays for `pi_iso` and `N_iso`
- published statistical error fields
- published continuum / volume / chiral systematics

So this profile is the honest coding surface that upgrades the hadron engine
from “a JSON skeleton exists” to “an external RHMC/HMC code now has an exact,
fillable production export contract.”
