# QCD Thomson Backend Bundle Packet

Connect hadron production systematics to the hadronic spectral object needed by
the Ward-projected Thomson endpoint.

Scope update: GitHub issues #153 and #157 are closed as out-of-scope /
computationally blocked, not solved. The current local and Chrome-worker
environment cannot execute or promote this backend. Reopen only when a working
OPH hadron backend such as GLORB/Echosahedron emits the production export and
systematics.

Return `backend_contract` only after a production artifact for `rho_had(s;P)` or
the more primitive OPH spectral measure exists, plus quadrature/matching bounds
and continuum/volume/chiral/statistical budgets. Until then return
`out_of_scope` with the missing OPH backend primitive.

Surrogate hadron artifacts are not endpoint inputs.

## Worker Result

`constructive_spectral_measure_contract_emitted`.

The first worker pass found that the current stable-channel hadron backend is
not the QCD object needed by the Ward-projected Thomson endpoint.  That is not
accepted as a stopping point; it has been converted into a constructive
production schema:

```text
code/particles/hadron/ward_projected_spectral_measure.schema.json
```

The stable-channel backend remains a production/export surface for stable
hadron rows such as `pi_iso` and `N_iso`, while the Thomson endpoint now has a
separate Ward-projected electromagnetic current-current spectral-measure target.

The future endpoint backend must emit either `rho_had(s;P)` directly or a more
primitive production spectral measure whose pushforward defines it, with all
limiting and numerical budgets surfaced.

## Required Production Artifact

The missing artifact is:

```text
oph_qcd_ward_projected_hadronic_spectral_measure
```

It must be non-surrogate and unquenched, and it must include:

- publication-complete backend manifest provenance;
- Ward-projected `U(1)_Q` current normalization;
- finite-volume energy levels and electromagnetic current residues, or an
  equivalent primitive spectral measure;
- a defined pushforward to `rho_had(s;P)`;
- statistical, continuum, finite-volume, chiral, current-renormalization or
  matching, quadrature, and endpoint-remainder budgets.

## Claim Boundary

- Stable `pi_iso`/`N_iso` backend contract: selected-class production contract,
  still requiring real arrays and budgets for public rows.
- Stable-channel surrogate bridge: diagnostic only.
- `rho`/finite-volume resonance and Ward-projected spectral readout:
  missing production primitive.
- Ward-projected Thomson endpoint: theorem path blocked on backend spectral data.
- Imported Thomson/CODATA endpoint: compare-only metadata.
- Issues #153/#157: closed out-of-scope / computationally blocked, not solved.

## Smallest Missing Primitive

`production_ward_projected_hadronic_spectral_measure_export`.

Stable masses, surrogate correlators, or the existing stable-channel backend
cannot be promoted into this role without silently replacing a spectral-measure
problem by a stable-mass problem.

The next in-scope local patch may improve the schema, loader, or endpoint
contract. Backend execution itself is out of scope until OPH hadron hardware is
available; Chrome workers are not an execution substitute.
