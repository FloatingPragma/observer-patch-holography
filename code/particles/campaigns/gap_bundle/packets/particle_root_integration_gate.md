# Particle Root Integration Gate Packet

Integrate returned bundle packets and decide whether promoted prediction roots may
change.

Return `promote` only if the endpoint, RG/matching, interval, hadronic, and
spectrum-source gates all close. Otherwise return `keep_candidate` and list the
next coupled bundle revision.

Default decision: `keep_candidate`.

## Integration Result

`keep_candidate_with_constructive_next_artifacts`.

All three first-wave packets failed to close the full theorem gates, but the
campaign now treats obstruction-only output as unacceptable.  The useful result
is the constructive artifact set:

- `code/P_derivation/runtime/thomson_endpoint_contract_current.json`
- `code/particles/hadron/ward_projected_spectral_measure.schema.json`
- `code/particles/runs/hadron/ward_projected_spectral_measure_contract.json`

The compressed `P` trunk remains candidate/audit metadata and must not be
promoted into certified particle builders until those constructive artifacts are
populated and certified.

Next local target: populate or load the production spectral-measure contract and
wire it into the endpoint certificate interface before asking Chrome Pro for
another proof/audit pass.
