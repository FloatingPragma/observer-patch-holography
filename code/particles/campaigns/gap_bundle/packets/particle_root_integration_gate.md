# Particle Root Integration Gate Packet

Integrate returned bundle packets and decide whether live prediction roots may
change.

Return `promote` only if the endpoint, RG/matching, interval, hadronic, and
spectrum-source gates all close. Otherwise return `keep_candidate` and list the
next coupled bundle revision.

Default decision: `keep_candidate`.

## Integration Result

`keep_candidate`.

All three first-wave packets returned claim-safe blockers rather than closure
proofs.  The compressed `P` trunk therefore remains candidate/audit metadata and
must not be promoted into live particle builders.

Blocking results:

- Electroweak root: no theorem-grade `Delta_Th(P)` yet; the missing pieces are
  `rho_had(s;P)`, OPH-derived RG/matching/scheme conversion, a certified
  `Delta_EW^src(P)` treatment, and interval-level endpoint certification.
- Spectrum source: only a reusable source-normalized trace-lift schema is
  available; charged normalization remains open, quarks remain selected-class,
  and neutrino PMNS comparison tension remains visible.
- QCD/Thomson backend: the current stable-channel hadron backend is not the
  Ward-projected electromagnetic spectral-measure export needed by the endpoint.

Next local target: implement the missing production spectral-measure contract
and the endpoint certificate interfaces before asking Chrome Pro for another
proof/audit pass.
