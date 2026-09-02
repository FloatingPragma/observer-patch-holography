# Local-domain evidence archive

This package is a byte-exact mirror of the OPH-FPE local-domain receipt
family used by the spacetime, particle, and cosmology papers. The manifest
binds the compressed stage-one arrays; stages one through four; and the
defect-sector, clock-unit, classical-realization, matter-attachment, and
source-gap receipts.

Run the producer-independent archive check from the repository root:

    python3 evidence/local_domain/verify_local_domain_archive.py

The checker imports no simulator code. It verifies every manifest digest,
strict JSON decoding and schemas, the 2,304-event finite-domain identity,
semantic replay, and the recorded nonpromotion boundaries. It verifies
archive custody and theorem-level fields; it does not independently regenerate
the simulator's arrays or producers.
