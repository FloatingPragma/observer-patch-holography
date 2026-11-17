# Method and reproducibility

The complete frozen experiment is `headline.py --source`; its exact-only
replay is `headline.py`. Both use PREREG.json and reject changed bound inputs.
The exhaustive all-seam run is intentionally retained. The short reviewer
command is:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B openwave/shared/magnetic_source_join/verify.py
```

Run from the oph-meta root, with NumPy and SciPy available. The exact matrix
checks use int64 with entries and products far below its range. The Y spin
component is represented by i times a real integer matrix, so no floating
zero test enters the algebra. The native source replay instead uses the
existing binary64 six-level compiler; its Frobenius/intertwiner errors remain
numerical diagnostics with the preregistered bounds.

`spin` constructs the field bilinear generators. `branches` calls the existing
full-domain repair and independently cross-checks all binary occupancies
against the pinned integer native endpoint arithmetic. `exact_audit` tests
complete effects, unitality, su(2), Casimir, every branch commutator and every
site's conditional spin continuity. `execute_source` uses the existing actual
source dilation, then executes source fermionic swaps and complete pointer
copies. It never substitutes an ideal unitary for the executed source map.

The source arrays retain amplitudes on every complete joint address, including
numerical leakage and otherwise empty/padded outputs. In an output address,
bits0..3 are field occupancy, bit4 is coin, bits5..8 are raw, and bits9..12
are the four pointer occupancy bits. Columns are all16 initial Fock masks.
No outcome is conditioned away, clipped, selected or renormalized. The ideal
compiler's complete proper-code exports contain all complement branches;
their effect sum is checked. This is an ideal controlled-source construction,
not a physical hardware demonstration or an independent calibration of hbar.

The short verifier independently recalculates all complete two-site and
crossing audits and representative twelve-port seams, checks the exhaustive
receipt's complete case inventory, and recomputes all four pointer laws from
the retained amplitudes. The general relabel-covariance proof is in PROOF.md;
the exhaustive all132 seam/law cases are reproduced by the longer command.

The changed-law control preserves the occupation/source/raw object and changes
only transported spin by a source-executed flavor swap. It has identical raw
effects and reverses the field-spin read. This supports the conservation claim
and apparatus sensitivity, not a magnetic field or a gyromagnetic prediction.
