# Independent spatial winding review

This review reconstructs the carrier surface from the seam graph using a planar rotation embedding and half-edge face traversal. It does not call the producer's surface, path or phase-meter functions. The dual port boundaries are recovered as directed pentagonal walks between independently reconstructed triangular faces. General separating contours are directed Euler boundary walks, including disconnected and holed regions.

The actual native repair is the only transition law used. Each history begins with its actual occupation G and a zero integer current tape. The review retains that G without resetting it to a convenient background. For all 4,094 nonempty proper port subsets S, it checks

    sum_{v in S}(n_v − G_v) = −degree(boundary S).

Boundary coefficients come from actual dual contours rather than the producer incidence matrix. Complex phase samples are integrated independently on selected complete contours; the sample density is eleven subdivisions per absolute turn plus eleven. This avoids branch-cut aliasing for the tested integer tape values. Native noops retain their unchanged state and tape. A mutation that forgets a nonzero current is rejected. Complete one-load instruments retain both fair-coin outcomes. Separately, every one of the producer’s 456 complete prescribed-schedule histories is checked against actual source transitions, exact coin probabilities, its archived tape and all 4,094 separating-region identities; every archived port phase degree is recomputed with independent contour sampling. The review also realizes an actual degree-two current and its subsequent split into two integer units; the meter does not clip all nonzero values to a sign.

The exact general identity follows locally: a source current j changes occupation by −j at the first endpoint and +j at the second, while the oriented dual boundary degree changes by +j and −j respectively. Summing over an arbitrary region cancels every internal edge. This is a state-dependent **present relative charge**, measured from current n and retained initial G, rather than a claim that winding equals absolute occupation n on every state. The explicit zero-tape histories establish more than an arbitrary initial choice of tape.

The punctured sphere has a dual graph deformation retract. With all twelve port faces punctured its first homology has rank eleven; filling all but two removes ten independent boundary cycles and leaves rank one. This is the topology of the native carrier's spatial dual, rather than a disconnected phase loop assigned to a port. Integer tape produces an S1 map with a common endpoint phase at every dual vertex. A half-integer tape modification fails that gluing. No independently supplied texture or charge–texture code restriction is used.

For the declared integer native process and the declared spatial S1 readout, an elementary nonzero relative defect has an integer topological charge and cannot split into two half-integer closed-loop defects. This is a substantive constructive realization of the literal mathematical row mechanism. It does not identify an absolute occupation with relative charge, normalize physical electromagnetic charge, or establish that the interpolated diagnostic phase field is an independently prepared dynamical physical field. Official acceptance remains a separate platform result; this review changes no official row status. Mass, spin, force and architectural uniqueness are not additional premises of this finite topological statement.

`SPATIAL_INDEPENDENT_REVIEW.json` contains exact source pins, measured errors, native event counts, complete instrument checks and two explicit zero-tape histories. The producer and this reviewer are separate Codex workers; this is not an opposing-team Fable review.

Reproduce from the workspace root:

```sh
OPENBLAS_NUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/ROW1_NATIVE/reviews/review_spatial_independent.py
```
