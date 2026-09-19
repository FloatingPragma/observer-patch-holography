# Canonical support-wiring diagnostic (#776)

The frozen [specification](SPECIFICATION.md) states all inputs, the event
convention, interval selection and limits. The complete receipt and trace live
in [`evidence/support_wiring_776/`](../../evidence/support_wiring_776/README.md).
No unpublished simulator module or #777 implementation is required for replay.
The [post-run control addendum](CONTROL_ADDENDUM.md) retains the q13 family's
existing one- and two-dimensional scientific controls beside the primary
comparison. The `readouts.py q13` command executes all three families; full
verification and fresh reproduction require every family and its complete tape.

Run the independent checks from the repository root:

```sh
python -m pip install -r requirements.txt
python -m pytest -q code/support_wiring/test_support_wiring.py
python code/support_wiring/verify.py
```

The verifier authenticates all files, reconstructs current writers, checks
every exact mean and refinement copy, computes every writer readback again,
and derives the order from these reads. It does not import either producer.
The record-metric tapes require exact integer payloads, array shapes and every
read offset; fractional values and unused appended reads cannot be silently
accepted. Historical comparisons require their complete law and result fields.
It counts comparable pairs using reverse descendant bitsets; the producer
uses forward ancestor bitsets. Exhaustive small-DAG tests check both against
a third, Boolean-matrix oracle. Kernel verification uses the full-field
`(T^n B)^T (T^n B)` Gram rather than the producer's local `T^(2n)` return.

Fresh reproduction (choose an empty output directory):

```sh
python code/support_wiring/experiment.py trace --output temp/support-wiring-reproduction
python code/support_wiring/experiment.py controls --output temp/support-wiring-reproduction
python code/support_wiring/readouts.py provenance --output temp/support-wiring-reproduction
python code/support_wiring/readouts.py q13 --output temp/support-wiring-reproduction
python code/support_wiring/assemble.py --output temp/support-wiring-reproduction --no-mirror
python code/support_wiring/verify.py --archive temp/support-wiring-reproduction
python code/support_wiring/check_reproduction.py temp/support-wiring-reproduction
```

Optional geometry recapture requires the public simulator revision recorded in
`geometry/geometry.json`:

```sh
python code/support_wiring/capture.py /path/to/oph-physics-sim
```

Capture checks the actual source bytes against their Git blobs, then reproduces
the original L3 fixture. Replay checks the complete neighbour incidence, port
slot uniqueness and each cell's geometric assignment objective. The endpoint
assignment can have tied optima; the pinned capture fixes their actual labels.
Antipodal label matches and global port-direction dot products are separate
statistics. Antipodal gluing is measured, not imposed.

## Why the confluence statement is limited

Let an attempted seam join scalar loads `a,b`. Its exact mean preserves their
sum and decreases the squared Euclidean norm by `(a-b)^2/2`. It is the
orthogonal projection onto the hyperplane where its endpoint loads agree.
The intersection of these hyperplanes is the space constant on each connected
component. Fixing each component total therefore specifies one common fixed
point. The verifier checks connectedness and conservation separately.

For a fixed finite sweep containing every seam, let `A` be the product of
these orthogonal projections. On the component-mean-zero subspace, every
nonzero vector has `||Ax|| < ||x||`: equality at the end would require equality
at every projection, hence agreement across every seam and thus a zero vector.
Compactness of the unit sphere gives `||A|| < 1` on this finite-dimensional
subspace. Repeated identical sweeps converge geometrically to the component
mean. This argument supplies no regulator-independent contraction bound.
It covers the declared forward and reversed cyclic schedules, and does not
assert convergence for unfair schedules, finite exact settlement or equality
of finite schedule endpoints. Four-sweep residuals and their endpoint
disagreement are measured, not substituted for a convergence theorem.

Refinement copying preserves the declared area-weighted conditional expectation
and the copied coordinate value. It does **not** preserve the unweighted sum
when one parent is replaced by four children. Conservation is asserted within
each repair level, not across a change of population.

The cumulative signed load starts with preparation and adds each exact repair
increment; telescoping makes it equal to the current port load. Rank-three
positions are floating approximations to that exact dyadic load readback.
They are snapshots after a simultaneous matching, not additional causal reads.

The pair of writes in a mean action is a declared event granularity: both have
the two old writers as parents and neither precedes the other. Parent IDs are
checked against actual current versions before they can enter the order. The
SHA256 chain authenticates the retained execution bytes and ordering under
the repository's custody; it is not an externally signed physical observation.
