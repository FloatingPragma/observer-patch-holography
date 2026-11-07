# Record capacity in the executable observer

For the native interpreter's `export_oldest` policy, changing the admitted
FIFO capacity from any integer 1 through 8 to another leaves the complete
repair/hold instrument unchanged after reversible record resegmentation. The
carrier, input, parameters, commands and evaluator limits are held fixed.
Every exported record is included in the output. The result holds with an
arbitrary reference, for adaptive commands and every retained zero block.

This derives a capacity dependence from the actual source: FIFO size changes
the placement and resource accounting of records, with zero error in the
complete native experimental channel. It supplies no relation between FIFO
slots and physical global capacity, nor between that capacity and a field
Hamiltonian.

## Operator proof

Write the complete chronological tape as `h = archive + memory`. At capacity
`K >= 1`, its storage is uniquely the suffix of length `min(K, len(h))` in
memory and the remaining prefix in the archive. Empty tape means empty
memory. For every nonempty tape, `memory[-1] = h[-1]`.

In `native_program.interpreter._choices`, uniform and fixed selectors are
independent of memory. The `last_current` and mixed observer selectors depend
only on whether the current of `memory[-1]` is positive. They therefore give
the same edge labels and weights at every admitted capacity. A hold record
has no positive current; it resets this particular feedback flag identically
at every capacity. The theorem does not replace the last record by the last
repair record.

The raw law, fair or biased coin, label decoder, native map and all its
weights in `_outcomes` depend on the common command, parameters, carrier,
repetition index and selector output. None depends on capacity. With export
enabled, the FIFO halt test never fires. `_commit` appends the same record to
`h` and makes precisely the new suffix/prefix split. Consequently, for each
raw outcome, the quantum linear map is identical at the two capacities and
the record difference is exactly resegmentation. This statement holds on
each matrix unit, so it holds on arbitrary operators and after tensoring an
arbitrary reference identity. It does not require diagonal states, normalized
nonzero branches or a choice of C/E coherence convention.

Induction over all command repetitions proves equality of the full recorded
instruments. Their output algebras are related by a permutation of classical
record labels, with identity on the quantum carrier. The resulting channel
distance is zero, including its half-diamond distance. An adaptive selector
in the declared grammar cannot distinguish the two representations.

The evaluator also preserves this equivalence. At each common prefix, the
number of raw outcomes depends on the same selector list; the step and branch
limits therefore stop both evaluations at the same prefix. The retained
`RESOURCE_STOP` record and unchanged quantum blocks coincide. The unevaluated
carrier guard is independent of FIFO size and has the same outcome on both
sides. The mathematical induction over the one-step kernel does not require
the expanded evaluator to materialize an arbitrary long tree.

## Boundaries that matter for a particle calculation

The literal program digest and resource counters are different. The result
does not assert equal physical costs for storing or exporting a record.
It compares the declared experimental outputs, preserving all tape content.
A specification read returning the FIFO size distinguishes the programs
exactly. Such reads are allowed; this equivalence concerns the repair/hold
instrument after omitting storage-placement metadata, not every observable
of the readable implementations.
An observer restricted to the FIFO, without the exported archive, has a
different observation map; that restriction is not covered. A scheduler
allowed to inspect older retained records would require another proof.

The `halt` policy gives an exact countercontrol: two hold attempts with
capacity one return `FIFO_TIMEOUT`, while capacity two completes. The two
status events have disjoint support, so their complete classical outputs
have total variation distance one. Capacity independence cannot be inferred
without specifying the memory policy.

Increasing a federation level, adding carrier factors, changing seam weights
or modifying a prepared state does not meet the fixed-program hypotheses.
Those are separate source-family questions. In particular, this result does
not establish the regulator limit, a relativistic field action, or equality
of particle masses at different physical global capacities.

## Executable verification

`native_capacity.py` reuses the admitted interpreter and recomputes resource
declarations for the capacity change. Its retained fixture compares every
raw output block across all eight capacities, including a reference-entangled
input, a hold reset and adaptive repair. Tests cover both native shells,
record laws, selector types, evaluator stops and the halt countercontrol.
These checks exercise the implementation; the all-operator argument is the
induction above.

```sh
python -m pytest model/particle_identification/test_native_capacity.py
python -m model.particle_identification.completion_certificate --check
```
