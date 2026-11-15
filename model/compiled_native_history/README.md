# Compact exact native histories

This package compiles an admitted native repair program into one exact scalar
weight and partial source-address map per complete history. It reuses
`native_program.interpreter.execute` for control choices, FIFO export,
intrinsic timeout and host resource stops. All raw histories remain present,
including zero outcomes. The result can then act on any admitted exact
source/reference operator without recompiling the program.

This reduces host storage for repeated channel evaluations. It adds no
physical source operation, source-word control, particle interpretation or
clock calibration, and changes none of the declared physical resource counts.

## Why diagonal tagged inputs suffice here

Every native raw C/E Kraus operator has the form `sqrt(w) U P`: its domain is
a subset of occupation addresses and `U` is injective on that domain. Its
amplitude is the same positive scalar on every surviving address. Once a
complete history is fixed, its record-only controller choices and weights
are fixed too. Composition, holds and terminal identities preserve this form.

For a history with partial injection `f`, domain `D` and squared scalar `w`,
the complete reference-valued action is exactly

```text
|i,a><j,b|  ->  w 1[i in D] 1[j in D] |f(i),a><f(j),b|.
```

The tagged diagonal operator `sum_i |i,i><i,i|` reveals `D`, `f` and `w`
without materializing the off-diagonal Choi entries. This argument depends on
the native branch form. Tagged diagonal data do not determine an arbitrary
quantum channel: identity and dephasing agree on every diagonal
input. A coherent service inserted between native attempts would require a
different compiler and is not admitted by this package.

The underlying partial-injection and signed-difference results are
documented in [complete record recovery](../observer/RECORD_RECOVERY.md).
This package supplies an execution representation of those results; it
adds no protection or particle theorem. Exact real and imaginary operator
parts are applied separately, preserving arbitrary untouched references.

## Explicit host bounds

`compile_history` enumerates the whole source carrier, up to 8,192 columns.
It supports the local shells, bounded independent products and evaluated
fixed-number federations. Oversized carriers raise `ValueError`; they remain
members of the native declaration grammar. Near full binary occupation, the
compiler enumerates holes to avoid expanding every occupied port list.

The default map-entry budget is `2**20`, with hard maximum `2**24`. It bounds
the branch cap conservatively by
`min(max_branches, max_map_entries // source_columns)`. The returned
`effective_max_branches` is the cap actually passed to the shared interpreter.
If this cap or `max_steps` stops execution, the ordinary identity
`RESOURCE_STOP` outcome is retained. Direct comparisons must use that same
cap. The number of distinct complete histories remains exponential; this
optimization does not compress away classical record distinctions.

The compiled maps and rational payloads are tuples, the admitted program is
an immutable canonical string, and returned histories own fresh mutable
data. There is no global compiled-object cache. These objects are ordinary
in-process evaluation results, not an external serialization or admission
certificate.

## Validation

The tests compare every compiled branch with full Choi execution on H66 and
H78, and cover complex reference coherence, record feedback, holds, products,
federation sectors, FIFO timeout and resource stops. A two-attempt H66
feedback witness retains all 625 histories using 1,056 map entries instead of
38,788 materialized Choi coefficients. Those counts concern host
representation, not a physical memory gain or a change to source costs.

```bash
python -m pytest model/compiled_native_history/test_compiler.py -q
```
