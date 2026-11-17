# Coverage amendment to the record-only admission verifier

This is a separately frozen verification repair of known evidence, not a new
native experiment, new scientific candidate or prospective discovery. The
ancestor is `A-native-record-only-admission-v2`, seal
`ba4fc9e36bf41592b272feee134c4bdd54d67c51416c8f992416535d536c9f2f`.
Its source/controller and all recorded evidence remain unchanged.

Root review identified that the old independent reader checked each supplied
consumer column, but not the exact catalog. It counted histories, native
successes and zero outcomes without requiring each registered input/stream
combination exactly once. A correct strict producer replay and the mathematical
source proof do not excuse this independent-checker gap. The initial separately
frozen mutation probe demonstrates four accepted malformed records: a missing
column, a duplicate replacing another column, a duplicated same-input stream
history and a duplicate replacing a different total. An uncompensated dropped
row already fails. Those failures of the old verifier are retained.

The amended checker requires the exact 30 unique `(position, edge, coin)` keys
for every one of the 864 successful native histories, and 25,920 columns in
total. It independently reconstructs the initial core from census currents,
checks that vector against the registered preparation formula and vacuum
ancillas, and identifies each public source stream by comparing the complete
**observed coin tape** with the two publicly declared streams. It then requires
each `(B, M, decoded total, observed stream, public law setting)` exactly once:
B=2/3/4, M=4/8, all totals from −12B to 12B, both streams and both laws,
1,752 cases. It also requires the exact twelve signed complete-prefix keys.

The stream identity check is a separate finite-assay provenance audit. It does
not give a seed or total to the controller, does not use a seed to infer a
hidden load and does not claim randomness from pseudorandom tapes. Both same
public streams are used for every input. The quantum consumer and current-only
controller remain byte-identical to the ancestor. The semantic checker ignores
`seed_audit_only`, `total_audit_only`, filenames and row order. Public geometry
and the declared source/control law setting are experimental inputs. The
full initial-state audit vector is still required to agree with the state
reconstructed from currents; it is not treated as an irrelevant label.

Before amended execution, freeze this design, amended checker, public graph
copy and validation source. Gates: the complete saved source passes normal
and optimized Python; all five catalog/domain mutations reject under both
modes at their intended conditions; a shuffled record with false seed/total/
filename/order labels passes; all saved evidence bytes, sizes, mtimes and
inodes remain unchanged. Every check is unconditional. The checker never
writes; the validation program defaults to saved-result comparison and only
`--record-to NEW_FILE` creates a new result exclusively.

This amendment earns verification coverage only. It does not increase the
input domain, number of native source runs, 12 executed complete prefixes or
the physical scope of the ancestor. No external G6 or official award is added.
