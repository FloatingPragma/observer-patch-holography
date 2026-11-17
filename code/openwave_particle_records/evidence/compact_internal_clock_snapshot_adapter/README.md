# Read-only document snapshot binding

The sealed compact_internal_clock packet pinned two theorem documents that
are expected to receive future additions. This adapter substitutes exact
byte-identical local snapshots only for those two read-source bindings.
All scientific code, provider/source inputs, old preregistrations, evidence
and receipt comparisons are unchanged. No sealed file is written. The
original strict replay still rejects changed live theorem files; this
adapter explicitly follows their original source bytes. Run replay.py.
This registration precedes the adapter replay and changes no scientific gate.
