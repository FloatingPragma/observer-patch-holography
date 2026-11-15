# Replay and evidence custody

The frozen design is bound by PREREG.txt. SOURCE_PINS.json binds 152 source
and theorem files. The source run emits RECEIPT.json and its evidence hashes;
headline.py independently checks those bindings and replays the compact
source branches, source matrices and detector laws. The independent
acceptance_review.py audit recomputes every detector probability array from
the raw source matrices and effects, checks all reset tapes and exact energy
intervals, and composes the complete preparation/readout error bound.

From this workspace root:

    python openwave/shared/localized_rest_join/headline.py
    python openwave/shared/localized_rest_join/acceptance_review.py

The first replay takes about 0.2 seconds; the full stored-matrix audit takes
about 8.3 seconds on the development machine. Both retain the unresolved
pointer outcomes and auxiliary-reset failure certificates. The post-reset
pointer arrays are identified as such in README.md.

For a complete source regeneration, use a disposable copy, move its evidence
directory aside, and run:

    python openwave/shared/localized_rest_join/run.py

The source run took about 11.6 seconds and 332,712 native calls. It refuses
to overwrite evidence. Timing is not a hardware guarantee. Numerical execution
uses Python, NumPy and SciPy; declared proper-code services and precision are
the scientific apparatus contract, not inferred laboratory services.

The standalone snapshot preserves the pinned source files and the exact
packet in their original relative layout. Its root is `standalone/`.
Execute the same two replay commands from that root. This protects the
evidence from subsequent edits to workspace theorem/status files. Its
SNAPSHOT.sha256 binds every copied file; the outer MANIFEST.sha256 binds the
complete delivered packet. No upstream publication or column-pin change is
performed by either replay.
