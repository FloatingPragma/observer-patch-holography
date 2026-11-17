# Frozen-input replay adapter

This reader changes only resolution of the two documentation dependencies
`openwave/row_05_magnetic_moment/THEOREMS.md` and
`openwave/row_06_spin_half/THEOREMS.md`. The exact bytes at the primary freeze
are copied into frozen_inputs and their original source path/hash pairs are
retained in FROZEN_INPUT_MAP.json. This permits maintained row documentation
to include results without changing their scientific provenance.

All other primary source pins, the four primary frozen files, the complete
scientific receipt and every retained array are checked by the same rules.
No scientific function or source input used to generate data changes. The
adapter overrides only check_pins in memory and calls the original main
with no --record argument. It refuses command-line options. The primary
probe, preregistration and original pins remain byte-identical.

Freeze this task, replay.py, map and both copies before running the adapter.
Success means canonical equality with the original receipt, exact equality
with every array, no evidence write, and elapsed time below60seconds.
