# Continuous geometry and finite matter sources

These executable constructions test how finite observer records can support
continuous geometry and particle descriptions. They establish finite regional
entropy, a certified interacting wall vacuum, colour and spin bounds, and a
complete conditional charge instrument. The [proof](PROOF.md) states their
domains, the [source register](SOURCE_REGISTER.json) separates the physical
claims, and the [export manifest](EXPORT_MANIFEST.json) binds their code and
retained counterexamples. Physical particle predictions and an Einstein
realization require source, state and clock identifications beyond these
component results.

From the RER root:

~~~sh
python3 model/openwave_obstruction_realizations/verify.py
python3 model/openwave_obstruction_realizations/verify.py --replay quick
~~~

Custody verification uses the Python standard library. Replay requires
NumPy, SciPy, SymPy and mpmath. The trace replay recomputes the independent
1024-holonomy interval trace; the charge replay executes the fixed-policy
source and pointer instrument. The all replay runs all three groups.
Replays use a temporary copy of the frozen inputs and retain every outcome.
They print the location of their logs and receipts.

The additions replay checks the complete neutral-pole intervals and
current records, and the held-out gravity-fit failure. The all replay
includes this group as well. Successful reproduction of a failed
physical comparison preserves its negative conclusion.

The inputs tree preserves workspace-relative paths and exact source bytes.
Its imported records include construction history and failed attempts as
scientific provenance. The workspace export tool is not required to verify
or run this standalone package.
