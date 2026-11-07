# Weighted repair across cells

The same weighted endpoint rule can operate on the retained connected sphere
graph. Its diagonal stationary state has exact local restrictions, and a cell
with two loads has the native H66 reference. Population equilibrium is unique
within each fixed-number sector; coherent-state uniqueness and a uniform
preparation bound are unproved. Complete record bounds quantify exchange with
adjacent cells. The connected weight assignment is an explicit extension;
selecting it as the physical source requires an additional argument. See the
[proof](PROOF.md) and [receipt](results.json).

Run `python -m model.weighted_source_family.packet --check` and
`python -m pytest -q -W error model/weighted_source_family` from the repository root.
The evaluator takes exact rational real coefficients. Complex operators use
separate calls for their real and imaginary parts, as in the native interpreter.
