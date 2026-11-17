# Clarifications after the first run

The frozen proof's parenthetical count “(B+1)^12 signed initial vectors” is
a typographical cardinality error: the declared domain[-B,B]^12 contains
**(2B+1)^12** vectors. Its actual domain, theorem, code, inputs, gates and
probability bound use[-B,B] and are unchanged. The frozen proof is retained.

At the extraction checkpoint center0 is empty and the core total is Q−sign(Q).
During the quantum prefix center0 belongs to the target triangle. The fixed
reservoir spectator is therefore the core **excluding center0**. That reservoir
retains Q−sign(Q) throughout; the selected target may occupy center0 as part
of its actual source dynamics. No target/reservoir connecting seam outside
the triangle is addressed during the prefix.

The first independent checker used first-row equality to bind a selected prefix
representative. A complete representative may occur anywhere in an unordered
evidence set. Before the shuffled-metadata test, the checker was amended to
require exact membership of the full recorded-output set. The preserved earlier
checker, registry and append-only ledger show this order. All source dynamics
and first-run evidence are unchanged.
