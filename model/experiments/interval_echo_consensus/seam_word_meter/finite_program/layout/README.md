# How many sharp registers fit on a native federation?

A connected federation of $C$ complete twelve-port owners with three distinct external gateways per owner supports exactly **$6C-6$ data dimers** for the finite-program apparatus. A program with $s$ logical registers, $k$ coefficient axes per register and $W\ge1$ binary places per axis fits if and only if $s k W\le6C-6$. One material midpoint at a time can use the native writing seam and two separate pointer stations. The twelve other ports hold the writing endpoints, the two heads and their four standards, and the two publication dimers. Their distinctness is part of this physical interface. [Prospective protocol](PROTOCOL.md) · [Mount and router](mount.py) · [Native material operations](../machine.py)

| Board | Owners | Available data dimers | Evidence |
| --- | ---: | ---: | --- |
| Level 0 | 20 | 114 | Complete packing, routed registers and both return reads |
| Level 1 | 80 | 474 | Complete packing, routed registers and both return reads |
| Level 2 | 320 | 1,914 | Complete packing, labelled route checks |
| Level 3 | 1,280 | 7,674 | Complete packing and labelled routes on the SHA-verified deposit |
| Levels 4–7 | 5,120 through 327,680 | 30,714 through 1,966,074 | Formula conditional on the supplied registered array and its connected owner graph |

On the fully packed level-0 board, 19 logical registers with two three-place coefficients occupied all 114 data dimers. The two component pairs had inputs $(2,6)$ and $(0,4)$. Both endpoint observers read the resulting means 4 and 2 at both destinations, the untouched sentinel 7, and the high publication flag. On a level-1 enrollment with nine registers, three axes and thirteen places, the two operated axes gave 6,142 and 2 at both stations, and the unchanged sentinel read 8,191. Both boards refused the attempted half-integral write and preserved every load through the refusal. These are deterministic native sharp effects on the specified contacts, routes and input words. [Raw per-place reads, stores and compact graph records](records.json) · [Level-3 certificate](level3.json)

## Constructing the routes

Each registered internal cell is an icosahedral graph with 12 labelled ports. Every external seam pairs one gateway with another owner. Build a spanning tree of the owner graph, recording its actual site pair per edge. There are $3C/2$ external seams and only $C-1$ tree edges, so at least one seam lies outside the tree. Choose the first eligible seam by its smallest endpoint as the writer. Its endpoints belong to two different owners. On each such owner, mount one head and two adjacent standards wholly away from the three gateway ports. The other two gateways remain available for all tree traffic.

For each of the twelve possible writer ports, there are $\binom{11}{2}$ choices of its other two gateways. On all **660 labelled cases**, the mount returns a head with two neighbors as standards, plus four disjoint internal data dimers covering the eight other ports. The remaining eight ports induce a connected graph. Each of the four service ports has a neighbor in that graph. These are exact checks of the finite shared cell template. All other owners have an internal perfect matching of six dimers. Reserve two of the resulting dimers for the material publication flags. The count is $4+4+6(C-2)-2=6C-6$. With 12 service and flag ports and two distinct ports per data bit, the port count forbids another bit on the same board. [Independent graph and native service checks](test_layout.py) · [Lean counting and exchange identities](Routing.lean)

To exchange two endpoint loads, ascend and descend their owner-tree path through recorded external seams. Within each owner, use a path avoiding every other writer, station and reference port. An addressed endpoint itself may be a head or a writing port. No transit path uses the active writing seam. For a simple physical path, traverse every adjacent seam forward, then return through all but the last seam. Induction by conjugating the shorter endpoint swap with the first seam swap proves that the resulting permutation exchanges only the two endpoints. Occupied shelves and publication flags can be intermediate ports; each is restored at the end of that exchange word. The second operand transfer and both station reads run serially under the original shutter model. Data from different read stations need no independent routing aisles. A disconnected control graph retains all its local cell and gateway incidences and is rejected.

The apparatus keeps one three-entry gateway table and a spanning tree over owners, and generates data shelves and physical paths as needed. Numbering logical register $r$ and axis $j$ from zero, its scalar shelf is $kr+j$ with places in increasing binary order. The tree route, writer and two heads are supplied services. An example uses the existing physical machine:

~~~python
from model.executable_model.carrier import federation
from model.experiments.interval_echo_consensus.seam_word_meter.finite_program.layout.mount import Mount

bench = Mount(federation(0)).bench(registers=19, axes=2, width=3)
bench.start([2, 0] + [7] * 34 + [6, 4])
assert bench.mean_axes([(0, 36), (1, 37)]) == "accepted"
assert bench.read(0, 0) == bench.read(36, 1) == 4
~~~

## Keeping the historical addresses

A separate matching check keeps the original 240-port bench's writer, four standards, two read heads and both publication dimers fixed. Its 228 other ports admit 114 data dimers. If all 84 historical data dimers must also keep their addresses, the 60 unused ports have connected component sizes 9, 10, 12, 12 and 17. Each additional dimer consumes two ports in a single component, so at most 29 can be added. The explicit matching attains this bound, giving **113 data dimers while preserving all 84 original dimers**. This is a combinatorial check on the fixed addresses, separate from the prospective movable-station experiment. [Both complete matching witnesses](original.json) · [Verifier](legacy.py)

The bundled loader supplies levels 0–2. The optional refinement reader accepts deposited levels 0–7 by byte identity. Level 3 was checked against its deposited SHA-256; the other four high-level arrays were unavailable for this replay. On any supplied board the mount checks connectivity and the gateway matching before allocating shelves. An unbounded sequence of material capacities requires an unbounded family of admitted connected boards. The finite deposit menu has a largest possible capacity. Width for exact dyadic computations, decision margins on pooled nonsharp inputs and physical derivation of the controller are separate requirements. [Refinement reader](../../../../../executable_model/refinement.py) · [Program domain](../README.md)

From the repository root the following commands use the bundled boards and records.

~~~sh
python -m model.experiments.interval_echo_consensus.seam_word_meter.finite_program.layout.measure
python -m model.experiments.interval_echo_consensus.seam_word_meter.finite_program.layout.legacy
python -m pytest -q model/experiments/interval_echo_consensus/seam_word_meter/finite_program/layout
(cd Lean && lake env lean -j1 ../model/experiments/interval_echo_consensus/seam_word_meter/finite_program/layout/Routing.lean)
~~~

To replay a deposit, add its file and number to the measurement command, for example `--deposit /path/to/FEDERATION_LEVEL3.npz --level 3`.
