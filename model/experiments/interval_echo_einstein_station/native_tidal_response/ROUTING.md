# A shared routing signal and its source

Changing the internal route of one cell changes the outgoing flights of both probe species by the same amount. In freely evolving Interval Echo the blank or written state of each routing seam is conserved. A separate exact test determines which sharp stationary material sources can produce the shared response through local energy coupling.

## Matched responses

The central cell's first perfect matching is replaced by its second one. Both preparations include the same parked resident and a sharp 32-unit coin-odd source on a guarded seam. The stored-energy gap has the same value for every endpoint-load word. Holding the guards in place while removing that 32-unit gap changes only a global phase in these launches. Changing their locations changes which rails exchange. Both tag kinds retain the same kinetic operator on each prepared routing pattern.

The table gives terminal counts as percentages for rays released from central port 2 and measured at two native time units. The receiver at port 3 was connected to the sender by a guarded bond on the original route and a blank bond on the replacement route. The exterior receiver counts all cells outside the central cell.

| Cells | Tag | Port 3, original | Port 3, replacement | Exterior, original | Exterior, replacement |
|---:|---:|---:|---:|---:|---:|
| 20 | 1 | 5.861% | 8.275% | 7.479% | 12.513% |
| 20 | 2 | 5.862% | 8.274% | 7.478% | 12.510% |
| 80 | 1 | 5.859% | 8.274% | 7.482% | 12.520% |
| 80 | 2 | 5.861% | 8.274% | 7.484% | 12.523% |
| 320 | 1 | 5.862% | 8.275% | 7.485% | 12.520% |
| 320 | 2 | 5.861% | 8.274% | 7.482% | 12.522% |

The two exact outgoing distributions coincide pointwise on each route. The 120 simulated frequencies at five appointments give simultaneous 99% half-widths of 0.01939 percentage points for differences. Observers can read a common change in routing from freely released rays. Here the experimenter supplied the change by preparing different sets of written seams.

## Complete stationary contact census

Let the resident load be r. Vary the neighbouring mobile load q between 0, 1 and 2. For each q, embed the entire 11-dimensional memory space at the endpoint word (r,q) in the local 99-dimensional Hamiltonian. For each possible energy E, solve the full-column equation (H − E)Jv = 0. The intersection across the three q values consists of all pure memories v that stay fixed up to phase and keep those three load words sharp. The calculation enumerates the complete eigenspaces over the rationals, including coherent memory superpositions. Reversing the resident and probe ports gives identical dimensions and shifts.

| Resident | Dimensions of admissible memory eigenspaces, summed | Possible conditional energy shifts for (q=1, q=2) |
|---:|---:|---|
| 0 | 8 | (0, 0) |
| 1 | 6 | (0, 0), (−32, 0), (−33, 0) |
| 2 | 8 | (0, 0) |

Every conditional shift of a sharp double tag vanishes in this census. None of the sharp stationary sources gives both tags a common nonzero potential. The finite proof relies on a fixed product memory ray and sharp resident and endpoint words for all three calibrations. Mobile matter, collective occupied backgrounds and memories that entangle with those loads require different evidence.

## What free evolution conserves

The blank memory projector commutes with the entire 99-state seam Hamiltonian. It also commutes with the other seams' terms and with the site energy, so its expectation and sector weights are conserved under free evolution of the complete board. An addressed write or a change of preparation can move weight between these sectors.

The exact memory-to-memory block for flight from (2,0) to (0,2) is −|blank⟩⟨blank|/2. A fixed memory ray supporting a coherent double-tag hop with the same memory can provide the blank amplitude −1/2 or zero. Moreover, on each written seam with endpoint words (0,0), (2,0) or (0,2), the native Hamiltonian equals the bare record generator on the entire ten-dimensional written memory. On otherwise empty material this gives exact factorization within any fixed blank/written pattern. The double tag sees the active blank graph and its native port diagonal, while the written memories evolve independently of its position. The result applies even to entangled, nonstationary written memories in that sector. Superposed routing sectors give an incoherent mixture of tag-position distributions when the orthogonal memory flags are unread. Thus storing energy within a fixed pattern on empty material cannot produce a new double-tag trajectory in this channel.

This command replays the [protocol](ROUTING_PROTOCOL.md) and checks the [exact census and terminal counts](routing.json).

```sh
python -m model.experiments.interval_echo_einstein_station.native_tidal_response.routing
```
