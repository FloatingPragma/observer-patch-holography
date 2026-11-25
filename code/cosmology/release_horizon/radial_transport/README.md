# Finite shell currents and a three-axis record census

A source sends 432 sharp packets through each of two rectilinear arrays of 125 owners on the registered 320-cell Interval Echo carrier. One source cycle visits three principal directions and all pairs of twelve transverse slopes once. Material arrival records give the same total current across four successive shells. The number of crossing edges is 3, 12, 27 and 48. Every edge of its shell receives respectively 144, 36, 16 and 9 packets. The current per injection and per edge is exactly 1/(3n²) for shell index n from one through four.

## Apparatus and normalization

The [protocol](protocol.json) fixes the finite graph, launch cycle, clocks and interventions. A native π exchange swaps the complete loads on a seam. A routed hop swaps two specified endpoints and restores every borrowed port in reverse order. The compiler checks the permutation with distinct port labels. A sharp echo read against an owner-local blank returns the load and restores its input. Both placements share the registered carrier. They compare distinct endpoint and transit assignments under the same supplied routing rule. The octant is a graph of services with reflecting lower coordinate planes. Its corner and its layer-zero injection owner are distinct. Elementary carrier distances vary from link to link.

Each box has 300 positive coordinate links and one injection link. A fresh load challenges each link and both endpoints are read. A gate parked on one crossing reduces the third cut from 27 links to 26. All offered links, including those that receive zero traffic in a cycle, enter the flux analysis. A ray with principal axis speed one has transverse speeds (2i+1)/24 and (2j+1)/24. It advances on the next integer-plane crossing, with axis order breaking transverse ties. At the nth principal crossing the two face addresses are the integer parts of n times the two speeds. Each bin receives 12/n settings per transverse axis on the four measured shells. The local balance counts departures minus arrivals over a complete injection cycle. All 64 grid stations through layer four balance exactly, the injection owner supplies 432 loads and layer five collects them. Repeating the cycle leaves the transport register in its starting state. Selecting a single principal face concentrates the flux. Absorbing at (1,0,0) gives 144 more arrivals than departures (net outflow −144).

The joint shell-two and shell-four records are a secondary analysis of the frozen cycle. All 432 packets reach a fine edge whose transverse address, divided by two and rounded down, recovers their coarse edge on the same face. Pairing independent launches across the complete product of measured marginals gives this match with frequency 1/12. This is directional memory of one prepared packet across two radii. A primordial curvature covariance and the sky transfer require a physical source identification.

## Volume and time on the same board

A range probe leaves the corner, meets a persistent sharp reflector and returns. Each of the three reflectors is transported along its own axis between epochs. The longest link requires 119 elementary seam exchanges. Each route receives the same allocation, using parked exposures after its final exchange. A separate dimer on owner 273, outside every transit path, ticks in each active or parked slot. Its observed round trips are 238, 476 and 952 changes at widths one, two and four. Parking the last outbound link yields a missing reflection. These distances are in calibrated units of the routed apparatus. The local crest clock runs on the conformal service schedule, independently of this range clock.

On each of 36 conformal layers the corner and three widths define the half-open census box. All owners carry sharp clocks. Only a changed marker inside the box contributes a publication. On the two arrays the counts are 1, 16 and 256 per box per layer. The pulse train on the same carrier gives successive received-to-emitted periods of 2 and 4 on the two expanded plateaus. The angular-only control gives 1, 8 and 64 per base face. At cadence four, unequal measured widths 2, 3 and 4 give 96 publications. Parking one interior marker for the four strokes gives 252 with all widths four. The two controls distinguish record yield from measured extent. The range transactions preserve the shared crest and census registers at their service boundaries.

This volume is the product of three routed axis lengths. The expansion profile, direction ensemble and measurement schedules are construction choices. Identifying this volume with the physical volume of the thermal source and calibrating the radiation horizon require additional measurements. The finite equal-current result concerns a cubical cut with absorbing collectors beyond its declared depth. A mass or charge coupling and a spherical physical metric require separate evidence.

## How far apart are the correlated heat probes?

A separate range calibration uses the fixed thermal contrasts on owners (12, 6) and (35, 33). Their energy data come from the [36-layer clock archive](../background_clock/records.json). The [calibration](thermal_span.py) independently repeats all 320 single-owner impulse interventions for the two sweeps between release-marker readings 52 and 60. It then addresses a reflector at each named target and counts a clocked round trip from each possible input. For a coincident source and target the return takes zero strokes.

All four cross-contrast endpoint pairs give 714 range-clock changes, or three routed distance units. The largest measured input distances for the four endpoints are respectively 3, 2, 3 and 4. The two sets of delivered input addresses are disjoint. Their radial bounds overlap. Every cross-pair separation is at most the larger of its two measured input radii. The positive energy covariance reconstructed from all sixteen energy-word counts is 6.41 × 10⁻⁷ in squared native gaps, with a simultaneous 99% half-width of 2.25 × 10⁻⁷. It measures shared preparation in this circuit. These data do not establish primordial coherence at a separation beyond the reach of radiation. The energy and ruler preparations are separate, matched by carrier owner. The pulse channel and physical radiation require independent identification.

## Replay and record layout

From the repository root:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
.venv/bin/python code/cosmology/release_horizon/radial_transport/measure.py
.venv/bin/python -m pytest -q code/cosmology/release_horizon/radial_transport
(cd Lean && LEAN_NUM_THREADS=1 lake build +ObserverPatchHolography.Provenance.RoutedOctant:olean)
```

The [records](records.json) use local owner index 25x+5y+z. Source index −1 maps to owner 250 or 251. Grid offsets are 0 and 125. Edge current vectors use the ordered pairs in `edge_order`. Vertex outflow vectors start at −1. Sparse joint rows contain coarse edge index, fine edge index and joint count. Absent joint entries have zero count. The [experiment](measure.py), [seam compiler](carrier_routes.py), [independent record checks](test_radial_transport.py) and [Lean proofs](../../../../Lean/ObserverPatchHolography/Provenance/RoutedOctant.lean) specify the readouts and their finite scope.
