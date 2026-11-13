# Staggered finite Ward identities

The complete linearized curvature uses compatible scalar, vector and tensor
locations. It preserves exact Bianchi/gauge identities, reproduces the native
nearest-neighbour scalar Laplacian and removes the collocated stencil's
checkerboard ambiguity. All halos remain present.

`python -m model.openwave_staggered_ward --write` produces the receipt;
`--check` independently recomputes it. See [the proof](PROOF.md) for the
remaining observer co-location and physical metric/stress requirements.
