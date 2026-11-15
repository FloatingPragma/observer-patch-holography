# Three-dimensional native propagation

Twelve-port cells connected on a periodic cubic grid support an isotropic
scalar dispersion branch. Its exact coefficient is `1/8` in every direction,
with a controlled fourth-order error. The same native endpoint rules give
complete position records on 96, 324 and 768 ports. A declared refinement
clock gives dimensionless wave speed `1/2`; rational spectral intervals
verify convergence on larger cubes without enumerating their full carriers.

The routing, coherent event services and clock are declared architecture.
Physical light speed, interacting field dynamics and particle masses require
an observation and unit map beyond this scalar dispersion theorem.

[Proof](PROOF.md) · [Source](source.py) · [Results](results.json)

Run `python3 -m model.native_spatial_refinement --write`, followed by
`python3 -m model.native_spatial_refinement --check`.
