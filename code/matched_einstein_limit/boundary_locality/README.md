# How far can the exterior change a regional energy reading?

Trace out the ports outside a region of a thermal spin system. Its entropy response depends on the logarithm of the resulting density matrix. Couplings across the cut can change that response several contacts inside. At each finite order in inverse temperature, the exterior contribution must follow a connected set of contacts back to the cut. On the binary Interval Echo cell the first correction reaching a port two steps from the exterior occurs at cubic order. The proof also bounds the summed tail at fixed temperature on each finite graph. The bound depends on the spectral width of that graph. [Proof and fixed-temperature bound](PROOF.md)

The exact calculation retains all 30 contacts of the icosahedron, traces out port 0 and reads the pair at ports 1 and 8. The quadratic correction vanishes. At cubic order the two transverse Pauli coefficients are each 145/96 and the longitudinal coefficient is 17/96. Halving the native contrast changes these to 41/96 and 9/96. On the two contacts along the path, a rational enclosure also proves that the transverse correction is positive at inverse temperature 1/16. [Exact calculation](coefficients.py)

Run from the repository root:

```sh
PYTHONPATH=code:. python -m matched_einstein_limit.boundary_locality.coefficients
PYTHONPATH=code:. python -m matched_einstein_limit.boundary_locality.transport_cmi
python -m pytest -q code/matched_einstein_limit/boundary_locality
(cd Lean && lake env lean ObserverPatchHolography/EinsteinBranch/MeanForceBoundary.lean)
```
