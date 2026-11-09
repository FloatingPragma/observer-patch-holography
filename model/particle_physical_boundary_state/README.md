# Charge at the boundary of a physical region

The selected spatial state is gauge neutral as a whole. A region cut from
it can carry charge, balanced by the outside. Actual Wilson-loop states
prove that the prepared state contains ten nonzero boundary-charge
sectors and the vacuum sector. A complete finite charge measurement
retains every outcome. The exact state preparation and this boundary
readout use the [same source and clock](../particle_common_contract/PROOF.md).

The [proof](PROOF.md) gives the positive weights, the distinction between
the outer boundary and a regional cut, and the complete119-state
coordinate-link control. Compatibility under simultaneous spatial
refinement and cooling requires a separate limit theorem.

Run with the pinned project interpreter:

```
python -m model.particle_physical_boundary_state
python -m pytest -q -W error model/particle_physical_boundary_state
```
