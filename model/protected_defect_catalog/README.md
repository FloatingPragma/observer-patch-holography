# Protected objects before particle names

The [catalogue](CATALOG.md) distinguishes source memory charges, conserved
spin sectors and prescribed topological textures. It records how each is
prepared, protected and measured, together with its exact classification
domain. The producer replays complete bare-source and ordered-cover charge
classifications and exhausts all 65,536 Pauli-product charges of the selected
spin cube. The same source also has an exact integer spin charge; a local
creation followed by one native local pointer prepares two coherent charge
channels with certified weights. A separate fine readout audits their total
charge probabilities. No Standard Model particle or physical mass is identified.

The same native parity acts on the reconstructed canonical field. Its odd
Stone-field excitations have the inherited positive-energy spectrum. A
complete odd/even bounded preparation has a uniformly nonzero probability
in the source scaling, including its filter cost. The continuous internal
charge acts on the transverse configuration-field subnet. Its conjugate
charge excitations have equal nonzero energy-momentum measures on the same
state and clock. A sharp particle shell and physical spin require further
source estimates.

The one-hole sector has a complete source position detector that retains
internal spin and reference information within each recorded position
block. A uniform global-attempt schedule bounds total source/archive
disturbance relative to identical propagation with identity repair by `4/n`
under its stated propagation and resource conditions.
Its preparation and clock are supplied; general Hubbard evolution can leave
this sector, so the entry identifies no physical particle or mass.

The [proof](PROOF.md) gives the graph classification and the remaining
source-to-particle requirements. A complete catalogue over unspecified source
architectures is not claimed.

The [kinetic-node map](../source_family_compatibility/KINETIC_NODE_MAP.md)
adds all eight momentum-space windings and their conditional free Weyl-field
dictionary, including helicity and conjugate holes. It does not assign them
to measured species. On the separate prepared collective field,
[bounded charged records](../spin_joint_field_limit/ARCTAN_CHARGED_BANDS.md)
retain a uniform positive probability in explicit energy--momentum bands.
The same-state dispersion and scattering determine whether those records
contain a particle branch.

On the interacting centered source, [actual nodal symmetries](../source_family_compatibility/NODAL_CHANNELS.md)
map all eight original valley probes to one even scalar nodal measure,
using the actual preparation's inversion symmetry. The spectral measure
is defined at every momentum, but its
low-energy scaling remains unresolved. The [complete spectral instrument](../source_family_compatibility/SPECTRAL_WINDOW_ACCESS.md)
provides a measurable score with explicit errors and a ground-pressure
certificate. On the separate prepared collective family, [directed caps](../spin_joint_field_limit/DIRECTED_CHARGED_BANDS.md)
retain positive original weight in each charge and momentum direction.
These entries explain specific source-to-observer responses; none supplies
a measured mass or completes the Standard Model identification.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.protected_defect_catalog --check
python -m pytest -q -W error model/protected_defect_catalog
```
