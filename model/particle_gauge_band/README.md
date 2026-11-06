# Gauge fields and compact particle bands

Classical unitary links preserve a separated compact band on the complete
343-site source. Including the pairing in its spectral projector gives 180
canonical annihilation modes, a gauge-covariant quartic interaction, and an
explicit polar controller transporting the occupation grading between records.
The proof uses exact norm bounds; finite matrix execution checks covariance
and all thirteen raw family branches. The background links and transport
controller are supplied operations, without a physical mass identification.

[The classical-link proof](GAUGE_BAND_PROOF.md) retains the entire lattice,
pairing and complementary band. [The quantum-link proof](QUANTUM_LINK_PROOF.md)
constructs all sixty-state coordinate matrices under a compact-group rotor
attachment, checks 1,416 covariance identities and gives their exact norms.
The rotor attachment changes the native raising-only dynamics. Covariance
permits five independent link-energy differences; it selects none of them.

Quantum-link projectors require a separate fermion-algebra construction. The
[exact native witness](NATIVE_LINK_CAR_WITNESS.md) uses the literal raising link
and a specified two-site Hamiltonian: its gauge-invariant projector has a
noninteger naive occupation even on a Gauss singlet. Classical-band occupation
therefore cannot be transferred to quantum links merely by replacing matrix
entries with operators.

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_gauge_band.verify
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m model.particle_gauge_band.verify --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q model/particle_gauge_band
```

The check independently reruns every producer in a temporary directory,
compares exact identities and link coordinates, and checks numerical invariant
agreement. It also verifies stored source and binary-witness hashes. The
[gauge receipt](gauge_band.json), [coordinate receipt](quantum_link_completion.json)
and [native CAR receipt](native_link_projector_witness.json) distinguish exact
rational/algebraic results from numerical diagnostics. The latter are not
outward interval certificates or physical uncertainty estimates. No source
selection of a background, vacuum, link energies or laboratory clock is
asserted. The charged coefficient architecture retains its target-informed
provenance; these results prove conditional finite compatibility.
