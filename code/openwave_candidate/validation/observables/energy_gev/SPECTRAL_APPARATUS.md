# Complete spectral/native-event apparatus with work and reference ledger

`spectral_apparatus.run()` instantiates the energy-shell isometry method used by `model/coherent_defect_wave/mass_instrument/model.py::complete_experiment` on the actual modular uniform-fiber triangle Hamiltonian and native C/E Kraus maps. It does not import that program's CAR reservoir, filled sea or mass-two result.

The admitted initial source is I/d with retained conjugate purification vec(I)/sqrt(d). The declared experiment measures initial H energy, executes one source-seam native commit with every letter retained, then measures final H energy. Every (i,r,j) record, including zeros, has probability |<j|K_r|i>|²/d. The initial energy label is retained on the reference; its conditioning is explicit and differs from an untouched purification. Thus this is a real sequential joint law, distinct from `ENERGY_GEV.json`'s separate marginals.

A declared battery carries change Ei−Ej in every record. On its admitted initial support, the shell map |i,B0>→Σjr <j|K_r|i>|j,B0+Ei−Ej,r> conserves source+battery energy and is an isometry by ΣK_r†K_r=I. The code verifies that completeness identity and every energy balance. This does not claim a finite repeatedly reusable autonomous battery or a native implementation of the spectral apparatus. Preparation and spectral access are supplied resources; no unaccounted claim that they cost zero is made.

The source-neutral `finite_spectral_adapter.py` additionally evaluates the existing original-transfer/final-excitation coupling theorem WITHOUT an initial energy measurement. Its original weights are not normalized by probe success. This mathematical coupling is kept distinct from the sequential apparatus, so neither reference disturbance nor failure records are hidden.

All energies and clocks retain source units. No vacuum, species or physical mass identity is added. The full criterion remains false. Tiny tests verify complete records, negative transfers, shell isometry and work balance; no reserved acceptance or expensive source rerun is involved.
