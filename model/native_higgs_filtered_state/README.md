# A shared filtered fermion state

An elementary fermion and a quantum-Higgs composite prepare almost the
same normalized state when their energy is restricted to a certified
window of the same finite source. The bound includes the actual finite
preparation error and retains every pump and filter failure. Any fixed
subsequent complete detector history inherits the state comparison.
The scalar remains a quantum operator; no constant Higgs background or
observed particle mass is inserted.

From the repository root run
`python -B -m model.native_higgs_filtered_state --write`,
`python -B -m model.native_higgs_filtered_state --check`, and
`python -m pytest -q model/native_higgs_filtered_state` with the scientific
environment. The [proof](PROOF.md) states the spectral-service premise,
source clock and the separate free-spinor consequence of the refined
shared energy. The calculation supplies no Standard Model scattering
identification.
