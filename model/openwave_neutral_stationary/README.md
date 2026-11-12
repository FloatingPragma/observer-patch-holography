# A stationary neutral excitation

A declared eight-cell source with periodic fermions and link penalty two has
an isolated neutral band, an explicit finite spectral filter and complete
charge, position and occupation records. Its zero-momentum energy above the
actual ground is 0.089555330476 source units. An electron-energy calibration
maps this to 66.1638659224 keV. The band width is below 4×10⁻¹⁴ source units;
positive propagation and a physical particle identification are not derived.

The accepted preparation is exactly stationary. Starting from the completely
mixed full carrier makes its success probability exceptionally small, of
order 10⁻¹²⁶. The polynomial filter has degree 167608; this counts applications
of the specified Hamiltonian, not primitive gates. The full carrier is not
materialized. The [proof](PROOF.md) distinguishes the exact operator bounds
from the numerical replay of the native 16-state local block and retains every
preparation, energy and position complement in the [receipt](results.json).

Run `python3 -m model.openwave_neutral_stationary --check` and
`python3 -m pytest -q model/openwave_neutral_stationary` from the repository
root. The physical interpretation requires a resolved dispersion, source
clock and length calibration, and observable interactions on this member.
