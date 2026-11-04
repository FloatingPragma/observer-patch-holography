# Observable charged excitations of finite cores

The finite interacting sources have charged addition and removal poles
with a proved positive spectral weight. This connects their energy levels
to actual fermionic two-point functions and native local detectors.
One source has a specified mass texture; the other is the original open
Weyl cube with no imposed mass field.

For the 27-site cubes at U=1/65536, the charged gap is U/16 with an
explicit full-Fock error. Both zero-mode channels have low-pole weight
between 1/2-3t and 1/2, where t=w/(Delta_lower-2w) and w=27U/4.
The other spectral weight stays above a proved threshold. The
[proof](PROOF.md) derives the actual residue matrices, native number and
energy currents, and the complete Poisson-record comparison.

All repair events and failed quiet windows are retained. The Hamiltonian
poles do not imply autonomous repair protection. An explicitly specified
positive seam rate gives a continuing window with complete-output error
below 1/1000 for the exact ground input. The finite background
and energy units require physical identification; no observed particle
or common interacting field/gravity model is assigned.

~~~sh
python -W error -m pytest -q model/defect_core_response
python -W error -m model.defect_core_response --write
python -W error -m model.defect_core_response --check
~~~
