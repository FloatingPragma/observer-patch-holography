# Finite charge motion with a retained field

A bounded Peierls hopping term moves a declared fermionic charge on the
finite cone while retaining its fixed compensator, exact Gauss constraints
and finite transverse field. The actual neutral ground stays exact. The
same Hamiltonian gives charged energy windows, neutral field spectra and
strictly positive complete arrival records, including arrivals with
nonzero transverse occupation.

The [proof](PROOF.md) states the changed Hamiltonian, all-order error bound
and complete failure laws. The source is distinct from the native 119-state
link program. It proves no Coulomb phase, continuum scattering or measured
particle assignment.

From the repository root with `model/requirements.txt` installed:

```sh
python -m model.affine_gauss_mobile --write
python -m model.affine_gauss_mobile --check
python -m pytest -q model/affine_gauss_mobile
```

[Provenance](PROVENANCE.md) identifies the source replay and independent
scientific comparison. The finite matrices are declared algebraic services;
no efficient native compilation or unrestricted occupation limit is claimed.
