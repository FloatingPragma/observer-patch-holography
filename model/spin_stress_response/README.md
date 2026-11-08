# Energy flow and entropy response

The prepared spin source has a measurable entropy response on the same clock
as its energy transport. Exact weighted-energy identities give the current and
work of a local source pulse. A complete native experiment prepares two spins,
applies a bond pulse and records both outcomes of the local modular-energy
measurement, with all failure records retained. Its regional first law has an
explicit quadratic error bound. Relating these quantities to a geometric
modular boost and a conserved Lorentz stress tensor requires a separate field
and geometry construction. See the [proof](PROOF.md) and [receipt](results.json).

```
python -m model.spin_stress_response --check model/spin_stress_response/results.json
python -m pytest -q model/spin_stress_response
```
