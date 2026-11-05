# All quadratic conserved charges of the interacting source

The complete Weyl/contact source has exactly four independent traceless
self-adjoint quadratic conserved charges under its original Hamiltonian,
and six under its centered Hamiltonian. The classification includes every
pair-creation and pair-annihilation term. These are total number, the three
hidden-spin charges and, for the centered action, two staggered eta-pair
charges. A scalar identity is additional in each case.

The [proof](PROOF.md) applies to every finite connected open rectangular
box with positive hopping and contact. Exact full-algebra ranks cover the
three axis dimers, a three-site chain, a plaquette and the eight-site cube;
independent occupation-matrix ranks check the complete plaquette. No
nonconstant conserved quadratic is confined to a proper site region.
This closes the anomalous-pairing alternative to the finite quadratic
momentum classification. Higher-degree charges, a continuum stress tensor
and physical particle identification require separate constructions.

```sh
.venv/bin/python -W error -m pytest -q model/quadratic_conservation
.venv/bin/python -W error -m model.quadratic_conservation --check
```

The deterministic [receipt](results.json) retains all number sectors and
binds the proof and executable source dependencies. It establishes no
physical completion gate or additional empirical comparison.
