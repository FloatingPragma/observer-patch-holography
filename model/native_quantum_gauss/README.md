# Local charge conservation on a finite native gauge link

The actual fifteen-state native matter currents admit a chosen 60-state
link register with exact local non-Abelian Gauss generators. Matter can move
between endpoints within each multiplet while the link carries compensating charge.
All binary position readouts preserve Gauss, and their energy exchange is
computed from that same Hamiltonian. The link design, coherent control and
static dual reference are explicit inputs; the [proof](PROOF.md) does not
identify this register with a gauge dressing of the compact defect or a
chiral continuum field.

```sh
python -m model.native_quantum_gauss --write
python -m model.native_quantum_gauss --verify
python -m pytest -q model/native_quantum_gauss/test_algebra.py
```
