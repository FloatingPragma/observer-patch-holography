# Quantum weights reconstructed from native histories

Closed native repair histories determine a positive matrix functional on a
reconstructed finite algebra. Its trace form follows from the chosen linear
extension of classical endpoint probabilities, and native C and E evolution
converges to the uniform state without preparing it. The complete-attempt transfer
and twelve occupation fields generate all 4356 matrix coordinates; the
Hamiltonian is obtained from that same transfer. The [proof](PROOF.md)
separates this count-derived quantum description from implementing arbitrary
quantum measurements or identifying the reconstructed time with physical time.

```sh
python -m model.native_quantum_reconstruction.certificate --check
python -m pytest -q -W error model/native_quantum_reconstruction
```
