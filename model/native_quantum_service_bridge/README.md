# Native quantum service bridge

This package executes two finite field services through the existing twelve-port response controls and an explicitly supplied coherent code-transfer interface.

* **Native Yukawa singlet sector:** exact `HV=Vh` in the complete 320-state parent, with all twelve gauge charges annihilating the two-state code. One simultaneous source pulse implements its continuous evolution, and six existing source pulses connect it to a complete native C/E pointer.
* **Compensated U(1) link:** all 491,520 CAR/Gauss columns are checked. Fifteen source-X pulses reproduce the parent's complete 32,768-amplitude Cayley state.
* **Finite matrix compiler:** exact Pauli decomposition, emitted source words, retained padding and an explicit product-formula error/cost bound. The executable checks the native five-state block in an eight-state memory and a complex state preparation.

The [proof and service ledger](PROOF.md) separates exact identities, rational precision bounds and binary64 diagnostics. Coherent load/export, signed analog port control, memory routing and physical service accuracy remain explicit architecture premises. Real Observer repair alone and the passive classical simulator do not implement this program. A general nonabelian spatial program is not executed here.

The finite-precision replay retains all export failures and feeds the unnormalized success branch to the native pointer. Native reset successes and timeouts retain their tapes and final work; an explicit finite recipe represents the complete record tree.

From the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 -m model.native_quantum_service_bridge --check
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 python3 -m pytest -q -p no:cacheprovider model/native_quantum_service_bridge
```

`--write /new/directory` emits fresh receipts without overwriting an existing directory. `pins.json` authenticates the source/data inputs before replay. `results.json` contains the native services, exact checks, all emitted pulse requests and complete record semantics; `finite_compiler.json` contains the generalized word and error certificate.
