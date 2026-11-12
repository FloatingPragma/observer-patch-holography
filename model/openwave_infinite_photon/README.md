# Infinite-photon dressed particle

The full three-mode photon Fock space has an isolated dressed spin-half
band at the fixed complete spatial regulator. Finite trial vectors certify
its energy and photon tails using fourteen padded oscillator moments and a
global square-root remainder. The finite cap is removed without identifying
a truncated square root with a compression of the untruncated operator.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_infinite_photon --write
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m model.openwave_infinite_photon --check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q model/openwave_infinite_photon
```

`PROOF.md` states the source, domain, exact rotation argument, arithmetic
contract and scope. `results.json` contains the complete normalized integer
witnesses and source hashes. The GeV interval retains the same electron unit
as the finite-cap source. It is a convergence check, not a physical mass fit.
