# Charge, energy and finite transport

The same interacting source supplies spatial charge records, positive excitation energy and continued source states on the declared event clock. Two executed regional detectors enclose a clock interval; a grouped energy measurement follows. The complete output includes every clock count class, detector timeout, export failure and zero energy bin. A positive charged excitation is continued as an additional conditional witness, with its full poststate retained.

The calculation measures finite regional transport. It supplies no asymptotic stability theorem or particle-species identification. The energy display uses one supplied 511 keV calibration.

Reproduce from the workspace root:

```sh
OPENBLAS_NUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/CHAIN/particle/particle_chain.py
OPENBLAS_NUM_THREADS=1 reverse-engineering-reality/.venv/bin/python openwave-tuning/next/CHAIN/particle/verify.py
```

[The contract](PARTICLE_CONTRACT.json) defines the source, clock, complete outcomes and state reconstruction. [The result](RESULT.json) contains the measured law. [The proof](PROOF.md) states the measurement order and apparatus premises.
