# Particle measurement audit

This independent audit checks the complete finite experiment: native regional charge readouts, the shared source clock, energy outcomes, and the continuation of an actual charged excitation. It replays source pulses and all failure maps, then compares every retained compact energy state.

The experiment supports ordered charge and transport records on the declared source. Its single GeV scale is supplied calibration; physical particle identity, stability and a physical clock require further attachment.

Run from the workspace root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 openwave-tuning/next/CHAIN/audit_particle/audit_particle.py
```

The [proof](PROOF.md) states the mathematical checks and numerical tolerances. [AUDIT_CONTRACT.json](AUDIT_CONTRACT.json) supplies exact consumer pins and the accepted scope. The script shares declared primitive source definitions, while using independent pulse contractions, clock generators and poststate comparisons.
