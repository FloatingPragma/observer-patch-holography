# Native oscillator field

Each of the 190 photon modes and three freely spreading harmonic modes has
its own finite number register. Actual routed currents couple the charge
to those registers through explicit native controls. A probe measures a
Fourier component of the physical electric field on any of the 288 edges.

The [proof](PROOF.md) gives the cutoff, spectral and compiler bounds.
[The executable](run.py) checks complete small native words and emits
finite full-program resource contracts. Its source domain is the size-two
one-charge construction; interacting spatial refinement and a complete
joint sharp-edge apparatus require separate results.

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 run.py` from this
directory. `providers.py` imports only the snapshots bound by `INPUTS.json`.
The full field iterator is deliberately not materialized by the run.
