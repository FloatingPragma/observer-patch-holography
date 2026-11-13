# Energy-dependent full-Haar Schur certificate

The producer evaluates the complete shifted Gram correction at nonzero
energies and certifies full-source first-excited eigenvalue brackets by finite
inertia. Every omitted Haar representation remains in the analytic remainder.
The complete coarse and refined integer frames are checked against the exact
source; the inherited ground bound supplies the sharper gap interval.

Three spatial patches also quantify why a global fixed-electric-energy
projection cannot have a volume-uniform off-block norm. The existing local
propagation theorem supplies a separate path for local refinement.

Run `python -W error -m model.openwave_energy_schur --write`, then run
`python -W error -m model.openwave_energy_schur --check` in a separate process.
Focused checks use `python -m pytest model/openwave_energy_schur/test_schur.py`.
`prepare_witnesses.py` proposes proof coordinates; it is not an acceptance
check. See `PROOF.md` for domains, complete-complement bounds and scope.
