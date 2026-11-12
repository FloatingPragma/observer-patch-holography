# Magnetic moment under controlled refinement

At fixed matter and field volumes, the static bound-energy magnetic
factor rises from 1.47517 to 1.92738 as the lattice is refined. A full
precession experiment on the finest matter grid has a finite value near
1.92095 under its distinct doublet convention. Exact rational source
residuals and a global self-field bound control all unprepared modes.
The small self-coupling is an accuracy input; continuum precision and a
physical electron identification are separate questions.

Run `python3 -m model.openwave_magnetic_refinement --write` or `--check`
to rebuild and verify the operators against the stored numerical
witnesses. `python3 -m model.openwave_magnetic_refinement.regenerate`
regenerates those witnesses by solving the sparse eigenproblems.
See the [proof](PROOF.md) for the finite bounds and fit limitations.
