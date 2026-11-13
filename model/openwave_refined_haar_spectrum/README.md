# Common-clock coarse and refined Haar spectra

This calculation certifies full-source gaps and local field variances on a
coarse cube and its actual spatial refinement. Every omitted Haar channel enters
the Schur and trial residual bounds. An exact double-commutator identity yields
Wilson first spectral moments, including the eight-link coarse field embedded
in the refined graph. Positive spectral measures then give Euclidean correlation
intervals without assuming the two models agree.

The declared pair is (a,g²)=(1,64) and(1/2,16). Values are in the retained source
clock, not GeV. This does not establish a derived coupling flow or a4D continuum.
Run `python -m model.openwave_refined_haar_spectrum --write`, then a fresh
`python -m model.openwave_refined_haar_spectrum --check`. `prepare_witnesses` only
proposes fixed integer proof coordinates; the producer rechecks them completely.
