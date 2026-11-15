# Integer face repair and the real Green minimum

`audit_integer_relaxation.py` checks the actual oriented icosahedral chain
complex, with the previously proved primitive integral cycle lattice. For a
unit dipole whose source and sink have graph distance d=1,2,3, the minimum
integer residual energy is respectively 1/2,1,3/2. The corresponding exact
real Green minimum is 11/60,7/30,1/4.

The integer lower bound follows from ||r||² >= ||r||₁ >= d for an integer
unit flow. A shortest path attains d. Real face coordinates instead give the
orthogonal Hodge projection and its lower energy. For the adjacent dipole,
face-coordinate denominator 60 already attains the real minimizer exactly.
Nearest rounding to a 1/k face grid has energy error at most 57/(4k²), as
verified against the exact rational minimizers.

An integer random field can have this same real mean field. Its mean energy
also contains its variance. The receipt explicitly constructs independent
floor/ceiling rounding of the real face coordinates: its mean field has the
real Green energy, while its mean instantaneous energy is larger and obeys
the integer lower bound.

This does not refute the real-valued Hodge lemma in `RF1_BRIEF.md`, which
already labels that lemma's domain as real and separately requests a native
integer process and time averages. It prevents silently transferring that
real minimum to the instantaneous energy of a single integer native state.
A time-averaged, noise-subtracted, or release observable needs its own
process/readout calculation. No RF1 solver, release estimator, or reserved
acceptance case was run in this mathematical premise audit.
