# Clarifications frozen before execution

This append-only clarification corrects the numerical-error wording of the
frozen proof without changing its preregistered bytes.

1. The unitarity and determinant residual tests in the colour preparation
   audit are geometry checks. They verify that each reported matrix is close
   to SU(3); they do **not** certify the forward arithmetic error relative to
   the exact Box-Muller/QR decoder. A different unitary could pass both tests.
   The preparation convergence theorem is an exact-arithmetic statement.
   The finite NumPy matrices are numerical diagnostics, with forward error
   certification separate and unprovided. No source energy or full triality
   lower bound relies on these floating-point diagnostics.

2. The frozen executable imports the pinned native function from the local
   theory path and checks its exact hash before execution. `inputs/native.py`
   supplies byte-identical standalone source custody for an independent
   replayer. `SOURCE_PINS.json` binds every copied theorem/policy input to
   its original path and exact content hash. No new Lean build is claimed.

3. Finite-B colour preparations have exact centre symmetry after the stated
   three-preparation mixture. Full SU(3) Haar/gauge symmetry is a source-limit
   statement. The complete centre charge-sector bound is already uniform
   at finite counter level and over every SU(3) representation; no uniform
   spectral gap is inferred from a finite list of sampled matrix moments.
