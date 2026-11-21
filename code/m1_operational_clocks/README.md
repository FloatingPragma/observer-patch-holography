# Massive operational clocks

An explicit mass coupling of the charged tetrahedral walk supplies massive
Dirac dynamics, a local interference clock and quantitative control of
resolved reads in one native time coordinate. The maximum microscopic record
speed remains c, three times the continuum matter speed. Every smooth
quasienergy band has speed at most `(c/3) cos(m tau)`; the fixed-spacing
ballistic distribution obeys the same bound. No species or high band is
discarded to obtain these results.

The central result concerns actual operations: compact H2 smear functions
give local CAR number reads and phase encodings whose influence outside the
Dirac cone vanishes as `O(sqrt(a))` at fixed duration and resolution,
uniformly over Fock states and spectator ancillas. A delta-site preparation
retains its exact speed-c record. The proof specifies the operation class
instead of assuming a low-energy state or a successful M1 read law.

Two normalizable mass packets with a controlled finite preparation cutoff
produce a positive local Ramsey readout with a bounded-error comparison to
the beat law `Delta m/gamma`. Its reference fringe and finite-lattice error
are derived. A finite parameter example covers an entire moving-clock cycle
with reference fringe coefficient above 0.7457, probability error below
0.035153 and **actual probability swing above 0.6753** after every error.
The finite signal is not assumed exactly monochromatic. Its tiny spacing
is an analytic existence witness, not an executed huge lattice. The finite
executions are periodic waveguide packets at four coarser resolutions,
including all modes of their transverse-constant invariant sector.

Read [the contract](CONTRACT.md).
The quantum gate family, graph, native speed and masses are declared inputs;
A1--A3 source admission and physical energy calibration are not derived.
The operator proofs are analytic, not claimed Lean formalizations.

## Reproduction

From the repository root with `requirements.txt` installed:

```sh
python code/m1_operational_clocks/build.py
python code/m1_operational_clocks/verify.py
python -m pytest -q code/m1_operational_clocks
```

The producer assembles the walk and diagonalizes it. The checker reconstructs
individual gate entries, uses the closed eight-band spectrum, and evolves
packets by binary multiplication rather than diagonalization. Both calculate
every declared grid mode. The checker imports no producer functions. It
replays 12,672 eigenphase magnitudes, 960 full-channel dynamics cases,
96 all-band velocity slopes, eight lossless native records and twelve
two-mass packets, including actual reversed evolution. The receipt includes
the complete spectral values, not just a fitted slope or success flag.
Float comparisons use relative 2e-9 tolerance, with absolute 2e-10 for
numerical residuals and 1e-12 for other expected zeros; catalog shape, keys,
integer counts, JSON types and source custody are exact. Small residuals are
not certified exact zeros. Positive spacing, Fourier tails and error bounds
cannot be replaced by zero within tolerance.

The tests independently prove finite Clifford identities with symbolic
algebra, exercise three-dimensional momenta and all valleys, test wrong
flight adjoints, odd-tick carriers, retiming, phase erasure, missing bands,
malformed JSON and source tampering. Selected forgeries go through the
optimized Python CLI, so assert removal cannot bypass acceptance.
The analytic limits are proved in the note; numerical catalogs support
their finite implementation and do not stand in for an all-level proof.
