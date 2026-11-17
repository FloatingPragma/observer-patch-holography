# Full five-state memory complement

The inherited memory is a classical recorded reset chain, not a sterile
Hamiltonian. This package preserves that type and first tests its actual
complete channel. It does not call relaxation eigenvalues particle masses.

Fixed source: p=1/12; advance around the inherited fivefold port orbit
(1,5,11,10,7); otherwise reset to memory0; retain old memory in the output
record; feedback exp(−ih X tensor P_new). Use h=1/8,1/64,1/512, which are
existing producer/readout step counts, as a convergence/control ladder.
No step is selected by a mass or angle comparison. Retain all results.

Calculate all180 CQ poles, project the complete resolvent onto the stationary
memory marginal, and verify the exact Schur identity with its144-dimensional
complement. Derive the first two coefficients of the slow generator, keeping
both damping and its coherent part. Record their ancestry and compare the
leading Hamiltonian with the already fixed Y. The charged current is an
external probe of this complete response; a three-state unitary PMNS matrix
must not be manufactured by throwing away damped modes.

Only after that audit could an explicitly added coherent memory/sterile
service define a genuinely new Hamiltonian source. The old recorded reset
instrument itself does not supply one. No empirical mass scale, phase,
permutation, or counterterm is inserted here.
