# All-character Z6 energy and native-accessible excitation support

The same triangle at m=.23 with three CAR occupations and Z6 records has120 full states and six20-dimensional character blocks. This reuses `character_reference.blocks`, hence the existing Fable chord-character and Peierls APIs, the same signed Hamiltonian and all native C/E/MIX phases. No potential or new interaction is introduced.

The unique global ground is in character3 (holonomy pi), energy −4.065399873582047. The full fixed-number spectral gap is .05474076633844. The global lowest **signed** one-load energy also prefers character3. Only the minimum **positive-band** one-load edge prefers character0. Energetic claims must preserve this branch distinction.

The global spectral gap includes other characters. Every native branch preserves each character, so a prepared ground in character3 cannot be driven into a different-character low-lying state by this instrument. The receipt therefore also evaluates the complete joint native-event / full-H energy-eigenspace law starting from the actual global ground. Every distinct energy and its degeneracy is retained, including zero-probability groups. Within a degenerate eigenspace probabilities are aggregated, avoiding arbitrary rank-one basis dependence. The energy PVM is supplied apparatus and is separate from the native local readout.

The first numerically nonzero native excitation is selected with an explicit probability threshold1e−12 and energy grouping tolerance1e−9. Other-character exclusion is algebraic from the verified invariant projectors; finer zeroes inside the character are numerical diagnostics. The full spectral gap must not be substituted for the native-accessible support or called a physical mass. All branch probabilities, event/energy marginals and mean monitoring work remain available for replay. Character preservation also provides no cooling or automatic preparation of the energetically preferred sector.

CAR statistics, transport phase, finite topology, finite modulus and initial ground preparation remain explicit. No physical particle, spin, binding, physical vacuum selection mechanism or calibrated mass is established. This is a scoped experiment separate from the candidate runner and official assessment.

For portable replay, runtime dependencies use the byte-identical providers in `validation/bridges/providers/`, with source hashes, licence and provenance. All receipt pins are relative to the public bundle root; no external workspace source is imported.
