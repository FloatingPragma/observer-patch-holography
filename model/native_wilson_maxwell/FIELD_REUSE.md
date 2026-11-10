# Scope of reused Maxwell mathematics

`model/affine_gauss_source` supplies the finite total-occupation cutoff,
positive boundary completion and affine Gauss embedding pattern.
`model/affine_gauss_mobile` supplies the literal finite Peierls exponential
and the total-occupation enumerator. The cube calculation executes the
actual cube incidences, Green matrix, mode projectors and all56 photon
states. Cone eigenvalues, cone charged energies and cone current
probabilities are not transferred to this geometry.

The dimension-independent rational inverse in
`code/electromagnetism/discrete_coulomb_green.py` constructs the actual
cube Green matrix. Its twelve-port wrapper is not used. The finite
Gauss and minimum-energy identities are replayed on the cube matrices.

`WhitneyNormalModeConstruction.lean` proves existence of a positive normal
frame when the supplied constrained stiffness is invariant and strictly
positive relative to its mass form. Here both hypotheses and a complete
frame are executed for K on range PT. `WhitneyQuantumBridge.lean` uses
the polynomial finite-particle domain of unbounded oscillators; its
canonical commutators are not full identities on F_3. The finite-cap
completion and the explicit Q-commutator boundary replace that domain
assumption in this source.

`PositionSpaceMaxwellAction`, `TemporalMaxwellEvolution` and
`DiscreteCoulombGreen` have committed twelve-port/face domains. Their
finite algebraic mechanism is applicable after the cube's own incidence
and exact-sequence identities are checked. Their particular matrices,
time recursions and field normalizations do not define this Hamiltonian.

`SeamCurrentFreePhotonLift` gives a rank-two transverse fibre for nonzero
Cartesian momentum and a conditional scalar oscillator action.
`SeamMaxwellContinuum` bounds the specified FZ-12 symbol; it is not a
dispersion theorem for this cubic complex. `ModalMaxwellFactorizationBoundary`
states a modal curl factorization without a local gauge or source join.
Those statements alone supply neither this finite carrier nor the
coupled Dirac-current detector law. `SeamCurrentAuxiliaryOscillatorLift`
is an algebraic first-order lift without a physical trajectory or clock.

The Whitney interacting quantum and Gaussian-state executables use real
configuration spaces and unbounded local Hilbert spaces. Their finite
mesh does not make their quantum carrier finite dimensional. No Gaussian
or continuum source is substituted for the literal 56-state field.

The common member fixes its own source clock and edge/charge scale.
Any growing-box massless-field or helicity theorem requires an explicit
family and record-limit proof. Such a theorem is separate from the
nonzero finite charged-current/field record on the present coupled H.
