# Four-cell full-CAR dynamical wall source

Use L4 spatial cells, S2 regulator slices, one physical spin and all16
CAR modes. The h=1/8 Wilson/fifth/boundary operators and compact U1
electric coefficient1/4 are unchanged. Keep both bare-mass cases, zero
and the actual432 electron coefficient. No holonomy pinning, vacuum
subtraction, critical mass or target-selected wall length is added.

Every N0..16 sector belongs to the source. Static references are −S at
each vertex and an additional −(N−LS) at vertex zero. Exact Gauss fields
are `E_x=ell+sum_(v<=x)n_v+S(L−x−1)−N`. Tree hops and the closing
holonomy act on every allowed occupation. The field energy is the sum
of all four E_x squared divided by4.

Use matrix-free sector actions at flux cutoffs K2 andK4. Resource ceilings
are16 CAR modes, K<=4 and250000 states per individual sector. At K4 the
full direct sum has589824 states and its largest sector has115830.
All sectors are retained analytically even when a uniform lower bound
excludes them from the ground search.

At each holonomy the full matter spectrum includes every wall band and
all four spatial momenta. Positive electric energy yields an exact
sector lower bound from its complete matter minimum. An interval
holonomy grid and a global derivative bound enclose that minimum.
For lower-bound verification also retain every Barta test
`f(theta)=exp(a cos(theta))`, a in {0,±1/4,±1/2,±1,±2}. These functions
alter no source parameter. Vector diamagnetism bounds the complete
electric energy below by the scalar kinetic energy of the wavefunction
norm, supplying the Barta inequality on the same source.
Compare every sector bound to an N8 variational energy before claiming
the global vacuum number. Sparse low eigenvectors alone do not certify
neutral-sector spectral order.

A warm Krylov seed can lie in a neutral symmetry block. The retained
warm-seed control has inconsistent second Ritz values across nearby
bare masses. Every subsequent warm seed is supplemented by a normalized
generic vector with rational amplitude1/16, orthogonalized against the
warm vector. This is a numerical verification choice and changes no
Hamiltonian parameter. Sparse neutral spectral order stays uncertified.

For the actual finite variational preparation, solve both complete
charged resolvents for the two local boundary Dirac probes at vertex
zero. If the sector lower bounds put the reference energy below both
charged spectra, use the positive gap to bound the infinite-flux
resolvent error from its full interior and exterior residual. Preserve
all preparation amplitudes, occupation/flux probabilities and CAR
weights. This controlled prepared-state response is distinguished from
an exact vacuum propagator and from a physical lepton pole.
