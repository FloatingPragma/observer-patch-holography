# An autonomous finite collision clock and its phase-onset boundary

## Native operators and complete carrier

Take N copies of the declared longitudinal matter source, its quantum Gauss
reference, and the neutral fermionic memory of the native rest-clock
construction. On copy i the bounded even operator A_i raises the actual
chosen charged branch from the source vacuum and creates the memory fermion.
It includes the original Gauss dressing. Its complete source identity is
`[H_i,A_i]=m A_i`; its inactive Fock and constraint complements are retained.
Write n_i=A_i A_i*, n=sum n_i and J+=sum A_i. Operators on distinct copies
commute because each A_i is even. Each active pair is a spin one half, while
the two inactive displayed pair occupations and all original source
complements are spin-zero blocks. No collective-spin projection is placed
in the action.

Declare the bounded native-current function

`a=J- (N-n+1)^(-1/2)`.

The denominator is positive on every retained occupation. On the symmetric
active sector reached from the native vacuum, `a|j>=sqrt(j)|j-1>` for
j=0,...,N. This is a representation consequence of the actual collective
currents. Other spin sectors remain part of the operator. Their ladder
coefficients satisfy `J+J- <= n(N-n+1)`, so `||a||²<=N` and
`||a²||²<=N(N-1)` on the complete carrier. The source energy differs from
m n by an operator commuting with all the collision operators. Therefore
`[H,a]=-m a` globally, including the inactive complements.

The normalization function and the pump/loss law below are declared
architecture, rather than an OPH uniqueness claim or a particle fit.

## A complete energy-conserving collision

Let gamma=1/16, kappa=2 gamma/N and delta=2^(-16). Define

`S=gamma a a* + kappa a*² a²`,

`K0=sqrt(I-delta S)`, `K+=sqrt(delta gamma) a*`,
`K-=sqrt(delta kappa) a²`.

For each executed N, `delta [gamma N+kappa N(N-1)]<1/4` is checked. Hence
the square root is defined on the entire carrier and
`K0* K0+K+* K+ +K-* K-=I`. Each collision is followed by the same free
native evolution exp(-i H delta). This specifies every outcome and its
complement; there is no successful-jump postselection.

A fresh three-state battery has energies m,0,3m and starts in the first
state. Label its three outgoing records 0,+,-. The isometry

`V psi = K0 psi tensor|0> + K+ psi tensor|+> + K- psi tensor|->`

intertwines total energy with H+m. It extends within each energy eigenspace
to an energy-conserving unitary, retaining the unused battery inputs.
The mean source work of one collision, on any input rho, is exactly

`delta m [gamma Tr(rho a a*) - 2 kappa Tr(rho a*² a²)]`.

The battery loses this energy. A sequence uses a fresh nonpassive battery
per step and may change their entropy even when the net stationary mean
work vanishes. It is not an equilibrium or resource-free clock. The source
parameter is the native collision/event count; conversion of m to GeV and
seconds uses the inherited empirical electron anchor and SI definitions.

## Exact frequency and finite linewidth

On the symmetric sector let s_j=gamma(j+1) for j<N, zero pump rate at j=N,
plus kappa j(j-1). The k-th coherence band rho_(j+k,j) evolves under a real
nonnegative matrix B_k times the sole complex phase exp(-i k m delta).
Conjugating B_k by `D_j=sqrt((j+k)!/j!)` leaves diagonal entries
`sqrt[(1-delta s_(j+k))(1-delta s_j)]`; its pump and two-loss off-diagonals
are respectively `delta gamma j` and
`delta kappa(j+k+2)(j+k+1)`. All these entries are computed from the source.

For a positive dyadic witness v, exact rational square-root enclosures give
the Collatz bounds

`min_j (B_lower v)_j/v_j <= spectral_radius(B_k)
 <= max_j (B_upper v)_j/v_j`.

The positive dyadic witnesses are fixed certificate inputs, not accepted
solver answers. Every replay verifies them against the actual exact
collision entries. A hostile witness that cannot certify the bound fails.
Changing a LAPACK frame therefore cannot change the rational receipt.
The executable verifies a strict interval inside (0,1). The dominant
coherence therefore has exactly the source frequency k m, with a nonzero
finite linewidth. For the first band and N>=4 the positive directed graph
is irreducible and has positive diagonal, so its leading Perron root is
real, positive and simple. The fixed collision step resolves all three
inherited charged masses without phase aliasing. Bounds
`1-r <= -log r <= (1-r)/r` give a rational outward linewidth interval.
No floating logarithm or fitted sinusoid supplies the certificate.

The phase-insensitive rates change this linewidth and the populations, but
cannot change the separated factor exp(-i m delta). A differential
Hamiltonian drive proportional to n changes the source gap itself; a
coherent drive proportional to a+a* breaks this covariance. Neither is
asserted to preserve the original frequency.

## Onset, shared records and the mass domain

The vacuum population has a strictly positive pump probability delta gamma
at the first collision. It evolves autonomously to a populated stationary
state without a supplied phase or a time-dependent driving coefficient.
The entire channel is covariant under exp(i theta n). Starting
from the phase-invariant vacuum, every unconditional state remains phase
invariant, so `<a>=0` at every step. This is a statement about the ensemble,
not a requirement that observer records must have a nonzero ensemble mean.
An ensemble with uniformly distributed phases can have oscillatory shared
conditional records. The actual acceptance questions are their autonomous
onset, frequency and drive response, with the apparatus and battery
resources included.

In the executed active symmetric sectors the nonzero linewidth excludes a
permanent oscillating one-point mode. A vanishing-linewidth thermodynamic
limit, spontaneous phase selection, and its order of limits would require
additional proof; that limit is not silently made an additional criterion.

For N=1 the same complete source autonomously prepares the actual single
charged branch from vacuum, with its inherited rest-energy gap m. The
two-loss operator vanishes, and the exact coherence multiplier is
`sqrt(1-delta gamma) exp(-i m delta)`. Its phase therefore uses the same
single-particle source gap, while its finite linewidth and conditional
readout preparation remain explicit. No phase was supplied to the pump.

For N>1 the carrier is an array of native charged branches. Its first-band
frequency m is the constituent rest gap, not the total array energy
m times its occupation. This narrowing-array result cannot by itself be
identified with an internal oscillator of one particle. Such a transfer
requires an explicit one-particle internal-state construction. The present
source is not an equilibrium ground-state time crystal; its complete
conditional records, rather than the zero ensemble mean alone, determine
what it establishes about an autonomous observed clock.
