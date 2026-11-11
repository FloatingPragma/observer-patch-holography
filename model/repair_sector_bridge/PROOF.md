# Finite repair-to-quotient contract

## Native source quotient

`bridge.native_quotient` executes every branch returned by
[`coherent_continuation.branches`](../observer/coherent_continuation.py)
on the 66 binary two-load states. There are 30 seams and two coin values,
so each branch carries probability `1/60`. Self-transitions are retained.
A declared classical diagnostic reads the graph distance between the
occupied ports from the full source state. The three fibres have
cardinalities 30, 30 and 6. This diagnostic is not a compiled measurement
from the fixed observer's native current-only records.

For every pair of source states in the same fibre, the checker compares the
sum of their outgoing probabilities into each target fibre. Equality proves
strong lumpability. The exact native kernel has uniform stationary measure;
its quotient therefore has weights `(30,30,6)/66`. Detailed balance,
irreducibility and positive self-loops are checked on the actual native
kernel and its quotient. The finite chain consequently converges to that
unique law from every initial distribution.

In distance order `(1,2,3)`, the quotient is

    K = [[14/15, 1/15, 0], [1/15, 9/10, 1/30], [0, 1/6, 5/6]].

Exact multiplication gives `K**2(i,j) >= (11/450) pi(j)`. The usual
minorization decomposition therefore bounds the distance-record law after
`t` native attempts by

    TV(mu K**t, pi) <= (439/450)**floor(t/2).

This is an explicit statistical error for the declared diagnostic on the
native attempt clock. It is not a bound for every microscopic observable,
does not supply its current-only measurement implementation, and does not
identify that attempt clock with physical time.

This directly supplies the representative-independence and measure
pushforward premises of the
[finite-quotient ensemble theorem](../../paper/tex_fragments/FINITE_QUOTIENT_ENSEMBLE_THEOREMS.tex).
No independent equilibrium simulation or fitted Gibbs parameter is needed.
The distance fibres are classical pair configurations, not the irreducible
electric sectors of a gauge register. Their counting measure carries no
automatic transfer to the one-sided gauge-edge trace.

## Sector proposal measure

Let `q` be a finite stochastic proposal reversible for a positive measure
`m`, let `e_i` be nonnegative integer energies, and choose rational
`0 < y <= 1`. Accept a proposed `i -> j` with

    a_ij = min(1, y**(e_j-e_i)).

The off-diagonal accepted kernel obeys

    m_i y**e_i q_ij a_ij
      = m_i q_ij min(y**e_i,y**e_j)
      = m_j y**e_j q_ji a_ji.

Rejected probability is recorded separately and returned to `i`. It follows
that `pi_i` proportional to `m_i y**e_i` is stationary. Connectivity of the
positive proposal support gives irreducibility, because all accepted edge
probabilities are positive. Positive self-loops then suffice for convergence.
Stationarity alone proves neither uniqueness nor convergence.

The SU(2) controls use labels `n=2j=0,1,2,3`, dimensions `d=n+1`, and
`e=12 C2=3n(n+2)`. At the cutoff, excluded proposals stay at the source.
Neighbour proposals proportional to `1`, `d_j`, and `d_j/d_i` are reversible
for `1`, `d_i`, and `d_i**2`, respectively. Their accepted laws have the
corresponding dimension powers. Only the one-sided control leaves
`d_i y**e_i` stationary. To use a different target with the same proposal,
the acceptance needs the full Hastings ratio. An identity proposal leaves
every distribution stationary and cannot select any of them.

## Scope of the bridge

The exact checker rejects floats, negative probabilities, non-stochastic
proposals and partitions that are not lumpable. It distinguishes complete
records from accepted records and reversible equilibrium from unique
selection. Its tests mutate these premises, including an invalid proposal
whose excess probability would be hidden by acceptance.

The sampler's energy penalty and temperature parameter are declared inputs.
Realizing a Casimir-weighted stationary law does not derive the group
Laplacian as the physical Hamiltonian, a modular clock, or a gauge coupling.
It also does not prove a native operation that changes electric sectors.
A source-to-gauge transfer needs the actual sector projector map, proposal
measure, sector-changing events and physical clock on one admitted source.
