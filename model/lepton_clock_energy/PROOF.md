# Lepton cycles on a complete native event timeline

The declared source is the tensor product of the 66-dimensional native
clock carrier and a three-dimensional lepton register with fresh pairs of
60-address cells. It uses the literal C and E native clock channels from
`clocked_common_observer`, and the complete weak and mass instruments from
`lepton_cycle`. The shared register identification is explicitly
`weak j=k-1 <-> mass character k`. It is an architectural wire map.

## Source, references and compression

The actual address partition has sizes 50, 5 and 5. Let B contain its three
normalized indicator vectors. Then B*B=I. Every prepared source column,
native kick and return projection used here preserves
`range(B tensor B)`. Pulling the literal 3,600-address source columns
through this isometry gives the nine-coordinate columns stored in the
receipt. Applying `B tensor B` reconstructs every original source vector;
no return environment is traced out. This equivalence holds on the stated
source preparation, not on arbitrary 3,600-dimensional inputs.

For kind t and return r, write the pulled column for family k as v[t,r,k].
The branch map sends `|k> -> v[t,r,k] tensor |k>`. Orthogonal return tests
give `sum_r ||v[t,r,k]||²=1`. With fresh source pairs, successive archive
columns are their tensor products, conditioned on the same k. This is an
isometry into the direct sum of complete return histories. It extends by
the identity on any quantum reference and retains cross-family coherence
inside each history. A final Fourier-basis measurement retains its outcome
and the resulting source archive. It is a complete three-outcome PVM.

The proof uses the finite Kraus completeness and composition argument in
`Lean/Dynamics/ChoiCPTP.lean`; a literal square-matrix instantiation needs
the finite tape embedding. The numeric array replay checks the stated
finite source identity but is not an interval certificate for hardware
errors or arbitrary floating density-matrix execution.

## Retained stopping protocol

In each controller slot, an independent activity coin chooses an identity
idle with probability 1-q, or executes the actual native instrument with
probability q. Idle, zero current, negative current and positive current
are distinct records. The first nonzero current executes one weak cycle;
the second executes one mass cycle. Both signs and all intervening states
are retained. The controller stops after the second event or at cap N.
Every incomplete path is a timeout, without conditioning on success.

At N=3 there are 20 two-event clock histories and 32 timeout histories.
Appending the complete return instruments and final detector gives 1,272
records. The six declared runs retain 7,632 records, including zero weights
at q=1. Native source output matrices are stored exactly as rational
entries, independently checked against execution of the original channel.
Archive columns and these native matrices specify the complete joint
output for the product preparation in the receipt. Arbitrary-input and
reference completeness follows separately from composition of the maps.

The zero-current effect satisfies `T_zero*(I)<=13I/15`. Conditional on any
past history with positive weight, the next event probability is therefore
at least delta=2q/15. Coupling each Bernoulli trial to an independent uniform
variable gives stochastic domination by Bernoulli(delta), even with an
adaptively changing quantum state. Consequently

`Pr(fewer than two events by N) <= (1-delta)^N + N delta (1-delta)^(N-1)`.

The bound holds with untouched references because the effect inequality
tensors with their identity. It is a bound on unresolved mass, not an
execution of an infinite history. At q=1 the executed cap-three timeout
is exactly 15806/16875. The native attempt count, event count, controller
slot count and coherent service operations are separate resources.

## Waiting and reset controls

The event-only policy performs no register evolution between events.
Inserting one deterministic idle before the entire protocol leaves its
complete output probabilities unchanged while incrementing its slot count.
The slot-work policy additionally applies `exp(-i h J)` after every slot,
where h=1/7 and J=diag(-1,0,1). These diagonal operations commute with the
controlled source columns, so their accumulated phase can be applied at
the end without changing the map. The inserted idle changes the final
Fourier detector law. This is a constructive distinction between declared
service schedules, not a proof that time is an event count.

Fresh pairs are supplied preparation services. The countercontrol instead
passes the full correlated source state after the weak return into the
mass kick. Its sixteen joint return probabilities differ from the fresh
pair law by total variation approximately 0.00234192393460. All source
vectors and failures survive in both laws. Thus a fresh reset cannot be
silently replaced by reuse of an entangled source, or by erasure.

## Energy and calibration

The limiting weak generator is diag(-2/9,0,2/9), while the mass generator
is diag(e_0,e_1,e_2), with the rational effect intervals proved by
`lepton_mass_compiler`. For all six family permutations (a,b,c), the
interval for `e_a-2e_b+e_c` excludes zero. Hence no affine identification
of these generators is possible. Distinct generators are compatible with
the explicit common timeline above. Finite cycles are dephasing channels;
the complete reference-inclusive approximation bounds from
`lepton_cycle/composition.py` apply before any Hamiltonian-limit reading.

For any positive c, replacing a dimensionful Hamiltonian H by cH and its
time unit t by t/c preserves `exp(-iHt)` exactly. Every dimensionless native
map, service angle, complete history and recorded count in this protocol
is unchanged. This exact units freedom excludes determination of an
absolute GeV or SI conversion from these records alone. It does not
exclude a declared internal relation between model coordinates. An
empirical electron-mass anchor fixes one energy scale; it neither identifies
the interacting particle nor measures the duration of a native attempt.
No physical mass, coupling, laboratory clock or autonomous preparation
device is inferred from this finite join.
