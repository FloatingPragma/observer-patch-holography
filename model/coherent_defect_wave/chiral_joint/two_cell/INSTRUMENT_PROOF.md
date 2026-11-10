# Complete compensated charged instrument

This instrument measures finite source energy transitions using the same
two-cell Hamiltonian and its complete physical ground. It keeps the
preparation reference, boundary charges, unsuccessful branches and original
probabilities. The reservoir operation is driven; no apparatus work or
relativistic rest mass is inferred from it.

## Physical preparation and its retained reference

The complete neutral physical Gauss space has dimension 1053. Choose an
orthonormal basis p_j with definite source L, 2Sz and fermion parity. The
neutral graded preparation ancilla has matching parity and opposite L and
2Sz for its corresponding state j. The purification

`Omega_ref=(1/sqrt(1053)) sum_j p_j tensor |j>_anc`

has total Gauss charge zero, total L and 2Sz zero and even total parity.
Its source marginal is exactly `P/1053`. A specified finite unitary takes
the empty source and a neutral even ancilla state to this purification,
acting within their common conserved sector and extending as the identity
outside that sector. The basis can include the empty source vector so the
initial ancilla state belongs to the declared register. No efficiency or
native-idle preparation claim is made.

The full spectral certificate proves that the same physical Hamiltonian has
a unique ground in `(L,2Sz)=(2,0)`. The complete two-outcome preparation
instrument is `{P_ground,I-P_ground}`. Starting from `P/1053`, its original
ground probability is `1/1053`, and its retained failure probability is
`1052/1053`. The source success state is the actual unique ground. Its
retained ancilla has L minus two, 2Sz zero and even parity, so preparation
has not borrowed a conserved quantity from an unrecorded reservoir.
The global energy projector commutes with local Gauss, L, 2Sz and parity.

## Full apparatus carrier and actual permutation

Choose source e^c mode `m=0,3,6,9`, covering both cells and both spin
components. Append reservoir CAR mode 12 and compensator CAR mode 13.
The full fermion dimension is 16384; including the six link and two
five-state scalar factors gives 2457600. These factors are spectators
during the swap. The preparation reference remains retained separately.

At the chosen source cell, r carries the same hypercharge plus one, L minus
one and 2Sz as the chosen e^c mode. The static compensator b has opposite
hypercharge, L and 2Sz. Both are weak/color singlets. The two recorded initial
settings are `(n_r,n_b)=(1,1)` and `(0,0)`, each with probability one half.
Both contribute zero net gauge charge, L and 2Sz and have even parity.

The driven unitary is the exterior lift of the one-body permutation swapping
m and r and fixing every other mode, including b. Every one of its 16384
fermion columns is computed with the exact permutation sign. The code checks
its inverse, occupation and parity conservation, and all fourteen
annihilator intertwiners on all columns. Thus doubly occupied, compensator
occupied and arbitrary source occupation sectors are included.

The actual fourteen-mode one-body Gauss generator contains the source
current at each endpoint and, at the probed endpoint, reservoir and opposite
compensator currents. The permutation commutes with all twenty-four local
generators, checked using the actual native matrices. Its exterior lift
therefore commutes with their full second-quantized currents. The scalar
and quantum-link currents commute because the unitary acts as identity on
those factors. The same argument verifies total L and 2Sz. No equality
between a charged source creator and a source-only gauge-invariant operator
is assumed.

## Complete branches and boundary sectors

After the swap, read r and retain the unchanged b. At the chosen endpoint,
the apparatus boundary hypercharge is

`q_boundary=n_r,out-n_b=n_r,out-initial_setting`.

The add-success branch has boundary minus one, its source L changes by
minus one, and source 2Sz increases by the chosen mode's spin sign. The
remove-success branch has boundary plus one, source L increases by one and
source 2Sz decreases by that sign. Both unsuccessful branches retain
boundary zero and the original source L and 2Sz. These apparatus charges
are present during the subsequent source energy measurement. The correct
physical block is therefore selected by the recorded boundary, rather than
projecting every charged branch into the neutral source kernel.

For the even-parity neutral ground, the reduced Kraus operators in the
fixed source-before-reservoir-before-compensator ordering are

`a^dagger/sqrt(2), -n/sqrt(2), -a/sqrt(2), (I-n)/sqrt(2)`

in add-success, add-failure, remove-success, remove-failure order. Their
signs are derived by executing the actual permutation on every even source
basis state, not selected after comparing probabilities. Their effects sum
to the identity. The all-parity instrument is also checked directly on all
4096 source columns: each retained branch is an injective signed partial
permutation, and each source column has two unit-weight alternatives before
the common setting factor one half. The displayed simple Kraus signs are
not asserted as an unqualified cross-parity formula; the full permutation
defines that case.

The source energy effects are the complete spectral projectors of the same
Hamiltonian in the appropriate recorded boundary sector. They commute with
the total Gauss constraint, L, 2Sz and parity. For a conditional ground
reference rho_0, an original swap/energy record has probability
`Tr(P_E K rho_0 K^dagger)` with the `1/sqrt(2)` already in K. To obtain the
unconditional preparation-plus-readout law, multiply these successful
preparation records by `1/1053` and retain the `1052/1053` preparation-failure
record. No success-conditioned probability is silently presented as an
original full-protocol probability.

The spectral module represents physical vectors by its phased embedding
`D B`, with `D=i^(N_left)`. Probe matrices must use that same returned
embedding on both ends. Discarding D would substitute a different operator
while leaving its reported energy labels unchanged.

This finite instrument provides complete gauge-compatible charged energy
records. Its conservation laws do not imply conservation of source energy
during the driven swap. Relativistic dispersion, removal of apparatus and
boundary effects, and a physical energy convention remain separate tasks.
