# Full fermionic feedback for the coherent charged-family action

The coherent analyzer in
[`KOIDE_SOURCE.md`](../standard_model_source_map/KOIDE_SOURCE.md) supplies a
unitary on eight flag coordinates and three family coordinates. Its
blank-input, flag-zero block is `K0`, and `Y=K0* K0=C_ex^2/3`.
The construction here lifts that **entire unitary** to fermion Fock space
before compressing any coefficient. It gives a complete weak-feedback
limit for the charged Higgs pairing on every physical occupation sector.
The added modes, their fresh vacuum preparation, the controlled Hamiltonian
and their quantum archive are supplied apparatus. There is no physical
vacuum, particle pole, phase-selection or laboratory-clock conclusion.

## The one-particle dilation and its exterior lift

Use the configured physical space `h=C2_spin tensor M15 tensor C3_family`
and its ninety CAR modes. Extend the family coordinate to
`C8_flag tensor C3_family` in every internal/spin component. The enlarged
space has 720 modes. The physical modes occupy flag zero, with all 630
other modes initially vacant. The vacuum embedding

    J: F(h) -> F(h_extended),   P=J J*,

preserves the full physical occupation space, rather than selecting a
one-particle input. A flag coordinate is gauge neutral. Its extra fermionic
copies carry the same internal representation as the original modes;
they are not electrically neutral particles.

Let `U` be the complete 24-dimensional analyzer, and let
`Q=|flag0><flag0| tensor I3`. Then

    W=U* Q U,   W^2=W,   rank W=3,
    E* W E=K0* K0=Y,                                  (1)

where E embeds the original three family modes into flag zero. The code
reconstructs the actual source's flag-zero row, verifies that it is a
coisometry on all eight input flags, and compares (1) with its exact Y.
It does not rerun or replace the source's native analyzer execution.

On the entire lepton doublet use the family/flag map U; on the singlet
conjugate use `conjugate(U)`. Use identity on all other M15 components and
on the spatial spin coordinate. Call this unitary S. Both lepton projectors
commute with all twelve native currents, so S and its exterior lift
`Gamma(S)` commute with their full-Fock charges. The doublet's neutral
component must share its family frame. It cannot be transformed separately
without changing the weak currents.

For creation pairing `H_Delta=1/2 sum Delta_ij a_i* a_j* + h.c.`, the CAR
functor gives

    Gamma(S)* H_Delta Gamma(S)=H_(S* Delta conjugate(S)). (2)

If the oriented doublet/singlet block of the kick is Q, its transformed
block is `U* Q U=W`; its reverse is `W^T`. Using U on both legs instead
would give `U* Q conjugate(U)`, which need not equal W. The finite complex
control explicitly rejects that substitution. Equation (2) uses the
unitary S, never an identity between `Gamma(sum A_r)` and `sum Gamma(A_r)`.

## The actual Higgs pairing and compression on every occupation

Let D be the native directed charged-lepton Higgs tensor for unit neutral
Higgs norm. It has rank one and disjoint doublet/singlet supports. In the
source's all-left convention the kick is

    Delta_kick=epsilon_spin tensor
               [conjugate(D) tensor Q + conjugate(D)^T tensor Q],
    V=Gamma(S)* H_(Delta_kick) Gamma(S).                (3)

The physical CAR map and the transpose/conjugation conventions are those of
[`FERMION_ATTACHMENT.md`](../standard_model_source_map/FERMION_ATTACHMENT.md)
and `native_matter_attachment.kinetic_attachment`. Vacuum compression kills
each monomial containing an ancillary creation or annihilation operator.
The remaining terms are exactly

    J* V J=H_Y,
    Delta_Y=epsilon_spin tensor
            [conjugate(D) tensor Y + conjugate(D)^T tensor Y^T]. (4)

This is the source's charged Higgs pairing with family coefficient Y;
it is an operator identity on all `2^90` physical states. It is not an
identity only on a one-particle wavefunction. Pairing changes total number
by two and preserves total fermion parity. At neutral fixed Higgs it also
preserves the unbroken electric charge. The two Higgs components transform
covariantly under all twelve currents; a fixed neutral Higgs does not
commute with all the broken generators. Ancillary vacuum and P are gauge
invariant, so compression respects this covariance.

The kick in (3) is a sum of six disjoint two-mode pair operators: two
spatial-spin orientations and three selected family modes. Each pair has
eigenvalues `-1,0,0,1`. Their even operators commute and their maximum
eigenvalues can be attained simultaneously. Hence

    ||V||=6 |v_H|=:L,
    ||H_Y||=2 |v_H| Tr Y=5 |v_H|/9.                    (5)

Unitary conjugation preserves the first norm; diagonalizing the positive
Y gives the second. The neutral rank-one tensor normalization is checked
against both actual source Higgs components. No independent normalization
of the charged-lepton channel is introduced.

## Complete reset archive and the weak limit

Apply `exp(-ih V)` to the physical input and fresh ancillary vacuum. Record
whether the ancillary output is vacuum. Retain the entire quantum state
of every nonvacuum output, including its charges and coherences. Transfer
that output to an archive and use fresh vacant modes for the next round;
this is not erasure or deterministic reset of an unknown charged state.
The vacuum/nonvacuum measurement commutes with the gauge action and parity.
A finer occupation-basis readout is complete but need not preserve all
nonabelian currents; it is not used for the gauge-covariant conclusion.

The vacuum block and its complementary isometry are

    A_h=J* exp(-ih V) J,
    B_h=(I-P) exp(-ih V) J,
    A_h* A_h+B_h* B_h=I.                               (6)

For arbitrary physical input and untouched reference,

    ||A_h-exp(-ih H_Y)|| <= h^2 L^2,
    ||B_h|| <= |h| L.                                 (7)

The first bound follows by subtracting the first-order expansions: their
linear coefficients agree by (4), and each integral second remainder is
at most `h^2 L^2/2`. The second bound is Duhamel with `(I-P)J=0`.
Both are uniform over every occupation sector. For N rounds, `h=t/N`,
contraction telescoping bounds the all-vacuum block by `L^2 t^2/N`.
The union bound gives total original probability of any nonvacuum record
at most `L^2 t^2/N`, even after arbitrary prior failures and correlations.
All those nonvacuum quantum histories remain in the comparison output.

The ideal comparison applies `exp(-it H_Y)` and retains an all-vacuum
archive. The zero-record density error is at most `2L^2 t^2/N` in
unhalved trace norm. Other classical history blocks are positive and have
total trace at most `L^2 t^2/N`. Therefore the complete instrument obeys

    ||actual-ideal||_1 <= 3 L^2 t^2/N
                       =108 |v_H|^2 t^2/N.            (8)

This includes arbitrary references, all reset flags and the full quantum
complement, with no success conditioning. Before recording the vacuum
flags, the coherent archive has the weaker isometry norm bound
`sqrt(2) L |t|/sqrt(N)`, hence unhalved output error at most
`2sqrt(2) L |t|/sqrt(N)`. Classical flag commitment is an explicit operation
in (8), not a claim that the coherent dilation loses its information.

A common Hamiltonian H0 can be included if it reduces the physical-vacuum
subspace and its physical restriction is the desired target H0. In the
interaction picture, `P V(s) P` equals the interaction-picture H_Y.
The same first-order cancellation and integral bounds prove (7)-(8)
without a factor `||H0||`. For an analyzer implemented by conjugating a
whole pulse, the extended H0 must also commute with S, or its compensation
must be supplied. Family/flag-independent free kinetic copies obey this
condition; an arbitrary interacting H0 is not assumed to do so.

Certified per-round complete dilation error epsilon adds `2N epsilon`
to (8) and `N epsilon` to the nonvacuum probability bound. A supplied
input-state trace error adds directly. The model time t counts weak-kick
duration. The analyzer, uncomputation, fresh modes, vacuum readout and
archive transfer have their own cost. No equality to an elapsed native
clock, unbounded-rate service or volume-uniform limit is inferred.

## Exact finite controls and interpretation

`car_koide.py` checks both source phase controls, the complete analyzer row,
all twelve current support identities, both Higgs tensor normalizations,
and the exact coefficient (4). A separate eight-mode control has a complex
off-diagonal Y, four physical and four ancillary modes. It checks all 256
Fock states through an exact generator polynomial and full clock unitarity,
all sixteen physical compression columns and all sixteen ancillary output
labels. Its fine occupation instrument is an algebraic control, not an
assertion of nonabelian-invariant fine charge measurement. The production
bound uses the complete invariant vacuum/nonvacuum instrument instead.

The construction realizes the configured charged Higgs action dynamically
on full CAR at fixed finite volume. It does not alter the four distinct
current-gated feedback actions, select the coherent phase, suppress mirrors,
prepare an interacting vacuum, or identify an action singular value with
a measured lepton pole. The positive-chamber Koide equality of Y retains
those physical matching requirements.
