# Native internal currents on the configured fermion operators

The finite same-current attachment is supplied by
[`native_matter_attachment.action_intertwiner`](../common_source/native_matter_attachment.py):
its permutation `P15` maps the selected auxiliary occupation-mask order to
the existing kinetic/Higgs action. The separate field architecture has

    h = C2_spin tensor M15 tensor C3_family,
    mode(spin,internal,family) = 45 spin + 3 internal + family,
    P90 = I2 tensor P15 tensor I3.                     (1)

The fifteen auxiliary masks have even occupation degree in a five-mode
register. The ninety propagating CAR modes are a separate configured
architecture. Their Fock space is `F(h)=direct_sum_r wedge^r h`, of dimension
`2^90`. Equation (1) is not an embedding of this space into the thirty-two
states of the auxiliary register. Neither the odd field statistics nor the
three-family multiplicity is inferred from the even auxiliary parity.

The owning `kinetic_attachment()` verifies the same `P90` against all three
supplied Weyl principal matrices, both Higgs components and twenty-seven
arbitrary complex family coefficients. The
[`native_gauge_matter_transport`](../common_source/native_gauge_matter_transport.py)
implementation compiles the same internal currents and background-linked
hopping on separate CAR modes, including full local Fock factors and routing
signs. These results are reused. The operator formulas below do not execute
another native tape or identify these modes with the default two-flavour
interacting centered-ground source.

## The complete exterior map and its current action

For a unitary one-particle map `P`, define

    Gamma(P) = direct_sum_r wedge^r P.

If `P e_j=e_pi(j)`, its action on an ascending occupation basis is

    Gamma(P)|j1<...<jr>
      = sign(pi(j1),...,pi(jr)) |sort(pi(j1),...,pi(jr))>. (2)

The sign is the parity of inversions in the displayed image sequence.
This is a bijection of complete orthonormal bases, preserves occupation
number and parity, and includes both the vacuum and fully occupied state.
Wedge multiplication gives, for every field vector `f`,

    Gamma(P) a^dagger(f) Gamma(P)^* = a^dagger(P f).     (3)

No one-particle or dilute-state assumption enters (3). In particular,
changing an exterior sign to a bosonic permutation preserves the
one-particle map but can fail this identity on occupied inputs.

For a Hermitian one-body matrix `T`, set

    dGamma(T) = sum_ij T_ij a_i^dagger a_j.

The full CAR relation gives the exact operator identity

    [a_i^dagger a_j,a_k^dagger] = delta_jk a_i^dagger,
    [dGamma(T),a_k^dagger] = sum_i T_ik a_i^dagger,
    [dGamma(T),a_k] = -sum_j T_kj a_j.                 (4)

The degree-three terms cancel using anticommutation; they are not discarded
by a particle-number truncation. The same formulas hold for a general
complex coefficient matrix. The native port currents are anti-Hermitian
`G_A`; their physical Hermitian convention is `T_A=iG_A`, as in the existing
compiler. With `T_A^action P90=P90 T_A^native`, equations (3)-(4) imply

    Gamma(P90) dGamma(T_A^native) Gamma(P90)^*
      = dGamma(T_A^action),                           (5)

for all twelve currents on all `2^90` occupation states. The all-left
singlet-conjugate convention is inherited from the source representation.
The dual transformation of annihilation operators does not by itself
prepare an antiparticle; a vacuum/reference and its particle-hole map are
separate inputs.

The executable checks all twelve source intertwiners, exhausts all 32768
occupation columns of `Gamma(P15)` and all 491520 corresponding creation
columns, and checks ninety-mode vacuum, filled and spectator-rich columns.
It also applies every actual current directly to four occupation patterns
of the fifteen-mode CAR algebra, giving 720 commutator-column controls.
These finite controls test the conventions in the universal algebraic
proof; their count is not an exhaustive enumeration of `2^90` states.
Independent tensor-Pauli tests cover every primitive three-mode identity,
including occupied spectators and zero creation columns.

## Charge labels of the field operators

The existing source-weight calculation fixes `q=6Y`, `T3`, `Q=T3+Y`, two
color Cartan weights and the color/weak Casimirs before inserting target
names. Equation (4) transfers each diagonal weight to the corresponding
creation operator on every occupation sector. For example,

    [dGamma(Q),a^dagger_spin,internal,family]
       = Q_internal a^dagger_spin,internal,family.      (6)

The two spatial spin indices and three declared family indices copy the
same internal representation; they do not create additional internal
multiplet types or identify electron, muon and tau mass eigenstates.
`field_operator_charge_records` retains all fifteen source masks, their
six field indices each, the action-order indices, and the existing seven
internal comparison-component labels.

Casimir labels refer to this operator representation. For a normalized
current basis `T_a`, they are established through

    sum_a [dGamma(T_a),[dGamma(T_a),a^dagger(f)]]
       = a^dagger((sum_a T_a^2)f).                    (7)

They are not an identity between `dGamma(sum T_a^2)` and the many-body
Casimir `sum dGamma(T_a)^2`. The latter contains cross-particle terms.
The color basis in the literal canonical-current order has six off-diagonal
matrices of twice the usual normalization and two simple-root Cartan
matrices. Writing these Hermitian matrices as `H_0,...,H_7`, its Casimir is

    sum_(a=0)^5 H_a^2/4 + (H_6^2+H_7^2+H_6 H_7)/3.

The weak matrices `H_8,H_9,H_10` have the usual one-half normalization.
The code compares both Casimirs with the independently computed auxiliary
source records exactly.

For any chosen mode the ideal addition/removal instrument
`{a_j^dagger,a_j}` is complete, since its effects are `1-n_j,n_j`.
Its two original probabilities are `1-omega(n_j),omega(n_j)`. They need not
be equal, and neither is normalized to one. This algebraic availability
does not constitute a ground-preparation or native instrument-execution
certificate for the configured field.

The attachment supplies an explicit current and fermion-operator dictionary
on an existing configured architecture. It selects no Yukawa matrix, family
eigenbasis, interacting vacuum, particle pole, physical species or mass.
An experimentally identified current and propagating sector require the
same state, dynamics and observation dictionary; agreement of the finite
operator representation alone does not supply them.
