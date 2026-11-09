# A complete coherent source for a balanced family generator

This is a configured alternative to the four current-gated memory duties.
Its coherent flag preparation, controlled family rotation, phase reference,
and uncomputation are part of the action. The default phase is zero; a
separate exact rational unit-circle control tests a nonzero phase without
an empirical input. Neither phase is selected by classical repair records.
The construction supplies a finite logical family generator and its complete
record law. Its physical CAR/Higgs interpretation is specified separately.

## Native weights and a coherent analyzer

The literal C/E instruments in `native_repair_mass.native_records` have 128
labels, including 104 zero-weight labels. Write `A+` for C coin0/raw1,
`A-` for C coin1/raw2, and `B` for E coin0/raw4 or5. The disjoint events
`B only`, `A+ only`, `A- only`, their joint event, and the complement have
probabilities

    (10,5,5,2,50)/72.                                 (1)

All 128 original labels are checked before constructing a coherent service.
This service prepares square roots of (1) on five orthogonal flags, padded
to eight. It is not the classical coarse graining of a measured tape.
Ordinary classical commitment would remove the required interference.

The actual native half-turns of the two seams have product `gE gC` of
order three. Its action on the native geometric family band is a real
orthogonal matrix `R` with trace zero. The code constructs an orthogonal
frame `O` with `R O=O P`, where `P` is the cyclic permutation of three
coordinates. Gauge currents act as identity on this separate family factor.
All formulas below use that cyclic frame and a fourth padding coordinate.

For an explicitly supplied unit phase `z`, the five work operators are

    A0=sqrt(10/72) I,
    A+=sqrt(5/72) z P,       A-=sqrt(5/72) conjugate(z) P^2,
    AJ=sqrt(2/72) I,         AF=sqrt(50/72) I.          (2)

Their effects sum to identity. A three-point Fourier readout of the first
three flags, with the others fixed, gives a complete unitary analyzer `U`.
Its blank-input Kraus maps are

    K_m=(A0+omega^(-m) A+ +omega^(-2m) A-)/sqrt(3), m=0,1,2,
    K_3=AJ, K_4=AF, K_5=K_6=K_7=0.                   (3)

For a tracial family input, the original eight probabilities are
`(5/54,5/54,5/54,1/36,25/36,0,0,0)`. The three active probabilities sum to
`5/18`. Each individual effect `K_m^dagger K_m` depends on the input.
An output unitary cannot make all three Kraus maps equal, because it cannot
change these effects. No such replacement is used.

Let

    C_ex=A0+A++A-,     Y=K_0^dagger K_0=C_ex^2/3.      (4)

Here `C_ex` is Hermitian. Its exact traces are
`Tr(C_ex^2)=5/6`, `Tr(C_ex)=sqrt(5)/2`, and `Tr(Y)=5/18`.
Writing `z=exp(i theta)`, its eigenvalues are
`sqrt(10/72)+2sqrt(5/72)cos(theta+2pi k/3)`.
In the positive chamber `|theta|<=pi/12` modulo `2pi/3`, the action
coefficients `m_k=eigenvalue_k(Y)` satisfy

    sum m_k/(sum sqrt(m_k))^2 = 2/3.                  (5)

The executed exact controls are `z=1` and `z=(99+20i)/101`; positivity is
checked by exact principal minors. The zero-phase spectrum has a repeated
lower coefficient. Equation (5) does not select distinct lepton ratios.
Outside the positive chamber the signed-root identity alone is insufficient.

## Unconditioned weak feedback from the same analyzer

Let `Q` project onto output flag0, and define `P_U=U^dagger Q U`. It is a
projection on the complete flag/family register. Introduce a logical Nambu
qubit with Pauli `X`. One round prepares the blank flag, applies `U`, applies
`exp(-ih X tensor Q)`, uncomputes `U`, and reads and resets every flag. Its
complete isometry is the blank-input restriction of

    V_h=exp(-ih X tensor P_U).                        (6)

The vacuum-flag compression of `P_U` is exactly `Y` in (4). The retained
zero-return Kraus operator and the sum of all other effects are

    A0(h)=I+(cos h-1) I tensor Y-i sin h X tensor Y,
    sum_(r!=0) A_r(h)^dagger A_r(h)
       =4sin^2(h/2) I tensor(Y-Y^2) <= h^2 I/4.       (7)

These identities hold on all logical target inputs and external references.
Every other flag is retained, including padding flags of zero weight.
There is no conditioning on the return flag or division by `5/18`.
Using `Q` equal to all three active flags would instead compress to
`(5/18)I`, losing the nontrivial family generator.

The first derivative of (7) is `-i X tensor Y`. The second-order integral
remainders of the compression in (6) and of `exp(-ih X tensor Y)` each
have norm at most `h^2/2`, because both generators have norm at most one.
Hence

    ||A0(h)-exp(-ih X tensor Y)|| <= h^2.              (8)

For `N` fresh blank/reset rounds and `h=t/N`, contraction telescoping gives
`||A0(h)^N-exp(-it X tensor Y)||<=t^2/N`. The probability of any nonzero
return history is at most `t^2/(4N)`. The no-return density operator differs
from the ideal unitary output by at most `2t^2/N` in unhalved trace norm,
including arbitrary untouched references. Adding all other positive
history blocks gives the complete recorded-output bound

    ||actual - ideal unitary with all-zero archive||_1
       <= 9t^2/(4N).                                 (9)

This is a deterministic complete-instrument comparison; no independence
assumption about target states or selected successful histories is used.
Fresh blank preparation and reset are explicit services. An independently
certified per-round dilation error `epsilon` adds at most `2N epsilon`
to (9), by an isometry/channel telescope, and adds at most `N epsilon`
to the nonzero-history probability bound. A small floating probability
residual is not a certificate of `epsilon`.

## Native execution and the physical action boundary

The analyzer uses the existing native `compile_real_orthogonal`, diagonal
phase and tensor-Pauli controls. Controlled cyclic permutations are exact
products of conditional one-bit transpositions. The three-point Fourier
matrix has a real orthogonal eigenbasis and eigenvalues `1,-1,-i`, so it
uses the same real-basis/diagonal compiler. The weak kick is a product of
eight commuting tensor-Pauli rotations, since

    X tensor Q = X tensor product_(j=1)^3 (I+Z_j)/8.

`native_execution` runs the emitted words through the actual M6 generators
of `native_matter_attachment.execute_word`. It checks all 32 analyzer input
columns and all 64 weak-round input columns, with the family and flag
padding retained. Every primitive load/export has its indexed complement
failure outcomes; their weights are exactly zero on the invariant proper
code. They are not removed from the instrument's definition.

The executable compares floating native contractions with exact ideal
Kraus matrices at the declared admission. This is a numerical execution
control, not an interval certificate for accumulated roundoff. The exact
ideal formulas (1)-(9) do not use that numerical residual as a proof.
The finite control takes `h=2atan(1/64)` and retains each of its eight
original return probabilities. The analogue native controls, coherent
code transfers, blank preparation and final readout are specified apparatus.
No experiment on a physical detector or primitive default-source selection
is asserted. The time `t=Nh` in (9) is the accumulated weak-kick angle.
The displayed word uses 9,388 native pulses per round, including analyzer
and uncomputation. Their duration does not disappear when `h` shrinks.
No equality between this angle clock and global native-attempt time or
physical time is supplied; no physical continuum speed follows from (9).

The family coefficient `Y` can be inserted into the verified same-M15
Higgs invariant through `native_matter_attachment.kinetic_attachment`,
with canonical kinetics and the declared physical CAR functor. That map
accepts arbitrary complex family matrices and acts on the same two spinor
components, fifteen internal states and three neutral family copies.
Its quadratic action masses are a common Higgs factor times the singular
values of `Y`, so (5) applies on the positive chamber.

This handoff is a configured quadratic action. The execution above proves
the logical Nambu/family generator `X tensor Y`; it does not compile a
weak-kick history on the full ninety-mode physical CAR algebra. In
particular, `Gamma(sum a_r R_r)` is not `sum a_r Gamma(R_r)`. Applying
controlled second-quantized rotations to every occupation sector would
therefore require a separate interaction and limit calculation. Neither
this construction nor the arbitrary-Yukawa intertwiner supplies that
missing dynamical equality, an interacting pole-mass map, a running-to-pole
correction or a selected physical vacuum.

The configured coherent alternative realizes finite Koide balance without
replacing the four failed fixed-action tests in
[`KOIDE_ATTACHMENT.md`](KOIDE_ATTACHMENT.md). It supplies an explicit
source and clock for its logical generator. The phase and the physical
interpretation retain their stated premises.
