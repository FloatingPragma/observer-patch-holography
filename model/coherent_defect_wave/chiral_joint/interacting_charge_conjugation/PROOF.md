# Antiunitary charge and reflection symmetry of the native source

The source is the complete periodic Peter-Weyl Hamiltonian on
`F(C^(30 N^3)) tensor C^5^(tensor N^3) tensor C^119^(tensor 3N^3)`, with
an even side N, the actual antiperiodic kinetic coefficients, real native
Yukawa coefficients, and the scalar and link excitation penalties. No
reference system is included in this Hamiltonian. Its clock is
`exp(-it H_lambda)`. Scaling the kinetic term by a real lambda preserves
the identities below. The ground and charged-response consequences use
exactly the frozen N=2, lambda=2^-34 member.

## The antiunitary on the complete carrier

Write x'=-x modulo N and epsilon_x=(-1)^(sum_j floor(-x_j/N)). Define
Theta by its action on the native CAR operators,

`Theta c_(x,s,a) Theta^-1 = i epsilon_x c_(x',s,a)^dagger`.

On the scalar register it exchanges states 1 with 3 and 2 with 4,
fixing state 0, and moves the register from x to x'. A link at (x,j)
moves to (-x-e_j,j). Its complete label transforms as
`(sign,a,b) -> (sign,b,a)`, with the blank state fixed. Coefficients are
complex conjugated. The executable uses the full exterior lift, including
all occupied-mode permutation signs and antiperiodic seam phases.

Factor Theta=P A. The local antiunitary A sends c to i c^dagger and
uses the scalar exchange. Its link action is the antiunitary C K, where
`C|sign,a,b>=|-sign,a,b>`. For every actual compressed coordinate U_ab
and every endpoint current L_e^A, the complete 119-dimensional matrices
obey

`C conjugate(U_ab) C = U_ab^dagger`,
`C conjugate(L_e^A) C = -L_e^A`.

These identities are verified for all 59 native coordinates and all 24
endpoint currents. They cover occupied-to-occupied transitions. Inversion
uses the existing complete-link involution J, with
`J U_ab J = U_ba^dagger` and `J L_0 J=L_1`. Thus J C is precisely the
same-sign transpose permutation used above.

The complete local CAR/scalar polynomial satisfies A H_local A^-1=H_local
for each of the three Yukawa terms separately and for the scalar penalty.
Every local gauge and spin current is negated. For an off-site kinetic
word `k c_x^dagger sigma_j W c_y`, antiunitarity and CAR reordering give
minus its Hermitian adjoint. The spin-index transposition introduces no
extra sign because each Pauli matrix is Hermitian. The path becomes its
adjoint through the complete coordinate identity. Consequently
`A H_kin A^-1=-H_kin`.

Spatial inversion P leaves spin indices unchanged, reverses every path,
and uses J on every reversed edge. The scalar kinetic coefficient is
purely imaginary. The exact seam identity in the native periodic kinetic
formula therefore gives `P H_kin P^-1=-H_kin`; onsite terms and excitation
penalties are preserved. This is the same path calculation as the native
proper-rotation proof with the axial spin minus removed. The executable
checks every directed path word at N=2 and N=4, including a failing
control when the seam phases are omitted. The formula applies to every
even N and every full-carrier basis state.

It follows that

`Theta H_lambda Theta^-1=H_lambda`,
`Theta G_(x,A) Theta^-1=-G_(x',A)`.

The physical Gauss kernel is preserved. Theta squared is source fermion
parity. The explicit CAR formula has this square on every occupation
sector; the scalar and link permutations square to the identity. The
inversion exterior lift commutes with A: its real one-particle determinant
is one because each spatial mode occurs thirty times. Full Hamiltonian
and Gauss column controls additionally exercise filled and partially
filled Fock states, all five scalar values, and nonblank link values.

## Time parameter and physical interpretation

Antiunitarity conjugates i. Although Theta preserves H_lambda, its action
on the source clock is

`Theta exp(-it H_lambda) Theta^-1=exp(+it H_lambda)=U(-t)`.

The symmetry therefore reverses the evolution parameter. The exact claim
is a finite antiunitary charge/reflection symmetry and its equal spectral
measures. It is not a proof of a unitary physical CP symmetry, physical
CPT, Lorentz invariance or spin-statistics. In particular, this
antiunitary identity does not establish absence of measured weak CP
violation. Complex family couplings, CKM dynamics and experimental weak
CP observables are separate constructions. The local polynomial proof
here explicitly uses the declared real native Yukawa coefficients.

## Actual source addition and removal spectra

For the frozen N=2 member the physical ground g is unique and even.
The inherited ground theorem and the full symmetry imply that its
projector is Theta invariant. At the origin epsilon_0=1 and x'=0.
For every real spectral interval I and every native mode,

`||1_I(H_lambda-E_g)c_(0,s,a)g||^2`
` = ||1_I(H_lambda-E_g)c_(0,s,a)^dagger g||^2`.

The equality includes every actual eigenvalue and the complete spectral
complement. It is an equality of source measures in opposite boundary
charge sectors. It also gives exact origin-mode occupation 1/2, hence
exact pump-failure probability 1/2 for each normalized native family
pump on the true interacting ground with its reference initially empty.
No local atomic energy is substituted for a full interacting energy.

For e^c removal, the source Gauss charge at the origin is -Y_ec=-1;
addition has charge +1. Theta pairs these sectors at identical source
energies. Their centered raw lepton numbers are opposite, since
`Theta L Theta^-1=2 N^3-L` and the neutral ground has raw L=N^3.
Scalar and link contributions remain in the exact paired eigenvectors.
The frozen positive e^c bin therefore has an exactly paired addition
measure. This proves equal source energies without a mass calibration.
It does not identify either finite excitation with an empirical particle.

## The actual finite rotation type

The unique full interacting ground has a nonzero overlap with the native
cube trial vector. On the three-dimensional polynomial subspace for the
roots of `x^3+20x^2-16x-576`, the code applies the actual quarter turn about
z to the twisted cube vectors. The calculation includes the permutation
sign of eight odd local ground doublets, their antiperiodic seam phases,
and their inherited spin-half matrices. The action is +I. All three
roots therefore have quarter-turn character +1.

A one-dimensional representation of the proper cube group is determined
by this character: quarter turns about coordinate axes are conjugate and
generate the group. The spin center acts +1 on the even ground. Thus the
trial and, by positive overlap, the unique full interacting ground have
trivial proper-cubic character. At the fixed origin the two native e^c
annihilators transform in the irreducible binary-cube type G1. Every
nonzero energy-projected pair therefore has scalar Gram matrix
`<g,c_s^dagger P_E c_t g>=w_E delta_st`. Its normalized two-dimensional
range is an invariant odd G1 eigenspace, even when the full energy
eigenspace is larger.

This is a finite-group statement. G1 also occurs in the restriction of
continuous spin 7/2: the proper-cube angle classes have counts
`(1,6,9,8)`, G1 characters `(2,sqrt(2),0,1)` and spin-7/2 characters
`(8,0,0,1)`, giving inner product one. No continuous physical spin or
Lorentz little-group action on the full interacting charged sector has
been proved. Assigning an abstract SU(2) action within a degenerate pair
would be an extra symmetry choice and would not establish spatial spin.

## Reference and boundary energies

The frozen charged-response apparatus has a full thirty-mode static CAR
reference with `H_ref=4 N_ref`. Its standard particle-hole map sends
`N_ref -> 30-N_ref`, hence `H_ref -> 120-H_ref`. The vacuum and the filled
reference have energies 0 and 120. This reference Hamiltonian is not
charge-conjugation invariant. Therefore the source symmetry does not
by itself prove equality of total source-plus-reference energies or
of complete instruments with that unchanged reference preparation.

A source-only removal sector is compensated by a boundary charge +1;
its conjugate addition sector needs boundary charge -1. A boundary
extension may give these two retained states equal declared reference
energies, in which case the two corresponding total sector spectra are
equal. Both new states, their grading and their energies must be declared.
Keeping the old full reference complement does not make it globally
charge-conjugation invariant. A full reference symmetry would additionally
need an energy-preserving dual reference carrier, for example a direct
sum with the conjugate full Fock representation and an antiunitary
exchange between the summands. No such enlarged reference is substituted
into the frozen charged-response package here.

The public receipt freshly replays these identities, compares every inherited
scientific field with the immutable research receipt, retains all three
frozen prerequisite receipts and validates
their inherited source pins. The current result adds an exact symmetry
and spectral pairing; it does not change the source, add transport,
provide a continuum particle dictionary, or promote an empirical gate.
