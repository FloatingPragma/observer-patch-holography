# The complete physical local ground and gap

## Source and domain

The carrier is the native thirty-mode fermion Fock space tensor the
five-state scalar, of dimension `5*2^30`. The Hamiltonian is the sum of
all three native quantum Yukawa channels and the scalar excitation
projector, with all four coefficients one. The domain is the simultaneous
kernel of all twelve local gauge generators. There is no boundary charge
or link compensating a nonzero local gauge representation. All raw baryon
and lepton sectors, fermion parities and spin multiplets in that physical
domain are compared.

The source compiler is checked against the full normal-ordered quantum
Yukawa polynomial, rather than only selected occupation columns. The
full gauge and spin Ward identities are inherited from the pinned
[native composite proof](../higgs_composites/PROOF.md) and recomputed as
operator polynomial identities in `reference.py`. The local clock is
`exp(-itH)` with hbar one. No intersite or composite hopping is included.

## Exact particle-hole transformation

Exchange creation and annihilation on the fourteen native weak-singlet
modes, whose internal indices are `0,1,4,9,12,13,14` in each spin copy.
This is a finite canonical anticommutation relation automorphism and
therefore has a unitary Fock implementation. Exact normal ordering gives

`H' = N_H + sum_(r=1)^4 (dGamma(M_r) tensor |r><0| + adjoint)`.

All four rational thirty-dimensional matrices are extracted from the
actual native tensors; rebuilding the entire operator polynomial from
them gives exact equality. They obey `M_3=M_1^transpose` and
`M_4=M_2^transpose`. The one-body graph consists of six quark components
of dimension four and two lepton components of dimension three. Exact
commutator closure has dimension eighteen, with commuting quark and
lepton parts of dimensions ten and eight.

On every quark component the first two transition matrices and their
adjoints obey the two-generator one-body Clifford relations. Write

`gamma_(2i)=(F_i+F_i^dagger)/2`,
`gamma_(2i+1)=(F_i-F_i^dagger)/(2i)`.

The four Hermitian gamma matrices satisfy
`{gamma_a,gamma_b}=delta_ab I/2`. Together with the six Hermitian
`T_ab=i[gamma_a,gamma_b]`, they form ten trace-orthonormal generators
of the compact algebra sp(4), equivalently so(5). The code verifies
the root spaces, fundamental weights and the universal tensor identity

`sum_T ([T,M_r] tensor T + T tensor [T,M_r])=0`.

Consequently the collective quadratic Casimir commutes with every
transition on all exterior powers. Its fundamental spinor value is
5/2. These are one-body Clifford identities: collective second
quantizations of the F matrices are not canonical fermion annihilators.

Let `A_q` be the column of the four collective quark transitions from
the scalar blank state. Pairing each F with its adjoint gives

`A_q^dagger A_q = 2 sum_a dGamma(gamma_a)^2 <= 2 C_q`.

The omitted six squared Hermitian rotation generators are positive.
For the leptons, the full sixty-four-dimensional Fock matrix is built
by exact canonical signs. Its squared-transfer spectra at transformed
occupations zero through six have maxima `0,2,5,6,5,2,0`; thus
`A_l^dagger A_l <= 6 I` on the complete lepton space.

## Complete quark representation inventory and the complement bound

In doubled orthogonal weight coordinates the fundamental four weights
are `(1,1),(1,-1),(-1,1),(-1,-1)`, and the positive roots are
`(2,-2),(0,2),(2,0),(2,2)`. The full occupation character is the exact
product of the twenty-four factors `1+t e^weight`. Multiplication by
the Weyl denominator extracts the dominant highest-weight coefficients.
The code reconstructs the entire alternating polynomial at every weight,
including all nondominant weights, and checks the dimension at occupation
n against `binomial(24,n)`. The total is `2^24`.

For doubled highest weight `(a,b)`, `a>=b>=0`, the dimension and Casimir are

`dim = (a+3)(b+1)(a-b+2)(a+b+4)/24`,
`C = (a(a+6)+b(b+2))/4`.

These are the Weyl character and quadratic Casimir formulas in the
normalization verified above; see the primary derivation in
[Etingof, lecture 26](https://ocw.mit.edu/courses/18-745-lie-groups-and-lie-algebras-i-fall-2020/mit18_745_f20_lec26.pdf).
The finite polynomial reconstruction and dimension checks are part of
the replay, so no unexamined representation table supplies the bound.

The inventory contains twenty-eight distinct collective irreducibles.
The unique maximal Casimir is 54, at transformed quark number twelve,
highest weight `(12,0)`, dimension 140 and multiplicity one in the
entire quark Fock space. The transformed charges are checked as full
operator identities:

`N_q'=12+B3`, `N_l'=2+L`.

The two transformed color Cartans assign weights `(1,0),(-1,1),(0,-1)`,
each to eight quark modes. A color singlet has equal occupations in the
three colors; hence its transformed quark number is divisible by three.
Among all such occupations, every module other than the unique C=54
module has C at most 42. This uses only a necessary color-singlet
condition, so it includes an upper-bound superset of every excluded
physical sector.

For those excluded sectors the complete scalar transfer is `A=A_q+A_l`.
The operator triangle inequality gives

`||A||^2 <= (sqrt(84)+sqrt(6))^2 = 90+12 sqrt(14) < 135`.

The strict last inequality is checked by rational squaring. Since
`H=[[0,A^dagger],[A,I]]`, every energy is bounded below by
`(1-sqrt(1+4||A||^2))/2`. Exact arithmetic gives
`(-1113/100)*(-1213/100)>135`, and therefore

`H > -1113/100 = -11.13`

on the entire excluded physical complement. This bound retains all
lepton occupations, all scalar states, all spins and both parities.
Without the color-singlet restriction the Casimir 95/2 is admitted and
this exclusion at 135 does not follow.

## An equivalent maximal module with all physical states

All eight native color generators and the three physical spin generators
commute with the collective quark algebra. Their transformed quark
constants are zero. The maximal irreducible appears with multiplicity
one in the full quark Fock space. Each commuting connected semisimple
color or spin group therefore acts through a one-dimensional
representation on its multiplicity space, and that action is trivial.
The complete 140-dimensional quark module is a color and physical-spin
singlet.

For a compact exact realization, take the C=4 five-dimensional part of
the exterior square of the four-dimensional spinor, then its sixth
symmetric power. This space has dimension 210 with positive metric
`prod_i n_i! g_i^n_i` in its monomial basis. Its C=54 eigenspace has
dimension 140 and highest weight `(12,0)`, so it is unitary-equivalent
to the unique actual Fock module. Every transition and current is
restricted through an exact invariant intertwiner. Thus this is an
equivalent representation of the actual operator, rather than a
variational truncation.

The actual transformed weak and hypercharge matrices are expressed in
the same ten generators on all twenty-four native quark modes. The
native quark hypercharge contains the scalar term `N_q'/6-2`; it
vanishes exactly at `N_q'=12`. This fixes the hypercharge origin of the
abstract realization. The full sixty-four lepton states and five scalar
states retain their actual native current matrices.

At each lepton number and highest physical spin, the code takes the
complete weight-zero basis for weak T3 and total hypercharge, imposes
the actual weak raising current into its entire adjacent weight space,
and imposes the physical-spin raising current. Exact kernel dimension
agrees with the independent weak/spin character subtraction. The
Hamiltonian preserves each kernel exactly. Its restricted Gram matrix
is the pullback of a positive metric, and `G H=H^dagger G` holds.
Including all magnetic partners gives exactly 104 physical states.

The scalar blank and excited parts of each physical block are retained
separately. Their H diagonals are zero and identity; the full squared
transfer is nonsingular on each nonempty blank part. Let

`p_v(x)=x^3-313x^2+31738x-1041088`,
`p_o(x)=x^3-313x^2+31708x-1041468`.

The complete highest-spin characteristic polynomials are:

| Raw L | Spin | Dimension | Characteristic polynomial |
|---|---|---:|---|
| -1, 3 | 1/2 | 2 each | `(E-1)^2` |
| 0, 2 | 0 | 11 each | `(E-1)^5 p_v(E(E-1))` |
| 0, 2 | 1 | 5 each | `(E-1)^3 (E(E-1)-103)` |
| 1 | 1/2 | 16 | `(E-1)^8 (E(E-1)-105) p_o(E(E-1))` |
| 1 | 3/2 | 3 | `(E-1) (E(E-1)-108)` |

No other lepton occupation or spin contributes a physical state in this
maximal module. Magnetic multiplicities are `2s+1`. The total dark
energy-one rank is 56. The remaining forty-eight states give eighteen
distinct nontrivial energy levels.

## Global physical ground, first gap and reference

Exact rational root isolation and outward square-root intervals order
all nineteen maximal-module energy levels strictly. The lowest is

`E_g=(1-sqrt(1+4 xi_o))/2`,

where xi_o is the largest root of p_o. It has rank two, raw L=1,
B3=0 and spin one-half. The next energy uses the largest root of p_v
and has rank two, comprising one spin singlet in each of L=0 and L=2.
Both complete levels lie below -11.13. The excluded-sector bound
therefore proves that these are the ground and first excited levels
of the entire local physical source, including states invisible to
the vacuum and neutral cyclic preparations.

The exact ground interval is

`[-1580564367072699/140737488355328, -3161128734145397/281474976710656]`.

The full physical gap is enclosed by

`[25207584560067/281474976710656, 25207584560069/281474976710656]`.

The chosen reference is `rho=P_g/2`. It is stationary, gauge-invariant
and invariant under the actual local spin rotations. This states a
reference-selection rule; it does not construct a preparation from an
arbitrary source state. The reference has odd fermion parity because
the transformed total occupation is fifteen and the particle-hole
map exchanges fourteen modes.

Independently of the fourteen-mode transformation, the full antiunitary
`A c_i A^(-1)=i c_i^dagger`, with scalar conjugation, preserves H and
sends gauge and spin generators to their negatives. For every native
irrep of internal dimension d it sends its occupation `N_r` to
`2d-N_r`. These are exact full polynomial identities. The ground
projector is A invariant, so its mean occupation in each native irrep
is d: `(N_uc,N_Q,N_ec,N_L,N_dc)=(3,6,1,2,3)`, with total fifteen.
Each occupation commutes with physical spin; Schur's lemma on the
irreducible ground doublet gives the same means in every pure ground
spin state. This is a filled interacting reference. Raw L is a source
charge, and its centered value L-1 is zero on the doublet.

## Scope and replay

The theorem covers the full zero-boundary-charge one-cell physical
space. It does not give the unconstrained or boundary-compensated
local ground, a complete spectrum of the excluded quark modules, a
relativistic particle dictionary, or the vacuum of an interacting
many-cell Hamiltonian. Nonblank gauge links can compensate charged
local sectors; their energies require bounds on those sectors.

The receipt rebuilds every exact source matrix, full quark character,
physical kernel and spectral interval. It hashes all declared inputs
and their repository-local import closure and fails on missing files.
The checks reject omitted native Yukawa words, altered fundamental
characters, a missing maximal module, omitted current images, false
scalar penalties, incomplete spectral ranks and overlapping level
claims. No floating eigensolver or persisted numerical eigenbasis
enters this certificate.
