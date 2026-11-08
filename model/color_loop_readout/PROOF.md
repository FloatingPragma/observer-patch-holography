# Electric and magnetic coefficients from complete color-loop records

The construction uses the complete 438- and 1,596-dimensional Gauss spaces
of the [finite Higgs join](../interacting_extension/higgs_join/PROOF.md).
Every Higgs, lepton and native-link state in those declared sectors is
retained. Six finite-time, unconditioned record distributions determine
the positive bare electric coefficient relative to the magnitude of the
magnetic coefficient, with exact rational error bounds. This is a finite
coefficient measurement on that specified interacting Hamiltonian.

## Source, input and measuring clock

There are three family copies, two native 60-state links, one scalar doublet
at vertex 0 with occupation at most one, and the charged-lepton Higgs vertex.
Colored fermion occupation is zero. The up and down pair vertices are
disabled. All three family matrices are the identity, all phases vanish,
and the fermion and external-reference energy offsets vanish. These are
literal defaults of `build_join`, rather than parameters fitted to a
measured strong coupling. The independently variable rational controls are
the electric strength e, magnetic strength m, hopping strength h, Higgs
strength g and scalar energy s. Their default values are all one.

The source Hamiltonian and its one clock are

    H = e E - m (L+L*) + h T + g S + s N_H.                    (1)

Here T is the source hopping operator including its minus sign, and
S is the sum of its three diagonal-family Higgs creation operators and
their adjoints. The exact operator numerator convention is

    <v_i,H v_j> = A_ij,       ||v_i||² = n_i,
    H_ij in raw coordinates = A_ij/n_i.                       (2)

All inner products below use this Gram metric. Normalized matrix entries
are A_ij/sqrt(n_i n_j). `source_replay` rebuilds every basis vector and
every loop, diagonal-family hopping and Higgs numerator used by (1) from
the actual native source. The four-dimensional flux compression is not
used to evolve the state.

The Q*=0 input chi_0 is the empty matter/scalar/link vector. The Q*=1
input chi_1 is exactly the image of chi_0 under the source's
reference-neutral addition operator for family 0 and spin up. The odd
external reference is included in that source sector. Both vectors have
unit norm. They are a specified preparation pair, rather than interacting
ground states. Every time setting evolves this same input under (1).

The detector records (N_H,r_0,r_1), the scalar occupation and the irreducible
blocks of both links, with -1 denoting an empty link. Its orthogonal
projectors P_z exhaust the full sector. They commute with every Gauss
generator: the native generators preserve each link representation and
scalar occupation. The complete instrument is

    K_z(t) = P_z exp(-itH),      sum_z K_z(t)* K_z(t) = I.       (3)

Its definition includes the output state and any untouched quantum
reference, not only a list of probabilities. The source declares finite
Hamiltonian control, the input and this gauge-invariant measuring
instrument. This package supplies no compiler of these flux measurements
from bare occupation-repair records and no identification of this clock
with physical time. Those interfaces are explicit premises of the finite
experiment; a numerical probability calculation does not discharge them.

## The cancellation on the full interacting spaces

The actual three colored link blocks have

| r | six times hypercharge | dimension d_r | one-link native Casimir c_r |
|---|---:|---:|---:|
| 0 | -4 | 3 | 4 |
| 2 | 1 | 6 | 25/8 |
| 3 | 2 | 3 | 3 |

The trace convention of the source is

    C_native = 2 C_color + C_weak/2 + (6Y)²/12.                (4)

For each q=0,1 and r, let ell_r be the normalized loop vector with the
same matter/reference input as chi_q. The two oriented links are contracted
as (a,b) and (b,a). The source operator gives

    L chi_q = sum_r d_r ell_r,    sum_r d_r² = 54,
    E ell_r = 2 c_r ell_r.                                  (5)

The factor 2 in the energy is the number of occupied links. It is not a
conversion between e/|m| and a renormalized gauge coupling.

Write P_r for the coarse sum of all P_z with r_0=r_1=r, retaining both
Higgs occupations. The source identities needed below are

    [S,L+L*] = 0,       [S,P_r] = 0,       T P_r = 0,
    ||S chi_q||² = 12-2q,       ||T chi_q||² = q.              (6)

The first identity follows because the Higgs operator changes only scalar
and matter factors, while the loop changes only link factors. The hopping
operator vanishes on a colored link pair because the sector contains only
leptonic matter; neither endpoint index can match a colored annihilation,
and link creation requires an empty link. A single hopping from chi_1
creates one noncolored link, so both L and L* annihilate that image.
These statements can also be checked in the full exact matrices.

Let u_r=S ell_r and a_q=12-2q. Directly applying (1) gives

    P_r H chi_q = -m d_r ell_r,
    P_r H² chi_q = -m d_r (2 e c_r ell_r + 2g u_r),
    <ell_r,H³ chi_q> = -m d_r [4e²c_r²+54m²+3g²a_q+q h²].   (7)

The states ell_r and u_r have different scalar occupation, and
||u_r||²=a_q. The scalar energy s cannot contribute to the matrix element
in the third line: a single Higgs vertex changes occupation, while two
vertices plus a scalar-energy insertion leave no step for loop creation.
Two hoppings contribute only through their return to the seed. These
identities are checked for both complete sectors and several nonunit
rational controls, with no restriction to the flux compression.

All matrices, input coordinates and projectors are real, so p_r(t) is
even. Expanding the complete record probability gives

    p_r(t) = ||P_r exp(-itH) chi_q||²
           = A_r t² + B_r t⁴ + O(t⁶),
    A_r = ||P_r H chi_q||²,
    B_r = ||P_r H² chi_q||²/4
          - <P_r H chi_q,P_r H³ chi_q>/3.

Inserting (7) proves

    A_r/d_r² = m²,
    B_r/d_r² = -m²[e²c_r²/3 + 18m² + q h²/3].               (8)

The positive g²a_q term in the norm in the second line cancels its
contribution to the third-order amplitude. Discarding the Higgs outcomes
would remove that positive probability and yield a different coefficient.
Thus (8) relies on complete readout of the interacting system.

The interaction does change the color response. At the default controls,
turning off g changes the exact sixth-order coefficients. For the two
dimension-three blocks, B_6(g=1)-B_6(g=0) is -81/5 in Q*=0 and -66/5 in
Q*=1; the dimension-six block has four times these differences. Setting
h=0 removes this sixth-order Higgs dependence. These are exact full-source
operator moments, not inferred continuum masses.

Taking the difference of r=0 and r=3 eliminates both magnetic and hopping
terms in (8):

    M = A_0/9 = m²,
    D = B_0/9 - B_3/9 = -(7/3)e²m²,
    e² = -3D/(7M),       (e/|m|)² = -3D/(7M²).              (9)

This inversion requires m nonzero. The positive square root reports the
positive electric-action convention e>0. The record probabilities do not
determine the sign of m. The source default has e/|m|=1; the coefficient of
the color Casimir in (4) is 2e, so its ratio to |m| is 2. Neither number is
alpha_s. Multiplying the entire Hamiltonian by a positive clock conversion
multiplies A by its square and B by its fourth power, leaving (9) invariant.

For these complete records, sign blindness also holds beyond the quartic
expansion. The diagonal involution J equals -1 on every colored link pair
and +1 elsewhere. It commutes with E,T,S,N_H and each P_z, anticommutes
with L+L*, and fixes both input seeds. Therefore J H(m) J=H(-m), and every
record probability from these seeds is invariant under m changing sign.
The full-matrix certificate checks these commutation signs on every
nonzero generator entry.

## Finite times and exact errors

A source-computed norm bound is obtained by rounding every absolute
normalized matrix entry outward to a dyadic rational and taking the largest
row sum. Hermiticity implies ||H||<=B. This uses every row in the full
sector, including states with zero probability in the selected preparation.

For T_n(t)=sum_{k=0}^n (-itH)^k/k!, Taylor's integral remainder along real
time and unitarity give

    ||exp(-itH)-T_n(t)|| <= R_n(t)
                         = (B |t|)^(n+1)/(n+1)!.             (10)

The PVM dilation Dv=sum_z |z> P_z v is an isometry. Hence the exact and
polynomial complete instruments differ in dilation norm by at most R_n.
Their output trace-norm difference, for any normalized input and any
untouched reference, is at most 2R_n+R_n². Exact rational weighted norms
of T_n(t)chi_q therefore give rigorous probability intervals for every
record, with no branch division. All intervals retain their joint unit
mass enclosure. Reported interval endpoints are rounded outward again
to 128 binary fractional digits for compact storage.

For K positive settings t_j=j tau, use the exact inverse of the matrix
V_jk=j^(2k), j,k=1,...,K. The first two rows recover the t² and t⁴
coefficients of every even polynomial of degree at most 2K. Since
||ad_H^(2K+2)(P_r)|| <= (2B)^(2K+2), Taylor's real-time integral bound gives

    |p_r(t)-sum_{k=1}^K c_2k t^(2k)|
       <= (2B |t|)^(2K+2)/(2K+2)!.                          (11)

Interval arithmetic applies these weights to complete probability intervals
widened by (11). The t² interval has a strictly positive lower endpoint in
the retained experiment, so interval division and the positive square root
give rigorous finite-time bounds in (9). The executed protocol uses K=6,
tau=1/1024 and n=28. There is no extrapolated or floating-point remainder.

An additional absolute empirical error delta per coarse probability enters
the same interval calculation. With independent repetitions at each setting,
the Hoeffding union bound for all three colored probabilities, two sectors
and six settings is at most 72 exp(-2N delta²). `sample_budget` gives an
integer sufficient N using ln(x)<=ceil(log2(x)). Its large recorded value is
a finite resource bound, not a practical experimental proposal. A clock
setting error epsilon can also be included by adding 2B epsilon per
probability. State-preparation trace distance contributes its corresponding
probability error. Such errors must be supplied rather than set to zero
when the specified ideal interfaces are replaced by physical operations.

## Premise consumption and physical boundary

The result consumes the declared Higgs-join action, its exact Gauss sector,
identity family controls, prepared addition pair, Hamiltonian clock and
complete gauge-invariant flux/scalar instrument. It derives the coefficient
inversion, Higgs cancellation, interacting sixth-order response and all
finite-time numerical errors on that same finite system. It supplies no
bare-repair selection of these controls, compiler of the flux instrument,
growing spatial family, chiral colored matter, Wilson-link continuum,
renormalization scheme or correspondence to a laboratory strong coupling.
No claim of an interacting 3+1-dimensional continuum or particle spectrum
follows from this finite readout. The measurement isolates a coefficient
that a source-to-field construction can carry without replacing it by a
dimension count or a measured QCD input.
