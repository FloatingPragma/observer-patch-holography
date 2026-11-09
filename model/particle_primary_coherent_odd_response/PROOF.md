# Populated fixed spatial analyzers on the original selected source

## 1. One physical source and its actual six odd fields

Fix member zero of the [primary contract](../particle_common_contract/PROOF.md):
spacing a=1, eight cells x in {0,1}³, twelve internal nearest-neighbor links,
the interacting Wilson Hamiltonian H, product link-word cutoff, actual
Gauss space, selected stationary polynomial density rho, idle quantum
reference, and the same complete source-energy clock. The state satisfies

    Tr rho(H-E0) <= epsilon = 1/16.                         (1)

No state preparation, source action, reference or time scale is changed.
The six actual local gauge-singlet degree-three creation polynomials
F_(x,j), j=1,...,6, are those of the
[selected spin response](../particle_primary_selected_spin_response/PROOF.md).
Two transform as the G1 probe multiplet and four as H. Their native
occupation increment is Delta=(0,2,0,0,1), in native order
(-4,-3,1,2,6). They have zero total gauge charge. Their local full-CAR
calculation includes all occupations of six active modes, with every
other native mode retained as a spectator.

The proved local effect bound and selected occupation bound give

    sum_j F_(x,j)* F_(x,j) >= 6I - 3N_x,
    w = 1 - <D_L+D_e>/32 > 0.84402.                        (2)

Here sum_x <N_x> = <D_L+D_e>/2. This is the literal source conjugation
identity, not an assumption of spatially independent or empty sites.
The finite cube representation labels describe the fields/probes;
they are not final-state total spin, Lorentz spin or a particle name.

## 2. A complete signed apparatus for every fixed spatial phase

Let f_x have |f_x|²=1/8, with arbitrary complex phases, and define

    A_(f,j) = sum_x f_x F_(x,j).                            (3)

At distinct sites the odd local polynomials anticommute, including
with each other's adjoints. Thus

    {A_(f,j), A_(f,j)*}
        = sum_x |f_x|² {F_(x,j), F_(x,j)*}.                 (4)

Each normalized F has F²=0 and norm at most one, as checked on the
full64 local CAR states. A square-zero contraction obeys
{F,F*}<=I: F*F and FF* have orthogonal supports and each is bounded by I.
The receipt also calculates the exact anticommutator spectra for every
actual G1/H component. Hence the thirteen Kraus operators

    K_(j,+)=A_(f,j)/sqrt(6),
    K_(j,-)=A_(f,j)*/sqrt(6),
    K_fail=[I-(1/6)sum_j{A_(f,j),A_(f,j)*}]^(1/2)           (5)

form a complete instrument. The failure is a positive square root,
not a projector. All signs, components and failures remain as outputs.
Equations (2) and (4), and flat |f_x|, give the original total success
probability

    p_++p_- >= w.                                         (6)

This holds for every chosen spatial phase without a success normalization
or an average over phase settings. Interference terms cancel only in the
sum of both signs. A creation-only lower bound is obtained below using
energy and the selected state, not by assuming its spatial correlations.

An apparatus memory has loaded and vacant physical singlet slots and
supplies -Delta for creation and +Delta for annihilation; source and memory together preserve native occupations
and CAR parity. All source, memory, reference and failure correlations
are retained. The finite isometry in (5) admits a unitary completion in
physical source-plus-memory blocks. This is a declared collective
apparatus construction, not a claim that autonomous repair implements
it efficiently. The source input remains rho before intervention.

## 3. The actual Wilson residual and cutoff locality

The full source commutator is

    [H,F_(x,j)] = 9 F_(x,j) + R_(x,j).                      (7)

Onsite terms give9, and electric and magnetic terms commute with the
local gauge-singlet F. The residual comes only from the three incident
Wilson hops. Its support comprises matter at x and its three neighbors,
and the three incident links. It is odd and gauge invariant.

Locality is proved before physical restriction. Let P_link be the tensor
product of the actual link-word projectors. It preserves the support of
each incident-link hopping residual. Because F preserves the cutoff,

    [P_link H_parent P_link,F]-9F
        = P_link [H_parent,F]_residual P_link.              (8)

Identity factors on all other links can be omitted when this is viewed
as an operator on the retained tensor-product carrier. Therefore each
compressed residual has its original star support. The Gauss
projector commutes with each F and residual because they are gauge
invariant. Restricting these operators and their identities to the
physical space does not change their products or graded commutation.
A generic global compression would not justify this step; it is the
actual product cutoff followed by a reducing Gauss restriction that does.

For each site, the six Wilson hopping images are orthogonal to the six
original modes and to one another, with Gram z_x/2 times I6, where z_x=3.
The full12-mode CAR calculation in the selected spin proof gives

    sum_(j in G1) R_(x,j)*R_(x,j) <= (19/3)(z_x/2) I,
    sum_(j in H)  R_(x,j)*R_(x,j) <= 12(z_x/2) I.          (9)

No hopping-image occupation is assumed empty. Particle-hole conjugation
on all twelve CAR modes sends the actual real cubic residual polynomials
to their adjoints, up to a common sign. The code verifies the exact
integer-matrix conjugation of the summed Grams. Thus (9) also holds
with RR* in place of R*R. Compression preserves these upper bounds by
(P R P)*(P R P)<=P R*R P. Combining both multiplets and signs gives

    sum_j {R_(x,j),R_(x,j)*} <= 55 I.                     (10)

## 4. Four disjoint-star colors give a spatial moment bound

Color the eight cube vertices into the four opposite-corner pairs
{x,1-x}, with x_1=0. In every pair, the two residual stars have disjoint
matter supports and disjoint incident links. Equation (8) therefore
preserves their graded anticommutators exactly.

For X_j=sum_x f_x R_(x,j), let X_(c,j) be the sum in one color. Then

    {X_(c,j),X_(c,j)*}
        = sum_(x in c)|f_x|² {R_(x,j),R_(x,j)*}.

The operator Cauchy inequality for four colors and its adjoint imply

    {X_j,X_j*} <= 4 sum_x |f_x|² {R_(x,j),R_(x,j)*}.        (11)

For (5), creation is centered at transfer+9 and annihilation at-9,
since [H,A*]+9A*=-X*. The exact stationary spectral identity for the
complete two-energy experiment gives the summed signed second moment

    M = sum_j integral (t-9)² dmu_(j,+)(t)
              + integral (t+9)² dmu_(j,-)(t)
      = (1/6) sum_j Tr rho {X_j,X_j*}
      <= 4*55/6 = 110/3.                                 (12)

All measures are original joint probabilities. The same actual source
Hamiltonian and clock enter both energy records. The reference is kept;
only source stationarity is needed for this spectral identity.

## 5. An original positive-energy event for every fixed phase

Set r=17/2. Markov's inequality applied to (12) bounds the union of
out-of-window signed records by

    M/r² <= 440/867.

Thus the creation window J_+=[1/2,35/2] together with the annihilation
window J_-=[-35/2,-1/2] has original probability at least w-440/867.
For any complete instrument, a transfer t<=-1/2 requires the initial
energy to be at least E0+1/2. Summing its entire Kraus family retains an
effect at most I on every initial energy state. The finite ground-tail
[theorem](../../Lean/QFT/GroundRecordTail.lean)
`OPH.QFT.GroundRecordTail.downwardWeight_le_excess_div` and (1) therefore
bound all these downward transfers by

    Pr{t<=-1/2} <= epsilon/(1/2) = 1/8.                    (13)

Its arrays are the actual eigenvalues E_i, the diagonal state weights p_i,
and w_(i,k)=p_i sum_b |<k|K_b|i>|² for the retained success rows b.
The spectral lower bound supplies hfloor, squared amplitudes supply hw,
and completeness supplies hrow; (1) supplies hexcess. Degenerate energies
can be resolved in any joint eigenbasis of H and its polynomial rho.
There is no factor for the twelve success operators. Removing that whole
negative event proves the creation event on the original state:

    mu_+(J_+) >= w - 440/867 - 1/8 > 0.2115 > 1/5.        (14)

This is a positive response at each specified phase, not a phase chosen
after a successful record. The same complete two-clock plan as the spin
proof has per-clock resolution1/8 and total event error at most1/4096.
It reports the widened interval [1/4,71/4] with weight

    > w - 440/867 - 1/8 - 1/4096 > 1/5.                  (15)

The receipt binds the actual clock plan and exact rational lower bound.
No division by success probability, original occupation probability or
source-reference conditional weight appears.

The multiplets can be separated without paying the whole downward bound
for each. Their signed subinstrument effects obey E_lambda<=d_lambda I/6
by (4). The same finite theorem has row bound (d_lambda/6)p_i
and excess bound (d_lambda/6)epsilon. It therefore gives downward weight
at most d_lambda/48. Their moment bounds are M_G1=38/3 and M_H=24.
Consequently, at every selected flat phase,

    mu_G1(J_+) >= w/3 - 152/867 - 1/24,
    mu_H (J_+) >= 2w/3 - 96/289 - 1/12.                  (16)

After the same complete-clock event allowance1/4096, these give reported
G1 weight greater than1/16 and H weight greater than9/64. They sum over
components inside each multiplet, and do not assume the spatial phase
setting is invariant under source rotations.

## 6. One complete apparatus with all eight spatial rows

A separately selected normalized phase instrument (5) must not be
confused with eight disjoint rows of one experiment. For the latter,
prepare a coherent uniform address over the eight cells, apply each
complete local signed six-component apparatus, and retain its eight
address failures. Before measuring successful addresses, rotate each
success register with the real Walsh matrix

    W_(q,x)=(-1)^(q dot x)/sqrt(8),   q,x in {0,1}³.

The resulting96 success Kraus operators are

    K_(q,j,+)=A_(f_q,j)/sqrt(48),
    K_(q,j,-)=A_(f_q,j)*/sqrt(48),                         (17)

with eight failures K_(x,fail)=[I-sum_j{F_(x,j),F_(x,j)*}/6]^(1/2)/sqrt(8).
The exact unitary identity W*W=I proves completeness. Every spatial row,
sign, internal component, address failure and clock record remains.
The apparatus has104 records before the source-energy clocks and optional
native-number meters. The address is not measured before this rotation.

For a given q, its implemented creation-event effect is precisely1/8 of
that of (5) with f=f_q, because the source clock is the same and every
success Kraus operator is scaled by1/sqrt(8). Therefore each joint
positive spatial row retains weight greater than1/40. Even the coarser
bound that subtracts the entire global clock error after division by8
exceeds1/40 for this source. Eight such disjoint lower bounds sum to1/5,
not8/5. The same exact1/8 effect scaling gives, at every positive spatial
row, G1 weight greater than1/128 and H weight greater than9/512. For
complex Fourier rows, use the conjugate pointer rotation on
the negative branch to obtain A_f*; the same complex rotation would give
A_(conjugate f)*. Walsh rows are real.

## 7. Quantum labels and scope

Each creation/annihilation operator has native increment +/-Delta.
The [five coherent radix-three meters](../particle_primary_selected_spin_response/PROOF.md#6-conserved-native-current-readout-without-initial-charge-measurement)
therefore attach that signed increment without first measuring unknown
absolute occupations. The complete failure commutes with native numbers.
Literal source charge conjugation transports the full apparatus to the
corresponding negative-beta hole apparatus with identical complete
classical source, probe and clock record probabilities. Equality of the
retained source-reference output states is not asserted.
This is not an assertion that A* itself is the charge-conjugated field.

G1/H remain actual internal probe multiplets. Both are populated at every
spatial row, with their distinct stated lower bounds. These bounds sum
over their components and do not prove a positive event for each individual
internal component.
The fixed Walsh/Fourier label is a finite collective spatial analyzer,
not an eigenvalue of a conserved translation generator. The construction
is certified at member zero. It proves neither a fixed-energy continuum
residue nor relativistic dispersion, an isolated pole, a rest frame,
a physical species assignment or a GeV mass. It gives an actual populated
joint spatial, native-current and source-energy record without changing
the source state to a trial vector.


## 8. Exact four-measure reduction of the spatial response

The actual proper-cube symmetry reduces each internally summed G1 or H
energy-transfer Gram to four real measures, possibly signed off diagonal,
one for each Hamming distance between source cells. The complete Walsh rows are four exact
Krawtchouk combinations with multiplicities1,3,3,1. The
[finite reduction proof](KERNEL_PROOF.md) identifies these measures on the
same source and apparatus and supplies the tested inverse transformation.

This reduces the spatial-response calculation to four measures per
multiplet; it does not compute their values. The complex time transforms,
full component kernels and cross-multiplet kernels are not reduced by
this scalar argument. A difference of unnormalized spatial responses
need not imply a difference of their normalized energy shapes or a
dispersion relation. The original positive event bounds above are unchanged.
