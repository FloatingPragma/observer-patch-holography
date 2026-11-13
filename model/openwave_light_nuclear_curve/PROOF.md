# Complete finite nuclear bands

## Common source and calibration

The source is the same native light-quark action used for the common hadron
spectrum and the degree-corrected infinite-chain rest band. On a finite path
of three or four cells, each cell has the complete 12-mode light-quark
exterior algebra. Every occupation and native colour representation is
retained. Antiquark and extra-family vacua are exact reducing spectators
of this number-preserving strong sector. The full-source number floor below
also covers populated sectors of the declared extension; the nuclear
spectral calculation is restricted to the stated pure-light-quark sector.

The interior degree terms are

    −lambda S_x² − eta I_x² − t N_x + delta_U C_x² + delta_m N_x.

Their exact coefficients are imported from the hadron refinement source.
They vanish identically on two cells. The additional interaction is

    B sum_x P[N_x=N_(x+1)=N_(x+2)=3]
      ((S_x+S_(x+1)+S_(x+2))²−7/4),
    B=865442239271/68719476736.

Number projectors commute with every spin in the same expression. This is
one range-three native spin interaction. It vanishes as an operator for
total quark number at most six, so all seven two-cell light calibration
domains are preserved. Its coefficient is chosen using the stated triton
rest energy. The helion and alpha energies are held-out comparisons. The
constant 7/4 specifies the interaction, not a change of the source vacuum
or an independently chosen energy origin for each nucleus.

On every occupied triple window, nine spin-half quarks have total spin at
least one half. The interaction is therefore bounded below by −B per
occupied window. Each such window can be charged to a cell with N_x=3;
there are no more windows than such cells. The count-resolved full-Fock
inequality in `openwave_hadron_refinement_flow.positivity` reserves the
negative nuclear contribution only on these cells. The charge-conjugate
number correction is paid only by populated quark/antiquark sectors.
`joint_stability(B)` proves a strictly positive number floor with both
terms and the common degree correction. Coarse unrelated linear penalties
are not substituted for this occupancy-dependent bound.

## Exact local reduction

For one cell, take the exterior character of

    C³_colour tensor C²_spin tensor C²_isospin.

Multiplying its twelve factors and applying the colour and two spin Weyl
denominators gives 97 irreducible rows
(n,p,q,2S,2I,multiplicity). The producer checks that their dimensions sum
to binomial(12,n) separately for every n and to 4096 in total. No local
representation is omitted by a spin or occupation cutoff. Colour tensor
products are computed from exact weight characters. Each product checks
nonnegative multiplicities and equality of its full representation
dimension.

The path is in the physical Gauss carrier with trivial exterior electric
flux. Sequential colour fusion enumerates all intermediate edge irreps.
The complete singlet dimension at N=3A on A cells follows independently
from exterior Schur--Weyl duality: it is the dimension of the GL(4A)
representation with Young diagram (3 repeated A times). The hook formula
gives 572572 for three cells and 176729280 for four cells. A second
three-cell tensor-character count checks the first value directly.

Let P be the product of the local n=3, colour-singlet, S=I=1/2 quartet.
Its dimensions are 64 and 256. On this reducing local representation the
native spin and isospin currents are Pauli matrices divided by two. Thus
`projected_source` is the literal PHP matrix, not an independently assigned
few-state Hamiltonian. The full hopping creates occupations 2 and 4 on an
edge, so its off-diagonal compression in P is zero. Its original positive
Laplacian contributes the actual diagonal number term. Every spin and
isospin multiplet of PHP is retained.

## Complete complement bound

Write Q=1−P. The exact dynamic program retains total quark number, the
outgoing electric irrep, whether any local state is outside the quartet,
and the last two local occupations. At a cell with irrep (p,q), it adds

    A n²/9 + (U−kappa/2) C2(p,q)
      + Js S(S+1) + Jw I(I+1),

and the displayed interior degree terms. Each interior edge adds
kappa C2(edge). An occupied triple subtracts the rigorous lower bound B.
Positive spin-gradient and isospin bond squares and the original positive
kinetic Laplacian are discarded. Boundary flux, total N=3A, and the Q flag
are imposed exactly. The final minimum, plus the common mass term and
−delta I_z, is an operator lower bound C on the complete QHQ. The minimum
need not be an eigenvalue or an attainable joint state.

A scalar upper bound E on PHP follows from the edge estimates
(S_x−S_y)²<=3 and (I_x+I_y)²<=2 and the triple estimate
(S_x+S_y+S_z)²−7/4<=2 on P. The producer verifies C>E in each box.
Only the hopping connects P to Q. Its complete fixed-number norm obeys

    ||QHP||² <= (9 t sqrt(2))² = 81/32     (three cells),
    ||QHP||² <= (12 times 2t)² = 9        (four cells).

These are bounds on the full second-quantized path adjacency, not on a
selected transition. The finite block spectral theorem then places every
PHP eigenvalue lambda_j in the full-source interval

    [lambda_j−b²/(C−E), lambda_j].

Exactly dim(P) eigenvalues lie below C; all remaining eigenvalues are at
least C. The separated spectral subspace theorem also gives, for arbitrary
normalized input supported in P, high-complement probability at most
b²/(C−E)². These bounds retain the complete Q sector without enumerating
its occupation basis.

## Rational levels and observer records

The finite matrices are reconstructed as rational matrices from exact
quarter-integer spin entries. At fixed S_z=S and the stated I_z, exact
nullspaces isolate each total-spin and total-isospin representation. The
generalized matrices A=P*HP and G=P*P are not normalized by a floating
square root. Deterministic rational inertia bisection certifies all their
eigenvalue intervals. The sum of multiplicities times (2S+1)(2I+1) is
exactly 4^A. Every other S_z has the same energies; another I_z shifts them
exactly by −delta times its change. The complement bound transfers under
the same shift. No higher magnetic states are discarded.

The original input is I/d on the full physical N=3A singlet carrier. Each
charge block has an opposite-charge reference. The P filter has original
success weight 4^A/d and retains its full failure. The spectral apparatus
records every low interval and the full higher complement, conditioned on
the charge block. The receipt bounds both joint probabilities before
conditioning on preparation success. The complete preparation filter is
unital, so its average work on I/d is zero. The complete spectral readout
commutes with H and has zero energy work. All operations preserve charge
and admit arbitrary external quantum references. Inputs outside the
stated number sector receive the identity failure outcome.

## Physical comparison

In electron-calibrated GeV units, the complete three-cell triton-like
S=I=1/2 ground band is enclosed by approximately
[2.8089192030,2.8089230704]. It contains its calibration target
2.80892113668. Every other three-cell spin/isospin ground interval lies
above it. The helion-like partner is approximately
[2.8076258705,2.8076297379] GeV and is disjoint from the held-out target
2.80839161112. The four-cell S=I=0 alpha-like interval is approximately
[3.7391239822,3.7391394985] GeV and is disjoint from 3.7273794118.
Exact rational endpoints, rather than these printed approximations,
are authoritative.

These are complete finite-box source spectra. Different box boundaries
are explicit parts of the operator. No asymptotic breakup threshold,
physical halo, observed binding curve or infinite-volume nuclear mass
is inferred from these boxes. The two helium failures constrain this
particular common action despite the triton calibration.
