# Charged regional spinor endpoints with original energy records

## 1. One physical source and one regional charged factor

Use member zero of the [primary contract](../particle_common_contract/PROOF.md):
the eight-cell open cube, all sixty fermion modes per cell, the native
Wilson action with its product link cutoff, full Gauss support, selected
polynomial state rho, idle canonical reference and source Hamiltonian H.
The source is unchanged. A detector acts on it through a declared complete
instrument. Its work need not vanish.

For adjacent interior sites x and y, take R={x}. The
[induced regional cut](../particle_primary_geometric_cut/PROOF.md)
places the crossing link and every mode at y in the complement. This
is a restriction of the original source, with no added edge carrier.
The regional matrix algebra and its completely positive restriction use
the finite [regional-net theorem](../../Lean/QFT/CPRestrictionNet.lean).
Only gauge-invariant regional matrices inject as physical observables.
A charged regional field is a map between charge sectors; it is not a
physical operator on the full Gauss kernel by itself.

Let p_(x,a), m_(y,b), a,b=1,2, be the actual positive- and negative-beta
spin modes of the one-dimensional native weak/color singlet with Q6=6.
Write chi_xy for its compensating oriented link character and P0_xy for
the Haar-constant projection. The
[charged-string construction](../particle_primary_charged_string/PROOF.md)
supplies the physical operators

    B_xy,ab = c*(p_x,a) c(m_y,b) chi_xy P0_xy.                 (1)

Their endpoint matter charge increments are (+6,−6), while every
global native occupation is conserved. The inside factor is odd in
the canonical anticommutation relations and transforms as a two-component
spinor. The outside factor is odd too; (1) is even and gauge invariant.
The link carries the compensating flux. The outside endpoint and link
are retained in the experiment and in its quantum correlations with the
reference. They are not a classical background or a replacement state.

Summing Gauss constraints over the full open source gives total native
charge zero. A single operator with nonzero total charge therefore has
P_G A P_G=0. The zero-charge ghost endpoints cannot absorb that charge.
Replacing a ghost by an external charged memory would require a declared
different Gauss domain and its Hamiltonian. Equation (1) avoids that
change by retaining the opposite charge in the original complement.
The finite [seam Gauss theorem](../../Lean/Screen/DiscreteGauss.lean)
expresses a related neutrality condition on its rational seam graph;
it does not supply the present source Hilbert space or its spin dynamics.
Here Gauss matching follows directly from the native character and CAR.

## 2. The exact spin coefficients and complete apparatus

In the actual beta chart the proper binary cube lift restricts to the
same two-dimensional unitary g on the p and m rows. Its character is G1,
and the central 2pi element is −I. This is a finite source transformation
law, without an identification with Lorentz spin. The outside annihilator
c(m_b) transforms by the conjugate representation. For

    t_0=I_2, t_i=−i p* gamma_(i+1) gamma_(i+2) p,
    B_xy,a = (1/sqrt(2)) sum_(b,c) (t_a)_bc B_xy,bc,          (2)

the four t_a are the source identity and Pauli matrices. The exact
source spin action is

    g t_0 g* = t_0,
    g t_i g* = sum_j R_ji t_j.                              (3)

Thus the internal probe coefficients carry A1 plus T1. A source rotation
also rotates x, y and the link. These are probe transformation labels;
the mixed original source is not assigned a final total-spin eigenstate.
In particular, the full even pair has central sign +1. The two-dimensional
regional charged endpoint carries G1, without an H component in this
single native spinor fiber. Orbital or composite additions require their
own coefficient spaces.

On all sixteen occupations of the four active endpoint modes, each
matter polynomial A_a in (2) is a contraction and satisfies

    A_a* A_a >= I − D_xy/2,
    D_xy = sum_b n_p(x,b) + sum_b (1−n_m(y,b)).               (4)

These are the exact full-Fock identities from the
[pair coefficient construction](../particle_primary_pair_channels/PROOF.md),
also executed here. They hold for every unitary t divided by sqrt(2).
They do not replace the original state by the filled-beta anchor.

The joint menu has 24 directed bonds and four coefficients:

    K_xy,a = B_xy,a/sqrt(96),
    K_fail = sqrt(I − sum_(xy,a) K_xy,a* K_xy,a).             (5)

It has 97 outcomes. Alternatively, a fixed coefficient a gives a complete
25-outcome apparatus with K_xy=B_xy,a/sqrt(24) and its own failure.
The successful history probabilities in (5) for a fixed a are exactly
one quarter of those in that fixed-coefficient experiment. This scaling
also holds with identical initial and final clocks and any reference.
It does not assert equality of the two failure maps.

Every success commutes with all gauge constraints, global native
occupations and total fermion parity. So does the positive square-root
failure. The isometry sum_r K_r tensor |r> extends to a unitary separately
in the finite conserved blocks, using a neutral outcome register. All
zero-weight outcomes and failures are retained.

The exact 48-dimensional local CAR/character control has a five-dimensional
Gauss kernel. The full spin coefficients send its neutral anchor to four
orthogonal charged spin-pair states. Each endpoint's spin density in this
anchor control is I_2/2. This illustrates why a spin coefficient label is
not a definite pre-existing endpoint spin. It is an algebra control, not
the primary selected histogram.

## 3. An original probability for each spin coefficient

Because P0 commutes with the matter operators, (4) implies

    I − B_xy,a* B_xy,a <= D_xy/2 + (I−P0_xy).                (6)

Average (6) over all 24 directed bonds. Each native-six positive component
appears three times as a target, and each negative component three times
as a source. Every link appears twice. Therefore

    1 − (1/24)sum_xy Tr(rho B_xy,a* B_xy,a)
      <= Tr rho [D6/16 + (5/18)sum_link Casimir].             (7)

The charged-string full-carrier occupation/electric allocation gives
the same upper bound on the right for every a. Its authenticated selected
source mean-energy and variational estimates imply

    w_a >= w > 0.79446654.                                  (8)

No covariance assumption on rho is needed for (7). No selected-spin
occupations, constant links or normalized success state replace rho.

## 4. Full-source energy moments for the spin coefficients

The onsite and selected-link electric contribution has the exact center

    [H_onsite+H_E,B_xy,a] = (87/4) B_xy,a.                   (9)

This follows from two beta defects costing six and the character Casimir
3/10 contributing 63/4. It is a pair transfer center, not one endpoint's
rest energy. Character degree one and the constant character lie in the
actual cutoff. B commutes with that product projection, so the compressed
commutator is the compression of the parent commutator.

For the selected bond, Haar orthogonality and the complete native hopping
block give ||T P0||<=sqrt(40) and the same bound on the character image.
The two terms TB and BT have orthogonal input and output supports. Since
||B||<=1, their commutator has norm at most sqrt(40), as in the
charged-string proof. The two magnetic plaquettes contribute at most4/7.

For the remaining incident bonds, the actual Wilson columns map each
two-dimensional p or m spin space isometrically into its neighboring
modes. The two endpoint image spaces are disjoint on this eight-cell
graph and are orthogonal to the original endpoint modes. After the
unitary coefficient t is absorbed into the source and image spin bases,
the remaining commutator is the sum of four disjoint normalized raising
pairs, each with coefficient1/sqrt(2). On a fixed block, a pair that is
empty or doubly occupied contributes zero; a singly occupied pair is a
spin-one-half raising matrix. For r active pairs the maximum squared
norm is

    r=0,1,2,3,4:  0, 1/2, 1, 2, 3.                       (10)

Equation (10) is checked by exact matrices in every active block. It is
a full-occupation norm bound. Using sqrt(2), the value for a single
matched-spin polynomial, would omit active terms in (2).

Consequently, on the complete physical carrier,

    ||[H,B_xy,a]−(87/4)B_xy,a||
       <= sqrt(40)+sqrt(3)+4/7 < sqrt(75).                  (11)

For each fixed coefficient's original two-energy transfer measure,
stationarity of rho turns (11) into

    integral (omega−87/4)^2 dmu_a(omega) <=75.               (12)

Chebyshev at radius13 and the unchanged source clocks give

    mu_a([35/4,139/4]) >= w−75/169,
    P_a(recorded transfer in [69/8,279/8])
       >= w−75/169−1/32 >5/16.                             (13)

The clocks have resolution1/16, tail allowance1/128 per clock and
implementation allowance1/128 per clock. Every fine initial-energy,
success/failure and final-energy record is retained. The joint97-outcome
probe menu therefore has original weight >5/64 in each chosen coefficient,
>5/64 in the scalar coefficient, and >15/64 in the three vector coefficients,
in the same reported transfer window. These disjoint events have total
lower bound5/16, not four times5/16.

## 5. Regional charge, spin readout and exchange domains

For every coefficient at a directed bond,

    [Q6(z),B_xy,a]=6(delta_zx−delta_zy)B_xy,a.                (14)

The [two-coupling charge meter](../particle_primary_charged_readout/PROOF.md)
uses one radix-five register per cell. Its initial unknown charge cancels
coherently; the proved alphabet {−6,0,+6} is decoded exactly. Failure has
zero increment. This retains the entire original source/reference history,
so (13) is an original joint charge-transfer, spin-coefficient and energy
record bound. The positive event includes a successful probe outcome.
An initial absolute-charge measurement is not inserted.

There is also a complete raw-spin apparatus using the four entries (1)
instead of the four coherent coefficients (2). For one chosen spin axis,
define the actual regional doubled-spin operators

    Z_p(x)=n_p(x,1)−n_p(x,2),
    Z_m(y)=n_m(y,1)−n_m(y,2).

The spin-projection operators in this finite source representation are
Z/2. The exact selection rules are

    [Z_p(x),B_xy,ab]=s_a B_xy,ab,
    [Z_m(y),B_xy,ab]=−s_b B_xy,ab,    s_1=1, s_2=−1.       (15)

Each regional spin change is therefore plus or minus one half in this
finite source convention. These gauge-invariant operators commute with
all local native charges and with each other. The raw-spin success sum
on a directed bond is

    sum_ab B_xy,ab* B_xy,ab
       =P0_xy (2−N_p(x)) N_m(y),                            (16)

so the complete square-root failure commutes with every Z_p and Z_m.
Its increments are zero. One radix-three pointer per beta sign and cell
reads the alphabet {−1,0,+1} by the same exact two-coupling cancellation.
The sixteen spin pointers commute with the eight radix-five charge
pointers. Neither unknown initial spin nor initial charge is measured.
All source/reference histories of the raw-spin experiment are preserved.
This reads actual charge and spin-projection changes during the probe,
not a final spin eigenvalue after the energy clock. Noncommuting spin
axes are different apparatus settings.

For a fixed raw entry a,b, its missing-success bound is

    I−B_xy,ab* B_xy,ab
      <=n_p(x,a)+(1−n_m(y,b))+(I−P0_xy).                    (17)

The selected rho is invariant under actual cube rotations. Summing a
local spin vector over all eight sites gives a vector whose cube average
is zero, since the G1 matrix twirl sends every traceless2×2 matrix to zero.
Thus each chosen spin row contributes half the corresponding global
positive-particle or negative-hole deficit. Average (17) over the24
directed bonds to recover exactly the right side of (7). Each fixed raw
entry has success at least w. Its other-bond residual is sqrt(2), since
only two disjoint raising pairs occur. The charged-string moment70 and
its unchanged complete clocks consequently give

    P_(a,b)(recorded transfer in [77/8,271/8])>1/4.          (18)

In the complete96-success raw-spin menu, each of its four spin entries
has weight >1/16 in that same band, summed over directed bonds, with
coherent native charge and doubled-spin increments. The four entries
are normalized by1/sqrt(96); no per-entry probability is multiplied by
four. The raw-spin and coherent A1/T1 readouts are two apparatus choices
with the same summed success effect. Their different spin-record bases
are not asserted to have identical joint energy distributions.

A terminal regional spin apparatus can resolve the positive-beta q6
doublet's occupations zero and two and, in occupation one, the two
eigenprojectors of n dot sigma for a chosen unit axis n. These four
projectors are gauge invariant and complete on that doublet, with the
identity on all other modes. They can follow the final energy clock.
Summing their records preserves every earlier classical energy/probe
probability, while their source/reference disturbance is retained. They
are not claimed to commute with H or to report a spin that persisted
through the final energy operation.
The positive bound (13) sums the terminal occupation outcomes; it supplies
no positive lower bound for occupation one in that same final-energy event.

The regional charged factors satisfy the actual CAR: disjoint creators
anticommute, and disjoint creation/annihilation factors anticommute.
Even disjoint Wilson factors do not change this sign. The global
compensated pair operators instead have even parity and commute at
disjoint supports. This distinguishes regional fermion exchange from
the algebra of the complete neutral instrument. It does not establish
a continuum spin-statistics theorem or a single-particle scattering state.

The [literal source charge conjugation](../particle_primary_charge_conjugation/PROOF.md)
preserves H, rho, Gauss support and the cutoff. Its conjugated complete
apparatus has exactly the same original classical records and reversed
endpoint charge increments. The idle reference is retained; the two full
outputs are related by the corresponding source/reference unitary, not
asserted identical. This conjugated creation experiment is distinct from
replacing B by B*, which reverses transfer frequency.

## 6. What this energy can identify

The construction supplies an actual gauge-complete charged regional
spinor, an explicit compensator, complete spin-coefficient records and
positive original energy weight on one fixed source. The recorded energy
includes the link and both endpoints. Its exact electric contribution
63/4 already prevents treating the bare onsite value six as the entire
pair response. Neither dividing the transfer by two nor reading87/4 as
a mass has a proved physical basis.

A single charged particle mass requires a controlled charged-sector
threshold or branch, a treatment of its dressing and compensator, a
dispersion relation and a physical unit reference. A chiral physical
dictionary also has to satisfy the source's native parity-pairing
constraint. The present finite G1 endpoint and native Q6=6 do not alone
identify an electron, positron or a measured GeV mass.
