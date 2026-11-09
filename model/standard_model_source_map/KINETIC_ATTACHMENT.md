# Same-current Weyl frames and the mirror boundary

The fifteen internal states in `native_matter_attachment.representation()`
have an exact current-preserving attachment to a configured free Weyl
kinetic law. The unitary spin/internal frame, its free Lorentz current
identities, and the distinction between mirror fields and antiparticles
are checked here. No interacting gauge vacuum, mirror removal, physical
mass or Standard Model particle identification is inferred.

## Source architectures and their shared current map

The auxiliary even exterior register supplies the internal matrices. The
propagating fields use separate odd CAR modes, ordered

    spin*45+3*matter_index+family,
    C^2_spin tensor M15 tensor C^3_family.              (1)

The source-selected occupation-mask ordering of M15 is mapped into the
configured chiral action by the permutation P15 from
`native_matter_attachment.action_intertwiner()`. Its twelve exact
anti-Hermitian matrices obey

    T_A^action P15=P15 T_A^native.                     (2)

The same package's `kinetic_attachment()` proves the P90=I2 tensor P15
tensor I3 kinetic and generic complex Yukawa identities. This result
does not execute that mass-parameter producer or change a Yukawa value.

The configured continuum-momentum menu in
`common_source/chiral_candidate/kinetic.py` has h(p)=-v sigma.p tensor I45,
with v=Omega/2=log(5)/2 in its declared coordinate/modular-time convention.
Its finite menu is not closed under all Lorentz boosts. The separate local
sine regulator has h(k)=-c sum_j sigma_j sin(k_j) tensor I45 and eight
nodes. Native Gaussian Fock words implement its quadratic factors with
full fermion routing signs. A chosen scaling c=v/a is an explicit source
parameter/coordinate prescription; the field analysis does not derive
it from FIFO memory capacity.

The gauge-linked staggered action in `native_gauge_continuum.py` is another
declared regulator. Its constant-background compact-band limit retains two
Dirac tastes per internal component. It must not be identified with the
eight-valley two-spinor regulator by equating their internal dimensions.
Operator-valued quantum links are also distinct from prescribed unitary
background links. A nilpotent quantum-link raising operator cannot be set
to the identity, or replaced by an identity vacuum expectation, to import
the free dispersion.

## Exact valley and free Lorentz current intertwining

At node k_b=pi b, b in {0,1}^3, the sine symbol has principal part

    h_b(p)=-v sum_j (-1)^b_j sigma_j p_j,
    chi_b=-(-1)^(b1+b2+b3).

Use the actual constant SU(2) lift U_b of
`weyl_geometric_modular.algebra.frame`, for which
U_b^dagger h_b U_b=chi_b v sigma.p. Then the map from native/source
coordinates to canonical/action coordinates is

    W_b=U_b^dagger tensor P15 tensor I3.               (3)

It is unitary, preserves the CAR metric, intertwines all three kinetic
axes, and satisfies

    (I2 tensor T_A^action tensor I3) W_b
      =W_b (I2 tensor T_A^native tensor I3).           (4)

The executable verifies all 96 current identities and 24 kinetic-axis
identities on the full ninety-row one-particle space. No internal component
or family is removed. The source mass/kinetic dictionaries are reused;
the default two-flavour interacting centered ground is not relabeled as
this ninety-mode architecture.

In coordinates x0=v*t, let chi=+/-1 and define

    sigma_chi^mu=(I,chi sigma_1,chi sigma_2,chi sigma_3),
    R_j=-i sigma_j/2, B_j=chi sigma_j/2.

These satisfy the Lorentz Lie algebra and preserve epsilon=i sigma_2 by
L^T epsilon+epsilon L=0. For a boost
S=cosh(u) I+chi sinh(u) sigma_j,

    S^dagger sigma_chi^mu S=Lambda^mu_nu sigma_chi^nu,  (5)

where Lambda is the boost of rapidity 2u in direction j. The checked
rational example cosh(u)=5/3,sinh(u)=4/3 has determinant one and preserves
diag(1,-1,-1,-1). The same formula follows for every real u by Pauli
multiplication. For each actual Hermitian internal charge Q_A=iT_A,

    j_A^mu=psi^dagger(sigma_chi^mu tensor Q_A)psi       (6)

therefore transforms as a Lorentz vector and obeys the free continuity
equation. Equation (4) identifies these currents with the retained source
ones. Physical model-time spatial flux includes the factor v.

A pointwise spinor boost is nonunitary: diag(3,1/3) changes the ordinary
two-spinor norm. A unitary boost on the CAR field representation also
transforms momentum/Cauchy data and uses the compatible free vacuum.
Equations (3)-(6) do not construct that transformation on a finite menu.
The existing complete free-field theorem provides it for its own admitted
source family, vacuum and clock. An interacting M15 gauge source requires
its own state and continuum attachment; a principal-symbol identity does
not supply that premise.

## Mirror fields are not the holes of a retained chiral field

With Y=q6Y/6 and Q=T3+Y, the actual one-generation internal hypercharges are

    (-2/3)^3, (-1/2)^2, (1/6)^6, (1/3)^3, (1)^1,

where exponents denote multiplicities. Direct evaluation gives

    Tr_M15 Y=Tr_M15 Y^3=0,
    Tr_M15 Y^5=5/9.                                   (7)

The hole/conjugate representation has Y_conj=-Y^T and fifth trace -5/9.
Unitary equivalence preserving hypercharge would preserve this trace, so
M15 and its conjugate are inequivalent. The fifth trace is a representation
invariant, not a four-dimensional gauge anomaly. The anomaly traces are
the established separate checks of the internal chiral representation.

Every positive-energy sine-regulator valley carries the same M15. A
positive-winding valley with M15 cannot be identified, preserving the
source currents, with the antiparticle of a negative-winding M15 field:
that antiparticle carries the conjugate representation. Opposite winding
therefore supplies genuine mirror field copies in this regulator. Nambu
doubling only represents the holes at opposite momentum of each existing
field; it adds no independent fields. The two spatial spinor components
also add no internal multiplet multiplicity to anomaly traces.

Counting all three family coordinates, the supplied continuum menu has
45 left-Weyl internal components and ninety Cauchy rows. Its 180 Nambu rows
are redundant particle/hole coordinates. The local sine regulator has
eight times this number of Weyl fields: 180 components of each chirality,
or 180 vectorlike Dirac pairings. These counts are internal/color/family
counts, not measured species assignments.

## Original mirror weight and failure of a local gauge restriction

For one internal/family component on an even periodic free lattice, use
the complete local instrument

    {a_(x,s)/sqrt(2), a_(x,s)^dagger/sqrt(2): s=0,1}.

Its four effects sum to I by the CAR, on the entire Fock space. In the
filled negative-band reference, addition and removal use complementary
spin projectors. Their spin-averaged combined momentum measure is uniform:

    (1/2) Tr_spin[P_+(k)+P_-(k)]=1.                   (8)

Removal momentum is the negative of the removed mode momentum. This
permutation preserves the uniform measure. For side L divisible by four,
the eight half-open pi-shift valley cells each contain L^3/8 grid modes.
Each cell has original combined probability 1/8. Excluding its one exact
zero mode leaves 1/8-1/L^3. Thus the four opposite-winding cells have
strictly positive-energy original weight 1/2-4/L^3. All 4L^3 spin,
charge and momentum records are retained. Translation-invariant choices of
zero-mode occupations do not change the combined law (8). These whole-cell
weights do not assert a shrinking physical-energy-window floor. The
calculation runs exact finite source
matrices and identities; no additional large native tape is claimed.

A valley filter can select a response, with its complement recorded.
It does not remove those fields from the complete source law. Let P be
one full valley-cell Fourier projector and X_x the projector
onto one lattice site. Since P_xx=1/8,

    ||[P,X_x]||_HS^2=2 P_xx(1-P_xx)=7/32.

The admitted local hypercharge generator is X_x tensor iY. The actual
source trace Tr(Y^2)=10/3 gives

    ||[P tensor I15,X_x tensor iY]||_HS^2=35/48>0       (9)

per spin and generation. Thus this valley restriction is not reducing for
all local gauge currents. It cannot be declared a local gauge-covariant
mirror-removal operation. Source packets can suppress unwanted valleys
in particular observations; this is not suppression of their original
weight in the full theory. The admitted finite-menu and staggered
constructions supply no regulator theorem deleting these mirrors. A
different selected regulator or an explicitly restricted observable/state
dictionary needs its own locality, gauge and continuum proof.
