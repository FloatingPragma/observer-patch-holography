# A declared fractional-charge Dirac/current source

This package constructs the previously missing kinematic and local-current
service as an explicit new source. It does not transfer the old Maxwell
vacuum, prove a baryon ground state, or predict a physical axial ratio.

## 1. Carrier, operators and common unit

At each cell the quark carrier is

    C²_Dirac tensor C²_physical_spin tensor C²_flavour tensor C³_colour.

The Dirac coordinate is independent of position, weak flavour and colour.
There are24 complex CAR modes per cell. With the Dirac coordinate first,

    beta=Z tensor I, alpha_i=X tensor sigma_i,
    gamma5=X tensor I, Sigma_i=I tensor sigma_i.

The existing autonomous Dirac Clifford convention is intertwined exactly
by swapping its two spinor factors. All four Clifford coordinates remain.
The vector, axial and left charged currents are the literal one-body
matrices

    V+=(I,alpha_i) tensor tau+ tensor I_colour,
    A+=(gamma5,Sigma_i) tensor tau+ tensor I_colour,
    J_L+=(V+−A+)/2.

Their second quantizations define the corresponding CAR bilinears. They
commute with every local colour generator and obey[Q,J+]=J+ with physical
Q=diag(2/3,−1/3), repeated over Dirac/spin/colour. The charge change is one
in the same electron unit. A lepton reference of charge−1 supplies the
opposite charged current for a neutral semileptonic contact.

The common bare mass is the exact rational midpoint of the inherited432
register light coefficient, about6.03008e−5 in its existing raw unit. It is
a declared common up/down mass, not a quark-mass fit. The spacing/clock is
h=1/8 and Hcell=h Hphys. Every auxiliary coefficient is expressed in the
same Hphys unit: g=1/8 and Delta=3kappa, so its cell-clock coefficient is
h times that number. The physical electron-anchored unit is inherited only
as a conditional unit map; no elapsed time in seconds is derived here.

## 2. Physical fractional charge without relabeling the electron

The already established group is

    G=(SU3 tensor R)/<(zI,2pi)> isomorphic to U3,
    W=exp(−i theta/3)S, z=exp(2pi i/3).

Use its actual representations

    R_down(W)=W,
    R_up(W)=det(W)^−1 W,
    R_electron(W)=det(W).

Along the same physical electromagnetic parameter theta with S=I, these
are exp(−i theta/3), exp(2i theta/3), and exp(−i theta), respectively. Thus
the colourless electron retains charge−1; no separate rescaling changes
its unit. The up/down representations close under the joint centre
identification because exp(2pi i q)z=1. Their conjugate hole
representations have opposite electric charge and conjugate colour.

This is a representation and charge/current construction. It is not an
isometry from the old noncompact193-mode Maxwell field to the new compact
link field. In particular, the old and new electric spectra and carriers
are not identified by changing a label from1 to3.

## 3. Gauge-covariant Wilson hopping and both bands

Let Wxy be a U3 link and R(Wxy)=Rup direct-sum Rdown. For each positive
oriented edge along axis i, the cell-clock one-body matrix contains

    Hxy=−[r beta+i alpha_i] tensor R(Wxy)/2,
    Hyx=Hxy†,

and adds r beta/2 at both endpoints. The mass is h m beta. Take r=1.
Under an arbitrary local gx in U3,

    Wxy -> gx Wxy gy†,
    psi_x -> (I_Dirac,spin tensor R(gx)) psi_x,

the complete one-body matrix transforms by the corresponding block
unitary. The code checks this on all192 one-particle coordinates of an
eight-cell periodic lattice with24 independent nontrivial links and
independent local transformations. Its largest covariance residual is
2.23e−15. The local weak current transforms with character det(gx)^−1,
which is precisely compensated by the electron's determinant character
in the semileptonic product.

At identity links the free Bloch symbol is

    Hcell(k)=beta[h m+r sum_i(1−cos k_i)]
                    +sum_i alpha_i sin k_i.

All positive and negative bands are retained. The L4 control keeps all64
momenta and all12 states per sign at each momentum. The r=0 control has
eight mass minima; r=1 has one. The two-cell periodic lattice alone is not
used to infer unaliased continuum dispersion.

There is also a uniform background-link statement. Writing H=beta M+K,
M=h m+r times the positive covariant graph Laplacian and beta K is
anti-Hermitian. Thus

    Re <v,beta H v> = <v,Mv> >= h m ||v||²,

which proves ||Hv||>=h m||v|| for every vector and every unitary link
background. The grading Y_Dirac tensor I anticommutes with H and commutes
with electric charge, so both charge sectors have paired positive/negative
spectra. No negative band is removed by this gap argument. Complete
band-frame transformations of every vector, axial and left current are
checked, including their positive/negative off-diagonal blocks.

## 4. Dynamical links, electric normalization, Gauss law and source clock

The dynamical link Hilbert space is full L²(U3) on every edge. No
Peter-Weyl representation is truncated. Gauge transformations act by the
left/right regular link actions, the matter representations above, and
the adjoint action on each local auxiliary. The physical domain is the
invariant subspace of their combined local action; this defines all colour
and central electric Gauss generators with their orientation fixed by
the displayed link transformation.

The field-energy completion is a declared new architecture. In Hphys
units it is

    H_g = sum_links [kappa C_SU3 + (a/2) Q_link²]
           +3b sum_plaquettes [3−ReTr W_plaquette],
    a=1/2, b=3, kappa in{1/2,1,2}.

Here Q_link generates translation in the same physical theta used for
the electron. Along a central plaquette Wp=exp(−iF/3)I, the magnetic term
is9b[1−cos(F/3)]=bF²/2+O(F^4), with absolute quadratic remainder at most
bF^4/216. The electric energy is(a/2)Q_link² in physical charge units.
Thus the central small-field normalization is explicit; simply counting
charges in thirds is not used to hide a coupling change. These coefficient
values and local expansion do not establish a photon mode, force law,
old-field state map, or a physical fine-structure constant. The compact
carrier and spatial graph are new. Field-energy coefficients were declared
as this source completion, not fitted to the background covariance test.

The full source is

    Hphys = H_g + dGamma(Hcell_Dirac)/h + sum_x Haux,x,
    Hcell_total = h Hphys.

The auxiliary quark current is now
Gai,x=psi_x†(Sigma_i/2 tensor I_flavour tensor Ta)psi_x. It couples to the
same25-state adjoint/vector auxiliary with its actual g and Delta. Each
term is gauge covariant/invariant under the combined action.

For the free filled-reference charge convention, introduce one colourless
reference of charge−2 per cell, represented by det(gx)^2. The free filled
negative bands carry+2 per cell, so the normal-ordered matter charge is
dGamma(Q)−2 per cell. This reference fixes the Gauss background; it is
not inferred from a numerical vacuum minimization. The interacting ground
has not been equated with the free filled state. Charged excitations require
their retained compensating reference or boundary charge.

At finite volume the positive electric Casimirs on the compact link group
have compact resolvent. Magnetic, matter and finite auxiliary terms are
bounded operators on that field carrier. The complete Hamiltonian is
therefore self-adjoint, bounded below, and has an attained ground in every
nonempty closed gauge sector. The invariant sector is nonempty: a local
three-up colour singlet has determinant character det(g)^−2, cancelled by
the prescribed reference. This is an existence statement, not a computed
vacuum, a mass gap, or a baryon-spectrum prediction.

## 5. Exact old auxiliary operator control and complete small-CAR execution

Embedding the old12 one-quark coordinates into the upper Dirac component
gives exactly its old colour-spin Gai and vector/spatial-axial matrices.
With kinetic terms off, the three-quark local algebra recovers the earlier
220×25 auxiliary source plus the scalar3m. This is a local operator
intertwiner. It does not transfer that source's state to the new gauge
background or identify its fixed-three-quark branch with the new vacuum.

The executed small CAR member is an open z bond with both Dirac
coordinates and both flavours, one fixed physical spin and one colour:
eight quark modes. Two lepton-reference modes make1024 full occupations.
This reduction is invariant for that identity-link one-particle z
Hamiltonian; currents that change the omitted physical spin are tested
instead on the complete24-coordinate local carrier.

The small control constructs every CAR annihilator, diagonalizes the
charge-resolved one-body matrix, and forms all positive-particle and
negative-hole operators. Their full CAR, Hamiltonian reconstruction and
charge reconstruction are checked on all1024 occupations. Particles have
charges(+2,−1)/3 and their holes(−2,+1)/3. The filled Slater state is
explicitly checked as a free reference. Six complete finite-time/contact
controls retain every occupation probability, including pair excitations
relative to that reference. None is called the full interacting vacuum.

All1024 occupation labels also retain left/right charge, electric flux,
and the quantum boundary-reference label. On the open U1 control the
Gauss map is E=3Qleft, R=3(Qleft+Qright). Every one of4096 hopping matrix
elements changes flux by its actual charge; holding flux fixed fails the
constraint. All semileptonic matrix elements conserve total charge, and
the512 entries moving charge between the remote lepton reference and the
other cell carry their required flux shift. This is a complete U1
kinematic dressing control; its fixed-colour slice is not presented as a
standalone colourless fractional-charge representation of U3.

The matter/flux map extends unitarily to all integer flux/reference labels
by controlled translations. It therefore defines the complementary
constraint outcomes as well. This small control uses no electric energy
and is not a simulation of H_g. Its full native implementation is not
claimed executed.

## 6. Resource and remaining calculation

At eight cells the full matter Fock dimension is2^192; all local
auxiliaries add25^8, and24 L²(U3) links retain infinitely many electric
representations. The192-dimensional one-particle code alone fits in eight
logical qubits with64 padding states, for which energy16 is declared in
cell units. That small code is not the many-body source.

Only the exact Clifford basis map, quotient charge representation, local
current/auxiliary intertwiners and finite CAR/gauge controls are transferred
here. An existing logical compiler can encode finite matrices, but no
full gauge cutoff certificate, native complete-source pulse word, physical
clock calibration, or vacuum preparation has been produced by that fact.
The next physical calculation must ground the full constrained source,
retain pair/auxiliary/field sectors, and measure its baryon vector and axial
currents. A positive-band three-particle compression would not meet that
requirement.
