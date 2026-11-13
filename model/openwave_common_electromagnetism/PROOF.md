# Native charge and longitudinal electric energy

## Matter source and exact charge

Let H0 be the entire 192-mode source in `openwave_common_mass_source`,
including its positive strong terms, all three charged and calibrated neutral
Dirac families, and the local semileptonic current with G=1/64. Its four strong
calibration targets, retrospective charged analyzer phase and neutral
oscillation inputs are unchanged. The exact native charge is T3+Y, with
values 0, ±1, ±1/3 and ±2/3. The sterile coordinate has charge zero.

The executable checks the complete formal bilinear polynomial commutator
[Q,H0]=0. Cross terms of the first two weak generators cancel before any
occupation restriction. It also checks qA=Aq, qB=-Bq and tr(q BB*)=0 for
the actual local Dirac transformation c=Aa+Ba*. Hence native local electric
charge equals sum q_j a_j* a_j on the lepton Fock space, without a reference
constant. Each neutral coordinate, conjugate, spin and family is retained.
The common source bound |V|<=18 N_colourful is inherited with its exact
64-dimensional characteristic-polynomial proof.

## Compact Gauss realization, including its reference

For V matter vertices add V full integer rotors, with E_x in Z and
[E_x,U_x]=U_x, and one full integer quantum reference R. Impose

    E_x - 3 V Q_x - R = 0,       sum_x E_x = 0.

All constraints commute. For each joint matter occupation, their unique
solution is R=-3 sum Q_x and E_x=3 V(Q_x-average Q). Because 3Q_x is
integral, the solution exists for every matter occupation. The map W taking
an occupation to that occupation together with these flux/reference labels
is unitary onto the complete physical constraint kernel. This proves
surjectivity and retains all matter Fock complements. The orthogonal
constraint-failure projector is an explicit zero-probability outcome for a
physical preparation.

For a matter matrix element changing charges by Delta Q_x, dressing multiplies
it by rotor shifts Delta E_x=3V Delta Q_x-3 sum Delta Q and a reference shift
Delta R=-3 sum Delta Q. Thus a charge-q hop y to x carries
U_x^(3Vq) U_y^(-3Vq). Cross-cell weak raising/lowering products require the
same dressing. Dressing only the hopping would fail. Local neutral
semileptonic terms require none. WH0W* is the exact restricted source;
equivalently every finite native polynomial term has the specified rotor
shift. No rotor truncation defines the theory. Gauge transformations generated
by each constraint leave the dressed action invariant. Global charge is
conserved, and the finite continuity equation is the operator identity
Qdot_x=i[H,Q_x]; hopping and cross-cell weak exchange supply its currents.
Summing this identity gives zero by the checked global commutator.

This reference implements a uniform compensating background. It is a physical
part of the declared finite construction, with an explicit recorded charge.
It is not an isolated-charge boundary condition at spatial infinity.

## Electric energy and the fixed self-energy subtraction

For the two-cell graph use L=[[1,-1],[-1,1]] and its Moore-Penrose inverse
G_L=[[1,-1],[-1,1]]/4. Declare kappa=1/8 and

    H_E = kappa E^T G_L E/(18 V^2).

Under W this is kappa Q^T G_L Q/2=kappa(Q0-Q1)^2/8. The electric divergence
identity L G_L Q=Q-average Q is exact. This is longitudinal energy with a
nonlocal Green metric on the rotor coordinates; there are no independent
transverse photon or magnetic-loop coordinates.

Define delta=kappa/8=1/64 and D=sum q_a^2 n_a, using native colourful
occupations and positive-Dirac lepton occupations. The declared source is

    H = W(H0 + delta[(Q0-Q1)^2-D])W*.

The subtraction is the diagonal self-energy of this exact Green matrix.
Its coefficient is fixed by the architecture, without an empirical mass
refit. The bare control omits -delta D and gives the exact one-charge rest
shift delta m_e=0.00000798435860453125 GeV. This control exposes the effect
of the chosen normal ordering. It is not a radiative renormalization proof
or a calculation of the physical fine-structure constant.

The electric square is positive and q^2<=4/9 on colourful modes, so

    H >= (m_min-18G-4delta/9) N_colourful
           +(1-delta) N_charged + H_neutral.

The first coefficient exceeds one. In the colourful vacuum, one charged
excitation has (Q0-Q1)^2=D=1, independent of its position superposition;
neutral excitations do not change this identity. Its exact energy is its
analyzer mass. Two or more charged excitations cost at least
2(1-delta)=63/32>1. Therefore the full physical charge ±1 threshold is
exactly the anchored electron mass. The positive neutral masses give a unique
finite physical vacuum and the inherited finite neutral gap. This is a
matter-plus-longitudinal model, so the statement does not assert a gap in a
photon-containing theory.

## Complete interacting baryon energy records

The source executes the complete 2024-dimensional three-quark sector and
its complete 120-dimensional colour-singlet reduction. The electric term
is diagonal in the native occupation basis and the intertwiner is checked.
All 120 strong-plus-electric eigenvalues are retained as floating diagnostics.
The baryon-plus-lepton Hamiltonian is larger; those eigenvalues are not its
interacting poles.

Use the exact integer neutron trial v and its literal native-current proton
partner from the common mass source, with the identical pure Dirac reference.
The local lepton charge annihilates that reference. Its weak means vanish
and its current covariance is 3I/4. Consequently the exact full-source moments
are

    mu = <v,(H_strong+C)v>/<v,v>,
    variance = ||(H_strong+C-mu)v||^2/<v,v> + G^2 kappa_current.

All sums in these formulas are rational and executed. The electric shifts
are derived, without refitting the four strong inputs. The trials are not
assumed to be stationary states of the changed source.

A source-native occupation seed followed by the trial projection has success
s=|v_seed|^2/||v||^2. The receipt retains failure 1-s, an energy window of
half-width m_e/2, and its complement. Chebyshev bounds give conditional
inside probability at least 1-variance/(1/2)^2 and original inside probability
at least s times this quantity. The probability bounds are not assigned exact
probabilities. The nonselective preparation work is evaluated exactly; the
spectral-window instrument commutes with H and has zero mean work.
Dressing every preparation and effect by W preserves their complete laws
and any external quantum reference. The induced channels are completely
positive and trace preserving on the physical space. Adding the
constraint-failure record makes the apparatus complete on its kinematic
extension. No observed nucleon pole mass, transverse photon or continuum
Lorentz symmetry follows from these finite identities.
