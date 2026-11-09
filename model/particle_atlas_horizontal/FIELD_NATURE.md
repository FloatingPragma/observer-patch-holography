# What the named fields do

The atlas follows source operations into a declared field action and then
into detector records. Its named-field panel explains the exact quadratic
consequences of the [conditional Standard Model action](../particle_candidate_atlas/CONDITIONAL_STANDARD_MODEL.md).
The action, couplings, preparation and apparatus may be architectural
choices. Their consequences must be calculated on that construction. A
coupling chosen using a measured mass is a calibration input, even when the
subsequent calculation is exact.

[`field_nature.py`](field_nature.py) reuses the owning action Hessians and
free-field Clifford calculation. It imports no numerical laboratory masses.
[`field_nature.json`](field_nature.json) contains eighteen named target rows,
their antiparticle meanings, diagram descriptions and references. The
formulas use canonical kinetic terms and natural units, except where the
speed parameter is written explicitly. A rest energy \(\mu\) corresponds
to inertial mass \(\mu/c^2\).

## From repair records to a field mode

A carrier record is a finite source event. The
[native current attachment](../standard_model_source_map/PROOF.md) constructs
an internal matter register whose exact gauge weights match one left-Weyl
Standard Model family. Its
[fermion attachment](../standard_model_source_map/FERMION_ATTACHMENT.md)
acts on every occupation sector, and its
[gauge attachment](../standard_model_source_map/GAUGE_ATTACHMENT.md)
couples the same currents to finite operator-valued links. These are
computed finite interfaces. A chosen field action adds a spatial and
temporal rule to them. A particle interpretation also uses the state,
propagating spectrum, transformations and detector response of that rule.

Protection and propagation answer different questions. A protected defect
can preserve a label as the source evolves. A propagating field excitation
has a dispersion relation, charge and detector signature. A map between
them must send the preserved source label and the actual source evolution
to that same field response. A topological invariant alone supplies neither
a particle name nor a rest energy.

The finite Maxwell construction provides a concrete field precursor. Seam
potentials \(A\) have face curvature \(CA\), and \(CD=0\) makes it
invariant under \(A\mapsto A+D\chi\). Under gauge invariance, strict
face-sharing Hessian locality and carrier rotation covariance, the quadratic
magnetic action has the form \(\beta\lVert CA\rVert^2/2\). With positive
electric and magnetic energy weights, Faraday kinematics, conservation and
the specified time reversal, its linear evolution is

\[
 \dot E=(\beta/\epsilon)C^{\mathsf T}B,\qquad \dot B=-CE.
\]

The [finite action classification](../maxwell_local_selection/PROOF.md)
and [weak-field estimates](../maxwell_weak_field/PROOF.md) state these
premises explicitly. The Lean
[face-curvature result](../../Lean/Screen/LocalFaceMaxwellAction.lean),
[Gauss propagation and wave law](../../Lean/Screen/TemporalMaxwellEvolution.lean)
prove the corresponding finite identities. Nineteen oscillator pairs on
this carrier are its finite mode count. They do not equal the two helicities
of each momentum mode of a continuum photon. The
[screen microphysics paper](../../paper/screen_microphysics_and_observer_synchronization.tex)
states these results in the propositions on local face curvature, local
Maxwell selection and temporal Maxwell evolution.

## Why the photon, W and Z differ

Choose a complex Higgs doublet with hypercharge \(Y_H=1/2\), canonical
kinetic term and a color- and charge-preserving reference phase represented
by \(H_0=(0,v/\sqrt2)\). This is a gauge-fixed representative; it is not
a claim that a gauge-charged expectation is a gauge-invariant observable.
The derivative is
\(D_\mu=\partial_\mu-i g W_\mu^a\sigma_a/2-i g'B_\mu/2\).
Differentiating \(\lvert D H_0\rvert^2\) gives

\[
 M^2_{W_1,W_2,W_3,B}=\frac{v^2}{4}
 \begin{pmatrix}
 g^2&0&0&0\\0&g^2&0&0\\0&0&g^2&-gg'\\0&0&-gg'&g'^2
 \end{pmatrix}.
\]

The neutral block has rank one for positive \(v,g,g'\). Its orthonormal
directions are

\[
 A=\frac{g'W_3+gB}{\sqrt{g^2+g'^2}},\qquad
 Z=\frac{gW_3-g'B}{\sqrt{g^2+g'^2}}.
\]

The electromagnetic direction annihilates the background:
\((T_3+Y_H)H_0=0\). Consequently the Higgs term contributes zero photon
mass. The orthogonal direction changes the background and has curvature
\(m_Z^2=(g^2+g'^2)v^2/4\). The charged combinations
\(W^\pm=(W_1\mp iW_2)/\sqrt2\) have
\(m_W^2=g^2v^2/4\). The source weak currents also obey
\([Q,T_\pm]=\pm T_\pm\), fixing their charged-current labels.
The coupling to the unbroken direction is
\(e=gg'/\sqrt{g^2+g'^2}\). None of these relations selects numerical
\(g,g',v\).

In intuitive terms, the photon is the neutral gauge direction that leaves
the selected background unchanged. The W and Z directions disturb that
background, so their quadratic equations contain a rest-frequency term.
The exact null direction concerns the quadratic action. Its physical
massless-photon interpretation also requires the electromagnetic Coulomb
regime and its state. A massive medium response is a different spectral
problem.

For the free Lorentz vacuum mode, \(E=c|p|\), giving group speed \(c\)
away from \(p=0\). A finite lattice ray can instead have
\(E_a(p)=c\sin(ap)/a\) for \(0<ap<\pi/2\), whose group speed is
\(c\cos(ap)\). It tends to \(c\) at fixed momentum as \(a\to0\).
The check therefore does not assign speed \(c\) to every carrier mode or
every wavelength. Absorber positions are detector outcomes; a universal
sharp photon position operator is not required by this explanation.

## What the Higgs excitation is

Take
\(V(H)=-\lambda v^2 H^\dagger H+\lambda(H^\dagger H)^2\),
with \(\lambda>0\), and write
\(H=(x_1+i x_2,v+h+i x_3)/\sqrt2\). The exact Hessian at zero
fluctuation is

\[
 D^2V=\operatorname{diag}(0,0,0,2\lambda v^2).
\]

The radial mode \(h\) changes the size of the scalar field and has
\(m_h^2=2\lambda v^2\). Its energy above the background is
\(\lambda v^2h^2+\lambda vh^3+\lambda h^4/4\), so the same action
also specifies scalar self-interactions. The three zero-curvature gauge
orbit directions furnish the longitudinal polarizations of W+, W− and Z
in the broken gauge theory. They do not add three observable scalar species.
The physical boson degree count is
\(2\times8+2+2\times3+3+1=28\).

The Higgs field gives masses through terms in the equations of motion.
This is not drag or energy loss during uniform motion. Massive modes have
nonzero rest energy; massless modes can carry energy at nonzero momentum.
Most of a hadron's mass is an interacting strong-field observable and is
not the sum of the quark Yukawa terms.

The [finite quantum-Higgs calculation](../particle_observer_outputs/HIGGS_RECORDS.md)
has a different domain. On its conserved 307-state block, the source vertex
couples the empty state to one bright pair state, retaining 305 dark states.
Its exact bright-level separation is
\(\Delta=\sqrt{d^2+4g_{\rm pair}^2}\), with
\(g_{\rm pair}^2=4\operatorname{Tr}(Y^\dagger Y)\). Complete sequential
detector records can read this separation on one declared clock. Here
\(g_{\rm pair}\) is a pair-vertex strength, not the electroweak gauge
coupling \(g\). Identifying this gap with the continuum radial Higgs mass
requires an embedding, state and scaling map for those same correlators
and instruments. No such identification is made by the symbolic Hessian.

## Why fermions have different masses

The Higgs background turns Yukawa interactions into chiral mass matrices:

\[
 M_e=vY_e/\sqrt2,\qquad M_u=vY_u/\sqrt2,\qquad M_d=vY_d/\sqrt2.
\]

Their nonnegative singular values are the tree action masses after
canonical normalization. Left and right frame changes diagonalize each
matrix. The mismatch between the up- and down-type left frames gives
\(V_{\rm CKM}=U_{uL}^\dagger U_{dL}\). A three-component family register
does not by itself fix these matrices or three numerical masses.

The [same-coefficient free-field calculation](../charged_lepton_response/RELATIVISTIC_MASS.md)
uses the source charged-Higgs tensor and checks the spin Clifford algebra.
For one diagonal mass channel its convention can be written
\(\alpha_j=\tau_3\otimes\sigma_j\),
\(\beta=\tau_1\otimes I\), with

\[
 H(p)=c\alpha\cdot p+\mu\beta,\qquad
 H(p)^2=(c^2|p|^2+\mu^2)I.
\]

The mass operator is off diagonal in the chiral blocks. At rest the
evolution is \(\cos(\mu t)I-i\sin(\mu t)\beta\) in units
\(\hbar=1\). A prepared chiral basis vector has opposite-chirality
probability \(\sin^2(\mu t)\); an energy eigenstate is a different
preparation. Thus chirality coupling creates a rest frequency without a
friction model. Physical charged masses require the interacting spectral
response, infrared treatment and specified energy units. The source's
finite sine regulator retains all eight valleys; a free Dirac limit does
not by itself supply an interacting chiral Standard Model continuum.

The gluon explanation uses actual color operators, rather than the count
eight alone. The eight independent native color currents have nonzero
commutators. Their declared spin lifts annihilate the Higgs doublet, so
the Higgs-generated color mass matrix vanishes. Their noncommutativity
enters Yang-Mills self-interactions. A zero gluon Lagrangian mass does not
mean that confinement supplies an isolated massless colored particle.
Quark masses require a stated scheme and scale; hadron and jet records
provide the comparison with experiment.

## Neutrinos and antiparticles

The minimal fifteen-state family contains a neutral weak-doublet component
but no neutrino mass interaction. A chosen Dirac completion adds neutral
singlets and \(M_\nu=vY_\nu/\sqrt2\). A Majorana effective completion
uses a symmetric coefficient with \(M_\nu=v^2\kappa/(2\Lambda)\).
Takagi factorization defines its nonnegative mass values. The weak flavor
channels are coherent combinations of mass eigenfields, with
\(U_{\rm PMNS}=U_{eL}^\dagger U_\nu\). Propagation phases and a complete
weak detector determine the oscillation records. Setting all neutrino
masses to zero does not reproduce oscillations. This atlas does not supply
absolute neutrino masses or choose Dirac versus Majorana behavior.

Conjugate fields carry opposite additive gauge charges. In the declared
CPT-compatible relativistic theory, a particle and its antiparticle have
the same mass. Charge conjugation alone is not a proof of that equality;
weak interactions need not respect charge conjugation. Electron, muon,
tau, each quark and W+ have distinct conjugate excitations. Photon, Z and
the neutral Higgs are self-conjugate. The gluon adjoint is a real multiplet,
and a color-basis label is not a distinct gauge-invariant species label.
Dirac neutrinos have distinct antiparticles; a Majorana mass field is
self-conjugate with different weak detection helicity channels.

The actual scalar vertex pairs spins with the invariant tensor
`epsilon=i sigma2`. The declared SU(2) spin generators therefore commute
with that Hamiltonian, its prepared neutral-sector ground is a singlet,
and its nonzero addition/removal probe vectors transform as doublets.
The [complete exterior attachment](../standard_model_source_map/FERMION_ATTACHMENT.md)
proves the CAR exchange signs and the current/creator commutator on every
occupation. This supplies finite transformation and statistics evidence.
The attachment to observer spatial rotations, unitary boosts and an
interacting relativistic spin-statistics theorem requires the corresponding
common field, positive-energy state and locality domain.

The atlas takes physical spin and statistics from the declared field action
and fermion operator algebra. Matching internal gauge weights alone does not
establish a physical spin-statistics theorem. Detector probabilities use
the same complete instrument contract described in
[READOUT_CONTRACT.md](../particle_observer_outputs/READOUT_CONTRACT.md).
Consensus compares the resulting records and their uncertainties; it does
not choose a desired particle identity or outcome.

## Reproduction

Run `python3 -m model.particle_atlas_horizontal.field_nature` to regenerate
the content and exact algebra receipt. Run
`python3 -m pytest -q model/particle_atlas_horizontal/test_field_nature.py`
to check the owner Hessians, neutral null direction, native color currents,
Dirac algebra, finite-lattice speed counterexample and all eighteen target
rows. The symbolic checks support the conditional field meanings, without
a claim of measured masses or a complete source-to-laboratory map.
