# Conditional Standard Model spectrum

A complete elementary-field benchmark follows from a supplied renormalized
Standard Model action, with a specified neutrino extension. Its algebraic
charges, quadratic masses and field counts are checked by
[`conditional_standard_model.py`](conditional_standard_model.py). This is a
conditional phenomenological completion. It neither derives its parameters
from OPH nor proves an interacting continuum construction. The atlas's
[measured reference values](reference_values.json) are external comparisons;
they are not all inverted simultaneously as tree-level Lagrangian parameters.

## Source attachment and physical assumptions

The [finite fermion attachment](../standard_model_source_map/FERMION_ATTACHMENT.md)
provides the same twelve currents on fifteen internal components, their full
fermion Fock action, three configured families and the exact kinetic/Higgs
intertwiner. These are positive source results. The conditional completion
requires the following additional choices and estimates:

1. A chiral gauge regulator and physical sector with precisely the stated
   light representations, controlled mirror removal and gauge-current Ward
   identities. Local anomaly cancellation alone does not construct its
   measure, dynamics or physical state. The retained naive sine regulator
   has eight valleys; its finite current map does not remove them.
2. A common state, physical clock and low-energy limit with positive,
   Lorentz-covariant canonical kinetic normalization, local unitary dynamics,
   and the assumed color-preserving, electromagnetism-preserving Higgs phase.
   The [finite quantum scalar](INTERACTING_CHIRAL_JOIN.md) supplies invariant
   vertices and Gauss identities. Its three-dimensional carrier is not a
   continuum Higgs doublet or a selected Higgs vacuum.
3. Renormalized gauge, scalar and Yukawa parameters at a stated scale and in
   a stated scheme, including mixing matrices and any allowed CP phases.
   A matching map relates them to specified physical observables with
   radiative, threshold and infrared errors. The parameters may be supplied
   or measured; neither permission makes them source-selected predictions.

Under these assumptions the following quadratic derivations are exact.
Using them as physical masses requires the matching in item 3. No theorem
here asserts nonperturbative existence of four-dimensional interacting
Standard Model dynamics.

## Representations, charges and statistics

Use hypercharge `Q=T3+Y`, with coupling `g_prime` to Y. This is half the
hypercharge convention `Q=T3+Y_PDG/2` in the
[PDG electroweak review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-standard-model.pdf).
It is not the grand-unified rescaling of that coupling. One family consists
of left-handed Weyl fields

| Field | SU(3) | SU(2) | Y | Electric weights |
|---|---|---|---|---|
| Q | 3 | 2 | 1/6 | 2/3, -1/3 |
| u^c | conjugate 3 | 1 | -2/3 | -2/3 |
| d^c | conjugate 3 | 1 | 1/3 | 1/3 |
| L | 1 | 2 | -1/2 | 0, -1 |
| e^c | 1 | 1 | 1 | 1 |

There are fifteen Weyl components per family, forty-five for three families,
and ninety spatial-spinor CAR coordinates. The two spinor rows are not two
families. In this convention a conjugate singlet supplies the right-handed
component of a charged Dirac field after conjugation. Nambu doubling is
particle-hole bookkeeping, not a second set of independently counted species.
The code checks the linear and cubic hypercharge sums, mixed color/weak
anomalies, the color cubic anomaly and the even number of weak doublets.
These checks are necessary representation identities, not a chiral-regulator
existence proof.

Fermions carry spin one-half and fermionic statistics; the gauge fields carry
spin one and bosonic statistics; the radial Higgs field carries spin zero.
The source CAR proves the finite fermionic statistics. Physical spin and the
spin-statistics connection require the assumed Lorentz/locality structure.
Antiparticles carry conjugate internal representations and opposite electric
charges. Equal physical particle/antiparticle masses additionally use the
local Lorentz-covariant quantum theory with a positive-energy unitary vacuum
representation and the usual spin-statistics/CPT hypotheses of the conditional
physical completion. Charge conjugacy by itself proves no such mass equality; the
weak interactions need not be charge-conjugation invariant.

## Supplied action and its exact quadratic spectrum

With canonical fields, the renormalized action has gauge terms
`-G^2/4-W^2/4-B^2/4`, the covariant fermion and Higgs kinetic terms, and

    V(H) = -mu_H^2 H*H + lambda_H (H*H)^2,
    L_Y = -bar(Q_L) tilde(H) Y_u u_R
          -bar(Q_L) H Y_d d_R -bar(L_L) H Y_e e_R + h.c.  (1)

Here `lambda_H>0`, `mu_H^2>0`, and the tree expansion uses
`H0=(0,v/sqrt(2))`, `v^2=mu_H^2/lambda_H`. The gauge-fixed representative H0
describes the chosen phase; it does not assert a nonzero gauge-invariant
expectation of a gauge-charged local operator. The parameters and the loop
counterterms have one specified renormalization convention. The
[PDG Higgs review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-higgs-boson.pdf)
states the one-doublet action and its coupling conventions.
A supplied QCD theta term and complex Yukawa phases are allowed parameters;
no source mechanism for strong-CP suppression is inferred. The theta term
does not change these tree-level mass Hessians. Nonzero charged masses
require nonsingular charged Yukawa matrices.

Differentiating `(D H0)*D H0` with respect to `(W1,W2,W3,B)` gives

    M_vector^2 = (v^2/4) diag(g^2,g^2,
                              [[g^2,-g g_prime],
                               [-g g_prime,g_prime^2]]). (2)

The neutral photon vector is proportional to `(g_prime,g)` in the `(W3,B)`
plane; the orthogonal Z vector is proportional to `(g,-g_prime)`. Thus

    m_W^2 = g^2 v^2/4,    m_Z^2 = (g^2+g_prime^2)v^2/4,
    m_gamma^2 = 0,        e = g g_prime/sqrt(g^2+g_prime^2),
    m_h^2 = 2 lambda_H v^2.                              (3)

The scalar Hessian has three zero Goldstone directions and one radial
eigenvalue in (3). The three Goldstone directions supply the longitudinal
massive-vector polarizations; they are not extra elementary particles.
Color is unbroken, so the eight gluon Lagrangian masses vanish.

For each `f=u,d,e`, biunitary diagonalization of
`M_f=v Y_f/sqrt(2)` gives nonnegative quadratic masses

    m_f,i = sigma_i(M_f).                                (4)

Family labels specify which diagonalized branches are matched to observed
flavors; a sorting convention does not make that identification. The left
rotations give `V_CKM=U_uL* U_dL`. With positive noncanonical Lorentz-covariant
kinetic matrices, the matrix to diagonalize is instead
`Z_L^(-1/2)(M_0+E) Z_R^(-1/2)`, under the assumptions of
[`MASS_RESPONSE_BRIDGE.md`](../charged_lepton_response/MASS_RESPONSE_BRIDGE.md).
The source's coefficient `v_H` agrees with `v/sqrt(2)` only after matching
the Higgs and kinetic normalization. The finite free theorem does not
compute E or those Z matrices in an interacting physical state.

## Neutrino extension and the complete elementary channels

The renormalizable forty-five-Weyl Standard Model has massless neutrinos.
It cannot reproduce the observed oscillations. Two explicit alternatives
for the conditional benchmark are:

- Three right-handed gauge singlets, with Dirac matrix
  `M_nu=v Y_nu/sqrt(2)`. This adds three Weyl components, giving forty-eight
  Weyl components and ninety-six spinor CAR coordinates.
- A dimension-five Weinberg operator with symmetric coefficient kappa. Fix
  its normalization by `M_nu=v^2 kappa/(2 Lambda)`. Takagi diagonalization
  gives three nonnegative Majorana masses without adding light sterile
  fields to the effective field content.

Both choices add data absent from the fifteen-component massless-neutrino
construction. In either case
`nu_alpha,L=sum_i (U_PMNS)_(alpha i) nu_i,L`,
`U_PMNS=U_eL* U_nu`, and the mass fields are `nu_1,nu_2,nu_3`.
Their absolute masses are symbolic here. Oscillations constrain squared-mass
differences and mixing, not the absolute mass scale or Dirac/Majorana choice.
The [PDG neutrino review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-neutrino-mixing.pdf)
gives these distinctions. No mass is assigned separately to each flavor
label in the atlas.

| Channels | Antiparticles | Spin and charge | Conditional mass/measurement meaning |
|---|---|---|---|
| e, mu, tau | positron, anti-muon, anti-tau | 1/2; -1 versus +1 | Three charged-lepton values in (4); electron charged infrared threshold, muon/tau resonance parameters |
| u,c,t | anti-u,anti-c,anti-t | 1/2; +2/3 versus -2/3; three colors | Up-sector values in (4), matched to an explicit quark scheme |
| d,s,b | anti-d,anti-s,anti-b | 1/2; -1/3 versus +1/3; three colors | Down-sector values in (4), matched to an explicit quark scheme |
| nu_1,nu_2,nu_3 | distinct Dirac conjugates or Majorana self-conjugate fields | 1/2; zero charge | Symbolic masses; nu_e,nu_mu,nu_tau are flavor superpositions |
| photon | self-conjugate | 1; zero charge; two helicities | Massless unbroken electromagnetic gauge field |
| eight gluons | real adjoint multiplet | 1; zero electric charge | Massless gauge coefficients, not eight isolated colored asymptotic particles |
| W+,W- | each other | 1; charges +1,-1 | Equal quadratic masses in (3), three polarizations each; unstable resonances |
| Z | self-conjugate | 1; zero charge | Neutral massive vector resonance |
| h | self-conjugate | 0; zero charge | Neutral scalar resonance |

This covers the atlas's eighteen channel labels, with W signs separate,
three neutrino flavor cards linked to their mass basis, and quark color
multiplicity retained. Antiparticles are included without duplicating the
field count. A graviton is outside this Standard Model benchmark. Hadrons
are composite channels, not additional elementary rows.

The quadratic bosonic degree count is `2(8+3+1)+4=28` before symmetry
breaking and `16+2+6+3+1=28` after it. Charged Dirac fields contribute
`4(3+6*3)=84` on-shell degrees, and three Weyl or Majorana neutrinos contribute
six, or three Dirac neutrinos twelve. These are free-field counts, not a
claim that every colored or unstable component exists as an asymptotic
one-particle state.

## Matching physical observables

In units with the matched limiting speed and hbar equal to one, free massive
branches satisfy `E(p)^2=|p|^2+m^2`. A stable massive asymptotic representation
has the same dispersion. The electron's infrared charged sector requires
dressing/inclusive observables and need not contain an isolated Fock pole.
Muon, tau, W, Z and Higgs masses require resonance/width conventions rather
than an exact stable-particle subspace. The
[electroweak mass and renormalization review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-standard-model.pdf)
distinguishes renormalized parameters, complex-pole observables and input
schemes. No exact muon or tau lifetime is computed by the quadratic action.

Quarks and gluons are confined colored fields. Running quark masses depend
on scale and scheme, and a quark pole-mass convention has a nonperturbative
ambiguity of order the strong scale. A top mass extracted with one
reconstruction convention cannot silently be equated to another convention.
These distinctions are stated in the
[PDG QCD review](https://pdg.lbl.gov/2026/reviews/rpp2026-rev-qcd.pdf).
The atlas preserves the reference source's scheme and uncertainty labels.

To reproduce numerical measurements conditionally, choose an independent
renormalized input set, specify the electroweak/QCD matching calculation and
its uncertainty, and compare the resulting observables in their stated
schemes. There is no assertion that every measured mass can be inserted
into (3)-(4) simultaneously while retaining an exact interacting prediction.
The calibrated charged-family example and the four P-law branches are
separate conditional source calculations. Neither exact Koide balance nor
its supplied deformation automatically survives this matching.

## Concrete source-to-benchmark obligations

| Required link | Finite result available | Calculation needed for the physical implication |
|---|---|---|
| Same state and clock | Complete finite records and reference-safe feedback | Select the joint state, compute its time-dependent correlators, control the growing-source limit and physical clock conversion |
| Kinetic and mass response | Exact free symbol and conditional whitening bounds | Extract temporal/spatial inverse two-point kernels and self-energy from those same correlators; bound residues, dispersion errors and matching corrections on the stated energy range |
| Chiral gauge dynamics | Full-CAR currents, finite Gauss identities, quantum-link controls | Construct a gauge-covariant chiral regulator with canonical CAR, current-preserving mirror control and a suitable physical state; a sign-function identity alone is insufficient |
| Higgs phase | Charged finite scalar vertices and common-Hamiltonian feedback | Establish the selected phase and scalar correlations, normalize its effective doublet and determine its potential and Yukawa coefficients |
| Defect-to-particle dictionary | Exact compact localized CAR/current map | Construct defect creation and motion with the same Gauss/charge operators, then compare its spectral measures; the fixed Pauli texture is not a scalar rest mass or a momentum shell |
| Neutrinos | Neutral massless components in M15 | Select and implement the Dirac-singlet or Weinberg interaction, with a mass matrix and flavor-current map |
| Measured masses and widths | Complete spectral scores and conditional error budgets | Specify charged infrared/resonance observables, renormalization scheme and inclusive records, then propagate actual source errors through the observable matching |

The compact-texture limitation follows from the actual kinetic matrices and
interior translation commutator in
[`DEFECT_FIELD_JOIN.md`](DEFECT_FIELD_JOIN.md#which-identification-is-excluded).
The [quantum chiral join](INTERACTING_CHIRAL_JOIN.md) separates the finite
Gauss construction from its continuum hypotheses. These restrictions define
the scope of their sources; they do not prohibit other explicitly derived
regulators, dynamical defects or physical completions.
