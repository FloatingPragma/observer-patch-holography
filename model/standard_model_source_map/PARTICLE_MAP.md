# Standard Model source correspondence

**Physically identified particle sources: none.**

The table compares exact internal source representations with Standard Model target channels.
It assigns no observed mass, topological defect or physical particle identity.
The configured field has three family coordinates. Their measured-species assignment is unproved.
The source action coefficients and four complete duty policies are exported in [SOURCE_MASSES.md](SOURCE_MASSES.md).
Neutrino rows are flavor channels, not identified mass eigenstates. Gluons retain eight color components.

| Target | Computed source interface | Physical source | Predicted mass |
| --- | --- | --- | --- |
| electron | charged_lepton_left: masks 23; charged_lepton_right_conjugate: masks 24 | Unidentified | Underived |
| muon | charged_lepton_left: masks 23; charged_lepton_right_conjugate: masks 24 | Unidentified | Underived |
| tau | charged_lepton_left: masks 23; charged_lepton_right_conjugate: masks 24 | Unidentified | Underived |
| up | up_left: masks 9,10,12; up_right_conjugate: masks 3,5,6 | Unidentified | Underived |
| charm | up_left: masks 9,10,12; up_right_conjugate: masks 3,5,6 | Unidentified | Underived |
| top | up_left: masks 9,10,12; up_right_conjugate: masks 3,5,6 | Unidentified | Underived |
| down | down_left: masks 17,18,20; down_right_conjugate: masks 27,29,30 | Unidentified | Underived |
| strange | down_left: masks 17,18,20; down_right_conjugate: masks 27,29,30 | Unidentified | Underived |
| bottom | down_left: masks 17,18,20; down_right_conjugate: masks 27,29,30 | Unidentified | Underived |
| electron_neutrino | neutrino_left: masks 15 | Unidentified | Underived |
| muon_neutrino | neutrino_left: masks 15 | Unidentified | Underived |
| tau_neutrino | neutrino_left: masks 15 | Unidentified | Underived |
| photon | unbroken Q direction conditional on the neutral Higgs vacuum | Unidentified | Underived |
| gluon | eight color-current directions; one color-adjoint field target | Unidentified | Underived |
| W_plus | weak raising current, [Q,T+]=T+ | Unidentified | Underived |
| W_minus | weak lowering current, [Q,T-]=-T- | Unidentified | Underived |
| Z | neutral span(T3,Y), with source-selected mixing required | Unidentified | Underived |
| Higgs | radial mode conditional on a selected neutral doublet vacuum | Unidentified | Underived |

The right-handed comparison components are charge-conjugated left-Weyl internal weights.
The even internal register has an exact current intertwiner into a separately configured odd CAR field.
Free Weyl propagation and finite dynamical gauge identities are checked; the interacting chiral continuum is unproved.
A conserved symmetry label, a topological invariant, an internal multiplet and a particle pole are distinct objects.

Target conventions: [Particle Data Group electroweak review](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf).
Source computation and admission requirements: [proof](PROOF.md), [machine export](particle_map.json).
