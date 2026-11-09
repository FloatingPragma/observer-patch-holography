# Candidate particle interfaces and the defect attachment

[`candidate_map.py`](candidate_map.py) connects all eighteen target channels
in the [internal source dictionary](../standard_model_source_map/PROOF.md)
to their constructive source interfaces. It retains two distinct kinds of
charged-family action: the calibrated coherent source and the conditional
pixel-dependent MCPR source. The result is a candidate map. Its physical
particle map is empty, and every physical mass, spin, lifetime, width,
topological invariant and defect-to-particle map is unassigned.

This boundary does not erase the finite constructions. The fifteen internal
weights, their ninety-mode CAR lift and same-current Higgs tensors are exact.
The coherent charged action has a complete full-CAR realization and a free
Dirac limit on a specified reference. The MCPR coefficient has a complete
native-memory realization. These statements concern defined operators,
records and source states. Assigning a measured particle to one of their
sectors needs the additional correspondence stated below.

## Meaning of a candidate row

The pure function
`summarize(parent_map, calibrated_summary, p_dependence_certificate, optional_defectinventory)`
checks the complete eighteen-target and fifteen-weight dictionaries before
projecting any data. The component match preserves every actual source mask,
color weight, T3, hypercharge and electric charge. Right-handed comparison
components use the charge-conjugated left-Weyl convention. A source charge
or representation label is not a topological invariant or a measured
spacetime spin.

The three charged names are an explicit hypothesis: ascending action mode
indices 0,1,2 are compared with electron, muon and tau. Their names do not
enter the source matrix or interval calculation. This sorting does not
supply a family-mixing, decay or defect-creation theorem. Neutrino rows are
flavor channels and carry no proposed sorted mass index. Quark rows retain
internal representation candidates without choosing one of the supplied
family actions. The gluon target denotes an eight-component color-adjoint
field, while W+ and W- are separate target channels.

| Target | Constructive source interface | Exact missing physical attachment |
|---|---|---|
| Electron | Ascending charged mode 0; charged left and conjugate-right CAR components | Interacting charged threshold, common clock and observable match; a defect map is unproved |
| Muon | Ascending charged mode 1 in the same chosen action | Interacting unstable-particle resonance, decay/width and observable match; a defect map is unproved |
| Tau | Ascending charged mode 2 in the same chosen action | Interacting unstable-particle resonance, decay/width and observable match; a defect map is unproved |
| Up | Up component of the color-triplet weak doublet and its conjugate-right component | A chosen family action, physical mixing and a gauge-invariant short-distance or hadronic mass convention |
| Charm | The same up-type internal interface, family coordinate unspecified | A family eigenmode and the interacting quark observable |
| Top | The same up-type internal interface, family coordinate unspecified | A family eigenmode, unstable-particle matching and the interacting quark observable |
| Down | Down component of the color-triplet weak doublet and its conjugate-right component | A chosen family action, physical mixing and a gauge-invariant short-distance or hadronic mass convention |
| Strange | The same down-type internal interface, family coordinate unspecified | A family eigenmode and the interacting quark observable |
| Bottom | The same down-type internal interface, family coordinate unspecified | A family eigenmode and the interacting quark observable |
| Electron-neutrino channel | Neutral component of the left weak doublet | Neutrino mass mechanism, flavor-to-mass mixing and oscillation response |
| Muon-neutrino channel | The same neutral internal component with an unspecified family coordinate | Neutrino mass mechanism and oscillation/mixing identification |
| Tau-neutrino channel | The same neutral internal component with an unspecified family coordinate | Neutrino mass mechanism and oscillation/mixing identification |
| Photon | Q=T3+Y direction, conditional on a neutral Higgs vacuum | Selected vacuum, physical transverse propagation, masslessness and observable correspondence |
| Gluon | Eight independent exact color-current directions | Dynamical gauge field and gauge-invariant interacting observables; a colored free pole is not assumed |
| W+ | Raising current T+ with [Q,T+]=T+ | Selected broken-phase gauge action, physical polarizations and charged resonance |
| W- | Lowering current T- with [Q,T-]=-T- | The same broken-phase construction and conjugate charged resonance |
| Z | Neutral current span(T3,Y) | Source-selected mixing, vacuum and neutral resonance |
| Higgs | Radial response conditional on a supplied neutral doublet direction | Selected vacuum, dynamical scalar kinetic response and resonance |

Every row also requires the claimed source/observable map, common physical
clock and controlled approximation errors. The map gives no physical spin
or lifetime by copying a target name. In particular, finite adjoint current
directions are not a proof of propagating spin-one particles. The minimal
M15 representation has no right-neutrino mode. Its neutral rows do not
specify a Dirac or Majorana mass mechanism.

## Two different charged-family actions

The calibrated coherent branch uses the specified positive circulant
`C(z)` and `Y=C(z)^2/3`. Its phase is chosen from the muon/electron ratio
alone by the [calibrated comparison](CALIBRATED_SANITY.md). The candidate
map copies exact rational intervals for the three eigenvalues of Y and
their ratios. Tau remains outside that phase choice. These are dimensionless
quadratic-action coefficients. The free energy uses `|v_H| sigma_j(Y)`;
this construction supplies no interacting physical energy unit.

The source objects for this branch are the three spectral subspaces of Y
and the corresponding same-current odd CAR creators. Its
[full-CAR theorem](CAR_KOIDE.md) and [free field theorem](RELATIVISTIC_MASS.md)
use a declared apparatus, filled-negative-energy Gaussian reference and
clock. All eight regulator valleys are retained. The three families are
not three lower thresholds of the entire electric-charge sector. Physical
interactions can mix channels or make higher excitations unstable.

The P-MCPR branch instead uses the unconditional C-duty coefficient
`Z_C`, with a selected gain g. Its constituent native projectors are
P2, P3 and P10. These projectors are not mutually orthogonal mass projectors
and are not assigned separately to electron, muon and tau. The mass-mode
objects are the spectral subspaces of the whole matrix `g Z_C`. The
[exact inverse memory theorem](../mcpr_native_memory_bridge/PROOF.md)
realizes the three MCPR coefficients with two admissible rate branches,
all disabled and failure records retained, and a correlated weak-pulse
error bound. A gauge-trivial family unitary preserves the current algebra
when diagonalizing this coefficient.

The [capacity dependence](P_DEPENDENCE.md) supplies four separately labelled
coordinate variants: the forward source-audit endpoint, the measured-alpha
endpoint, the supersymmetric-running proxy closure, and the minimal-census
color-scale closure. Each retains its source ancestry and both inverse
memory branches. Their phase, amplitude and squared-root shape intervals
are copied separately from the coherent action. For shape squares s_j,

    t_j=s_j [exp(3 kappa)/(2*6^14 product(s))]^(1/3),
    kappa=(u/50)/(1+u/31+u^2/310).

Thus s_j alone is not the absolute coefficient m_j/v. Its ratios equal the
action ratios, and the common normalization is explicitly stated. The
MCPR amplitude is `sqrt(2) exp(-chi)`, so its positive-root quotient is
strictly below 2/3 for the admitted positive u. It is not the exact balanced
coherent action. The MCPR architecture has target-informed provenance;
its source evaluation does not fit charged masses at runtime. A unique
conditional closure root does not select that architecture from primitive
repair. The [phase audit](SOURCE_PHASE_SELECTION.md) gives these distinctions.

No quark, neutrino or gauge-boson action mass is filled from these charged
calculations. A chosen action in another sector would require its own
operator and interval attachment. Numerical resemblance is not such an
attachment.

## A separate defect inventory

The optional input is the existing eleven-entry
[protected-object catalogue](../protected_defect_catalog/PROOF.md).
Its source domain, protected object and topology descriptions are retained
without pairing them to any of the eighteen target names. Its exclusions,
spin sectors, pair/exchange charges, centered-CAR responses and kinetic
nodes belong to different specified source problems. Their count and their
conserved labels do not establish a bijection with particle species.

Every candidate row therefore has null fields for the identified defect,
topological invariant, defect-creation map and physical spectral-sector map.
A positive defect identification would have to supply an explicit map from
one stated source defect and its admissible creation instrument into the
same particle or resonance sector, preserving the relevant observables,
charges, preparation and clock, with the approximation errors charged.
The construction here supplies no such map and does not claim an exhaustive
classification of all possible source defects.

## Executable scope

`certificate()` evaluates only the bounded calibrated interval projection,
the bounded capacity readout and the candidate factory. It does not repeat
the full native apparatus, full-Fock experiment or implicit pixel-root
producer. It binds the standard-model map, defect catalogue, MCPR memory
receipt, pixel closure receipt and calibration reference. No charged-owner
results file is read, so the exporter has no circular dependence on itself.

The tests use the actual fifteen-weight parent, complete defect inventory
and bounded four-branch evaluation. They reject missing or duplicate target
channels, changed component weights, physical admission, a fitted tau,
calibration relabelled as derivation, collapsed MCPR branches, false exact
balance and forced defect-species links. All returned maps are JSON-safe
and detached from their supplied input objects.
