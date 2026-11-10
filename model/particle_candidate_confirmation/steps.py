"""Physical confirmation tests attached to the existing joint dependency map."""
from copy import deepcopy


COMMON = (
    ('operator', 'Source operator and records', [3, 6],
     'Bind the actual source operator and accessible current/field records, conserved charges, conjugate channel and complete detector, including failure and boundary compensation.'),
    ('response', 'Joint energy and momentum response', [4, 5, 7],
     'Calculate the original populated joint spectral response using the operational translation generator and the same clock.'),
    ('transformation', 'Physical transformation and statistics', [5, 6],
     'Attach physical Lorentz spin, chirality, exchange and current transformations; finite cube labels alone do not establish these.'),
    ('interaction', 'Interactions and distinct species', [5, 6, 7],
     'Test the current vertices, distinguish flavors or polarizations, and retain extra modes, compensation and competing assignments.'),
    ('mass', 'Controlled mass observable', [4, 7],
     'Use the appropriate threshold, resonance or renormalized mass definition and bound volume, cutoff, capacity and numerical errors.'),
    ('comparison', 'Independent physical comparison', [8],
     'Attach one common unit convention, preserve calibration ancestry and compare complete records with experimental uncertainties.'),
)

FAMILIES = {
    'charged_leptons': {
        'targets': ['electron', 'muon', 'tau'],
        'source_evidence': 'A charged source spinor and its scalar dressing have complete compensated pair-energy records. Flavor labels are unassigned.',
        'response': 'Resolve a charged energy-momentum threshold while retaining the opposite endpoint and flux; their total energy cannot be divided by two to define a lepton mass.',
        'interaction': 'Identify electromagnetic and weak current responses and three distinct source flavor channels. Test conjugate records and allowed decays.',
        'mass': 'Specify the charged infrared threshold or resonance convention; distinguish stable electron and unstable muon/tau readouts.',
    },
    'active_neutrinos': {
        'targets': ['electron_neutrino', 'muon_neutrino', 'tau_neutrino'],
        'source_evidence': 'The internal weak-doublet neutral component exists. The tested neutral singlet trion has different unbroken weak-current transformations.',
        'response': 'Construct a populated neutral weak-current channel and its energy-momentum response; a singlet composite name does not supply it.',
        'interaction': 'Derive the production/detection flavor dictionary, mass eigenchannels and oscillation probabilities. Specify the neutrino operator and test its lepton-number consequences.',
        'mass': 'Extract mass eigenvalues or controlled splittings and the mixing matrix; flavor names are not mass eigenstates.',
    },
    'up_type_quarks': {
        'targets': ['up', 'charm', 'top'],
        'source_evidence': 'Exact up-type internal weights exist, with one local copy per native representation.',
        'response': 'Construct color-dressed, gauge-invariant production/current records in a specified effective dictionary.',
        'interaction': 'Test non-Abelian color and weak vertices, three distinct up-type flavors and charged-current mixing with the down-type channels.',
        'mass': 'Calculate renormalized quark masses in a declared scheme and scale with controlled physical matching; an isolated colored rest-energy pole is not required.',
    },
    'down_type_quarks': {
        'targets': ['down', 'strange', 'bottom'],
        'source_evidence': 'Exact down-type internal weights exist, with one local copy per native representation.',
        'response': 'Construct color-dressed, gauge-invariant production/current records in a specified effective dictionary.',
        'interaction': 'Test non-Abelian color and weak vertices, three distinct down-type flavors and charged-current mixing with the up-type channels.',
        'mass': 'Calculate renormalized quark masses in a declared scheme and scale with controlled physical matching; finite neutral sector bottoms are different observables.',
    },
    'photon': {
        'targets': ['photon'],
        'source_evidence': 'Finite polar and axial loop channels have controlled source transitions in a rare empty-matter sector.',
        'response': 'Identify an unbroken electromagnetic current and populated transverse response with a massless dispersion law that survives refinement.',
        'interaction': 'Identify the two physical helicities and electromagnetic vertices, with longitudinal gauge redundancy and charged infrared records treated explicitly.',
        'mass': 'Bound the invariant mass at zero in the same energy-momentum convention; a finite loop energy gap is not a photon mass.',
    },
    'gluon': {
        'targets': ['gluon'],
        'source_evidence': 'Eight color-current directions exist; their undressed local adjoint currents vanish under physical Gauss compression.',
        'response': 'Construct gauge-invariant color-field correlation and scattering records with an explicit field/gauge dictionary.',
        'interaction': 'Derive non-Abelian self-interactions and their coupling, plus the color responses of matter. Eight algebra directions are not eight free observed particles.',
        'mass': 'Control the massless short-distance gauge-field description and infrared gauge-invariant spectrum without assigning a free gluon pole mass.',
    },
    'W_plus': {
        'targets': ['W_plus'],
        'source_evidence': 'The native weak raising direction is defined; a physical massive charged vector is unassigned.',
        'response': 'Construct a charged vector resonance with its compensating charge and three physical massive-vector polarizations.',
        'interaction': 'Test charged-current couplings, the conjugate W channel, Higgs-phase relations and inclusive decay records.',
        'mass': 'Determine a charged-vector complex-pole or specified experimental mass/width convention with radiative and finite-detector errors.',
    },
    'W_minus': {
        'targets': ['W_minus'],
        'source_evidence': 'The native weak lowering direction is defined; a physical massive charged vector is unassigned.',
        'response': 'Construct a charged vector resonance with its compensating charge and three physical massive-vector polarizations.',
        'interaction': 'Test charged-current couplings, the conjugate W channel, Higgs-phase relations and inclusive decay records.',
        'mass': 'Determine a charged-vector complex-pole or specified experimental mass/width convention with radiative and finite-detector errors.',
    },
    'Z': {
        'targets': ['Z'],
        'source_evidence': 'A neutral weak/hypercharge span exists, without a physical mixing-angle or Z assignment.',
        'response': 'Construct the neutral massive vector resonance and separate it from the transverse massless electromagnetic channel.',
        'interaction': 'Calculate weak/hypercharge mixing and neutral-current vertices on the same Higgs-phase construction used for the charged vectors.',
        'mass': 'Determine a neutral-vector complex-pole or specified experimental mass/width convention with controlled corrections.',
    },
    'Higgs': {
        'targets': ['Higgs'],
        'source_evidence': 'The Higgs-conjugate composite supplies a neutral radial scalar with a capacity-uniform total detector bound and a separate nonzero original dynamical response. Total success includes elastic events.',
        'response': 'Measure a gauge-invariant scalar response and establish its populated branch in a demonstrated Higgs phase.',
        'interaction': 'Construct consistent scalar, vector and fermion vertices. The native-number obstruction to one sharp scalar serving every bare Yukawa channel requires an explicit resolution.',
        'mass': 'Extract a scalar resonance and width around the selected physical phase; neither a scalar operator nor a zero one-point function establishes or excludes that phase.',
    },
}


def candidate_steps(target):
    family = next(name for name, row in FAMILIES.items() if target in row['targets'])
    specifics = FAMILIES[family]
    return [
        {'id': key, 'title': title, 'joint_nodes': list(nodes),
         'requirement': specifics.get(key, requirement),
         'physical_confirmation_status': 'unproved',
         **({'finite_source_evidence': specifics['source_evidence']} if key == 'operator' else {})}
        for key, title, nodes, requirement in COMMON
    ]


def plan():
    return {
        'scope': 'Candidate tests within the existing joint dependency map; no additional physical gates.',
        'families': deepcopy(FAMILIES),
        'common_sequence': [{'id': key, 'title': title, 'joint_nodes': list(nodes),
                             'requirement': requirement}
                            for key, title, nodes, requirement in COMMON],
        'neutral_composite_requirement': (
            'Use the six signed conserved-increment labels, calculate populated joint energy/momentum responses, '
            'attach physical currents and spin, and test interactions before any neutrino or baryon assignment. '
            'Retain the explicit multiple-excitation alias of the seven-state counter.'),
    }


def render(value):
    data = value['confirmation_plan']
    lines = ['# Confirming candidate particles', '', data['scope'], '',
             'The six tests use one source, state, reference and clock. A finite source ',
             'result supports a test without establishing its physical interpretation.', '']
    for i, step in enumerate(data['common_sequence'], 1):
        lines += [f"{i}. **{step['title']}.** {step['requirement']}"]
    for family, row in data['families'].items():
        lines += ['', '## ' + ', '.join(row['targets']), '', row['source_evidence'], '',
                  '| Test | Required calculation |', '| --- | --- |']
        for step in candidate_steps(row['targets'][0]):
            if step['id'] in ('response', 'interaction', 'mass'):
                lines.append('| ' + step['title'] + ' | ' + step['requirement'] + ' |')
    lines += ['', '## Neutral composite hypotheses', '', data['neutral_composite_requirement'], '',
              '[Source calculations and their boundaries](RESULTS.md). ',
              '[Common field-and-particle dependency map](../particle_joint_reconstruction/DEPENDENCIES.md).', '']
    return '\n'.join(line.rstrip() for line in lines)
