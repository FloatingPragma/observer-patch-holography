"""Export a source correspondence without assigning unproved particles."""
from pathlib import Path
import hashlib
import json
from functools import lru_cache

from tools import strict_json
from tools.build_completion_program import integration_source_paths
from . import (algebra, fermion_attachment, kinetic_attachment, gauge_attachment,
              feedback_limit, source_mass_readout, koide_source)

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'particle_map.json'
MARKDOWN = HERE/'PARTICLE_MAP.md'
MASS_MARKDOWN = HERE/'SOURCE_MASSES.md'
PARENTS = ('model/common_source/native_matter_attachment_results.json',
           'model/particle_mass_readout_gate/table.json',
           'model/source_strong_coupling/results.json',
           'model/native_flavor_orbit/results.json',
           'model/native_quantum_gauss/results.json')
REFERENCE = 'https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    for name in PARENTS:
        parent = strict_json.load(ROOT/name)
        pins = parent.get('source_bindings', parent.get('source_sha256'))
        if not isinstance(pins, dict) or not pins:
            raise ValueError('unbound source parent')
        for path, expected in pins.items():
            target = (ROOT/path).resolve()
            if not target.is_relative_to(ROOT) or hashlib.sha256(target.read_bytes()).hexdigest() != expected:
                raise ValueError('stale source parent: '+path)
    paths = integration_source_paths(ROOT, {'domain_and_proof': 'model/standard_model_source_map/PROOF.md',
                                          'source_receipts': list(PARENTS),
                                          'independent_checks': [
                                              'model/standard_model_source_map/'+name for name in
                                              ('README.md', 'FERMION_ATTACHMENT.md', 'KINETIC_ATTACHMENT.md',
                                               'GAUGE_ATTACHMENT.md', 'FEEDBACK_LIMIT.md', 'MASS_READOUT.md',
                                               'KOIDE_ATTACHMENT.md', 'KOIDE_SOURCE.md')]})
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


@lru_cache(maxsize=1)
def _science_json(source_key):
    """Cache immutable bytes by the entire authenticated source content."""
    return canonical({
        'fermion_operator_attachment': fermion_attachment.certificate(),
        'kinetic_and_mirror_attachment': kinetic_attachment.certificate(),
        'dynamical_gauge_attachment': gauge_attachment.certificate(),
        'complete_feedback_limit': feedback_limit.certificate(),
        'source_action_masses': source_mass_readout.certificate(),
        'coherent_Koide_source': koide_source.certificate()})


def build():
    pins = source_bindings()
    admitted = strict_json.load(ROOT/PARENTS[1])
    color = strict_json.load(ROOT/PARENTS[2])
    if (admitted['physical_particle_mass_table'] != [] or
            color['scope']['physical_alpha_s_prediction'] is not None):
        raise ValueError('physical admission changed; the export needs a reviewed identification map')
    components = algebra.match_components()
    gauge = algebra.gauge_and_higgs()
    science = json.loads(_science_json(tuple(pins.items())))
    common = ['same-source interacting relativistic field and state',
              'source creation and complete detector map into a particle or resonance sector',
              'physical spin/statistics and current identification',
              'common clock, energy convention and controlled limiting errors',
              'independent experimental correspondence']
    rows = []
    groups = ((('electron', 'muon', 'tau'), ('charged_lepton_left', 'charged_lepton_right_conjugate')),
              (('up', 'charm', 'top'), ('up_left', 'up_right_conjugate')),
              (('down', 'strange', 'bottom'), ('down_left', 'down_right_conjugate')),
              (('electron_neutrino', 'muon_neutrino', 'tau_neutrino'), ('neutrino_left',)))
    for names, keys in groups:
        for name in names:
            rows.append({'target': name, 'target_kind': 'neutrino_flavor_channel' if 'neutrino' in name else 'fermion_species',
                         'source_interface': 'M15 internal weights attached to the configured 90-mode CAR field',
                         'internal_components': {k: components[k] for k in keys},
                         'source_family_index': None,
                         'family_or_mixing_identification': 'three family coordinates are configured; their identification with nature is unproved',
                         'missing': common+['physical family, mass/mixing and chiral regulator identification'],
                         'antiparticle_convention': 'conjugate internal weights; physical CPT/antiparticle attachment is unproved'})
    bosons = [('photon', 'unbroken Q direction conditional on the neutral Higgs vacuum'),
              ('gluon', 'eight color-current directions; one color-adjoint field target'),
              ('W_plus', 'weak raising current, [Q,T+]=T+'),
              ('W_minus', 'weak lowering current, [Q,T-]=-T-'),
              ('Z', 'neutral span(T3,Y), with source-selected mixing required'),
              ('Higgs', 'radial mode conditional on a selected neutral doublet vacuum')]
    for target, candidate in bosons:
        rows.append({'target': target, 'target_kind': 'boson_species',
                     'source_interface': candidate, 'internal_components': {},
                     'source_family_index': None,
                     'missing': common+['source-derived kinetic action and vacuum/mixing coefficients'],
                     'color_components': 8 if target == 'gluon' else 1,
                     'antiparticle_convention': 'W_plus and W_minus are listed separately; neutral fields require their physical field attachment'})
    for row in rows:
        row.update({'correspondence_status': 'internal_representation_candidate',
                    'identified_OPH_particle_source': None,
                    'identified_topological_defect': None,
                    'physical_particle_identified': False,
                    'predicted_mass': None,
                    'physical_pole_or_resonance_proved': False})
    result = {'schema': 'oph.standard_model_source_map.v1',
              'title': 'Standard Model target channels and computed source representations',
              'source_bindings': pins,
              'source_weights': list(algebra.source_components()),
              'complete_internal_component_match': components,
              'gauge_and_Higgs_readout': gauge,
              'target_dictionary_reference': REFERENCE,
              'target_dictionary_enters_source_weight_calculation': False,
              'empirical_mass_or_coupling_input_used': False,
              'source_family_copies_derived': False,
              'configured_family_copies': 3,
              'configured_family_factor_read_from_source': True,
              'finite_source_attachments': science,
              'even_internal_occupation_parity_identified_with_spacetime_statistics': False,
              'rows': rows, 'physical_particle_map': [],
              'physical_mass_prediction_count': 0,
              'physical_alpha_s_prediction': None,
              'parent_full_scientific_replays_performed': False,
              'parents_cold_scientifically_replayed': [],
              'parent_source_bindings_authenticated': True}
    if source_bindings() != pins:
        raise RuntimeError('source changed during export')
    return result


def validate(value):
    if canonical(value) != canonical(build()):
        raise ValueError('source correspondence, physical admission or provenance changed')


def render(value):
    lines = ['# Standard Model source correspondence', '',
             '**Physically identified particle sources: none.**', '',
             'The table compares exact internal source representations with Standard Model target channels.',
             'It assigns no observed mass, topological defect or physical particle identity.',
             'The configured field has three family coordinates. Their measured-species assignment is unproved.',
             'The source action coefficients and four complete duty policies are exported in [SOURCE_MASSES.md](SOURCE_MASSES.md).',
             'Neutrino rows are flavor channels, not identified mass eigenstates. Gluons retain eight color components.', '',
             '| Target | Computed source interface | Physical source | Predicted mass |',
             '| --- | --- | --- | --- |']
    for row in value['rows']:
        description = '; '.join(k+': masks '+','.join(map(str, v['source_masks']))
                               for k, v in row['internal_components'].items()) or row['source_interface']
        lines.append('| '+row['target']+' | '+description+' | Unidentified | Underived |')
    lines += ['', 'The right-handed comparison components are charge-conjugated left-Weyl internal weights.',
              'The even internal register has an exact current intertwiner into a separately configured odd CAR field.',
              'Free Weyl propagation and finite dynamical gauge identities are checked; the interacting chiral continuum is unproved.',
              'A conserved symmetry label, a topological invariant, an internal multiplet and a particle pole are distinct objects.',
              '', 'Target conventions: [Particle Data Group electroweak review]('+REFERENCE+').',
              'Source computation and admission requirements: [proof](PROOF.md), [machine export](particle_map.json).', '']
    return '\n'.join(lines)


def render_masses(value):
    """Display rounded action coefficients; the JSON retains exact interval bounds."""
    from fractions import Fraction
    def midpoint(pair):
        return float(sum(map(Fraction, pair))/2)
    data = value['finite_source_attachments']['source_action_masses']
    lines = ['# Source action mass coefficients', '',
             'These are the static mean-feedback coefficients with canonical kinetic form and unit neutral Higgs norm.',
             'They are not measured pole masses. No GeV conversion or physical family assignment is supplied.',
             'Every complete duty policy is retained; no policy is selected by comparison with experiment.', '',
             '| Duty | Internal channel | Ascending action coefficients | Ratios to smallest | Koide quotient |',
             '| --- | --- | --- | --- | --- |']
    for policy in data['policies']:
        for row in policy['channels']:
            masses = ', '.join(f'{midpoint(box):.10g}' for box in row['ascending_dimensionless_action_mass_intervals'])
            ratios = ('undefined (exact zero)' if row['ratios_to_smallest_interval'] is None else
                      ' : '.join(f'{midpoint(box):.10g}' for box in row['ratios_to_smallest_interval']))
            q = midpoint(row['Koide_diagnostic']['positive_root_Koide_ratio_interval'])
            lines.append(f"| {policy['policy']} | {row['channel']} | {masses} | {ratios} | {q:.10g} |")
    lines += ['', 'Exact interval endpoints and authenticated source provenance are in [particle_map.json](particle_map.json).',
              'All policies have three zero coefficients in the minimal neutrino action and zero charged-current CP invariant.',
              'These do not identify physical neutrino masses or reproduce observed quark CP violation.',
              'The [feedback theorem](FEEDBACK_LIMIT.md) controls the declared weak-pulse limit without changing these coefficients.',
              'The [Koide test](KOIDE_ATTACHMENT.md) excludes balanced positive-root Koide for all four charged-lepton rows.',
              'A [different coherent program](KOIDE_SOURCE.md) produces a balanced family generator with an unselected phase.',
              'See the [mass interpretation](MASS_READOUT.md) for the physical attachment and clock requirements.', '']
    return '\n'.join(lines)
