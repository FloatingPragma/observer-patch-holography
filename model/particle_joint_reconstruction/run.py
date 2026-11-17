"""Validate the CCF-LT-2 acceptance graph against its actual finite domains."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

from model.particle_joint_reconstruction import candidate
from tools import strict_json
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / 'results.json'
PROGRAM = 'CCF-LT-2'
NODE_IDS = (
    'fixed_program', 'physical_preparation', 'local_records', 'joint_refinement',
    'relativistic_field', 'chiral_species_attachment', 'particle_observables',
    'independent_comparison', 'pure_gauge_gap', 'geometric_gravity',
)
FINITE_CONSUMERS = {
    'fixed_program': ('source_edge', 'source_triangle', 'local_capacity'),
    'physical_preparation': ('reference',),
    'local_records': ('actual_local_capacity_join', 'phase_instrument'),
}
FINITE_DOMAINS = {'fixed_program': 'The finite CCF-LT-2 graph family, with 96 native CAR modes per vertex, integer electric links, and the independent per-vertex local-capacity response law.', 'physical_preparation': 'Declared preparations, retained charge references and the faithful beta=1/1024 Gibbs family on the complete Gauss-zero cap6 support; arbitrary native propagation is not asserted.', 'local_records': 'Full native polynomial and local-capacity basis columns, exact initial complete phase instrument, and separate complete finite field, Dirac and geometry diagnostics.'}
REQUIRED_EDGES = {
    node: ([] if index == 0 else ['relativistic_field'] if node in ('pure_gauge_gap', 'geometric_gravity')
           else [NODE_IDS[index-1]]) for index, node in enumerate(NODE_IDS)
}
COMPONENT_REQUIREMENTS = {
    'field_common_theory_ten': {str(n) for n in range(1, 11)},
    'native_full_spectrum_ten': {str(n) for n in range(1, 11)},
    'defect_identification_ten': {str(n) for n in range(1, 11)},
    'yang_mills_Y1_Y6': {f'Y{n}' for n in range(1, 7)},
}


def validate_dependencies(graph):
    if graph.get('schema') != 'oph.joint_reconstruction_dependencies.v1':
        raise ValueError('unknown dependency schema')
    if graph.get('primary_program') != PROGRAM:
        raise ValueError('dependency graph belongs to another program')
    rows = graph['nodes']
    if [row['id'] for row in rows] != list(NODE_IDS):
        raise ValueError('the unified graph must contain its ten distinct nodes')
    if [row['number'] for row in rows] != list(range(1, 11)):
        raise ValueError('dependency numbering differs')
    known = set(NODE_IDS)
    assumptions = {row['id'] for row in candidate.assumption_ledger()}
    covered = set()
    completed = set()
    for row in rows:
        dependencies = row['depends_on']
        if len(dependencies) != len(set(dependencies)) or not set(dependencies) <= completed:
            raise ValueError('dependency graph is cyclic or not in acceptance order')
        if dependencies != REQUIRED_EDGES[row['id']]:
            raise ValueError('required acceptance dependency was removed or changed')
        completed.add(row['id'])
        expected = ('finite_constructed' if row['id'] in FINITE_CONSUMERS else
                    'conditional_criterion' if row['id'] == 'pure_gauge_gap' else
                    'physical_proof_required')
        if row['status'] != expected:
            raise ValueError('a physical premise cannot be promoted by editing its status')
        if row['id'] in FINITE_DOMAINS and row['domain'] != FINITE_DOMAINS[row['id']]:
            raise ValueError('finite consumer domain differs from the actual experiment')
        if not row['domain'] or not row['acceptance'] or not row['evidence']:
            raise ValueError('every requirement needs a domain, acceptance conditions and evidence')
        for item in row['evidence']:
            path = (ROOT / item['path']).resolve()
            if not path.is_relative_to(ROOT) or not path.is_file() or not item['scope']:
                raise ValueError('evidence needs an existing repository file and its scope')
        listed = set(row['conditional_assumptions'])
        if not listed <= assumptions:
            raise ValueError('unregistered conditional assumption')
        covered |= listed
    if covered != assumptions:
        raise ValueError('conditional assumptions are omitted from the unified graph')
    extension = graph['conditional_candidate']
    if (extension['id'] != candidate.EXTENSION_ID
            or extension['discharges_primary_physical_dependencies']
            or extension['calibrated_tree_coefficients_are_independent_predictions']
            or extension['full_theory_positive_global_gap_claimed']):
        raise ValueError('conditional extension cannot discharge primary physics or a global photon gap')
    pure = rows[8]
    if (pure['required_for_full_particle_spectrum']
            or pure['implies_full_theory_positive_global_gap']
            or pure['required_shared_nodes'] != list(NODE_IDS[:5])):
        raise ValueError('pure Yang-Mills must instantiate its own source through reconstruction')
    gravity = rows[9]
    if gravity['required_for_full_particle_spectrum'] or gravity['proper_geometry_or_Einstein_limit_proved']:
        raise ValueError('finite stress cannot be promoted to a physical Einstein limit')
    maps = graph['requirement_maps']
    if set(maps) != set(COMPONENT_REQUIREMENTS):
        raise ValueError('component requirement coverage is incomplete')
    for name, mapping in maps.items():
        if set(mapping) != COMPONENT_REQUIREMENTS[name]:
            raise ValueError('component requirement coverage is incomplete')
        for targets in mapping.values():
            if not targets or any(target not in known and not target.startswith('outside_target:')
                                  for target in targets):
                raise ValueError('a component requirement has no unified destination')



def build():
    from tools.build_ccf_lt2_claims import generate
    graph = strict_json.load(HERE / 'dependencies.json')
    validate_dependencies(graph)
    official = generate()
    return {
        'schema': 'oph.joint_reconstruction.v2',
        'primary_program': PROGRAM,
        'official_model': official,
        'finite_primary_acceptance_count': len(FINITE_CONSUMERS),
        'acceptance': [{
            'id': row['id'], 'status': row['status'], 'domain': row['domain'],
            'depends_on': row['depends_on'],
            'finite_scope': list(FINITE_CONSUMERS.get(row['id'], ())),
        } for row in graph['nodes']],
        'full_native_state_propagation_claimed': False,
        'full_particle_spectrum_derived': False,
        'interacting_relativistic_primary_limit_proved': False,
        'proper_geometry_or_Einstein_limit_proved': False,
        'pure_Yang_Mills_gap_proved': False,
        'source_derived_particle_mass_table_GeV': [],
        'input_sha256': bound_pins({
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [HERE/'dependencies.json', HERE/'run.py',
                         ROOT/'claims/ccf_lt2.json']
        }),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    payload = json.dumps(build(), sort_keys=True, indent=2)+'\n'
    if args.write:
        OUTPUT.write_text(payload)
    elif OUTPUT.read_text() != payload:
        raise SystemExit('CCF-LT-2 dependency receipt differs; regenerate from its actual inputs')
    print('PASS_CCF_LT2_DEPENDENCY_DOMAINS')


if __name__ == '__main__':
    main()
