"""Complete source-bound replay for compact opposite-degree defects."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from tools import strict_json
from . import certification, spectral, topology

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def source_bindings():
    from model.common_source.chiral_candidate.action import SOURCE_FILES
    names = set(SOURCE_FILES) | {'model/requirements.txt', 'requirements.txt'}
    names.update(str(p.relative_to(ROOT)) for p in PACKAGE.iterdir()
                 if p.is_file() and p.suffix in ('.py', '.md'))
    names.update(('model/lattice_texture_selection/PROOF.md',
                  'model/defect_current_bridge/PROOF.md'))
    pending, visited = list(names), set()

    def add_module(name):
        if not name.startswith(('model.', 'tools.')):
            return
        stem = name.replace('.', '/')
        parts = stem.split('/')
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += ['/'.join(parts[:end])+'/__init__.py' for end in range(1, len(parts))]
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in names:
                names.add(candidate); pending.append(candidate)

    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removesuffix('.py').replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add_module(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add_module(target)
                for alias in node.names:
                    add_module(target+'.'+alias.name)
    if any(p.startswith('docs/') for p in names):
        raise ValueError('a prose note cannot be source evidence')
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(names)}


def scope():
    return {'compact_opposite_degree_texture_is_chosen_design': True,
            'exact_local_shell_degrees_and_constant_exterior': True,
            'uniform_finite_spin_energy_and_fermion_energy_difference': True,
            'four_dimensional_band_exhausted_in_finite_and_infinite_volume': True,
            'exact_full_spectrum_certificate_excludes_zero_modes': True,
            'same_native_currents_and_symmetric_pairing_on_isolated_band': True,
            'topological_charge_requires_exact_zero_mode': False,
            'mode_occupation_creates_texture': False,
            'autonomous_texture_formation_or_stability_proved': False,
            'pair_energy_is_sum_of_identified_particle_rest_masses': False,
            'observer_consensus_particle_identity_or_relativistic_dispersion': False,
            'interacting_continuum_or_Einstein_dynamics': False,
            'integrated_gate_accepted': False}


def produce():
    from fractions import Fraction as F
    pins = source_bindings()
    cases = [spectral.case(radius, mass, a) for radius, mass, a in
             ((3, 16, 1), (3, 64, 1), (4, 64, 1), (4, 64, 2))]
    certificate = certification.complete_certificate()
    finite_lower = F(certificate['certified_core_abs_energy_lower'])
    boundary_error = F(cases[1]['exact_bounds']['ordered_band_energy_infinite_volume_error_upper'])
    infinite_lower = finite_lower-boundary_error
    if infinite_lower <= 0:
        raise ArithmeticError('the infinite-volume nonzero-band certificate failed')
    energy_error = max(abs(a-b) for a, b in zip(cases[1]['isolated_numerical_energies'],
                                               cases[2]['isolated_numerical_energies'], strict=True))
    if energy_error > float(boundary_error)+1e-11:
        raise ArithmeticError('volume comparison violates analytic bound')
    result = {'schema': 'oph.compact_defect_pair.v1', 'all_pass': True,
              'issues': [1139, 1140],
              'design': {'texture': 'f=(x,y,z^2-a^2) inside radius a+1; north outside',
                         'zero_mass_cores': '(0,0,+a) and (0,0,-a)',
                         'clock_and_length': 'same source hopping v/a_lattice=1; chosen lattice spacing',
                         'spin_texture_core_fill': 'north; finite ordered-spin rotation orbit',
                         'native_internal_lift': 'same M15 times three declared families'},
              'exact_topology_and_energy': [topology.topology(a) for a in (1, 2)],
              'finite_spectra': cases,
              'exact_leading_schur_terms': [topology.leading_schur(a) for a in (1, 2)],
              'complete_exact_certificate': certificate,
              'infinite_volume': {'nonzero_isolated_band_dimension': 4,
                                  'core_abs_energy_certified_lower': str(infinite_lower),
                                  'boundary_error_at_L3': str(boundary_error),
                                  'numerical_L3_L4_max_difference': round(energy_error, 12),
                                  'continuum_limit_is_not_taken': True},
              'same_current_and_pairing_control': spectral.lifted_control(),
              'scope': scope(), 'canonical_physical_premises_discharged': [],
              'physical_comparison_count_added': 0,
              'source_bindings': pins}
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return result


def report():
    pins = source_bindings()
    with tempfile.TemporaryDirectory(prefix='oph-compact-pair-') as temporary:
        path = Path(temporary)/'worker.json'
        child = subprocess.run([sys.executable, '-W', 'error', '-m',
                                'model.compact_defect_pair.packet', '--worker', str(path)],
                               cwd=ROOT, env=dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1'),
                               check=False, capture_output=True, text=True, timeout=600)
        if child.returncode:
            raise RuntimeError('compact-pair worker failed:\n'+child.stderr)
        result = strict_json.load(path)
    if result['source_bindings'] != pins or source_bindings() != pins:
        raise ValueError('source changed across cold replay')
    return result


def validate(value, replay=None):
    if canonical(value) != canonical(report() if replay is None else replay):
        raise ValueError('compact-pair receipt differs from complete cold replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', nargs='?', const=str(OUTPUT))
    group.add_argument('--check', nargs='?', const=str(OUTPUT))
    group.add_argument('--worker')
    args = parser.parse_args()
    if args.worker:
        Path(args.worker).write_text(canonical(produce()))
    elif args.write:
        Path(args.write).write_text(canonical(report()))
        print('wrote compact opposite-degree texture and complete spectrum')
    else:
        validate(strict_json.load(Path(args.check)))
        print('compact topology, finite energy, complete spectrum and native currents verified')


if __name__ == '__main__':
    main()
