"""Recompute the native Gauss, transport and complete-record certificate."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

from tools import strict_json
from . import algebra

ROOT = Path(__file__).resolve().parents[2]
PACKAGE = Path(__file__).resolve().parent
OUTPUT = PACKAGE/'results.json'


def source_bindings():
    paths = {str(p.relative_to(ROOT)) for p in PACKAGE.iterdir()
             if p.suffix in ('.py', '.md')}
    paths.update(str(p.relative_to(ROOT)) for p in
                 (ROOT/'evidence/observer_dynamics_20260925/oph-physics-sim').rglob('*.py'))
    paths.update(('evidence/native_geometric_source_20260925/spec.json', 'model/requirements.txt'))
    pending = list(paths)
    seen = set()
    def add(module):
        if not module.startswith(('model.', 'tools.')):
            return
        parts = module.split('.')
        candidates = [module.replace('.', '/')+'.py', module.replace('.', '/')+'/__init__.py']
        candidates += ['/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for path in candidates:
            if (ROOT/path).is_file() and path not in paths:
                paths.add(path)
                pending.append(path)
    while pending:
        path = pending.pop()
        if path in seen or not path.endswith('.py'):
            continue
        seen.add(path)
        module = path[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/path).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                if node.level:
                    target = importlib.util.resolve_name(target, package)
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)}


def produce():
    pins = source_bindings()
    certificate = algebra.exact_certificate()
    histories = [{**block, 'preparation_weight': str(algebra.s.Rational(block['dimension'], 15)),
                  'history': algebra.complete_history(block['dimension'])}
                 for block in certificate['native_blocks']]
    if pins != source_bindings():
        raise ValueError('source changed during calculation')
    return {'schema': 'oph.native_quantum_gauss.v1', 'issues': [1140, 1141, 1142],
            'exact_native_current_certificate': certificate, 'complete_histories': histories,
            'scope': {
                'actual_native_fifteen_state_currents': True,
                'finite_nonabelian_local_operator_Gauss_law': True,
                'native_multiplets_preserved_by_link_transport': True,
                'dynamical_link_and_matter_transport_on_same_Hamiltonian': True,
                'all_site_readout_outcomes_retained': True,
                'detector_energy_exchange_accounted': True,
                'static_dual_reference_is_declared_input': True,
                'link_register_and_coherent_hopping_are_declared_design': True,
                'bare_CE_compilation': False,
                'Wilson_unitary_links': False,
                'same_compact_defect_band_gauge_dressing': False,
                'continuum_chirality_or_regulator_anomaly_cancellation': False,
                'physical_mass_or_coupling_calibration': False,
                'integrated_gate_accepted': False},
            'canonical_physical_premises_discharged': [], 'source_bindings': pins}


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def verify(path=OUTPUT):
    actual = strict_json.load(path)
    expected = produce()
    if actual != expected:
        raise ValueError('native Gauss receipt differs from exact recomputation')
    return True


def main():
    parser = argparse.ArgumentParser(__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(canonical(produce()))
        print(OUTPUT)
    else:
        verify()
        print('native quantum Gauss certificate OK')


if __name__ == '__main__':
    main()
