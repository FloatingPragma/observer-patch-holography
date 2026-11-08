"""Full same-source two-charge algebra, response and complete detector replay."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
from tools import strict_json
from . import algebra, response, continuum

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/'results.json'
PROOFS = ('model/native_ordered_ground/README.md',
          'model/protected_defect_catalog/PROOF.md',
          'model/native_spin_observer/README.md',
          'model/spin_operational_position/PROOF.md',
          'model/spin_operational_position/generic.py',
          'model/spin_operational_position/bounds.py',
          'model/spin_canonical_field/README.md',
          'model/spin_canonical_energy/README.md',
          'model/spin_charge_dispersion/PROOF.md',
          'model/scattering_packet_readout/PROOF.md',
          'model/charge_channel_lifetime/PROOF.md',
          'Lean/QFT/SourceOperatorGeneration.lean',
          'Lean/QFT/FiniteUnitaryScatteringNoGo.lean',
          'Lean/Dynamics/ProtectedCharge.lean',
          'model/requirements.txt',
          'code/source_selection_model/response.json',
          'code/m1_source_realization/receipt.json',
          'evidence/native_geometric_source_20260925/spec.json')


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n'


def bindings():
    legacy = 'evidence/observer_dynamics_20260925/oph-physics-sim/'
    names = set(PROOFS)|{'model/native_ordered_ground/witness.npz',
                        legacy+'oph_exact/carrier.py', legacy+'oph_exact/__init__.py'}
    names.update(path.relative_to(ROOT).as_posix() for path in HERE.iterdir()
                 if path.suffix in ('.py', '.md'))
    pending = list(names)
    visited = set()

    def add(module):
        parts = module.split('.')
        prefix = legacy if parts[0] in ('oph_exact', 'oph_fpe') else ''
        stem = prefix+'/'.join(parts)
        candidates = [stem+'.py', stem+'/__init__.py']
        candidates += [prefix+'/'.join(parts[:i])+'/__init__.py' for i in range(1, len(parts))]
        for name in candidates:
            if (ROOT/name).is_file() and name not in names:
                names.add(name)
                pending.append(name)
    while pending:
        name = pending.pop()
        if name in visited or not name.endswith('.py'):
            continue
        visited.add(name)
        module = name.removeprefix(legacy)[:-3].replace('/', '.')
        package = module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        add(module)
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add(alias.name)
            elif isinstance(node, ast.ImportFrom):
                target = '.'*node.level+(node.module or '')
                target = importlib.util.resolve_name(target, package) if node.level else target
                add(target)
                for alias in node.names:
                    add(target+'.'+alias.name)
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def build():
    before = bindings()
    algebra.data.cache_clear()
    value = {'schema':'oph.source_charge_scattering.v1',
             'source_bindings':before, 'actual_source_algebra':algebra.report(),
             'complete_response':response.report(),
             'separated_continuum_charge_pairs':continuum.report()}
    if bindings()!=before:
        raise RuntimeError('source changed during pair response certificate')
    return value


def validate(value):
    if canonical(value)!=canonical(build()):
        raise ValueError('pair receipt, full outcomes, source or physical scope differs')
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',type=Path,nargs='?',const=OUTPUT)
    group.add_argument('--check',type=Path,nargs='?',const=OUTPUT)
    args=parser.parse_args()
    if args.write:
        args.write.write_text(canonical(build()))
    else:
        validate(strict_json.load(args.check))
    print('PASS: complete source pair algebra, propagation and detector outcomes')
