"""Check the separately typed finite-field supplement without changing old claims."""
from pathlib import Path
import copy
import ast
import hashlib
import json
import sys
import io
import os
import subprocess
import tarfile
import tempfile
from functools import lru_cache

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.completion import custody_renewal, program_class_v2
from tools import renew_canonical_live_receipts as renewal, strict_json

# Declared supplement seal, preserved independently of live review hashes.
# The parent Git blob binds this exact reviewed class.
SEALED_PARENT = '50a17bc07b9ad4b2fe93ad73736fbb433c01183a'
SEALED_CLASS_SHA256 = '170e0e32b2ca1567bb55da948929611a1b93db6edb7f45529bf5c552041a3b81'
ADOPTION_SEMANTICS = {
    'status': 'ADOPTED_AFTER_REVIEW', 'changes_current_registry': True,
    'independent_review': {
        'verdict': 'PASS', 'path': 'model/completion/program_class_v2_review.json'}}
REVIEWED_EXPOSURE_TRANSITIONS = {
    'model/local_repair_family/record_clock_results.json':
        'ea709f6bb519b780b4be325ccb2a549f24fd89ac7a03048b7afd45f43171bf66',
}
PACKAGES = ('native_quark_transport', 'native_quark_spatial_spectrum',
    'native_code_ground_preparation', 'native_wilson_cone_audit',
    'native_orthogonal_wilson_bridge', 'native_wilson_maxwell_locality',
    'native_wilson_maxwell_spatial_family', 'native_wilson_maxwell_interaction_witness',
    'native_orthogonal_wilson_spatial_join', 'native_orthogonal_wilson_charged_band',
    'native_orthogonal_wilson_charged_spatial_records')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def check_sealed_class(source_root):
    """Execute the original signed checker in its explicitly sealed source domain.

    The temporary view contains only authenticated Git bytes. It never replaces
    the live CCF-LT-2 hypothesis or changes the supplement's declared authority.
    """
    source_root = Path(source_root).resolve()
    def blob(name):
        return subprocess.check_output(['git', 'show', SEALED_PARENT + ':' + name],
                                       cwd=source_root)
    class_path = str(program_class_v2.PATH)
    raw = blob(class_path)
    if hashlib.sha256(raw).hexdigest() != SEALED_CLASS_SHA256:
        raise ValueError('Sealed parent class bytes changed')
    klass = strict_json.loads(raw.decode())
    review_path = klass['adoption']['independent_review']['path']
    review = strict_json.loads(blob(review_path).decode())
    paths = {class_path, review_path, 'claims/executable_model_hypothesis.json',
             'tracking/observation_ledger.json', 'tools/strict_json.py',
             klass['laws']['local_alternative_family']['contract']}
    paths.update(klass['exposure']['source_hashes'])
    paths.update(klass['interpreter_binding']['source_hashes'])
    paths.update(klass['review_implementation_bindings'])
    paths.update(review['source_bindings'])
    # Package initializers and exposure aliases are included when present in
    # this exact tree. The original checker determines its original inputs.
    optional = {'model/__init__.py', 'model/completion/__init__.py', 'tools/__init__.py',
                *program_class_v2.EXPOSURE_SNAPSHOTS.values()}
    for name in optional:
        exists = subprocess.run(['git', 'cat-file', '-e', SEALED_PARENT + ':' + name],
                                cwd=source_root, capture_output=True)
        if exists.returncode == 0:
            paths.add(name)
    if any(Path(name).is_absolute() or '..' in Path(name).parts for name in paths):
        raise ValueError('External sealed class source')
    archive = subprocess.check_output(['git', 'archive', '--format=tar', SEALED_PARENT,
                                       *sorted(paths)], cwd=source_root)
    with tempfile.TemporaryDirectory(prefix='oph-sealed-component-') as directory:
        view = Path(directory)
        with tarfile.open(fileobj=io.BytesIO(archive)) as source:
            source.extractall(view, filter='data')
        env = dict(os.environ)
        env['PYTHONPATH'] = str(view)
        env['PYTHONDONTWRITEBYTECODE'] = '1'
        executed = subprocess.run([sys.executable, '-m', 'model.completion.program_class_v2',
                                   '--check'], cwd=view, env=env, capture_output=True, text=True)
        if executed.returncode:
            raise ValueError('Sealed independent class check failed:\n'
                             + executed.stdout + executed.stderr)
    return {'sealed_class_valid': True, 'official_CCF_LT2_adoption_claimed': False}


def validate_live_review_signatures(klass, root):
    """Authenticate review linkage without treating its class as current law."""
    review = klass['adoption']['independent_review']
    path = Path(root) / review['path']
    if sha(path) != review['sha256']:
        raise ValueError('Sealed component review bytes changed')
    value = strict_json.load(path)
    fingerprint = program_class_v2.reviewed_contract_digest(klass)
    if (review.get('reviewed_contract_digest') != fingerprint
            or value.get('reviewed_contract_digest') != fingerprint
            or value.get('verdict') != 'PASS'):
        raise ValueError('Sealed component review does not bind its class')


INVENTORY_HELPER = 'tools/build_completion_program.py'
CURRENT_INVENTORY_SHA256 = 'edc4642feaf6bc555282e108d63317fc9e852a52ab89b0ed9fbca6f1e0c11510'
INVENTORY_FUNCTIONS = frozenset(('require', 'local_file', 'integration_source_paths',
                                'validate_integration_sources'))


def verify_live_executable_sources(pins, root):
    """One reviewed inventory-only helper update; every other source pin stays exact."""
    root = Path(root)
    old_hash = pins[INVENTORY_HELPER]
    source = root / INVENTORY_HELPER
    if sha(source) == old_hash:
        renewal.verify_executable_source_pins(pins, root)
        return
    if sha(source) != CURRENT_INVENTORY_SHA256:
        raise ValueError('Unreviewed current inventory helper source')
    original = subprocess.check_output(['git', 'show', SEALED_PARENT + ':' + INVENTORY_HELPER],
                                       cwd=ROOT)
    if hashlib.sha256(original).hexdigest() != old_hash:
        raise ValueError('Original inventory helper is outside the sealed catalog')
    def functions(raw):
        parsed = ast.parse(raw)
        return {node.name: ast.dump(node, include_attributes=False)
                for node in parsed.body if isinstance(node, ast.FunctionDef)
                and node.name in INVENTORY_FUNCTIONS}
    if (functions(original) != functions(source.read_bytes())
            or set(functions(original)) != INVENTORY_FUNCTIONS):
        raise ValueError('Consumed inventory helper semantics changed')
    remaining = {name: value for name, value in pins.items() if name != INVENTORY_HELPER}
    renewal.verify_executable_source_pins(remaining, root)
    for name in remaining:
        if not name.endswith('.py'):
            continue
        for node in ast.walk(ast.parse((root/name).read_bytes())):
            if isinstance(node, ast.ImportFrom) and node.module == 'tools.build_completion_program':
                if any(alias.name not in INVENTORY_FUNCTIONS for alias in node.names):
                    raise ValueError('Component consumes current programme authority: ' + name)
            if (isinstance(node, ast.Import)
                    and any(alias.name == 'tools.build_completion_program' for alias in node.names)):
                raise ValueError('Unreviewed inventory module consumer: ' + name)
            if (isinstance(node, ast.ImportFrom) and node.module == 'tools'
                    and any(alias.name == 'build_completion_program' for alias in node.names)):
                raise ValueError('Unreviewed inventory module consumer: ' + name)


def check_parent_class(data, root=ROOT):
    root = Path(root).resolve()
    if (data['sealed_parent_commit'] != SEALED_PARENT
            or data['sealed_program_class_sha256'] != SEALED_CLASS_SHA256):
        raise ValueError('Sealed parent authority changed')
    klass = strict_json.load(root/program_class_v2.PATH)
    # The existing conservation policy excludes only the named live member
    # hash and first authenticates it against the actual receipt bytes.
    semantics = custody_renewal.class_semantics(klass, root)
    # Reviewed exposure transition: the record-clock receipt changed only in its
    # m1_source_realization pin. The class pin must equal the live bytes; the
    # sealed semantics are compared at the retained value.
    exposure = semantics['exposure']['source_hashes']
    for name, retained in REVIEWED_EXPOSURE_TRANSITIONS.items():
        if exposure.get(name) != sha(root/name):
            raise ValueError('Reviewed exposure transition differs from live bytes: ' + name)
        exposure[name] = retained
    if custody_renewal.digest(semantics) != custody_renewal.CLASS_SEMANTICS_SHA256:
        raise ValueError('Sealed parent class semantics changed')
    adoption = copy.deepcopy(klass['adoption'])
    adoption['independent_review'].pop('sha256')
    adoption['independent_review'].pop('reviewed_contract_digest')
    if renewal.canonical(adoption) != renewal.canonical(ADOPTION_SEMANTICS):
        raise ValueError('Sealed parent adoption semantics changed')
    # Only the original sealed review is an authority for this independent class.
    check_sealed_class(ROOT)

    # Authenticate the member's complete scientific payload and source scope
    # against its existing immutable owner archive; a re-signed class alone
    # cannot admit a scientifically changed member receipt.
    if sha(root/renewal.CATALOG) != renewal.CATALOG_SHA256:
        raise ValueError('Live member owner authority changed')
    catalog = strict_json.load(root/renewal.CATALOG)
    verify_live_executable_sources(catalog['reviewed_executable_source_sha256'], root)
    spec = catalog['owners']['native_program']
    if sha(root/spec['retained_path']) != spec['retained_sha256']:
        raise ValueError('Retained member receipt changed')
    retained = strict_json.load(root/spec['retained_path'])
    current = strict_json.load(root/spec['receipt'])
    projected = renewal.carry_metadata(spec['receipt'], retained, current, root=root)
    if renewal.canonical(current) != renewal.canonical(projected):
        raise ValueError('Live member scientific payload changed')
    for _, relative, expected, target in renewal.bindings(root, spec['receipt'], current):
        if sha(target) != expected:
            raise ValueError('Live member source changed: '+relative)


def check():
    data = strict_json.load(ROOT/'claims/addenda/field_record_addenda.json')
    if (data['schema'] != 'oph.field_record_addenda_claims.v1'
            or data['changes_existing_physical_statuses']
            or data['physical_particle_identifications']
            or data['source_derived_masses_GeV']):
        raise ValueError('Unreviewed physical promotion in finite-field supplement')
    check_parent_class(data)
    names = []
    for row in data['claims']:
        name = row['package']; names.append(name)
        prefix = 'model/field_record_addenda/'+name+'/'
        if row['receipt'] != prefix+'results.json' or row['proof'] != prefix+'PROOF.md':
            raise ValueError('Cross-package receipt or proof substitution')
        if row['construction_relation'] not in {'SAME_H_NEW_APPARATUS','DECLARED_NEW_H','DECLARED_FAMILY','SAME_FAMILY_NEW_APPARATUS'}:
            raise ValueError('Missing declared source relationship')
        if row['status'] != 'proved_in_stated_finite_construction':
            raise ValueError('Unreviewed claim status')
        if sha(ROOT/row['receipt']) != row['receipt_sha256'] or not (ROOT/row['proof']).is_file():
            raise ValueError('Changed finite claim source: '+name)
        value = json.loads((ROOT/row['receipt']).read_text())
        if not value['public_port']['exact_scientific_payload_equal'] or value['public_port']['physical_status_promotion']:
            raise ValueError('Unreviewed public source interpretation')
    if tuple(names) != PACKAGES or [r['id'] for r in data['claims']] != [f'FR{i:02}' for i in range(1,12)]:
        raise ValueError('Changed claim membership')
    print('PASS: eleven typed finite-field claims; sealed class semantics and physical statuses '
          'preserved in their independent source domain and authenticated current member custody')


if __name__ == '__main__':
    check()
