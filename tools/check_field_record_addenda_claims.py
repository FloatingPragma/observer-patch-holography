"""Check the separately typed finite-field supplement without changing old claims."""
from pathlib import Path
import copy
import hashlib
import json
import sys

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
PACKAGES = ('native_quark_transport', 'native_quark_spatial_spectrum',
    'native_code_ground_preparation', 'native_wilson_cone_audit',
    'native_orthogonal_wilson_bridge', 'native_wilson_maxwell_locality',
    'native_wilson_maxwell_spatial_family', 'native_wilson_maxwell_interaction_witness',
    'native_orthogonal_wilson_spatial_join', 'native_orthogonal_wilson_charged_band',
    'native_orthogonal_wilson_charged_spatial_records')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_parent_class(data, root=ROOT):
    root = Path(root).resolve()
    if (data['sealed_parent_commit'] != SEALED_PARENT
            or data['sealed_program_class_sha256'] != SEALED_CLASS_SHA256):
        raise ValueError('Sealed parent authority changed')
    klass = strict_json.load(root/program_class_v2.PATH)
    # The existing conservation policy excludes only the named live member
    # hash and first authenticates it against the actual receipt bytes.
    if (custody_renewal.digest(custody_renewal.class_semantics(klass, root))
            != custody_renewal.CLASS_SEMANTICS_SHA256):
        raise ValueError('Sealed parent class semantics changed')
    adoption = copy.deepcopy(klass['adoption'])
    adoption['independent_review'].pop('sha256')
    adoption['independent_review'].pop('reviewed_contract_digest')
    if renewal.canonical(adoption) != renewal.canonical(ADOPTION_SEMANTICS):
        raise ValueError('Sealed parent adoption semantics changed')
    program_class_v2.validate(klass, root)

    # Authenticate the member's complete scientific payload and source scope
    # against its existing immutable owner archive; a re-signed class alone
    # cannot admit a scientifically changed member receipt.
    if sha(root/renewal.CATALOG) != renewal.CATALOG_SHA256:
        raise ValueError('Live member owner authority changed')
    catalog = strict_json.load(root/renewal.CATALOG)
    renewal.verify_executable_source_pins(catalog['reviewed_executable_source_sha256'], root)
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
          'preserved under authenticated current custody')


if __name__ == '__main__':
    check()
