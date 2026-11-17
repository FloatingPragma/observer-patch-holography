"""Check the separately typed finite-field supplement without changing old claims."""
from pathlib import Path
import copy
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from model.completion import program_class_v2
from tools import strict_json

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


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def check_parent_class(data, root=ROOT):
    """The supplement names its sealed parent class; the live class keeps the
    adopted review semantics and passes the class contract checker."""
    root = Path(root).resolve()
    if (data['sealed_parent_commit'] != SEALED_PARENT
            or data['sealed_program_class_sha256'] != SEALED_CLASS_SHA256):
        raise ValueError('Sealed parent authority changed')
    klass = strict_json.load(root/program_class_v2.PATH)
    adoption = copy.deepcopy(klass['adoption'])
    adoption['independent_review'].pop('sha256')
    adoption['independent_review'].pop('reviewed_contract_digest')
    if canonical(adoption) != canonical(ADOPTION_SEMANTICS):
        raise ValueError('Sealed parent adoption semantics changed')
    program_class_v2.validate(klass, root)


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
        if not (ROOT/row['receipt']).is_file() or not (ROOT/row['proof']).is_file():
            raise ValueError('Missing finite claim receipt or proof: '+name)
        value = json.loads((ROOT/row['receipt']).read_text())
        if not value['public_port']['exact_scientific_payload_equal'] or value['public_port']['physical_status_promotion']:
            raise ValueError('Unreviewed public source interpretation')
    if tuple(names) != PACKAGES or [r['id'] for r in data['claims']] != [f'FR{i:02}' for i in range(1,12)]:
        raise ValueError('Changed claim membership')
    print('PASS: eleven typed finite-field claims; sealed class adoption and physical statuses preserved')


if __name__ == '__main__':
    check()
