"""Check the separately typed finite-field supplement without changing old claims."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ('native_quark_transport', 'native_quark_spatial_spectrum',
    'native_code_ground_preparation', 'native_wilson_cone_audit',
    'native_orthogonal_wilson_bridge', 'native_wilson_maxwell_locality',
    'native_wilson_maxwell_spatial_family', 'native_wilson_maxwell_interaction_witness',
    'native_orthogonal_wilson_spatial_join', 'native_orthogonal_wilson_charged_band',
    'native_orthogonal_wilson_charged_spatial_records')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check():
    data = json.loads((ROOT/'claims/addenda/field_record_addenda.json').read_text())
    if (data['schema'] != 'oph.field_record_addenda_claims.v1'
            or data['changes_existing_physical_statuses']
            or data['physical_particle_identifications']
            or data['source_derived_masses_GeV']):
        raise ValueError('Unreviewed physical promotion in finite-field supplement')
    if data['sealed_program_class_sha256'] != sha(ROOT/'model/completion/program_class_v2.json'):
        raise ValueError('Sealed parent class changed')
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
    print('PASS: eleven typed finite-field claims; sealed class and physical statuses preserved')


if __name__ == '__main__':
    check()
