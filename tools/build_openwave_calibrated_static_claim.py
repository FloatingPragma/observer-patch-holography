"""Bind audited static calibration and native q=0 readers to explicit evidence.

This source-bound export requires the sibling OpenWave workspace. It is not a
standalone simulation or another execution of the enormous field programs.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json

RER=Path(__file__).resolve().parents[1]
WORKSPACE=RER.parent
EXPECTED={
 'openwave-tuning/common_dynamics/calibrated_charge/MANIFEST.json':'66dbf120a4978101e9268c1a5ed841c6d2c87933a2d757666089ff2d65efc78c',
 'openwave-tuning/next/CHAIN6/calibrated_common_anchor/MANIFEST.json':'0046929344355b85501e1836a9145e03b86a6409b7054454de6d70278187d150',
 'openwave-tuning/next/CHAIN6/calibrated_anchor_native/MANIFEST.json':'ebbda3cc44d39803760f7e7f2c4cf845132a4f8e890030567d4d3be97ec51adf',
 'openwave-tuning/next/CHAIN6/calibrated_anchor_native_review/MANIFEST.json':'0ce008cefb6e4815fad307148bc493e0c3679434767dadedd39f7fbbec0c8716',
 'openwave-tuning/next/CHAIN7/calibrated_anchor_connection_review/MANIFEST.json':'ae2045648586972c2059b339556d04214ddc258458be6205bda0725614e1496c',
 'openwave-tuning/row_completion/reviews/CALIBRATED_ANCHOR_ROOT/MANIFEST.json':'2c33c877d9eab6aa50aa4f45c13f4289f6fd27ba308e50103fb1ae7306ef3be0',
 'openwave-tuning/next/CHAIN6/calibrated_anchor_delta_review/MANIFEST.json':'537120c2336dd688b5529fe099311e6928e715672d72b448a3397186410c3c3e',
 'openwave-tuning/row_completion/ROW_20/VALIDATION.json':'e543f927753c04f664d45d24a974ccf203f6930f528b0124828e47c10d4bc728',
 'openwave-tuning/row_completion/ROW_02/CALIBRATED_VALIDATION.json':'af1dda68e4431f7c08fe3ed215ed9fae6b897f5c722168cc76bab0a51ca87a13',
 'openwave-tuning/row_completion/ROW_19/CALIBRATED_VALIDATION.json':'4fbd0c30da8bcf32c63d454c81876cebf0b0ee0e045f40f461b3907a6c8c1ad8',
}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())


def build():
    verified={};resolutions={};cache={}
    def digest(p):
        stat=p.stat();key=(stat.st_ino,stat.st_mtime_ns,stat.st_size)
        if key not in cache:cache[key]=sha(p)
        return cache[key]
    for name,expected in EXPECTED.items():
        path=WORKSPACE/name;assert digest(path)==expected,name
        verified[name]=expected;manifest=read(path)
        for section in ('files','files_sha256','source_pins','consumed_files_sha256',
                        'frozen_inputs_sha256','snapshot_files_sha256','input_files_sha256',
                        'producer_files','evidence_files'):
            for leaf,h in manifest.get(section,{}).items():
                if not isinstance(h,str) or len(h)!=64:continue
                candidates=[path.parent/leaf,path.parent/'inputs'/leaf,path.parent/'inputs'/h,WORKSPACE/leaf,
                            WORKSPACE/'openwave-tuning/next/CHAIN6/source_bound_integration_v4/inputs'/h]
                chosen=next((p for p in candidates if p.is_file() and digest(p)==h),None)
                assert chosen is not None,(name,section,leaf,'exact frozen bytes absent')
                key=name+'::'+section+'::'+leaf
                verified[key]=h;resolutions[key]=str(chosen.relative_to(WORKSPACE))
    anchor=read(WORKSPACE/'openwave-tuning/next/CHAIN6/calibrated_common_anchor/RESULT.json')
    native=read(WORKSPACE/'openwave-tuning/next/CHAIN6/calibrated_anchor_native/RESULT.json')
    connection=read(WORKSPACE/'openwave-tuning/next/CHAIN7/calibrated_anchor_connection_review/RESULT.json')
    apparatus=read(WORKSPACE/'openwave-tuning/next/CHAIN6/calibrated_anchor_native_review/REVIEW.json')
    receipt=read(WORKSPACE/'openwave-tuning/row_completion/ROW_20/VALIDATION.json')
    assert native['source_native_pulses_executed']==282 and connection['checks_pass']
    assert apparatus['candidate_manifest_sha256']==EXPECTED['openwave-tuning/next/CHAIN6/calibrated_anchor_native/MANIFEST.json']
    assert apparatus['full_raw_M6_all96_source_columns'] and apparatus['native_control_certificates_rechecked']
    assert apparatus['all_first_source_failure_maps_checked']==1080 and apparatus['fine_CE_probabilities_checked']==7644
    assert not apparatus['success_normalization'] and not apparatus['full_criterion_met']
    assert receipt['full_criterion_met'] and not receipt['independent_strength_prediction_credit']
    assert len([r for r in receipt['independent_reviews'] if r['full_criterion_met'] and r['producer_independent']])>=2
    assert len(native['energy_readouts'])==2 and native['all8_source_columns_executed']
    assert native['common193_field_word_executed'] is False
    assert native['independent_alpha_prediction'] is False
    assert F(anchor['optical_rest_energy_GeV']['interval_exact'][0])>F(40,1000)
    compositions=[]
    for row in (2,19):
        path=f'openwave-tuning/row_completion/ROW_{row:02d}/CALIBRATED_VALIDATION.json'
        accepted=read(WORKSPACE/path)
        assert accepted['row']==row and accepted['full_criterion_met']
        assert accepted['candidate_manifest_sha256']=='606573167184e68200e186c34fa33269a92a31c290e5e0796017c65caa49a43e'
        assert accepted['row_count_increment']==0 and accepted['official_rows_changed']==[]
        assert not accepted['one_construction_all31_proved']
        reviews=[r for r in accepted['independent_reviews'] if r['producer_independent'] and r['full_criterion_met']]
        assert len(reviews)>=2 and len({r['reviewer'] for r in reviews})>=2
        compositions.append({'row':row,'receipt_path':path,'receipt_sha256':EXPECTED[path],
                             'literal_criterion':accepted['literal_criterion'],
                             'candidate_id':accepted['candidate_id'],
                             'full_source_specific_internal_criterion_met':True,
                             'independent_full_reviews':reviews})
    audits=[dict(path=n,sha256=h) for n,h in EXPECTED.items() if 'review' in n.lower()]
    claim={
      'schema':'oph.openwave.calibrated_static.source_bound.v1',
      'status':'audited_static_components_and_scoped_internal_compositions',
      'source_manifest_sha256':EXPECTED['openwave-tuning/common_dynamics/calibrated_charge/MANIFEST.json'],
      'calibrations':{'alpha_input_exact':'72973525643/10000000000000','isolated_energy_input_keV':'511',
                      'g_squared':'12*pi*alpha_input','independent_charge_strength_prediction':False},
      'source_identity':'J_g* H_g J_g = diag(g^2 C_1(q)) + H_T on the declared held Gauss fiber',
      'executed_reader_spacing_h':'1/2','cell_time':'T_g=hbar/U_g; tau0=T_g/h at the selected member',
      'fixed_tau0_refinement_reanchoring_claimed':False,
      'source_self_energy':anchor['self_energy_g_source'],
      'keV_per_cell_Hamiltonian_energy':anchor['keV_per_cell_Hamiltonian_energy'],
      'optical_rest_energy_GeV':anchor['optical_rest_energy_GeV'],
      'finite_half_pair_native_keV':[{k:r[k] for k in ('target','half_pair_native_keV_display','maximum_complete_TV')} for r in native['energy_readouts']],
      'fresh_source_native_pulses':282,'native_source_input_columns':96,
      'complete_fine_C_E_records':7644,'source_first_export_maps':1080,
      'full_spectator_field_words_executed':False,
      'finite_half_pairs_are_isolated_electron_energy':False,
      'stationary_moving_charged_particle_or_inertial_mass_proved':False,
      'full_calibrated_row2_award_by_this_export':False,
      'accepted_internal_compositions':compositions,
      'composition_scope':'Held calibrated field-energy anchor and neutral optical excitation in the declared source and preparation. No filled-reference moving charged pole, finite inertia or independent mass prediction is established.',
      'row20_accepted_receipt_sha256':EXPECTED['openwave-tuning/row_completion/ROW_20/VALIDATION.json'],
      'independent_audit_manifests':audits,
      'independence_scope':'Each new reader is opposed independently; unit-map and self-energy author roles remain disclosed in their separate reviews. Component flags are not rewritten as full-row awards.',
      'official_OpenWave_awards_added':0,'all31_complete':False,
      'public_source_closure':'requires the sibling oph-meta OpenWave workspace and its pinned snapshots; not standalone RER replay',
      'expected_evidence':EXPECTED,'verified_manifest_leaf_references':len(verified),
      'source_bindings_sha256':hashlib.sha256(json.dumps(verified,sort_keys=True).encode()).hexdigest(),
      'generator_sha256':sha(Path(__file__))}
    return claim


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true');args=p.parse_args()
    claim=build();target=RER/'claims/openwave_calibrated_static.json'
    if args.write:target.write_text(json.dumps(claim,indent=2,sort_keys=True)+'\n')
    else:assert read(target)==claim,'Claim differs from exact audited source bindings'
    print('CALIBRATED_STATIC_CLAIM_PASS',claim['verified_manifest_leaf_references'])


if __name__=='__main__':main()
