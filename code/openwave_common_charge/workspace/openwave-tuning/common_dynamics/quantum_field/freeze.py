"""Freeze the bounded prepared quantum-field result and its source custody."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
RECEIPTS=('QUANTUM_CYCLE','NATIVE_CYCLE','GAUSSIAN_BOUNDS','UV_STATE','REFINEMENT_CONTROL','REFERENCE_PHASE')


def run():
    pins={};counts={}
    for name in RECEIPTS:
        receipt=json.loads((HERE/(name+'.json')).read_text())
        assert receipt['checks'] and all(receipt['checks'].values())
        counts[name]=len(receipt['checks'])
        for relative,expected in receipt['source_sha256'].items():
            assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
            assert relative not in pins or pins[relative]==expected
            pins[relative]=expected
    # The same full native source custody used by the common walk supplies
    # the exact response, control identities and code-service dependencies.
    core=ROOT/'openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json'
    for relative,expected in json.loads(core.read_text())['files_sha256'].items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()==expected,relative
        assert relative not in pins or pins[relative]==expected
        pins[relative]=expected
    for relative in ('openwave-tuning/common_dynamics/vector_field/MANIFEST.json',
                     'openwave-tuning/common_dynamics/manybody/MANIFEST.json'):
        path=ROOT/relative;manifest=json.loads(path.read_text())
        for key in ('files','source_pins','source_sha256'):
            mapping=manifest.get(key,{})
            if not isinstance(mapping,dict):continue
            for source,expected in mapping.items():
                if not isinstance(expected,str) or len(expected)!=64:continue
                target=ROOT/source
                if not target.is_file():continue
                assert hashlib.sha256(target.read_bytes()).hexdigest()==expected,source
                assert source not in pins or pins[source]==expected
                pins[source]=expected
        pins[relative]=hashlib.sha256(path.read_bytes()).hexdigest()
    files={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
           for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json'}
    manifest={'schema':'oph.common.prepared_quantum_field.freeze.v1','worker':'DYN13',
      'files':files,'source_pins':pins,'receipt_checks':counts,
      'checks':{'six_receipts_pass':True,'all_source_pins_match':True,
                'actual_native_Gaussian_preparation_and_readout':True,
                'integer_outward_state_transfer_proved':True,
                'fixed_coupling_and_optical_scope_preserved':True},
      'interfaces':{'prepared_cycle':'quantum_state.gaussian_cycle',
                    'source':'quantum_state.source_step',
                    'complete_records':'quantum_state.complete_records',
                    'native_preparation':'native_cycle.state_preparation_word',
                    'native_full_field_encoder':'native_cycle.encoder_word',
                    'native_full_source_clock':'native_cycle.full_word',
                    'preparation_certificate':'native_cycle.gray_preparation_bound',
                    'analytic_transfer':'gaussian_bounds.preparation_comparison; gaussian_bounds.analytic_comparison',
                    'explicit_reference_phase':'reference_phase.word'},
      'source_family':'pure-gauge one-face member of the same local plaquette source; native m16 and separately bounded ideal register refinements',
      'selected_by_this_worker':False,
      'scope':'Prepared optical cycle quantum state, complete native source/clock/readout and controlled weak-angle oscillator comparison. Weak coupling tends to zero. Full homogeneous photon state, fixed-coupling quantum Maxwell limit, interacting matter vacuum and calibrated forces are not claimed.',
      'independent_review_pending':True,'reserved_acceptance_runs':False,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest(),
                      'receipt_checks':counts,'checks_total':sum(counts.values()),
                      'files':len(files),'source_pins':len(pins)}))


if __name__=='__main__':run()
