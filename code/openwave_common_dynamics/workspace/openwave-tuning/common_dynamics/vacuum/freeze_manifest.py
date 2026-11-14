"""Verify the reference package and its immutable producer custody."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    sources={};checks={}
    def admit(path,expected):
        key=str(path.relative_to(ROOT))
        assert path.is_file() and digest(path)==expected,(key,expected)
        if key in sources:assert sources[key]==expected
        sources[key]=expected
    for name in ('spin','vector_field','quantum_field','units'):
        path=HERE.parent/name/'MANIFEST.json';manifest=json.loads(path.read_text())
        admit(path,digest(path))
        for section in ('files','source_pins','source_sha256'):
            for relative,expected in manifest.get(section,{}).items():
                path2=ROOT/relative
                if not path2.is_file():path2=path.parent/relative
                admit(path2,expected)
    receipts=['REFERENCE.json','NATIVE_PREPARATION.json','EXCITATIONS.json','ENERGY_RECORDS.json','ELECTRIC_CODE.json']
    for name in receipts:
        receipt=json.loads((HERE/name).read_text())
        assert receipt['checks'] and all(receipt['checks'].values()),name
        checks[name]=len(receipt['checks'])
        for relative,expected in receipt['source_sha256'].items():admit(ROOT/relative,expected)
    files={str(p.relative_to(ROOT)):digest(p) for p in HERE.iterdir()
           if p.is_file() and p.name!='MANIFEST.json'}
    result={'schema':'oph.dirac.shared_reference.manifest.v1','worker':'DYN14',
      'files':files,'source_pins':sources,'receipt_checks':checks,
      'checks':{'all_five_receipts_pass':True,'all_frozen_provider_pins_match':True,
                'complete_native_preparation_and_code_readout_words':True},
      'interfaces':{'free_reference':'reference.complete_occupation_law; reference.reference_charge; reference.dressed_field',
                    'native_preparation':'native_preparation.prepare_blocks; native_preparation.word',
                    'field_code_test':'electric_code.encoder_word; ELECTRIC_CODE.json',
                    'certified_field_expectation':'energy_records.electric_interval'},
      'allocation':'32Fock+1152field+3scratch=1187source/work qubits; pointeraux1187/1188; no clock/probe',
      'scope':['Explicit CAR lift of the routed Dirac source and its selected principal-phase free reference.',
               'One sharedGauss field, full occupation laws, native Gaussian/Slater preparation and complete field-code measurement.',
               'Generic electric evolution leaves the uniform reference code; no interacting vacuum/species/pole or fullFock native evolution is certified.',
               'All low/optical/Floquet valleys and reference-charge domains retained.'],
      'independent_review_pending':True,'reserved_acceptance_runs':False,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':digest(HERE/'MANIFEST.json'),'producer_checks':sum(checks.values()),
                      'source_pins':len(sources),'files':len(files)}))


if __name__=='__main__':run()
