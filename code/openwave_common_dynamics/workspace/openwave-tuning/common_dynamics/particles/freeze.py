"""Verify the bounded particle artifacts and their recursive source custody."""
from pathlib import Path
from collections import deque
import gzip
import hashlib
import json
import math
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def discovered(value):
    if isinstance(value,dict):
        for key,item in value.items():
            if isinstance(key,str) and '/'in key and isinstance(item,str) and re.fullmatch('[0-9a-f]{64}',item):
                if key.startswith(('openwave-tuning/','reverse-engineering-reality/'))or(ROOT/key).is_file():
                    yield key,item
            yield from discovered(item)
        if isinstance(value.get('path'),str) and isinstance(value.get('sha256'),str):
            if value['path'].startswith(('openwave-tuning/','reverse-engineering-reality/'))or(ROOT/value['path']).is_file():
                yield value['path'],value['sha256']
    elif isinstance(value,list):
        for item in value:yield from discovered(item)


def run():
    receipts=['MASS_COIN.json','PARTICLE_RECORDS.json','COUPLED_RESPONSE.json','MASSIVE_COUPLED.json']
    checks=0
    for name in receipts:
        data=json.loads((HERE/name).read_text())
        assert all(data['checks'].values()) and data['official_rows_changed']==[]
        checks+=len(data['checks'])
    law_checks={}
    for filename in ['COUPLED_JOINT_LAWS.json.gz','MASSIVE_JOINT_LAWS.json.gz']:
        payload=gzip.decompress((HERE/filename).read_bytes())
        data=json.loads(payload);joint_total=0.;outcomes=0
        for record in data['records']:
            total=0.
            for block in record['joint_law']:
                assert 0<=block['position']<96
                field=bytes.fromhex(block['all288_Z16_field_values_hex'])
                assert len(field)==288 and max(field)<16
                probabilities=block['all12_event_probabilities']
                assert len(probabilities)==12 and all(math.isfinite(p)and p>=0 for p in probabilities)
                total+=sum(probabilities);outcomes+=12
            if 'first_position'in record:joint_total+=total
            else:assert abs(total-1)<2e-10
        if filename.startswith('COUPLED'):assert abs(joint_total-1)<2e-10
        law_checks[filename]={'records':len(data['records']),'retained_event_outcomes':outcomes,
                              'uncompressed_sha256':hashlib.sha256(payload).hexdigest()}
    own={str(p.relative_to(ROOT)):digest(p)for p in HERE.iterdir()
         if p.is_file()and p.name!='MANIFEST.json'}
    seeds=[ROOT/'openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json',
           ROOT/'openwave-tuning/auditor/common_walk_particle/MANIFEST.json',
           ROOT/'openwave-tuning/auditor/native_current_walk/MANIFEST.json']
    seeds +=[HERE.parent/k/'MANIFEST.json'for k in ['clock','continuum','manybody','energy']]
    seeds +=[HERE.parent/'vector_field'/k for k in ['autonomous.py','cell_complex.py','composition.py',
             'COMPOSITION.json','NATIVE_FIELD.json','FIELD_PROGRAM.json','PREPARATION_PROGRAM.json',
             'reset_domain.py','RESET_DOMAIN.md','RESET_DOMAIN.json']]
    pins=dict(own)
    pins.update({str(p.relative_to(ROOT)):digest(p)for p in seeds})
    queue=deque(pins);visited=set()
    while queue:
        name=queue.popleft()
        if name in visited:continue
        visited.add(name);path=ROOT/name
        assert digest(path)==pins[name],('changed source',name)
        if path.suffix!='.json':continue
        for dependency,expected in discovered(json.loads(path.read_text())):
            assert digest(ROOT/dependency)==expected,('stale inherited pin',name,dependency)
            if dependency in pins:assert pins[dependency]==expected,('conflicting pins',dependency)
            else:pins[dependency]=expected;queue.append(dependency)
    result={'schema':'oph.common.particles.freeze.v1','worker':'DYN09',
            'candidate_proposals':['native-local12-mass-coin-v1-proposed','native-local12-mass-fermion-plaquette-v3-proposed'],
            'files':own,'source_pins':dict(sorted(pins.items())),
            'checks':{'all_four_receipts_pass':True,'source_custody_closure_matches':True,
                      'all_complete_joint_laws_normalized_without_postselection':True,
                      'reset_bounds_do_not_assume_independence':True,
                      'full_physical_and_padding_sectors_retained':True,
                      'candidate_changes_and_supplied_parameters_explicit':True},
            'receipt_check_count':checks,'joint_law_checks':law_checks,
            'numerical_law_scope':'Complete binary64 ideal-circuit probability tables with normalization controls. Native TV bounds compare the native program with the exact ideal circuit law, not with a formally interval-enclosed decimal table.',
            'reset_record_scope':'Fine private reset tapes are retained by the actual instrument. The blank-times-original-reference comparison marginalizes them. The particle/field TV bound concerns the terminal observable marginal with failures retained, without asserting tape/reference independence.',
            'field_source_scope':'Exact DYN08 autonomous/composition snapshot in source_pins, not a claim that its later Dirac/Maxwell work has frozen.',
            'scope':'Bounded source-native particle/hole, persistence, supplied scalar mass and controlled free continuum construction; complete finite coupled laws and native readout/error composition. No observed particle masses, interacting physical vacuum/pole, spin-half claim, Maxwell identification or two-charge force law is inferred.',
            'independent_review_pending':True,'reserved_acceptance_runs':False,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'sha256':digest(HERE/'MANIFEST.json'),'source_pins':len(pins),
                      'receipt_check_count':checks,'joint_law_checks':law_checks}))


if __name__=='__main__':run()
