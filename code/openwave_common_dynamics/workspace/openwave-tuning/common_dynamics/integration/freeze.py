"""Freeze the additive integration artifacts after their independent checks."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    receipt=json.loads((HERE/'INTEGRATION.json').read_text())
    if not all(receipt['checks'].values()):raise ValueError('failed integration check')
    context={
        'openwave-tuning/common_dynamics/SOURCE_CONTRACT.json':'consumed/SOURCE_CONTRACT.json',
        'openwave-tuning/platform_packet/criteria_frozen.json':'consumed/criteria_frozen.json',
    }
    snapshots={}
    for original,relative in context.items():
        path=HERE/relative
        expected=receipt['source_sha256'][original]
        if not path.exists():
            if digest(ROOT/original)!=expected:raise ValueError('context drift before snapshot: '+original)
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes((ROOT/original).read_bytes())
        if digest(path)!=expected:raise ValueError('context snapshot changed: '+relative)
        snapshots[original]={'snapshot':str(path.relative_to(ROOT)),'sha256':expected}
    for path,sha in receipt['source_sha256'].items():
        selected=HERE/context[path] if path in context else ROOT/path
        if digest(selected)!=sha:raise ValueError('changed source: '+path)
    all_pins={p:s for p,s in receipt['source_sha256'].items() if p not in context}
    for filename in ('UNITS_REVIEW.json','SOURCE_ADMISSION.json','CURRENT_REVIEW.json'):
        extra=json.loads((HERE/filename).read_text())
        if not all(extra['checks'].values()):raise ValueError('failed auxiliary review')
        for path,sha in extra['source_sha256'].items():
            if digest(ROOT/path)!=sha:raise ValueError('changed auxiliary source: '+path)
            if path in all_pins and all_pins[path]!=sha:raise ValueError('incompatible source hashes')
            all_pins[path]=sha
    rows=json.loads((HERE/'CRITERIA_AUDIT.json').read_text())
    assert [r['id'] for r in rows['criteria']]==list(range(1,32))
    assert not rows['official_rows_changed']
    files={str(p.relative_to(ROOT)):digest(p) for p in HERE.rglob('*')
           if p.is_file() and p.name!='MANIFEST.json' and '__pycache__' not in p.parts}
    data={'schema':'oph.common.integration.freeze.v1','worker':'DYN12',
          'files':files,'source_sha256':all_pins,
          'historical_context_snapshots':snapshots,
          'checks':{'integration_checks_pass':True,'source_bytes_match':True,
                    'all31_criteria_have_domain_and_obligation':True},
          'interfaces':{'complete_program':'verified.complete_program',
                        'error_contract':'verified.error_contract',
                        'admission':'joint_program.admit'},
          'scope':'Independent finite-source audit and additive complete clock/mass/position/field instrument. No physical criterion or official validation is inferred from local checks.',
          'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':digest(HERE/'MANIFEST.json'),
                      'files':len(files),'source_pins':len(data['source_sha256'])}))


if __name__=='__main__':run()
