"""Check receipt and frozen dependency custody, then freeze this package."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    pins={}
    def add(relative,digest):
        if relative in pins:assert pins[relative]==digest,relative
        assert sha(ROOT/relative)==digest,relative
        pins[relative]=digest
    for name in ('clock','continuum','manybody','energy'):
        path=HERE.parent/name/'MANIFEST.json'
        data=json.loads(path.read_text())
        for key in ('source_sha256','files','source_pins'):
            for relative,digest in data.get(key,{}).items():add(relative,digest)
        for relative,digest in data.get('files_sha256',{}).items():
            add(str((path.parent/relative).relative_to(ROOT)),digest)
        add(str(path.relative_to(ROOT)),sha(path))
    checks={}
    for name in ('ROTATION_EXCHANGE','DIRAC_SPIN','DIRAC_NATIVE','CYCLE_PROBE','ROUTE_CUTOFF'):
        data=json.loads((HERE/(name+'.json')).read_text())
        for relative,digest in data['source_sha256'].items():add(relative,digest)
        if 'checks' in data:
            assert all(value is True for value in data['checks'].values())
            checks[name]=len(data['checks'])
        assert data['official_rows_changed']==[]
    # DYN08 confirmed these individual files frozen before its final manifest.
    for name in ('autonomous.py','cell_complex.py','CELL_COMPLEX.json','CELL_COMPLEX_PROOF.md',
                 'linear_response.py','LINEAR_RESPONSE.json'):
        path=HERE.parent/'vector_field'/name
        add(str(path.relative_to(ROOT)),sha(path))
    for path in sorted(HERE.iterdir()):
        if path.is_file() and path.name!='MANIFEST.json':add(str(path.relative_to(ROOT)),sha(path))
    result=dict(schema='oph.common.local12.spin.freeze.v1',worker='DYN10',
                base_candidate='native-local12-current-v0-proposed',
                proposed_replacement='native-local12-dirac-routed-v3-proposed',
                checks={'all_producer_checks_pass':True,'all_frozen_source_pins_match':True},
                producer_check_counts=checks,source_sha256=dict(sorted(pins.items())),
                source_pin_count=len(pins),
                scope='Geometric spin obstruction for the scalar source; explicit Dirac architecture, free limit, routed modular current, native spin and cycle probes, and factorized native endpoint compiler. Full physical rotation apparatus, native full-field initialization/reset, locally Gauss-preserving gate-prefix implementation, aggregate physical errors, interacting canonical-field limit and species control remain obligations.',
                official_rows_changed=[])
    path=HERE/'MANIFEST.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'manifest_sha256':sha(path),'source_pins':len(pins),'producer_checks':checks}))


if __name__=='__main__':run()
