"""Verify transitive frozen inputs and seal the bounded energy package."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def run():
    visited=set()
    pins={}

    def check(path):
        path=path.resolve()
        if path in visited:return
        visited.add(path)
        if path.suffix!='.json':return
        # RER evidence blobs are pinned inputs of the explicit201-file
        # current-source custody. Their historical producer namespaces are
        # not reinterpreted as this package's workspace-relative manifests.
        if ROOT/'reverse-engineering-reality' in path.parents:return
        data=json.loads(path.read_text())
        if not isinstance(data,dict):return
        for field in ('source_sha256','source_pins','files_sha256','files'):
            entries=data.get(field,{})
            if not isinstance(entries,dict):continue
            for relative,expected in entries.items():
                if not isinstance(expected,str) or len(expected)!=64:continue
                source=(ROOT/relative if relative.startswith(('openwave-tuning/','reverse-engineering-reality/','openwave-push/'))
                        else path.parent/relative).resolve()
                actual=hashlib.sha256(source.read_bytes()).hexdigest()
                assert actual==expected,(str(path),relative)
                pins[str(source.relative_to(ROOT))]=actual
                check(source)

    sources=[HERE/'ENERGY_CONTRACT.json',ROOT/'openwave-tuning/auditor/native_walk_bridge/COMMON_SOURCE_CUSTODY.json',
             ROOT/'openwave-tuning/auditor/common_walk_particle/MANIFEST.json',
             ROOT/'openwave-tuning/auditor/native_current_walk/MANIFEST.json',
             HERE.parent/'clock/MANIFEST.json',HERE.parent/'manybody/MANIFEST.json']
    for p in sources:check(p)
    data=json.loads((HERE/'ENERGY_CONTRACT.json').read_text())
    assert all(data['checks'].values())
    own=['energy.py','ENERGY_CONTRACT.json','PROOF.md','README.md','freeze.py']
    for name in own:
        p=HERE/name;pins[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'schema':'oph.common.local12.energy.freeze.v1','worker':'DYN07',
              'candidate_id':data['candidate_id'],'checks':{'all13_producer_checks_pass':True,
                  'current_package_manifests_and201_source_custody_match':True,'source_and_clock_unchanged':True},
              'source_sha256':dict(sorted(pins.items())),
              'scope':'Native base-W quasienergy instruments, conserved D density/current, explicit logarithm/reference/alias controls. Current package manifests and the explicit201-file source custody are checked. Historical RER evidence blobs are pinned as inputs; their producer namespaces are not recursively reinterpreted. Manybody and autonomous field amendments retain separate source and apparatus obligations.',
              'independent_review_pending':True,'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'files':len(pins),'manifest_sha256':hashlib.sha256((HERE/'MANIFEST.json').read_bytes()).hexdigest()}))


if __name__=='__main__':run()
