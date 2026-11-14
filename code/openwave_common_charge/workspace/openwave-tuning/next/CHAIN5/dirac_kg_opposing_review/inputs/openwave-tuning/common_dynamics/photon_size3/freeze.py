"""Freeze completed size-three producer without changing provider bytes."""
from pathlib import Path
import json,sys
sys.dont_write_bytecode=True
import providers
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
def run():
    pins=providers.verify()
    native=json.loads((HERE/'NATIVE.json').read_text());cert=json.loads((HERE/'CERTIFICATE.json').read_text())
    assert all(native['checks'].values()) and all(cert['checks'].values())
    assert native['source_sha256']==pins and cert['source_sha256']==pins
    assert native['program_sha256']==providers.sha(HERE/'PROGRAM.json')
    assert cert['frame_sha256']==providers.sha(HERE/'FRAME.json')
    assert cert['geometry_sha256']==providers.sha(HERE/'GEOMETRY.json')
    names=['providers.py','frame.py','build_frame.py','walsh.py','native_member.py','freeze.py','PROOF.md','FRAME.json','GEOMETRY.json','CERTIFICATE.json','NATIVE.json','PROGRAM.json','NATIVE_RUN.log','RUN_ATTEMPTS.json']
    out={'schema':'oph.quadratic.size3.producer_freeze.v1','worker':'DYN21','candidate_id':native['candidate_id'],
        'files':{str((HERE/n).relative_to(ROOT)):providers.sha(HERE/n) for n in names},
        'source_pins':pins,'checks':native['checks'],'independent_review_required':True,
        'scope':'Actual size3 coexact free source, all646 modes, cap1 spectral code, both field quadratures and shared clock, complete native records. Physical field limit is a separate ROW29 packet.',
        'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')
    print(providers.sha(HERE/'MANIFEST.json'),flush=True)
if __name__=='__main__':run()
