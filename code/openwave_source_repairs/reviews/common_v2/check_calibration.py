"""Separate calibration map/domain/cost check; no producer imports."""
from fractions import Fraction as F
from pathlib import Path
import argparse,gzip,hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('source',type=Path);ap.add_argument('--record',action='store_true');a=ap.parse_args();p=a.source.resolve();I=p/'inputs'
    assert sha(p/'MANIFEST.sha256')==read(HERE/'CALIBRATION_PREREG.json')['producer_manifest_sha256']
    common=read(I/'common_v1/RECEIPT.json');certificate=read(I/'energy/inputs/certificate.json');transfer=read(I/'calibration/TRANSFER_RESULT.json');r=read(p/'RECEIPT.json')
    old=np.load(I/'calibration/inputs/first_cell_maps.npz');new=np.load(I/'common_v1/inputs/energy/first_cell_certified/maps.npz')
    assert np.array_equal(old['certified_Q_proposal'],new['certified_Q_proposal'])
    assert sha(I/'calibration/native_event.py')==sha(I/'energy/inputs/native_event.py')
    graph=read(I/'extraction/graph.json')['seams'];assert graph==read(I/'energy/inputs/graph.json')['seams']
    edges={frozenset(x)for x in graph};assert edges=={frozenset(x)for x in read(I/'calibration/inputs/graph.json')['seams']}
    records=json.loads(gzip.decompress((I/'calibration/inputs/bridge_RECORDS.json.gz').read_bytes()))
    valid=0
    for trial in records['trials']:
        if trial['selection_law']!='native':continue
        for selection in trial['selection']:
            for row in selection['reference']+[selection['body']]:assert frozenset(row[:2])in edges;valid+=1
        for row in trial['census']:assert frozenset(row[:2])in edges;valid+=1
    routes=read(I/'calibration/ROUTE_RECORDS.json')
    for route in routes['trials']+[routes['forced_timeout']]:
        for leg in route['legs']:
            for row in leg['records']:assert frozenset(row[:2])in edges;valid+=1
    raw=F(transfer['complete_native_error_upper'])+F(transfer['ideal_bad_weight_upper'])+F(common['complete_broad_loader_error_upper']);full=raw/(1-2*raw)
    rounded=F(-((-full.numerator*10**15)//full.denominator),10**15)
    assert rounded==F(r['single_calibration_join']['complete_preparation_and_load_upper'])
    assert full<F(certificate['common_full_preparation_and_apparatus_allowance'])==F(common['spatial_complete_error_upper'])
    result=dict(same_Q_exact=True,same_native_kernel_sha256=sha(I/'calibration/native_event.py'),same_graph=True,valid_recorded_seams_checked=valid,calibration_fullstate_preparation_load_upper=str(rounded),inherited_allowance=str(F(certificate['common_full_preparation_and_apparatus_allowance'])),all_passed=True,producer_imports=False)
    if a.record:
        with(HERE/'CALIBRATION_RESULT.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    else:assert read(HERE/'CALIBRATION_RESULT.json')==result
    print(json.dumps({k:v for k,v in result.items()if k!='inherited_allowance'},sort_keys=True))
if __name__=='__main__':main()
