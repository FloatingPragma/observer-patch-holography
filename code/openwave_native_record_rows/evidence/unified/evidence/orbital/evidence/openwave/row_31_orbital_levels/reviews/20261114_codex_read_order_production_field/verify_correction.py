#!/usr/bin/env python3
"""Keep original failure, verify corrected retry and exact native ancestry."""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
HERE=Path(__file__).resolve().parent
P=HERE.parent.parent/'work/2026-11-14_codex_production_field'

def main():
    old=json.loads(gzip.decompress((P/'failed_initial/RECORDS.json.gz').read_bytes()))
    new=json.loads(gzip.decompress((P/'RECORDS.json.gz').read_bytes()))
    assert hashlib.sha256((P/'failed_initial/RECEIPT.json').read_bytes()).hexdigest()=='7df9e4455e35d23643a11201ad87f755498662c83e3f3b71f3bed734e7696be8'
    assert hashlib.sha256((P/'failed_initial/RECORDS.json.gz').read_bytes()).hexdigest()=='342b68a0cb4a3bb15b7fec5621d8c960a8fe00f2e91041938302c042c6dde9b0'
    assert hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest()=='36d2829573ff06ac8145d0753917f65b48b9c5274aebcd663df933d270bc028a'
    assert hashlib.sha256((P/'RECORDS.json.gz').read_bytes()).hexdigest()=='f36f7e8e7b1c179632fb70fdc17de0571431ee423842e4e5f9705fc49121a001'
    corrected=json.loads((HERE/'CORRECT_BAND_READOUTS.json').read_text())['bands'];cursor=0;maxerr=0.
    for i,(orun,nrun) in enumerate(zip(old['runs'],new['runs'])):
        assert orun['field_strength_readouts']==nrun['field_strength_readouts']
        for oc,nc in zip(orun['channels'],nrun['channels']):
            maxerr=max(maxerr,float(np.max(abs(np.array(nc['band_readouts'])-corrected[cursor]))))
            assert maxerr<1e-14
            assert np.max(abs(np.array(oc['band_readouts'])-corrected[cursor]))>1e-5
            # Delete only the corrected derived readout; every native and all other derived field records stay identical.
            oc.pop('band_readouts');nc.pop('band_readouts');assert oc==nc;cursor+=1
    assert old==new and cursor==8
    o=json.loads((P/'failed_initial/RECEIPT.json').read_text());n=json.loads((P/'RECEIPT.json').read_text())
    o.pop('raw_uncompressed_sha256');n.pop('raw_uncompressed_sha256');assert o==n
    raw=gzip.decompress((P/'RECORDS.json.gz').read_bytes());nr=json.loads((P/'RECEIPT.json').read_text())
    assert hashlib.sha256(raw).hexdigest()==nr['raw_uncompressed_sha256']
    result={'pass':True,'same_group':True,'G6_pass':False,'initial_decoder_failure_retained':True,
        'all_eight_corrected_band_readouts_independently_verified':True,'maximum_corrected_readout_error':maxerr,
        'every_native_event_start_occupation_clock_censor_and_other_record_identical':True,
        'curl_and_population_unchanged':True,'original_native_full_replay_transfers_by_exact_identity':True,
        'blocking':0,'corrected_receipt_sha256':hashlib.sha256((P/'RECEIPT.json').read_bytes()).hexdigest(),
        'corrected_records_sha256':hashlib.sha256((P/'RECORDS.json.gz').read_bytes()).hexdigest()}
    (HERE/'CORRECTION_VERIFICATION.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
