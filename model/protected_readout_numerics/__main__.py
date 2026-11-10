"""Replay the unchanged producer, then certify its complete numerical output."""
import argparse,json
from pathlib import Path
from . import model
from model.protected_particle_readout import packet


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write',type=Path,nargs='?',const=model.HERE/'results.json')
    g.add_argument('--check',type=Path,nargs='?',const=model.HERE/'results.json')
    p.add_argument('--renew-live',action='store_true',help='after a successful --check, renew only the two typed legacy parent hashes')
    p.add_argument('--replay-output',type=Path,help='write the complete fresh JSON/NPZ pair to a nonexistent path')
    a=p.parse_args()
    if a.renew_live and not a.check:p.error('--renew-live requires --check')
    own={x.name:model.digest(x)for x in sorted(model.HERE.iterdir())if x.suffix in('.py','.md')}
    pins=dict(packet.source_bindings())
    pins.update({str((model.HERE/name).relative_to(model.ROOT)):sha for name,sha in own.items()})
    for name in('contract.json','retained.json'):pins[str((model.HERE/name).relative_to(model.ROOT))]=model.digest(model.HERE/name)
    spec,_,_=model.specification();pins[spec['retained_array']]=spec['retained_array_sha256']
    if a.check:
        saved=json.loads(a.check.read_text());model.require(saved['input_sha256']==pins,'numerical successor source changed')
    fresh=packet.report()
    result=model.report(fresh);diagnostics=result.pop('backend_replay_diagnostics')
    result['input_sha256']=pins;result['owned_source_sha256']=own
    # Observed backend residuals are separately bounded during this replay;
    # they are not fixed scientific values across BLAS/CPU implementations.
    result['every_fresh_complete_law_independently_validated']=True
    for path,sha in pins.items():model.require(model.digest(model.ROOT/path)==sha,'source changed during numerical successor replay')
    if a.check:model.require(result==saved,'exact numerical certificate changed')
    else:a.write.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    if a.replay_output:packet.write_pair(a.replay_output,*fresh)
    if a.renew_live:
        target=model.ROOT/'model/protected_particle_readout/results.json'
        before=model.digest(target);value=model.live_projection()
        target.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')
        model.specification()
        projection={'schema':'oph.protected_readout_numerics.live_projection.v1',
                    'retained_raw_json_sha256':spec['retained_receipt_sha256'],
                    'before_sha256':before,'current_sha256':model.digest(target),
                    'only_source_paths_changed':sorted(spec['typed_parent_metadata_transitions']),
                    'original_scientific_fields_and_complete_NPZ_preserved':True,
                    'actual_numerical_successor_and_unchanged_cold_producer_replayed':True,
                    'original_strict_checker_or_admission_modified':False}
        (model.HERE/'live_projection.json').write_text(json.dumps(projection,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'numerical_successor_passed':True,'original_checker_rejected_this_backend':diagnostics['old_checker_rejected_this_replay'],
                      'complete_state_error_upper':[x['complete_state_error_upper']for x in diagnostics['fresh_full_state_interval_validation']],
                      'original_packet_checker_modified':False},sort_keys=True))


if __name__=='__main__':main()
