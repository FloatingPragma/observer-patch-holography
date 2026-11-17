"""Verify each compiler table equals actual records on every admitted background."""
import gzip,json
from pathlib import Path
from fractions import Fraction as F
HERE=Path(__file__).resolve().parent

def require(ok,why):
    if not ok:raise ValueError(why)

def run():
    data=json.loads(gzip.decompress((HERE/'EVIDENCE.json.gz').read_bytes()));join=data['source_join'];n=0
    for pair in join['pairs']:
        for source in pair['sources']:
            obs=source['observer'];sign=obs['admission']['sign'];rep=join['representatives'][f"{pair['clock_law']}:{sign}"]
            require(rep['prep_columns']==obs['prep']['columns'] and rep['clock_columns']==obs['clock']['columns'],'representative detached from actual probe histories')
            require(source['public_seeds'][0] in (72365001,72365002,72365011,72365012),'public admission stream');n+=1
    for side in (0,1):
        for stream in (0,1):
            chosen=[r['sources'][side] for r in join['pairs'] if r['public_stream_pair']==stream]
            require(len({tuple(s['public_seeds']) for s in chosen})==1,'stream encodes hidden input/law')
    require(n==64 and F(join['per_source_timeout_upper'])==F(1816,2**64),'full continuing domain/timeout')
    return {'passed':True,'actual_sources':n,'observed_table_bindings':2*n,'public_streams_independent_of_hidden_input_and_law':True,'verification_only':True}

if __name__=='__main__':print(json.dumps(run(),indent=2))
