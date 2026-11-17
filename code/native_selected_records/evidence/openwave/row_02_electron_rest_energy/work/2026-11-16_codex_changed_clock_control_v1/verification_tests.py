"""Post-run verification diagnostics, separate from the alternative-law evidence."""
import gzip,hashlib,json,os,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def signature():
    return {p.name:(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns,p.stat().st_ino) for p in HERE.iterdir() if p.is_file() and p.name in ('HEADLINE.json.gz','EVIDENCE.json.gz','INTERVALS.json','HEADLINE_INTERVALS.json','PREREG.json')}

def main():
    started=time.monotonic();before=signature();rows=[]
    with tempfile.TemporaryDirectory(prefix='changed-clock-contract-') as temp:
        temp=Path(temp);saved=json.loads(gzip.decompress((HERE/'HEADLINE.json.gz').read_bytes()));saved['cases'][0]['visibility']+=.001
        damaged=temp/'damaged.json.gz';damaged.write_bytes(gzip.compress(json.dumps(saved).encode()))
        cases=[('optimized_read_only',[sys.executable,'-B','-O','experiment.py'],0,'"passed": true'),
          ('changed_saved_evidence',[sys.executable,'-B','-O','experiment.py','--verify',str(damaged)],1,'saved changed-clock evidence differs'),
          ('exclusive_first_record',[sys.executable,'-B','-O','experiment.py','--record-to',str(HERE/'HEADLINE.json.gz')],1,'refusing evidence overwrite'),
          ('actual_native_frame_sign_plant',[sys.executable,'-B','-O','-c',"import experiment as e; old=e.p.run_case\ndef changed(*a,**k):\n r=old(*a,**k);r['visibility']=-r['visibility'];return r\ne.p.run_case=changed;e.run([0])"],1,'same native frame sign/return gate'),
          ('physically_ineffective_changed_clock',[sys.executable,'-B','-O','-c',"import experiment as e; old=e.source_join.event\ne.source_join.event=lambda q,edge,coin,law:old(q,edge,coin,'native')\ne.run([0])"],1,'changed-clock detector did not reject native prediction')]
        for label,command,status,needle in cases:
            t=time.monotonic();r=subprocess.run(command,cwd=HERE,text=True,capture_output=True,timeout=90);passed=(r.returncode==0 if status==0 else r.returncode!=0) and needle in r.stdout+r.stderr
            rows.append({'name':label,'command':command,'exit_code':r.returncode,'expected':status,'passed':passed,'seconds':time.monotonic()-t,'stdout':r.stdout,'stderr':r.stderr})
            if not passed:raise RuntimeError(json.dumps(rows,indent=2))
        for label,code,needle in [('invalid_observed_coin',"d['source_join']['pairs'][0]['sources'][0]['events'][0][1][2]=2",'actual coin domain'),('dropped_joint_column',"d['source_join']['pairs'][0]['joint_columns'].pop()",'joint column count'),('changed_counter_bin',"d['clocks']['native:-1:0']['input']['counter_probabilities'][0]+=.01",'full initial counter bins')]:
            code="import independent_check as a,gzip,json;d=json.loads(gzip.decompress((a.HERE/'HEADLINE.json.gz').read_bytes()));"+code+";a.verify(d)"
            t=time.monotonic();r=subprocess.run([sys.executable,'-B','-O','-c',code],cwd=HERE,text=True,capture_output=True,timeout=60);passed=r.returncode!=0 and needle in r.stderr
            rows.append({'name':label,'exit_code':r.returncode,'passed':passed,'seconds':time.monotonic()-t,'stdout':r.stdout,'stderr':r.stderr})
            if not passed:raise RuntimeError(json.dumps(rows,indent=2))
    if signature()!=before:raise RuntimeError('Saved source/evidence mutated')
    return {'passed':True,'tests':rows,'saved_bytes_mtime_inode_unchanged':True,'science_scope':'Detector and execution diagnostics only; changed-clock science is the separately frozen source experiment.','seconds':time.monotonic()-started}

if __name__=='__main__':print(json.dumps(main(),indent=2))
