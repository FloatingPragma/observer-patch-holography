"""Post-run execution contract diagnostics for the outward companion."""
import hashlib,json,subprocess,sys,tempfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def snapshot():
    return {n:(hashlib.sha256((HERE/n).read_bytes()).hexdigest(),(HERE/n).stat().st_mtime_ns,(HERE/n).stat().st_ino) for n in ('INTERVALS.json','EVIDENCE.json.gz','PREREG.json')}

def main():
    before=snapshot();rows=[]
    with tempfile.TemporaryDirectory(prefix='changed-clock-interval-') as directory:
        directory=Path(directory);data=json.loads((HERE/'INTERVALS.json').read_text());data['rows'][0]['same_anchor']=False;(directory/'INTERVALS.json').write_text(json.dumps(data))
        altered=f"import interval_certificate as c;from pathlib import Path;c.HERE=Path({str(directory)!r});r=c.run(Path({str(HERE/'EVIDENCE.json.gz')!r}));s=c.json.loads((c.HERE/'INTERVALS.json').read_text());r.pop('seconds');s.pop('seconds');c.require(r==s,'saved interval differs')"
        sign="import interval_certificate as c;old=c.returned\ndef changed(*a,**k):\n v,p,w,t=old(*a,**k);return -v,p,w,t\nc.returned=changed;c.run(c.HERE/'HEADLINE.json.gz')"
        for label,args,needle,status in [('optimized_verify',['interval_certificate.py'],'"passed": true',0),('changed_saved_interval',['-c',altered],'saved interval differs',1),('actual_interval_sign_gate',['-c',sign],'unchanged frame sign/return',1),('exclusive_interval_first_record',['interval_certificate.py','--record-to','INTERVALS.json'],'refusing interval overwrite',1)]:
            start=time.monotonic();r=subprocess.run([sys.executable,'-B','-O',*args],cwd=HERE,text=True,capture_output=True,timeout=90);passed=(r.returncode==0 if status==0 else r.returncode!=0) and needle in r.stdout+r.stderr
            rows.append({'name':label,'exit_code':r.returncode,'passed':passed,'seconds':time.monotonic()-start,'stdout':r.stdout,'stderr':r.stderr})
            if not passed:raise ValueError(json.dumps(rows,indent=2))
    if before!=snapshot():raise ValueError('Evidence changed')
    return {'passed':True,'verification_only':True,'saved_bytes_mtime_inode_unchanged':True,'tests':rows}

if __name__=='__main__':print(json.dumps(main(),indent=2))
