"""Exercise independent scientific predicates with real retained case records."""
import argparse,copy,hashlib,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
import audit_records as a

HERE=Path(__file__).resolve().parent
CASES={
 'current': 'Marker literal current',
 'coherence': 'Full actual arm density/crossblock',
 'missing_flag': 'Complete target flags',
 'padding': 'Which-arm padding budget',
 'selected_cost': 'Full selected clock cost',
 'ideal_ray': 'No ideal ray substitution',
 'tiny_positive': 'Exact positive branch probability lost exact zero/positive',
}

def child(folder,kind):
    prep=a.load(folder/'PREPARATION.json.gz')['1']['native'];rho,success=a.prepare(prep)
    cl=a.clock(0,1,a.load(folder/'INPUT_s1_k0.json.gz'),rho)
    law='swap'if kind=='tiny_positive'else'native'
    record=a.load(folder/f'CASE_s1_k0_r1_t1_d1_{law}.json.gz')
    catalog=a.load(folder/'CATALOG.json.gz')['1'][law]
    if kind=='current':record['marker_records'][0]['histories'][0]['currents'][0]+=1
    elif kind=='coherence':
        f=a.array(record['target_first_failures'][0]['input_factor']);f[3:]*=1j
        record['target_first_failures'][0]['input_factor']={'real':f.real.tolist(),'imag':f.imag.tolist()}
    elif kind=='missing_flag':record['target_first_failures'][0]['probabilities'].pop()
    elif kind=='padding':record['target_first_failures'][0]['ref_controlled_idle_R_slots']+=1
    elif kind=='selected_cost':record['selected_clock_probability_given_entry']*=.5
    elif kind=='ideal_ray':record['rotation_input_factor']={'real':[[1.,0.],[0.,0.]],'imag':[[0.,0.],[0.,0.]]}
    elif kind=='tiny_positive':record['returned']['path_probabilities_plus_minus'][1]=0.
    elif kind!='clean':raise RuntimeError('Unknown mutation')
    a.audit_case(record,cl,catalog,success*cl['success'])
    print(json.dumps({'clean':True}))

def run(folder):
    index=a.load(folder/'INDEX.json')
    paths=[folder/'INDEX.json']+[folder/r['path']for r in index['files']]
    before={str(p):(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns,p.stat().st_ino)for p in paths}
    rows=[]
    for optimized in(False,True):
        for mutation in ['clean',*CASES]:
            command=[sys.executable,'-B']+(['-O']if optimized else[])+[str(Path(__file__).resolve()),str(folder.resolve()),'--child',mutation]
            result=subprocess.run(command,capture_output=True,text=True,env={**os.environ,'PYTHONOPTIMIZE':'0','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
            if mutation=='clean':a.need(result.returncode==0,'Clean scientific case rejected')
            else:a.need(result.returncode!=0 and CASES[mutation]in result.stderr,'Scientific mutant wrong rejection cause: '+mutation+'\n'+result.stderr)
            rows.append({'mutation':mutation,'optimized':optimized,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr,'expected_predicate':CASES.get(mutation)})
    after={str(p):(hashlib.sha256(p.read_bytes()).hexdigest(),p.stat().st_mtime_ns,p.stat().st_ino)for p in paths}
    a.need(before==after,'Saved producer bytes/metadata changed')
    return {'passed':True,'same_group':True,'external_G6':False,'scientific_mutants':len(CASES),'process_checks':len(rows),'producer_artifacts_unchanged':len(paths),'index_sha256':before[str(folder/'INDEX.json')][0],'rows':rows}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('folder',type=Path);p.add_argument('--child');p.add_argument('--record-to',type=Path);p.add_argument('--verify',type=Path,default=HERE/'MUTATIONS.json');args=p.parse_args()
    if args.child:child(args.folder,args.child)
    else:
        if args.record_to:a.need(not args.record_to.exists(),'Exclusive mutation receipt exists')
        start=time.perf_counter();result=run(args.folder)
        if args.record_to:
            with args.record_to.open('x')as f:json.dump(result,f,indent=2);f.write('\n')
        else:a.need(result==a.load(args.verify),'Saved mutation verdict differs')
        print(json.dumps({'passed':True,'processes':result['process_checks'],'seconds':time.perf_counter()-start}))
