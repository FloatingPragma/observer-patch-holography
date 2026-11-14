"""Freeze a full Maxwell candidate for independent criterion review."""
from pathlib import Path
import hashlib,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'openwave-tuning/common_dynamics/photon_size3'))
import providers

def run():
    pins=providers.verify()
    size=ROOT/'openwave-tuning/common_dynamics/photon_size3/MANIFEST.json'
    assert providers.sha(size)=='90178885fcaaaba405b950a6f132fd47d519d21a130a4f982ec09423db129ec2'
    sm=json.loads(size.read_text())
    for path,digest in sm['files'].items():
        assert providers.sha(ROOT/path)==digest,path
        assert path not in pins or pins[path]==digest
        pins[path]=digest
    pins[str(size.relative_to(ROOT))]=providers.sha(size)
    reports={name:json.loads((HERE/name).read_text()) for name in ('MAXWELL.json','RADIATIVE.json','SOURCE_CONTRACT.json')}
    for name,data in reports.items():assert all(data['checks'].values()),name
    assert reports['RADIATIVE.json']['program_sha256']==providers.sha(HERE/'RADIATIVE_PROGRAM.json')
    criterion=providers.sha(HERE/'CRITERION.json')
    candidate={'schema':'oph.openwave.full_criterion_candidate.v1','row':29,
        'candidate_id':'native-local12-dirac-noncompact-quadratic-v4-proposed',
        'literal_criterion':'Maxwell recovered; transverse waves at c','criterion_sha256':criterion,
        'producer_claim_full_criterion':True,'independent_review_status':'pending; no self-approval',
        'construction_inputs':['actual translated LOCAL12 oriented faces and incidence','noncompact quadratic Hamiltonian with a=1/2,b=3','spectral-code preparation and terminal native apparatus','common h/8 phase clock'],
        'derived':['exact oriented curl factorization with rho=2/3','Faraday and Ampere vacuum equations with physical field normalization','two transverse polarizations in every direction','controlled finite periodic physical-field and frequency limit','complete native field/quadrature and common clock records'],
        'ancestry':'Speed matching is an explicit coefficient choice; neither empirical coupling nor uniquely forced architecture is claimed. No continuum Maxwell stencil enters the microscopic face matrix.',
        'scope':'Vacuum radiative Maxwell fields and their first-moment instruments on the declared spectral-code source; fixed-native-precision TV floor retained.',
        'max_refinement_TV':[r['complete_native_source_TV_upper'] for r in reports['RADIATIVE.json']['refinements']],
        'official_rows_changed':[]}
    (HERE/'CANDIDATE.json').write_text(json.dumps(candidate,indent=2)+'\n')
    names=['CRITERION.json','TASK.md','maxwell.py','radiative_limit.py','source_contract.py','freeze.py','PROOF.md','MAXWELL.json','RADIATIVE.json','SOURCE_CONTRACT.json','RADIATIVE_PROGRAM.json','MAXWELL_RUN.log','RADIATIVE_RUN.log','SOURCE_RUN.log','ATTEMPTS.json','maxwell_attempt1.py','MAXWELL_ATTEMPT1.log','CANDIDATE.json']
    out={'schema':'oph.openwave.row29.producer_freeze.v1','worker':'DYN21','row':29,
        'candidate_id':candidate['candidate_id'],'criterion_sha256':criterion,
        'files':{str((HERE/n).relative_to(ROOT)):providers.sha(HERE/n) for n in names},
        'source_pins':pins,'checks':{name:data['checks'] for name,data in reports.items()},
        'producer_full_criterion_request':True,'independent_reviews':[],
        'acceptance':'Two producer-independent full-criterion reviews are required; this producer does not write VALIDATION.json.',
        'official_rows_changed':[]}
    (HERE/'MANIFEST.json').write_text(json.dumps(out,indent=2)+'\n')
    print(providers.sha(HERE/'MANIFEST.json'),flush=True)
if __name__=='__main__':run()
