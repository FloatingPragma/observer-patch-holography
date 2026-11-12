"""Independent native matrix, spatial rotation and complete probability audit."""
import hashlib,json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from model.common_source import native_matter_attachment as matter
from model.common_source.chiral_candidate import kinetic

ROOT=Path(__file__).resolve().parents[3]
def review():
    source=json.loads((ROOT/'model/openwave_lepton_mass_shell/results.json').read_text())
    spin=json.loads((ROOT/'model/openwave_spin_frame_records/results.json').read_text())
    masses=source['complete_lattice_spectrum']['three_distinct_positive_branch_minima_k0_k1_k2']
    rep=matter.representation();ix=list(rep['selected']);weak=rep['operators']['weak_C2'].extract(ix,ix)
    import sympy as s
    L=s.diag(*(int(weak[i,i]!=0)for i in range(15)))
    row=next(r for r in matter.invariant_higgs_lines()if r['neutral_rank']==2)
    D=np.asarray(L*row['components'][1]*(s.eye(15)-L),complex)
    K=np.kron(D,np.diag(masses));K+=K.T
    rng=np.random.default_rng(612);max_error=0.
    for p in rng.uniform(-1,1,(12,3)):
        vals=np.linalg.eigvalsh(kinetic.bdg_symbol(p,K,speed=1))
        en=np.r_[np.repeat(np.sqrt(p@p+np.array(masses)**2),4),np.repeat(np.linalg.norm(p),78)]
        max_error=max(max_error,float(np.max(abs(vals-np.sort(np.r_[-en,en])))))
    assert max_error<2e-13
    # Rebuild every spin record using dense axial source evolution, with no
    # producer packet, frame, field or probability helper.
    sig=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    tz=np.diag([1.,-1.]);tx=sig[0];seed=np.ones(4)/2
    maximum=0.;checked=0
    for exp in spin['complete_source_field_frame_experiments']:
        m=masses[exp['character']];time=exp['time'];mom=[]
        for k in [0,1,0,-1]:
            H=-k*np.kron(tz,sig[2])+m*np.kron(tx,np.eye(2))
            ev,U=np.linalg.eigh(H);P=U[:,ev>0]@U[:,ev>0].conj().T
            mom.append(expm(-1j*time*H)@P@seed/2)
        w=sum(np.vdot(v,v).real for v in mom)
        phi=sum(mom)/2
        for row in exp['records']:
            R=np.kron(np.eye(2),expm(-1j*row['physical_quarter_turns']*np.pi*sig[2]/4))
            v=R@phi;prob=np.vdot(v,v).real/w
            currents=np.array([np.vdot(v,np.kron(np.eye(2),q)@v).real/w for q in sig])
            a=[1-w]
            for axis in range(3):
                for u in [-1,0,1]:a.append(w*((1-prob)if u==0 else(prob+u*currents[axis])/2)/3)
            b=[(1-w)**2,(1-w)*w,w*(1-w)]
            for axis in range(3):
                for other in range(3):
                    for u in [-1,0,1]:b.append(w*w*((1-prob*prob)if u==0 else(prob*prob+u*currents[axis]*currents[other])/2)/9)
            for expected,key in [(a,'complete_frame_records'),(b,'complete_two_copy_tensor_records')]:
                maximum=max(maximum,float(np.max(abs(np.array(expected)-[x['probability']for x in row[key]]))))
                assert min(expected)>-1e-14 and abs(sum(expected)-1)<1e-13
                checked+=len(expected)
    assert maximum<1e-13
    # Independent Fourier symmetry, including all site/spin rows.
    pts=list(__import__('itertools').product(range(4),repeat=3));P=np.zeros((64,64))
    for j,(x,y,z)in enumerate(pts):P[pts.index((-y%4,x,z)),j]=1
    U=np.kron(P,expm(-1j*np.pi*sig[2]/4))
    assert np.max(abs(np.linalg.matrix_power(U,4)+np.eye(128)))<1e-14
    assert np.max(abs(np.linalg.matrix_power(U,8)-np.eye(128)))<1e-14
    refs=['model/openwave_lepton_mass_shell/'+f for f in ['spectrum.py','PROOF.md','results.json']]+['model/openwave_spin_frame_records/'+f for f in ['records.py','PROOF.md','results.json']]
    return {'verdict':'PASS','independent_random_full_180_row_spectra':12,'maximum_spectral_error':max_error,
        'independent_complete_record_probabilities':checked,'maximum_probability_error':maximum,
        'actual128_row_spatial_rotation_checked':True,'scope':'Specified free branches and two-copy observer tensor; no full charged gap, fundamental nematic medium or total Noether angular momentum.',
        'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in refs}}
if __name__=='__main__':
    value=review();Path(__file__).with_suffix('.json').write_text(json.dumps(value,indent=2,sort_keys=True)+'\n');print(value['verdict'])
