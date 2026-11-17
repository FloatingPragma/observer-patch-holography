"""Independent dense/high-precision checks of quantum geometry and wave records."""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import numpy as np
from scipy.linalg import eigh,expm

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parent
GEOMETRY=ROOT/'model/openwave_quantum_geometry/results.json'
WAVES=ROOT/'model/openwave_wave_deflection/results.json'

def value(x):
    v=F(x);return mp.mpf(v.numerator)/v.denominator

def mp_matrix_data(cap,gaps):
    n=cap+1;N=n*n
    A=mp.zeros(n);B=mp.zeros(n)
    for j in range(1,n):
        A[j-1,j]=A[j,j-1]=mp.sqrt(j)/200
        B[j-1,j]=B[j,j-1]=mp.sqrt(j)/(200*mp.root(3,4))
    ea,eb,ec=mp.expm(A),mp.expm(B),mp.expm(-B)
    Q=[mp.matrix([[ea[i//n,j//n]*E[i%n,j%n]for j in range(N)]for i in range(N)])for E in(eb,ec)]
    Hg=mp.diag([i+mp.sqrt(3)*j for i in range(n)for j in range(n)])
    return Q,Hg,[sum(map(value,p))/2 for p in gaps]

def native_geometry_review(packet):
    mp.mp.dps=65;rows=[]
    for group in packet['cap_comparisons']:
        cap=group['cap'];Q,Hg,energies=mp_matrix_data(cap,packet['native_gap_intervals'])
        blocks={};maximum_residual_ratio=0.;maximum_error_ratio=0.
        for cert in group['all64_conditional_spectral_certificates']:
            i,j=cert['native_levels'];H=Hg+energies[i]*Q[0]+energies[j]*Q[1];blocks[i,j]=H
            D=2**cert['dyadic_bits'];Oint=np.array(cert['eigenframe_integer_rows'],dtype=object)
            eigint=cert['eigenvalue_integers'];dim=len(eigint)
            gram=Oint.T@Oint-D*D*np.eye(dim,dtype=object)
            exact_frob_sq=F(sum(int(x)**2 for x in gram.flat),D**4)
            eta=F(cert['orthogonality_Frobenius_upper']);assert exact_frob_sq<=eta*eta
            O=mp.matrix([[mp.mpf(int(x))/D for x in line]for line in Oint])
            ev=mp.diag([mp.mpf(int(x))/D for x in eigint])
            residual=mp.norm(H*O-O*ev)
            bound=value(cert['exact_source_eigenframe_residual_Frobenius_upper'])
            assert residual<bound
            maximum_residual_ratio=max(maximum_residual_ratio,float(residual/bound))
            epsilon=F(cert['polar_unitary_Hamiltonian_error_upper'])
            assert epsilon==2*F(cert['exact_source_eigenframe_residual_Frobenius_upper'])+120*eta
            # A numerical polar audit is separate from the integer/interval proof.
            on=np.array(O.tolist(),dtype=float);hn=np.array(H.tolist(),dtype=float)
            u,_,v=np.linalg.svd(on);W=u@v
            actual=float(np.linalg.norm(hn-W@np.diag(np.array(eigint,dtype=float)/D)@W.T,2))
            assert actual<float(epsilon)+2e-14
            maximum_error_ratio=max(maximum_error_ratio,actual/float(epsilon))
        item=group['entanglement'];t=item['source_time'];vectors=[]
        for i,j in((0,0),(0,1),(1,0),(1,1)):
            e,U=mp.eigh(blocks[i,j]);vectors.append(U*mp.matrix([mp.exp(-mp.j*t*x)*U[0,k]for k,x in enumerate(e)]))
        rho=mp.matrix([[sum(vectors[a][n]*mp.conj(vectors[b][n])for n in range(len(vectors[0])))/4 for b in range(4)]for a in range(4)])
        w=mp.matrix([mp.mpc(a,b)/mp.sqrt(item['witness_squared_integer_norm'])for a,b in item['normalized_witness_integer_pairs']])
        pt=mp.matrix([[rho[2*(a//2)+b%2,2*(b//2)+a%2]for b in range(4)]for a in range(4)])
        expectation=mp.re((w.H*pt*w)[0]);lower,upper=map(value,item['exact_source_partial_transpose_witness_interval'])
        assert lower<expectation<upper<0
        for c in group['all49_positive_energy_source_pair_interactions']:
            assert F(c['mutual_ground_energy_interval'][1])<0
        for record in item['complete_certified_joint_instruments']:
            probs=record['all_active_clock_geometry_probability_intervals']
            assert sum(x[0]for line in probs for x in line)<=1<=sum(x[1]for line in probs for x in line)
            assert record['complete_outcome_count']==64*(cap+1)**2
        rows.append({'cap':cap,'all_native_blocks_reviewed':64,
                     'max_high_precision_residual_to_certificate_ratio':maximum_residual_ratio,
                     'max_dense_polar_error_to_certificate_ratio':maximum_error_ratio,
                     'independent_high_precision_witness':mp.nstr(expectation,35),
                     'certified_witness_interval':item['exact_source_partial_transpose_witness_interval'],
                     'all49_mutual_energies_strictly_negative':True,'all9_complete_instruments_retained':True})
    return rows

def independent_planar_matrix(q,gamma):
    n=len(q);K=np.zeros((n*n,n*n));S=np.exp((1+3*gamma)*q.ravel()/2)
    for x in range(n):
        for y in range(n):
            i=x*n+y
            for dx,dy in((1,0),(-1,0),(0,1),(0,-1)):
                a,b=x+dx,y+dy
                other=q[a,b]if 0<=a<n and 0<=b<n else 0.
                weight=np.exp((1-gamma)*(q[x,y]+other)/2)
                K[i,i]+=weight*S[i]*S[i]
                if 0<=a<n and 0<=b<n:K[i,a*n+b]-=weight*S[i]*S[a*n+b]
    assert np.max(abs(K-K.T))<1e-14
    return K

def dense_wave(q,gamma,t,offset=6):
    n=len(q);x,y=np.meshgrid(np.arange(n),np.arange(n),indexing='ij');c=n//2
    seed=np.exp(-((x-(c-10))**2+(y-(c+offset))**2)/18+.9j*x).ravel();seed/=np.linalg.norm(seed)
    H2=independent_planar_matrix(q,gamma);e,U=eigh(H2);assert e[0]>7/(n+1)**2 and e[-1]<9
    omega=np.sqrt(e);frame=U.T@seed;state=U@(np.exp(-1j*t*omega)*frame);hstate=U@(omega*np.exp(-1j*t*omega)*frame)
    return {'position':np.r_[.5,abs(state)**2/2],
            'velocity':2*np.vdot(state,y.ravel()*hstate).imag,
            'mean_y':float(abs(state)**2@y.ravel()),'norm':float(np.vdot(state,state).real)}

def planar_review(packet):
    rows=[];flat=None
    for variant in packet['variants']:
        if variant['side']!=31 or variant['transverse_preparation_offset']!=6:continue
        q=np.array(variant['complete_profile_slice']);gamma=variant['gamma']
        stored=next(r for r in variant['all_mode_rows']if r['source_time']==16)
        direct=dense_wave(q,gamma,16)
        if flat is None:flat=dense_wave(np.zeros_like(q),1,16)
        difference=float(np.max(abs(direct['position']-stored['complete_vacuum_and_all_site_records'])))
        assert difference<2e-13 and abs(direct['norm']-1)<3e-13
        shift=direct['mean_y']-flat['mean_y'];velocity=direct['velocity']-flat['velocity']
        assert abs(shift-stored['conditional_transverse_position_shift'])<1e-10
        assert abs(velocity-stored['transverse_velocity_shift'])<1e-11
        if gamma==-1:
            assert abs(shift)>1e-9
            assert all(v['first_order_transverse_position_shift']==0 for v in stored['comparison_first_order_rays'])
        rows.append({'gamma':gamma,'dense_all961_modes':True,'all_position_record_max_difference':difference,
                     'independent_position_shift':shift,'independent_velocity_shift':velocity,
                     'velocity_record_sum':sum(stored['complete_vacuum_and_velocity_POVM_records']),
                     'principal_ray_residual_retained':gamma==-1})
    assert len(rows)==3
    return rows

def run():
    geometry=json.loads(GEOMETRY.read_text());waves=json.loads(WAVES.read_text())
    inputs=[GEOMETRY,WAVES,Path(__file__)]
    for package in('openwave_quantum_geometry','openwave_wave_deflection'):
        inputs += list((ROOT/'model'/package).glob('*.py'))+[(ROOT/'model'/package/'PROOF.md')]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted(set(inputs))}
    return {'schema':'oph.independent_quantum_geometry_wave_review.v1','source_sha256':pins,
            'quantum_geometry':native_geometry_review(geometry),'planar_wave':planar_review(waves),
            'analytic_review':'polar bound 2R+120eta covers every64 native sector; exact rational witness uses partial-transpose operator norm<=1; finite wave gamma=-1 residual distinguishes the principal ray law',
            'arithmetic_scope':'exact integer Gram and rational bound reconstruction, independent65-digit source/witness diagnostic, independent full961-mode dense wave diagnostic; source interval proofs retain their stated scope',
            'verdict':'PASS_INDEPENDENT_NATIVE_QUANTUM_GEOMETRY_AND_PLANAR_WAVE_REVIEW'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();path=HERE/'geometry_review.json'
    if args.write:path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:
        from model.openwave_larmor_dynamics.certificate import equivalent
        assert equivalent(json.loads(path.read_text()),result)
    print(result['verdict'])
if __name__=='__main__':main()
