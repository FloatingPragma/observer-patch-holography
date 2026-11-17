"""Separate full-pair geometry, instrument and budget verification."""
from fractions import Fraction as F
import hashlib,json,sys,time
from pathlib import Path
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
SOURCE=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else HERE.parents[1]/'work/2026-11-16_codex_aligned_pair_spatial'
SEAL='db48ac7ee3c80f75a63ea991907f984c54135282316a8ef023878fd014535ba8'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def tv(x,y):return float(np.abs(x-y).sum()/2)
def exact_R(d):
    a,b,c=int(d[:,0]@d[:,0]),int(d[:,0]@d[:,1]),int(d[:,1]@d[:,1]);det=a*c-b*b
    return [[F(c*int(d[i,0])*int(d[j,0])-b*(int(d[i,0])*int(d[j,1])+int(d[i,1])*int(d[j,0]))+a*int(d[i,1])*int(d[j,1]),det)for j in range(12)]for i in range(12)]
def exact_points(rotate_second):
    # Independent coefficient representation in the basis 1,phi.
    pairs=[(-1,0),(1,0),(0,0),(0,1),(0,-1)]
    minus,plus,zero,phi,nphi=pairs
    local=[(minus,phi,zero),(plus,phi,zero),(minus,nphi,zero),(plus,nphi,zero),
           (zero,minus,phi),(zero,plus,phi),(zero,minus,nphi),(zero,plus,nphi),
           (phi,zero,minus),(phi,zero,plus),(nphi,zero,minus),(nphi,zero,plus)]
    out=[]
    for cell in range(20):
        for a,b,c in local:
            if rotate_second and cell==1:b,c=tuple(-x for x in c),b
            out.append((a,(b[0]+10*cell,b[1]),c))
    return out
def permutation(points):
    lookup={q:i for i,q in enumerate(points)}
    assert len(lookup)==240
    out=[lookup[(tuple(-v for v in a),b,tuple(-v for v in c))]for a,b,c in points]
    assert all(out[out[k]]==k for k in range(240))
    return out+list(range(240,256))
def superoperator_pump(Q,seed,dephase=False):
    escape=np.eye(12)-Q.conj().T@Q
    M=np.kron(Q,Q.conj())+np.outer(np.eye(12).reshape(144)/12,escape.T.reshape(144))
    if dephase:M*=np.eye(12).reshape(144)[:,None]
    v=np.zeros(144,complex);v[13*seed]=1
    return (np.linalg.matrix_power(M,512)@v).reshape(12,12)
def local_flag_state(columns,rho,E):
    A=columns[:,:,:12].reshape(512,12)
    block=(A@rho@A.conj().T).reshape(2,256,2,256)
    return np.einsum('ab,biaj->ij',E,block)
def all_addresses(rho,effects):
    result=[]
    for word in range(256):
        E=np.ones((1,1),complex)
        for bit in range(7,-1,-1):E=np.kron(E,effects[(word>>bit)&1])
        result.append(np.sum(E.T*rho).real)
    return np.asarray(result)

def main():
    start=time.perf_counter();I=SOURCE/'inputs';mfile=SOURCE/'MANIFEST.sha256'
    assert sha(mfile)==SEAL
    count=0
    for line in mfile.read_text().splitlines():
        digest,rel=line.split('  ',1);assert sha(SOURCE/rel)==digest;count+=1
    receipt=read(SOURCE/'RECEIPT.json');saved=np.load(SOURCE/'PAIR_MAPS.npz')
    common=read(I/'common_RECEIPT.json');native=np.load(I/'common_MAPS.npz')
    first=np.load(I/'first_maps.npz');second=np.load(I/'second_maps.npz')
    graph=read(I/'graph.json')['seams'];L=np.zeros((240,240),int)
    for a,b in graph:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    basis=[first['dark_integer_basis'],second['dark_integer_basis']]
    assert all(np.array_equal(L@d,6*d)for d in basis)
    r=[exact_R(basis[0][:12]),exact_R(basis[1][12:24])]
    points=exact_points(True);p=permutation(points);oldp=permutation(exact_points(False))
    assert receipt['coords_exact_a_plus_b_phi']==json.loads(json.dumps(points))
    assert np.array_equal(p,saved['global_Ry_address_permutation'])
    p0=p[:12];p1=[v-12 for v in p[12:24]]
    for a in range(12):
      for b in range(12):
       for c in range(12):
        for d in range(12):
         assert r[0][a][c]*r[1][b][d]==r[0][p0[a]][p0[c]]*r[1][p1[b]][p1[d]]
    oldtv=sum(abs(r[0][a][a]*r[1][b][b]-r[0][oldp[a]][oldp[a]]*r[1][oldp[12+b]-12][oldp[12+b]-12])for a in range(12)for b in range(12))/8
    assert oldtv==F(1,3)==F(receipt['unaligned_joint_ideal_TV'])
    assert p[0]==1 and p[1]==0 and np.linalg.det(np.array(receipt['global_rotation_matrix']))==1
    rho0=superoperator_pump(first['certified_Q_proposal'],0)
    rho1=superoperator_pump(second['certified_Q_proposal'],2)
    control=superoperator_pump(first['certified_Q_proposal'],0,True)
    state_error=max(np.max(abs(rho0-native['first_final_density'])),np.max(abs(rho1-native['second_final_density'])),np.max(abs(control-native['dephasing_final_density'])))
    assert state_error<1e-12
    assert np.max(abs(np.kron(rho0,rho1)-saved['prepared_full144_density']))<1e-12
    pointer=np.load(I/'pointer_witness.npz')['pointer_effects_C']
    def pure(k):a=np.zeros((12,12),complex);a[k,k]=1;return a
    cases={'native':(rho0,rho1),'dephasing':(control,rho1),
           'polar_plant':(pure(0),pure(2)),
           'rotated_polar_plant':(pure(p0[0]),pure(p1[2]))}
    law_error=0;recomputed={}
    for name,states in cases.items():
        local=[];flags=[];dirty=[];scratch=[]
        for i,state in enumerate(states):
            blocks=native[f'loader{i}_clean_columns']
            outs=[local_flag_state(blocks,state,E)for E in pointer]
            flags.append([float(np.trace(z).real)for z in outs])
            local.append(all_addresses(outs[1],pointer))
            dirty.append(float(np.trace(native[f'loader{i}_dirty_gram'][:12,:12]@state).real))
            scratch.append([float(np.trace(g[:12,:12]@state).real)for g in native[f'loader{i}_scratch_word_grams']])
        full=np.multiply.outer(*local);reference=saved[name+'_all65536_position_centres']
        err=float(np.max(abs(full-reference)));law_error=max(law_error,err);assert err<1e-12
        row=receipt['laws'][name]
        rejects=[flags[0][a]*flags[1][b]for a,b in [(0,0),(0,1),(1,0)]]
        residual=1-full.sum()-sum(rejects)-dirty[0]-sum(flags[0])*dirty[1]
        assert abs(residual-row['unresolved_processor_plus_rounding_residual'])<1e-12
        assert np.max(abs(np.asarray(scratch)-row['scratch_words_comparison_centres']))<1e-12
        assert np.max(abs(np.asarray(rejects)-row['region_reject_centres']))<1e-12
        fulltv=tv(full,full[np.ix_(p,p)]);marg=full.sum(axis=1)
        assert abs(fulltv-row['joint_halfturn_TV_centre'])<1e-12 and fulltv+1e-12>=tv(marg,marg[p])
        recomputed[name]=fulltv
    bit=F(35,2**31)+F(200000000001,20000000000000000000000)
    receipts=[read(I/'first_RECEIPT.json')['scientific'],read(I/'second_RECEIPT.json')['scientific']]
    flow=F(read(I/'SHARP_FLOW_CERTIFICATE.json')['512_complete_U4_halfdiamond_upper'])
    raw=[flow+F(x['512epoch_apparatus_error_upper'])for x in receipts]
    gains=[x/(1-2*x)for x in raw];bads=[F(x['512epoch_ideal_bad_weight_upper'])for x in receipts]
    joint_bad=1-(1-bads[0])*(1-bads[1]);assert F(receipt['ideal_pair_bad_weight_upper'])==joint_bad
    loading=F(common['complete_broad_loader_error_upper']);route=F(read(I/'route_RECEIPT.json')['per_trial_timeout_upper'])
    prior_S=sum(gains)+sum(bads)+loading+route+8*bit
    assert prior_S/(1-2*prior_S)==F(common['spatial_complete_error_upper'])
    S=prior_S+8*bit-bads[0]*bads[1]
    error=S/(1-2*S);assert error==F(receipt['full_joint_complete_error_upper'])
    assert 2*error==F(receipt['full_joint_halfturn_TV_upper'])
    x=raw[0]+512*8*bit;C=x/(1-2*x)+gains[1]+loading+route+16*bit
    fa=read(I/'FIELD_A_RESULT.json');delta=F(receipts[1]['certified_subblock_error_upper'])
    kappa=4*delta+2*delta*delta;gamma=F(128,2**53-128)
    second_comparison=(512*kappa+512*12*4096*gamma)/(1-512*kappa)+F(1,10**12)
    control_error=C/(1-2*C)+F(fa['comparison_512_error_upper'])+second_comparison
    lower=F(fa['dephasing_coordinate_halfturn_TV'])-2*control_error
    assert lower==F(receipt['active_dephasing_full_joint_TV_lower'])
    assert recomputed['native']<float(2*error) and recomputed['dephasing']>float(lower)>.03
    result=dict(producer_manifest_sha256=SEAL,verified_producer_files=count,
        exact_full_tensor_entries=12**4,complete_joint_bins_per_case=65536,cases=4,
        independent_matrix_power_state_error=float(state_error),independent_joint_instrument_max_error=law_error,
        retained_old_chart_TV=str(oldtv),new_joint_TV_upper=float(2*error),active_control_lower=float(lower),
        all_outcomes_and_costs_verified=True,producer_imports=False,review_group='Codex same-group',
        seconds=time.perf_counter()-start)
    (HERE/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
