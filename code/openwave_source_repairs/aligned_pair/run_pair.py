#!/usr/bin/env python3
"""Complete joint-position extension. Frozen comparison, no producer imports."""
from pathlib import Path
from fractions import Fraction as F
import argparse, gzip, hashlib, json, math, sys, time
import numpy as np
sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent; I=HERE/'inputs'
def read(p): return json.loads(p.read_text())
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x): return (json.dumps(x,sort_keys=True,separators=(',',':'))+'\n').encode()
def tv(a,b): return float(np.sum(abs(a-b))/2)
def sci(n): return read(I/n)['scientific']

def projector(d):
    # Exact rational projection from the native full-source integer basis.
    gram=d.T@d; det=int(gram[0,0]*gram[1,1]-gram[0,1]*gram[1,0])
    inv=[[F(int(gram[1,1]),det),F(-int(gram[0,1]),det)],
         [F(-int(gram[1,0]),det),F(int(gram[0,0]),det)]]
    return [[sum(F(int(d[a,j]))*inv[j][k]*int(d[b,k]) for j in range(2) for k in range(2)) for b in range(len(d))] for a in range(len(d))]

def coordinates(aligned):
    z=(0,0);p=(1,0);m=(-1,0);g=(0,1);h=(0,-1)
    body=[(m,g,z),(p,g,z),(m,h,z),(p,h,z),(z,m,g),(z,p,g),
          (z,m,h),(z,p,h),(g,z,m),(g,z,p),(h,z,m),(h,z,p)]
    points=[]
    for cell in range(20):
        for x,y,z in body:
            # C=Rx(pi/2): (x,y,z)->(x,-z,y), then center=(0,10cell,0).
            if aligned and cell==1:y,z=tuple(-v for v in z),y
            points.append((x,(y[0]+10*cell,y[1]),z))
    assert len(set(points))==240
    return points

def rotate(points):
    lookup={x:i for i,x in enumerate(points)}
    p=[lookup[(tuple(-a for a in x),y,tuple(-a for a in z))] for x,y,z in points]
    assert sorted(p)==list(range(240)) and all(p[p[j]]==j for j in range(240))
    return np.array(p+list(range(240,256)))

def pump(q,seed,dephase=False):
    x=np.zeros((12,12),complex);x[seed,seed]=1
    escape=np.eye(12)-q.conj().T@q
    for _ in range(512):
        x=q@x@q.conj().T+np.trace(escape@x)*np.eye(12)/12
        if dephase:x=np.diag(np.diag(x))
    return x

def loaded(blocks,rho,effect):
    k=blocks[:,:,:12]
    return sum(effect[i,j]*k[j]@rho@k[i].conj().T for i in range(2) for j in range(2))

def address_law(rho,pointer):
    # All eight native bit effects, preserving every address and coherence.
    args=[rho.reshape((2,)*16),list(range(16))]
    for bit in range(8):args.extend([pointer,[16+bit,8+bit,bit]])
    args.append(list(range(16,24)))
    return np.einsum(*args,optimize='greedy').real.reshape(256)

def run():
    origins=read(HERE/'ORIGINS.json');assert sha(I/'common_MANIFEST.sha256')==origins['common_manifest_sha256']
    manifest={n:h for h,n in(line.split('  ',1)for line in(I/'common_MANIFEST.sha256').read_text().splitlines())}
    for n,m in origins['inputs'].items():assert sha(HERE/n)==m['sha256']==manifest[m['common_relative']]
    common=read(I/'common_RECEIPT.json');m=np.load(I/'common_MAPS.npz')
    first=np.load(I/'first_maps.npz');second=np.load(I/'second_maps.npz')
    d0=first['dark_integer_basis'];d1=second['dark_integer_basis']
    L=np.zeros((240,240),int)
    for a,b in read(I/'graph.json')['seams']:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    assert np.array_equal(L@d0,6*d0) and np.array_equal(L@d1,6*d1)
    assert not np.any(d0[0]) and not np.any(d1[14])
    r0=projector(d0[:12]);r1=projector(d1[12:24])
    coords=coordinates(True);oldcoords=coordinates(False);p=rotate(coords);oldp=rotate(oldcoords)
    p0=p[:12];p1=p[12:24]-12;oldp1=oldp[12:24]-12
    exact0=all(r0[p0[a]][p0[b]]==r0[a][b] for a in range(12)for b in range(12))
    exact1=all(r1[p1[a]][p1[b]]==r1[a][b] for a in range(12)for b in range(12))
    assert exact0 and exact1
    target=np.kron(np.array(r0,float)/2,np.array(r1,float)/2)
    jointperm=np.array([12*a+b for a in p0 for b in p1])
    assert np.array_equal(target,target[np.ix_(jointperm,jointperm)])
    diag0=[r0[j][j]/2 for j in range(12)];diag1=[r1[j][j]/2 for j in range(12)]
    oldtv=sum(abs(diag0[a]*diag1[b]-diag0[oldp[a]]*diag1[oldp1[b]]) for a in range(12)for b in range(12))/2
    assert oldtv==F(1,3)
    # One global proper rotation, derived from the position dictionary.
    assert p[0]==1 and p[1]==0 and coords[1][0][0]-coords[0][0][0]==2
    # The inherited exact source/record join is bound rather than new native calls.
    assert common['gates']['all12_same_source_handoffs']
    starts={x['index']:x for x in read(I/'route_HANDOFFS.json')['handoffs']}
    routes=json.loads(gzip.decompress((I/'route_RECORDS.json.gz').read_bytes()))
    reconstructed=0;live=0
    for trial in routes['trials']:
        q=starts[trial['index']]['observed_source'][:]
        for move in trial['moves']:
            for aa,bb,coin,current in move['records']:
                before=q[aa]
                if trial['law']=='native':
                    total=q[aa]+q[bb];q[aa]=total//2+coin*(total%2);q[bb]=total-q[aa]
                assert before-q[aa]==current;reconstructed+=1
        assert q==trial['observed_source']==trial['source_audit_only']
        if trial['law']=='native':
            assert trial['success'] and trial['pair_code']==14 and trial['positive_frame']==[0,1]
            assert q[0]==1 and q[14]==-1 and sum(abs(x)for x in q)==2;live+=1
        else:assert not trial['success'] and trial['terminal']=='timeout'
    assert live==12 and not routes['forced_native_timeout']['success']
    q0=first['certified_Q_proposal'];q1=second['certified_Q_proposal']
    rho0=pump(q0,0);rho1=pump(q1,2);control=pump(q0,0,True)
    assert max(np.max(abs(rho0-m['first_final_density'])),np.max(abs(rho1-m['second_final_density'])),np.max(abs(control-m['dephasing_final_density'])))<1e-12
    # The same full input density is used before the broad loader and either assay.
    joint=np.kron(rho0,rho1);assert joint.shape==(144,144)
    pointer=np.load(I/'pointer_witness.npz')['pointer_effects_C']
    reset=F(35,2**31);ptr=F(200000000001,20000000000000000000000);bit=reset+ptr
    flow=F(read(I/'SHARP_FLOW_CERTIFICATE.json')['512_complete_U4_halfdiamond_upper'])
    receipts=[sci('first_RECEIPT.json'),sci('second_RECEIPT.json')]
    raw=[flow+F(c['512epoch_apparatus_error_upper']) for c in receipts]
    gain=[x/(1-2*x)for x in raw];bad=[F(c['512epoch_ideal_bad_weight_upper'])for c in receipts]
    broad=F(common['complete_broad_loader_error_upper']);routing=F(read(I/'route_RECEIPT.json')['per_trial_timeout_upper'])
    pairbad=bad[0]+bad[1]-bad[0]*bad[1]
    S=sum(gain)+pairbad+broad+routing+16*bit;error=S/(1-2*S)
    oldS=sum(gain)+sum(bad)+broad+routing+8*bit
    assert oldS/(1-2*oldS)==F(common['spatial_complete_error_upper'])
    fa=read(I/'FIELD_A_RESULT.json');comparison=F(fa['comparison_512_error_upper'])
    delta=F(receipts[1]['certified_subblock_error_upper']);kap=4*delta+2*delta*delta;u=F(1,2**53);gamma=128*u/(1-128*u)
    comparison1=512*kap/(1-512*kap)+512*12*4096*gamma/(1-512*kap)+F(1,10**12)
    x=raw[0]+512*8*bit;C=x/(1-2*x)+gain[1]+broad+routing+16*bit
    control_error=C/(1-2*C)+comparison+comparison1
    # Data processing lower bound is for the full joint law, without postselection.
    control_lower=F(fa['dephasing_coordinate_halfturn_TV'])-2*control_error
    arrays={'prepared_full144_density':joint,'selected_full144_density':target,
            'global_Ry_address_permutation':p,'unrotated_chart_Ry_permutation':oldp}
    def basis(k):
        out=np.zeros((12,12),complex);out[k,k]=1;return out
    laws={}
    cases=[('native',rho0,rho1),('dephasing',control,rho1),
           ('polar_plant',basis(0),basis(2)),('rotated_polar_plant',basis(int(p0[0])),basis(int(p1[2])))]
    for name,a,b in cases:
        flags=[];positions=[];scratch=[];dirty=[]
        for k,rho in enumerate((a,b)):
            block=m['loader%d_clean_columns'%k]
            states=[loaded(block,rho,E)for E in pointer]
            flags.append([float(np.trace(state).real)for state in states])
            positions.append(address_law(states[1],pointer))
            grams=m['loader%d_scratch_word_grams'%k]
            scratch.append([float(np.trace(g[:12,:12]@rho).real)for g in grams])
            dirty.append(float(np.trace(m['loader%d_dirty_gram'%k][:12,:12]@rho).real))
        law=np.outer(*positions)
        assert law.shape==(256,256)
        rejects=[flags[0][a]*flags[1][b]for a in range(2)for b in range(2)if(a,b)!=(1,1)]
        dirty0=dirty[0];dirty1=sum(flags[0])*dirty[1]
        residual=float(1-law.sum()-sum(rejects)-dirty0-dirty1)
        assert abs(residual)<1e-7
        whole_tv=tv(law,law[np.ix_(p,p)])
        marginal=law.sum(axis=1);marginal_tv=tv(marginal,marginal[p])
        assert whole_tv+1e-12>=marginal_tv
        arrays[name+'_all65536_position_centres']=law
        arrays[name+'_full144_density']=np.kron(a,b)
        laws[name]={'joint_halfturn_TV_centre':whole_tv,'positive_marginal_TV_centre':marginal_tv,
          'complete_position_mass_centre':float(law.sum()),'region_reject_centres':rejects,
          'scratch_words_comparison_centres':scratch,'dirty_first_centre':dirty0,'dirty_second_after_first_clean_centre':dirty1,
          'unresolved_processor_plus_rounding_residual':residual,
          'collision_address_centre':float(np.trace(law[:240,:240])),
          'invalid_address_centre':float(law.sum()-law[:240,:240].sum())}
    unequal=.75*arrays['polar_plant_all65536_position_centres']+.25*arrays['rotated_polar_plant_all65536_position_centres']
    arrays['unequal_plant_all65536_position_centres']=unequal
    unequal_tv=tv(unequal,unequal[np.ix_(p,p)])
    # Position dephasing commutes with each bijective position permutation.
    deph_covariant=all((a==b)==(p[a]==p[b]) for a in range(256) for b in range(256))
    gates={'injective_whole240_chart_and_global_proper_halfturn':len(set(coords))==240,
      'exact_full_selected_density_return':exact0 and exact1,
      'old_unaligned_whole_pair_failure_retained':oldtv==F(1,3),
      'actual_recorded_frame_reverses':p[0]==1 and p[1]==0,
      'all12_actual_same_source_continuations':live==12,
      'same_full_preparation_before_unfiltered_joint_reader':np.max(abs(joint-np.kron(m['first_final_density'],m['second_final_density'])))<1e-12,
      'complete_native_joint_halfturn_bound':2*error<F(11,10000) and laws['native']['joint_halfturn_TV_centre']<float(2*error),
      'same_chart_active_control_separates':control_lower>F(3,100) and laws['dephasing']['joint_halfturn_TV_centre']>float(control_lower),
      'positive_body_direction_and_shape_retained':F(1,2)-error>F(1,3) and F(1,3)-2*error>F(1,5),
      'dephasing_memory_operation_covariant':deph_covariant,
      'same_joint_reader_sees_polar_and_unequal_plants':laws['polar_plant']['joint_halfturn_TV_centre']>.99 and unequal_tv>.49,
      'complete_outcomes_no_renormalization':all(abs(v['unresolved_processor_plus_rounding_residual'])<1e-7 for v in laws.values())}
    gates={k:bool(v)for k,v in gates.items()}
    result={'candidate':'A-old240-aligned-pair-spatial-readout-v1','gates':gates,'all_gates_pass':all(gates.values()),
      'common_manifest_sha256':origins['common_manifest_sha256'],'coords_exact_a_plus_b_phi':coords,
      'global_rotation_matrix':[[-1,0,0],[0,1,0],[0,0,-1]],'second_cell_chart_matrix':[[1,0,0],[0,0,-1],[0,1,0]],
      'geometry_alignment_target_informed_and_supplied':True,'unaligned_joint_ideal_TV':str(oldtv),
      'full144_density_not_marginal_target':True,'joint_address_outcomes':65536,'terminal_native_bit_readers':16,
      'full_joint_complete_error_upper':str(error),'full_joint_halfturn_TV_upper':str(2*error),
      'ideal_pair_bad_weight_upper':str(pairbad),'unequal_joint_plant_TV_centre':unequal_tv,
      'active_dephasing_full_joint_TV_lower':str(control_lower),'complete_control_error_upper':str(control_error),
      'terminal16_error_upper':str(16*bit),'source_and_loader_costs_unchanged':True,'laws':laws,
      'source_energy_density_unchanged':True,'finite_exact_halfturn_claimed':False,'full_archive_rotation_invariance_claimed':False,
      'inherited_route_records_reconstructed':reconstructed,'actual_handoffs_checked':live,
      'literal_new_native_calls':0,'new_compiled_fullflow_executions':0,'native_bit_effect_compositions':16,
      'same_group_only':True,'other_group_G6':False}
    return result,arrays

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args();start=time.monotonic()
    for n,h in read(HERE/'PREREG.json')['files'].items():assert sha(HERE/n)==h,n
    result,arrays=run()
    if args.record:
        with(HERE/'RECEIPT.json').open('xb')as f:f.write(canonical(result))
        np.savez_compressed(HERE/'PAIR_MAPS.npz',**arrays)
    else:
        assert read(HERE/'RECEIPT.json')==json.loads(canonical(result))
        old=np.load(HERE/'PAIR_MAPS.npz')
        for k,v in arrays.items():assert np.array_equal(old[k],v),k
    print(json.dumps({'all_gates_pass':result['all_gates_pass'],'joint_upper':float(F(result['full_joint_halfturn_TV_upper'])),
       'control_lower':float(F(result['active_dephasing_full_joint_TV_lower'])),'native_centre':result['laws']['native']['joint_halfturn_TV_centre'],
       'control_centre':result['laws']['dephasing']['joint_halfturn_TV_centre'],'seconds':time.monotonic()-start},indent=2))
    assert result['all_gates_pass']
if __name__=='__main__':main()
