#!/usr/bin/env python3
"""Certified full-coordinate source comparison; not raw native execution."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,time
import numpy as np
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tv(p,q):return float(np.sum(np.abs(p-q))/2)
def td(a,b):return float(np.sum(abs(np.linalg.eigvalsh((a-b+(a-b).conj().T)/2)))/2)
def permute(a,p):return a[np.ix_(p,p)]
def exact_projector():
    D=[[F(0),F(0)] for _ in range(12)]
    for j,c in zip([4,5,6,7],[1,-1,-1,1]):D[j][0]=F(c)
    for j,c in zip([5,6,8,11],[-1,-1,1,1]):D[j][1]=F(c)
    inverse=[[F(1,3),F(-1,6)],[F(-1,6),F(1,3)]]
    R=[[sum(D[a][i]*inverse[i][j]*D[b][j] for i in range(2) for j in range(2)) for b in range(12)] for a in range(12)]
    assert all(sum(R[a][c]*R[c][b] for c in range(12))==R[a][b] for a in range(12) for b in range(12))
    assert sum(R[a][a] for a in range(12))==2 and all(R[a][0]==0 for a in range(12))
    return D,R

def run():
    freeze=json.loads((HERE/'PREREG_FIELD_A.json').read_text())
    for n,h in freeze['files'].items():assert sha(HERE/n)==h,n
    origin=json.loads((HERE/'inputs/first_cell_RECEIPT.json').read_text())
    assert sha(HERE/'inputs/first_cell_maps.npz')==origin['artifacts']['maps.npz']
    s=origin['scientific'];data=np.load(HERE/'inputs/first_cell_maps.npz')
    Q=data['certified_Q_proposal'];assert Q.shape==(12,12)
    D,exact=exact_projector();R=np.array(exact,dtype=float);target=R/2
    wholeD=data['dark_integer_basis'];assert np.array_equal(wholeD[:12],np.array(D,dtype=int)) and not np.any(wholeD[12:])
    graph=json.loads((HERE/'inputs/graph.json').read_text())
    L=np.zeros((240,240),dtype=int)
    for a,b in graph['seams']:L[a,a]+=1;L[b,b]+=1;L[a,b]-=1;L[b,a]-=1
    assert np.array_equal(L@wholeD,6*wholeD)
    # Coordinates in Q(phi), encoded as integer coefficient pairs. Sign flips
    # are exact spatial rotations, not changes to probability-bin definitions.
    z=(0,0);one=(1,0);minus=(-1,0);phi=(0,1);nphi=(0,-1)
    coords=[(minus,phi,z),(one,phi,z),(minus,nphi,z),(one,nphi,z),
            (z,minus,phi),(z,one,phi),(z,minus,nphi),(z,one,nphi),
            (phi,z,minus),(phi,z,one),(nphi,z,minus),(nphi,z,one)]
    def spatial(signs):
        return [coords.index(tuple((sg*x[0],sg*x[1]) for sg,x in zip(signs,v))) for v in coords]
    ry=spatial([-1,1,-1]);rx=spatial([1,-1,-1]);rz=spatial([-1,-1,1])
    for p in (rx,ry,rz):
        assert sorted(p)==list(range(12)) and all(p[p[i]]==i for i in range(12))
        edges={tuple(e) for e in graph['seams'] if max(e)<12}
        assert {tuple(sorted((p[a],p[b]))) for a,b in edges}==edges
    assert all(exact[ry[a]][ry[b]]==exact[a][b] for a in range(12) for b in range(12))
    assert [exact[j][j]/2 for j in range(12)]==[F(int(j in (4,5,6,7,8,11)),6) for j in range(12)]
    assert ry[0]==1 and ry[1]==0 # Actual recorded frame arrow reverses.
    route=json.loads((HERE/'ROUTE_RECORDS.json').read_text())
    for trial in route['trials']:
        assert trial['observed_final_source']==[1]+[0]*239 and trial['frame_native_event']==[0,1,0,1]
    sigma=np.eye(12)/12;effect=np.eye(12)-Q.conj().T@Q
    delta=F(s['certified_subblock_error_upper']);kappa=4*delta+2*delta*delta
    u=F(1,2**53);gamma=128*u/(1-128*u)
    roundguard=512*12*4096*gamma/(1-512*kappa)+F(1,10**12)
    comparison=512*kappa/(1-512*kappa)+roundguard
    reset=F(35,2**31);pointer=F(200000000001,20000000000000000000000)
    address_reader=8*(reset+pointer)
    control_extra=512*address_reader
    reserved_native=F(1,1000)
    source_uncertainty=comparison+reserved_native
    control_uncertainty=source_uncertainty+control_extra
    rho=np.zeros((12,12),complex);rho[0,0]=1
    deph=rho.copy();history=[];density=[];dephdensity=[];escapes=[0.,0.]
    for epoch in range(513):
        history.append([epoch,float(np.trace(rho).real),float(np.trace(deph).real),
                        tv(np.diag(rho).real,np.diag(permute(rho,ry)).real),
                        tv(np.diag(deph).real,np.diag(permute(deph,ry)).real),*escapes])
        density.append(rho.copy());dephdensity.append(deph.copy())
        if epoch==512:break
        for which,state in enumerate((rho,deph)):
            amount=float(np.trace(effect@state).real);escapes[which]+=amount
            assert amount>=-float(comparison)
            updated=Q@state@Q.conj().T+amount*sigma
            assert np.max(abs(updated))<2
            if which==0:rho=updated
            else:deph=np.diag(np.diag(updated))
    p=np.zeros(257);p[:12]=np.diag(rho).real
    c=np.zeros(257);c[:12]=np.diag(deph).real
    # The comparison is not renormalized or clipped. Its trace drift is kept.
    p[-1]=1-sum(p[:-1]);c[-1]=1-sum(c[:-1])
    assert max(abs(p[-1]),abs(c[-1]))<float(comparison)
    ytv=tv(np.diag(rho).real,np.diag(permute(rho,ry)).real)
    controltv=tv(np.diag(deph).real,np.diag(permute(deph,ry)).real)
    fulltd=td(rho,permute(rho,ry));distance=td(rho,target)
    bad=F(s['512epoch_ideal_bad_weight_upper']);assert bad==F(1,65536)
    final_halfturn_bound=2*(bad+reserved_native)
    native_isotropic=tv(np.diag(rho).real,np.diag(sigma))
    native_x=tv(np.diag(rho).real,np.diag(permute(rho,rx)).real)
    # Finite support comparison at finer physical angles uses the full spatial
    # coordinates, retaining new rotated bins; it does not round to vertices.
    phinum=(1+math.sqrt(5))/2
    xyz=np.array([[a+b*phinum for a,b in v] for v in coords]);support=[4,5,6,7,8,11]
    fine=[]
    for angle in (math.pi/4,math.pi/2,3*math.pi/4):
        C,S=math.cos(angle),math.sin(angle);matrix=np.array([[C,0,S],[0,1,0],[-S,0,C]])
        rotated=xyz[support]@matrix.T
        overlap=sum(any(np.max(abs(v-w))<1e-12 for w in xyz[support]) for v in rotated)
        value=1-overlap/6;fine.append(value)
        assert value>0
    polar=np.zeros((12,12));polar[0,0]=1
    unequal=np.zeros((12,12));unequal[0,0]=.75;unequal[1,1]=.25
    gates={
      'exact_source_mode_and_spatial_algebra':True,
      'native_halfturn_enclosed':ytv<=float(2*(bad+comparison)) and fulltd<=float(2*(bad+comparison)),
      'conditional_complete_halfturn_bound':final_halfturn_bound<F(1,250),
      'active_control_escape':escapes[1]>0,
      'dephasing_separation':controltv>max(.01,10*float(control_uncertainty)),
      'retained_shape':native_isotropic-float(source_uncertainty)>1/3 and native_x-2*float(source_uncertainty)>1/5,
      'polar_and_unequal_plant_seen':tv(np.diag(polar),np.diag(permute(polar,ry)))==1 and tv(np.diag(unequal),np.diag(permute(unequal,ry)))==.5,
      'same_selected_density':distance<=float(bad+comparison),
    }
    result={'candidate':'A-old240-compact-energy-apolar-join-v1','stage':'F-A certified comparison',
      'native_implementation_transfer_complete':False,'fresh_native_flow_executions':0,
      'gates':gates,'all_F_A_gates':all(gates.values()),'epochs':512,
      'native_coordinate_halfturn_TV':ytv,'native_density_halfturn_TD':fulltd,
      'native_distance_to_selected_density':distance,'dephasing_coordinate_halfturn_TV':controltv,
      'native_and_dephasing_expected_escape_counts':escapes,
      'native_isotropic_TV':native_isotropic,'native_Rx_TV':native_x,
      'fine_y_angles_TV':fine,'source_coordinate_law':p.tolist(),'dephasing_coordinate_law':c.tolist(),
      'exact_rotations':{'Ry':ry,'Rx':rx,'Rz':rz},
      'comparison_512_error_upper':str(comparison),'reserved_native_error_allowance_NOT_YET_CERTIFIED':str(reserved_native),
      'control_extra_readout_error_upper':str(control_extra),'terminal_8bit_error_upper':str(address_reader),
      'conditional_final_halfturn_upper':str(final_halfturn_bound),
      'ideal_limit_exact_halfturn':True,'finite_exact_halfturn_claimed':False,
      'full_archive_invariance_claimed':False,'readout_merges_opposite_positions':False,
      'operator_input_receipt_sha256':sha(HERE/'inputs/first_cell_RECEIPT.json')}
    arrays={'native_density_history':np.asarray(density),'dephasing_density_history':np.asarray(dephdensity),'epoch_records':np.array(history),'source_position_record':p,'dephasing_position_record':c}
    return result,arrays

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();started=time.perf_counter()
    result,arrays=run()
    if args.record:
        with (HERE/'FIELD_A_RESULT.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
        np.savez_compressed(HERE/'FIELD_A_RECORDS.npz',**arrays)
    else:
        assert json.loads((HERE/'FIELD_A_RESULT.json').read_text())==result
        previous=np.load(HERE/'FIELD_A_RECORDS.npz')
        for k,v in arrays.items():assert np.array_equal(previous[k],v),k
    print(json.dumps({'all_F_A_gates':result['all_F_A_gates'],'gates':result['gates'],'native_TV':result['native_coordinate_halfturn_TV'],'control_TV':result['dephasing_coordinate_halfturn_TV'],'comparison_error':float(F(result['comparison_512_error_upper'])),'seconds':time.perf_counter()-started}))
    if not result['all_F_A_gates']:raise SystemExit(1)
