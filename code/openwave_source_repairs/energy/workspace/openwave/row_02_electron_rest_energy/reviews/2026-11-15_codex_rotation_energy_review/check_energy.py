#!/usr/bin/env python3
"""Retrospective same-group review; no producer modules are imported."""
from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,math,time
import numpy as np
HERE=Path(__file__).resolve().parent
PACKET=HERE.parents[1]/'work/2026-11-15_codex_source_repair_v2'
OW=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def run():
    r=json.loads((PACKET/'coupled_energy/RECEIPT.json').read_text());s=r['scientific']
    a=np.load(PACKET/'coupled_energy/maps.npz');sector=np.load(PACKET/'poisson_energy/source_sector_maps.npz')
    source=np.load(OW/'shared/compact_logarithmic_flow/source_columns.npz')
    L=np.zeros((240,240),dtype=int)
    for i,j in source['seams']:L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    assert np.array_equal(L,1260*np.eye(240,dtype=int)-source['counts'])
    assert max(np.diag(L))==6
    pump_errors=[];effect_errors=[];sector_failures=[];matrices=[];comparison=F()
    for folder,prefix,offset,initial in [('first_cell_certified','first',0,0),('second_cell','second',12,2)]:
        m=np.load(PACKET/folder/'maps.npz');Q=m['certified_Q_proposal'];D=m['dark_integer_basis']
        certificate=json.loads((PACKET/folder/'RECEIPT.json').read_text())['scientific']
        delta=F(certificate['certified_subblock_error_upper']);kappa=4*delta+2*delta*delta
        comparison+=512*kappa/(1-512*kappa)+F(1,10**6)
        assert np.array_equal(L@D,6*D)
        d=D[offset:offset+12];P=d@np.linalg.inv(d.T@d)@d.T
        E=np.eye(12)-Q.conj().T@Q
        superoperator=np.kron(Q,Q.conj())+np.outer((np.eye(12)/12).ravel(),E.T.ravel())
        vec=np.zeros(144);vec[12*initial+initial]=1
        rho=(np.linalg.matrix_power(superoperator,512)@vec).reshape(12,12)
        error=float(np.max(abs(rho-a['source_pump_'+prefix])));assert error<1e-11;pump_errors.append(error)
        effects=sector[prefix+'_effects'];failed=sector[prefix+'_processor_failure']
        ownflags=np.array([np.trace(e@rho).real for e in effects]);assert np.max(abs(ownflags-a['source_flags_'+prefix]))<1e-11
        residual=float(np.linalg.norm(effects[1]-P,2));assert residual<1e-7;effect_errors.append(residual)
        fail=float(np.trace(failed@rho).real);sector_failures.append(fail);matrices.append(rho)
    p0,p1=a['source_flags_first'],a['source_flags_second']
    pair=np.kron(p1,p0);assert np.max(abs(pair-a['source_flags_pair']))<1e-12
    phase_errors={};acceptances={};decoded={};negative_centers={}
    for t in ('1/10','4'):
        for kind,E,flags in [('single',6,p0),('pair',12,pair)]:
            rec=s['complete_source_coupled_records'][t][kind]
            meta=s['native_flag_controlled_words'][t][kind]
            indices=meta['success_outputs'];fringe=[]
            for axis in ('X','Y'):
                law=np.array(rec[axis]['complete_probabilities']);assert len(law)==2**meta['bits']
                negative=float(-sum(law[law<0]));negative_centers[t+'_'+kind+'_'+axis]=negative
                assert negative<float(comparison) # Signed interval centers, not physical negative probability.
                probability=float(law[indices].sum());assert abs(probability-rec[axis]['energy_success_probability'])<1e-15
                fringe.append(float((law[indices[0]]-law[indices[1]])/probability))
                assert abs(law.sum()+rec[axis]['processor_failure']-sum(flags))<1e-8
            z=complex(*fringe);target=np.exp(-1j*float(F(t))*E)
            residual=abs(z-target);assert residual<1e-7
            phase_errors[t+'_'+kind]=residual;acceptances[t+'_'+kind]=rec['energy_success_probability']
            if t=='1/10':
                value=-math.atan2(z.imag,z.real)/.1;decoded[kind]=value;assert abs(value-E)<1e-6
    # No source-value winding is selected: source degree bounds give
    # phase ranges[0,1.2] and[0,2.4], both wholly inside the principal branch.
    assert 2*12*.1<math.pi
    assert abs(decoded['pair']/decoded['single']-2)<1e-6
    for name,rec in s['negative_source_records'].items():
        assert rec['energy_success_probability']<1e-12 and rec['complex_accepted_fringe'] is None
    assert len(s['vacuum_native_records'])==1260 and all(row[-1]==0 for row in s['vacuum_native_records'])
    assert s['vacuum_transfer_eigenvalue']=='1' and s['vacuum_generator']=='0'
    assert len(s['vacuum_complete_native_address_law'])==256 and s['vacuum_complete_native_address_law'][240]>1-1e-7
    old=s['fixed_count_control_records'];zo=[]
    for kind in ('single','pair'):
        x,y=old[kind]['complex_accepted_fringe'];zo.append(complex(x,y))
    mismatch=abs(zo[1]-zo[0]**2);assert mismatch>.05 and abs(mismatch-s['fixed_count_phase_nonadditivity'])<1e-14
    # Expose the currently omitted upstream source-processor outcomes, without
    # pretending the tiny numeric sizes replace the complete error theorem.
    f0,f1=sector_failures;pairfirst=f0*np.trace(matrices[1]).real+sum(p0)*f1
    result={'retrospective_same_group_review':True,'producer_imports':[],
      'full144_superoperator_reconstruction_errors':pump_errors,
      'native_sector_effect_residuals':effect_errors,'phase_record_errors':phase_errors,
      'accepted_masses':acceptances,'decoded_principal_branch':decoded,
      'signed_negative_comparison_centers':negative_centers,'source_comparison_guard':str(comparison),
      'fixed_count_additivity_fringe_failure':mismatch,
      'upstream_sector_processor_failures':[f0,f1],'pair_upstream_first_failure':float(pairfirst),
      'complete_error_review_pending':True,
      'bindings':{str(p.relative_to(PACKET)):sha(p) for p in [PACKET/'PROOF.md',PACKET/'run_coupled_energy.py',PACKET/'coupled_energy/RECEIPT.json',PACKET/'coupled_energy/maps.npz',PACKET/'poisson_energy/RECEIPT.json',PACKET/'poisson_energy/source_sector_maps.npz',PACKET/'first_cell_certified/maps.npz',PACKET/'second_cell/maps.npz']}}
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args();started=time.perf_counter();r=run()
    if args.record:
        with (HERE/'ROUND1_RESULT.json').open('x') as f:json.dump(r,f,indent=2,sort_keys=True);f.write('\n')
    else:assert json.loads((HERE/'ROUND1_RESULT.json').read_text())==r
    print(json.dumps({'independent_core_pass':True,'decoded':r['decoded_principal_branch'],'fixed_count_failure':r['fixed_count_additivity_fringe_failure'],'pending_complete_instrument_repairs':True,'seconds':time.perf_counter()-started}))
