"""Independent dressed-incidence and complete boundary/cut projector review."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
sys.path.insert(0,str(BASE/'next/C10'))
from charged_reservoir import BoundarySource


def execute():
    path=HERE/'BOUNDARY_CHARGE.json';saved=json.loads(path.read_text())
    for relative,pin in saved['source_sha256'].items():
        if hashlib.sha256((BASE/relative).read_bytes()).hexdigest()!=pin:raise ValueError('producer custody drift')
    source=BoundarySource(6);states=source.states
    occupations=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in states]);numbers=occupations.sum(axis=1);boundary=3-numbers
    incidence=np.zeros((6,6),int)
    for e,(a,b) in enumerate(source.construction.seams):
        incidence[a,2*e]=1;incidence[b,2*e+1]=1;incidence[3+e,2*e:2*e+2]=-1
    records=np.array([f for _,f in states]);background=occupations+records@incidence.T
    background[:,0]+=boundary
    assert np.all(background%6==np.array([1,1,1,0,0,0]))
    rng=np.random.default_rng(778);amplitudes=rng.normal(size=(384,2))+1j*rng.normal(size=(384,2))
    amplitudes[numbers%2==1]=0;amplitudes/=np.linalg.norm(amplitudes)
    comparisons=0;outcomes=0;brancherror=0.;maxwidth=0
    for reported in saved['complete_regions']:
        region=set(reported['region']);indicator=np.array([int(z in region) for z in range(6)])
        charge=occupations@indicator;cut=(records@(indicator@incidence))%6;g=int(np.array([1,1,1,0,0,0])@indicator)
        for row in reported['outcomes']:
            total=row['total_number'];b=row['boundary'];c=row['cut_residue'];assert b==3-total
            admitted=[q for q in range(7) if q<=len(region) and q<=total and 0<=total-q<=6-len(region)]
            maxwidth=max(maxwidth,max(admitted)-min(admitted))
            possibilities=[q for q in admitted if (q+c+b*int(0 in region)-g)%6==0]
            assert len(possibilities)<=1
            expected_charge=possibilities[0] if possibilities else None
            assert row['regional_charge']==expected_charge
            cut_mask=(numbers==total)&(cut==c)
            direct_mask=np.zeros(384,bool) if expected_charge is None else (numbers==total)&(charge==expected_charge)
            assert np.array_equal(cut_mask,direct_mask)
            assert row['basis_indices']==np.flatnonzero(cut_mask).tolist()
            assert row['zero_outcome_retained']==(not np.any(cut_mask))
            brancherror=max(brancherror,float(np.max(abs(amplitudes*cut_mask[:,None]-amplitudes*direct_mask[:,None]))))
            comparisons+=int(cut_mask.sum());outcomes+=1
    # Same-even-parity cross-number coherence, entangled with an idle reference.
    chosen=[next(i for i,(mask,_) in enumerate(states) if mask==m) for m in (0,6)]
    assert all(occupations[i,0]==0 for i in chosen) and [numbers[i]%2 for i in chosen]==[0,0]
    psi=np.zeros(4,complex);psi[0]=1/np.sqrt(2);psi[3]=-1j/np.sqrt(2);rho=np.outer(psi,psi.conj())
    # Order(source0/ref0,source0/ref1,source1/ref0,source1/ref1).
    p0=np.diag([1,1,0,0]);p1=np.eye(4)-p0;dephased=p0@rho@p0+p1@rho@p1
    distance=float(np.linalg.norm(rho-dephased,ord='nuc')/2)
    assert abs(distance-.5)<1e-12 and comparisons==24576 and outcomes==2688 and maxwidth==3
    paths=[Path(__file__),HERE/'boundary_charge.py',path,BASE/'next/C10/charged_reservoir.py']
    return {'schema':'oph.c11.boundary_charge_reference_review.v1','reviewer':'Second distinct Codex worker; no decoder-function import','outcome':'PASS_WITH_COMPLETE_BOUNDARY_INSTRUMENT_SCOPE','checks':{'independent_dressed_incidence_Gauss':True,'all_24576_basis_comparisons':comparisons==24576,'all_2688_outcomes_including_zeros':outcomes==2688,'Pauli_interval_width_at_most_three':maxwidth==3,'quantum_reference_branch_equality':brancherror==0,'Q_only_crossN_continuation_rejected':abs(distance-.5)<1e-12},'independent_method':'Construct six-vertex half-edge incidence; enumerate admitted regional occupancy using only cardinality/total constraints; solve modular Gauss by set intersection; compare full diagonal branch masks and tensor-idle-reference amplitudes. Producer decode/cut functions are not imported.','basis_comparisons':comparisons,'outcomes':outcomes,'external_reference_branch_error':brancherror,'cross_number_control':{'actual_source_basis_indices':chosen,'occupation_masks':[0,6],'same_even_fermion_parity':True,'idle_reference_dimension':2,'Q_only_vs_B_then_cut_trace_distance':distance},'scope':'Equality is between complete(B,cut) and complete(N,Q) instruments with retained labels, including zero outcomes and arbitrary idle references. Q-only measurement alone preserves extra cross-number coherence and is not a continuing substitute. Integer B is supplied; no physicalspecies/chargeunit identification follows.','full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
 r=execute();HERE.joinpath('BOUNDARY_CHARGE_REFERENCE_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
