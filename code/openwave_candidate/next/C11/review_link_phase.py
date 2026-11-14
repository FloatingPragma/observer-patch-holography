"""Independent CAR bit/hopping, dense evolution and external-phase audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.linalg import expm
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(BASE/'next/C10'));sys.path.insert(0,str(BASE/'next/C9'))
from charged_reservoir import BoundarySource
from audit_car_native_algebra import phase


def run():
    paths=[Path(__file__),HERE/'link_phase_probe.py',HERE/'LINK_PHASE_PROBE.json',HERE/'link_phase_states.npz',BASE/'next/C10/joint_source_states.npz',BASE/'next/C10/charged_reservoir.py']
    before={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    data=json.loads((HERE/'LINK_PHASE_PROBE.json').read_text());target=np.load(HERE/'link_phase_states.npz');joint=np.load(BASE/'next/C10/joint_source_states.npz');h0=joint['hamiltonian'];rho=joint['global_allN_reference'];m=BoundarySource(6)
    a=np.array([0,0,0,0,.17,0.]);ts=[np.zeros_like(h0) for _ in range(6)];occupation=np.array([[(mask>>z)&1 for z in range(6)] for mask,_ in m.states]);numbers=occupation.sum(axis=1)
    for col,(mask,f) in enumerate(m.states):
        u=occupation[col]
        for e,(x,y) in enumerate(m.construction.seams):
            for side,port in enumerate((x,y)):
                if u[port] and not u[3+e]:
                    record=list(f);record[2*e+side]=(record[2*e+side]+1)%6
                    row=m.index[(mask^(1<<port)^(1<<(3+e)),tuple(record))]
                    ts[2*e+side][row,col]=(1 if side==0 else -1)*phase(u,port,3+e)
    def Hamiltonian(a):return h0.astype(complex)+sum(((np.exp(-1j*x)-1)*t+(np.exp(1j*x)-1)*t.T for x,t in zip(a,ts)),start=np.zeros_like(h0,dtype=complex))
    h=Hamiltonian(a);currents=[1j*(np.exp(-1j*x)*t-np.exp(1j*x)*t.T) for x,t in zip(a,ts)]
    ids=np.flatnonzero(numbers==3);u=expm(-.37j*h[np.ix_(ids,ids)]);small=u@rho[np.ix_(ids,ids)]@u.conj().T;state=np.zeros_like(h);state[np.ix_(ids,ids)]=small
    herror=float(np.max(abs(h-target['intervened_H'])));stateerror=float(np.max(abs(state-target['coherent_after_probe'])))
    qvalues=occupation[:,0]+occupation[:,3];q=np.diag(qvalues);regional=1j*(h@q-q@h);currenterror=float(np.max(abs(regional-(currents[5]-currents[1]))))
    lambdas=np.array([.2,-.37,.11,.43,-.19,.07]);g=np.exp(1j*occupation@lambdas);ga=a+np.array([lambdas[x]-lambdas[3+e] for e,ends in enumerate(m.construction.seams) for x in ends]);covariance=float(np.max(abs(Hamiltonian(ga)-g[:,None]*h*g.conj()[None,:])))
    kerror=0.;ledgererror=0.;posterror=0.;deltaerror=0.;lawresults={}
    for law in ('C','E','MIX'):
        out=.71*state;records=[]
        for e,ops in enumerate(m.operators[law]):
            for meta,k0 in ops:
                d=meta['record_value'];j=d if meta['law']=='C' else (d-(2*meta['coin']-1)*abs(d))//2
                k=k0*np.exp(-1j*j*(a[2*e]-a[2*e+1]));kg=k0*np.exp(-1j*j*(ga[2*e]-ga[2*e+1]));kerror=max(kerror,float(np.max(abs(kg-g[:,None]*k*g.conj()[None,:]))))
                part=.29*k@state@k.conj().T/3;out+=part;prob=float(np.trace(part).real);charge=float(np.trace(q@part).real);energy=float(np.trace(h@part).real)
                beforecharge=float(.29/3*np.dot(qvalues*np.sum(abs(k)**2,axis=0),np.diag(state).real))
                x,y=m.construction.seams[e];expected_delta=(int(y in (0,3))-int(x in (0,3)))*j*prob
                deltaerror=max(deltaerror,abs(charge-beforecharge-expected_delta));records.append((prob,charge,energy))
        records.append((.71,.71*float(np.trace(q@state).real),.71*float(np.trace(h@state).real)))
        expected=data['native'][law]['complete_outcomes']
        for own,saved in zip(records,expected):ledgererror=max(ledgererror,*[abs(x-saved[key]) for x,key in zip(own,('probability','unnormalized_regional_charge','unnormalized_energy'))])
        posterror=max(posterror,float(np.max(abs(out-target['native_'+law]))));lawresults[law]={'mean_work':float(np.trace(h@(out-state)).real),'sum_probability':sum(x[0] for x in records)}
    def ground(a):
        matrix=Hamiltonian(a);choices=[]
        for n in range(7):
            ix=np.flatnonzero(numbers==n);e,v=np.linalg.eigh(matrix[np.ix_(ix,ix)]);choices.append((float(e[0]),n,ix,v[:,0]))
        return min(choices,key=lambda x:x[0])
    e,n,ix,vector=ground(a);current=float(np.vdot(vector,currents[4][np.ix_(ix,ix)]@vector).real)
    eps=1e-4;energies={}
    for delta in (-2,-1,1,2):
        shifted=a.copy();shifted[4]+=delta*eps;energies[delta]=ground(shifted)[0]
    derivative=(energies[-2]-8*energies[-1]+8*energies[1]-energies[2])/(12*eps);hferror=abs(current+derivative)
    from link_phase_probe import ExternalPhaseProbe
    from types import SimpleNamespace
    # Use the actual publicAPI only for the subnormalizedcomposition regression.
    proxy=SimpleNamespace(n=m.n,field=SimpleNamespace(halfhops=ts,operators_by_law=m.operators),source=m,hamiltonian=h0,states=m.states)
    api=ExternalPhaseProbe(proxy);sub,subledger=api.native_step(.4*state,a)
    suberror=abs(sum(row['probability'] for row in subledger['complete_outcomes'])-.4)
    checks={'independent_H_and_signed_halfhops':herror<1e-12,'Krylov_matches_independent_dense_evolution':stateerror<1e-12,'regional_current_sign_identity':currenterror<1e-12,'continuous_external_probe_gauge_family':covariance<1e-12,'every_native_K_gauge_covariance':kerror<1e-12,'every_native_probability_charge_energy_ledger':ledgererror<1e-12 and posterror<1e-12,'every_applied_current_conditional_charge_delta':deltaerror<1e-12,'static_Hellmann_Feynman_fivepoint':hferror<1e-8,'subnormalized_branch_no_commit_probability':suberror<1e-12 and abs(np.trace(sub).real-.4)<1e-12}
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    if before!={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}:raise ValueError('sourcechanged')
    return {'schema':'oph.c11.link_phase_independent_review.v1','checks':checks,'errors':{'H':herror,'state':stateerror,'current_sign':currenterror,'H_gauge':covariance,'K_gauge':kerror,'ledgers':ledgererror,'poststates':posterror,'charge_delta':deltaerror,'HF':hferror,'subnormalization':suberror},'static_current':current,'static_ground_number':n,'native':lawresults,'scope':'DeclaredexternalPeierlsprobe on frozenfullC10source, orbitalgaugefamily notinternalcontinuousGaussgenerator atmodwrap. Nativeperrecordphase cancelsdensitybranch. Idealstaticperturbedground preparation differs fromquenchedoriginalground; no physicalB/spin/g2/alpha/radiation/physicalclock claim.','source_sha256':before,'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('LINK_PHASE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'errors':r['errors'],'current':r['static_current']}))
