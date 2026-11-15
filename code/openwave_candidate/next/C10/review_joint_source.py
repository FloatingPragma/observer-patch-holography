"""Independent full-state interacting all-number H/reference/reservoir audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1]
sys.path.insert(0,str(HERE));sys.path.insert(0,str(BASE/'next/C9'))
from charged_reservoir import BoundarySource
from audit_car_native_algebra import phase


def run():
    files=[Path(__file__),HERE/'joint_source.py',HERE/'JOINT_SOURCE.json',HERE/'joint_source_states.npz',HERE/'charged_reservoir.py',HERE/'modular_field_adapter.py',BASE/'next/C9/fermionic_shared_records.py']
    before={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    saved=json.loads((HERE/'JOINT_SOURCE.json').read_text());archive=np.load(HERE/'joint_source_states.npz');m=BoundarySource(6);numbers=np.array([mask.bit_count() for mask,_ in m.states])
    hbare=np.zeros((m.n,m.n));w=np.zeros_like(hbare);electric=np.zeros(m.n)
    gauss=True
    for col,(mask,f) in enumerate(m.states):
        u=[(mask>>z)&1 for z in range(6)];hbare[col,col]=.23*sum((1 if z<3 else -1)*u[z] for z in range(6));g=u.copy();g[0]+=3-sum(u)
        electric[col]=.37/2*sum(min(x,6-x)**2 for x in f)
        for e,(a,b) in enumerate(m.construction.seams):
            g[a]+=f[2*e];g[b]+=f[2*e+1];g[3+e]-=f[2*e]+f[2*e+1]
            for side,x in enumerate((a,b)):
                for src,dst,delta in ((x,3+e,1),(3+e,x,-1)):
                    if u[src] and not u[dst]:
                        record=list(f);record[2*e+side]=(record[2*e+side]+delta)%6
                        hbare[m.index[(mask^(1<<src)^(1<<dst),tuple(record))],col]+=(1 if side==0 else -1)*phase(u,src,dst)
        gauss &= tuple(x%6 for x in g)==(1,1,1,0,0,0)
        record=tuple((x+d)%6 for x,d in zip(f,(1,-1,1,-1,1,-1)));w[m.index[(mask,record)],col]=1
    field=np.diag(electric)+.61*(2*np.eye(m.n)-w-w.T);h=hbare+field;herror=float(np.max(abs(h-archive['hamiltonian'])))
    spectra={};energies=[]
    for q in range(7):
        ids=np.flatnonzero(numbers==q);e,v=np.linalg.eigh(h[np.ix_(ids,ids)]);spectra[q]=(ids,e,v);energies.extend(e)
    globalenergy=float(min(energies));refs={};spectrumerror=0.
    for q,(_,e,_) in spectra.items():spectrumerror=max(spectrumerror,float(np.max(abs(e-saved['actual_all_number_spectra'][str(q)]['all_levels']))))
    # Direct7orbitalSWAP, including spectatorparity; no producerannihilator/transport.
    ks=[];boundary=True;fieldcomm=0.
    for initial in (0,1):
        for final in (0,1):
            k=np.zeros_like(h)
            for col,(mask,f) in enumerate(m.states):
                bit=mask&1
                if bit!=final:continue
                target=mask if initial==bit else mask^1
                sign=(-1 if initial else 1) if initial==bit else (-1)**((mask>>1).bit_count())
                row=m.index[(target,f)];k[row,col]=sign/np.sqrt(2)
                boundary &= initial-(3-mask.bit_count())==final-(3-target.bit_count())
            # Check commutation via field's diagonal+permutation action.
            fieldcomm=max(fieldcomm,float(np.max(abs(field@k-k@field))))
            ks.append((initial,final,k))
    allreferenceerror=0.;brancherror=0.;readouterror=0.;workerror=0.;nativeerror=0.;fieldworkerror=0.
    for name,selected_numbers in [('global_allN',range(7)),('fixed_N3',[3])]:
        refenergy=min(float(spectra[q][1][0]) for q in selected_numbers);columns=[]
        for q in selected_numbers:
            ids,e,v=spectra[q]
            for j in np.flatnonzero(abs(e-refenergy)<1e-10):
                column=np.zeros(m.n);column[ids]=v[:,j];columns.append(column)
        r=np.column_stack(columns)/np.sqrt(len(columns));rho=r@r.T
        allreferenceerror=max(allreferenceerror,float(np.max(abs(rho-archive[name+'_reference']))))
        responses=saved['references'][name]['complete_reservoir_response'];fieldafter=0.;totalwork=0.
        for (initial,final,k),record in zip(ks,responses):
            x=k@r;branch=x@x.T;prob=float(np.sum(x*x));energy=float(np.trace(x.T@h@x));fieldafter+=float(np.trace(x.T@field@x))
            brancherror=max(brancherror,float(np.max(abs(branch-archive[record['poststate_archive_key']]))))
            work=energy-refenergy*prob+.23*(final-initial)*prob;totalwork+=work
            workerror=max(workerror,abs(work-record['unnormalized_apparatus_work']))
            expected=record['complete_source_number_energy_records'];cursor=0
            for q,(ids,e,v) in spectra.items():
                weights=np.sum((v.T@x[ids,:])**2,axis=1);groups=[]
                for value,p in zip(e,weights):
                    if groups and abs(value-groups[-1][0])<1e-9:groups[-1][1]+=p
                    else:groups.append([value,p])
                for value,p in groups:
                    readouterror=max(readouterror,abs(p-expected[cursor]['probability']));cursor+=1
        fieldworkerror=max(fieldworkerror,abs(fieldafter-float(np.trace(r.T@field@r))))
        nativeworks={}
        for law in ('C','E','MIX'):
            energy=0.;prob=0.
            for ops in m.operators[law]:
                for _,k in ops:
                    x=k@r/np.sqrt(3);energy+=float(np.trace(x.T@h@x));prob+=float(np.sum(x*x))
            nativeworks[law]=energy-refenergy
            nativeerror=max(nativeerror,abs(nativeworks[law]-saved['references'][name]['native_instruments'][law]['work']),abs(prob-1))
        refs[name]={'energy':refenergy,'multiplicity':len(columns),'native_work':nativeworks,'reservoir_work':totalwork}
    conjugate_images=sum((63-mask,tuple((-x)%6 for x in f)) in m.index for mask,f in m.states)
    checks={'fixed_background_charge_conjugation_rejected':conjugate_images==0,'independent_interacting_H':herror<1e-12,'allnumber_spectra_reproduced':spectrumerror<1e-12,'allground_multiplets_recomputed':allreferenceerror<1e-12,'every_complete_reservoir_poststate':brancherror<1e-12,'allnumber_energy_outcome_probabilities':readouterror<1e-12,'same_reference_source_reservoir_work':workerror<1e-12,'all_native_C_E_MIX_work':nativeerror<1e-12,'dressed_Gauss_and_reservoir_integer_charge':bool(gauss and boundary),'reservoir_no_unconditional_field_work':fieldcomm<1e-12 and fieldworkerror<1e-12}
    checks={key:bool(value) for key,value in checks.items()}
    if not all(checks.values()):raise ValueError(checks)
    after={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    if before!=after:raise ValueError('sourcechanged duringreview')
    return {'schema':'oph.c10.joint_source_independent_review.v1','checks':checks,'references':refs,'globalenergy':globalenergy,'errors':{'H':herror,'spectrum':spectrumerror,'reference':allreferenceerror,'poststates':brancherror,'complete_readout':readouterror,'work':workerror,'native_work':nativeerror,'field_commutator':fieldcomm,'reservoir_field_work':fieldworkerror},'source_sha256':before,'scope':'Onefull384state interacting source, completeGaussbackground/zeroextraBenergy, fairreservoirlabels andactualreferences. No baregap/species/physicalvacuum/spin/Maxwell claim.','full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('JOINT_SOURCE_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({'checks':r['checks'],'references':r['references']}))
