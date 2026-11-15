"""Independent bit-parity reservoir swap and complete source/reference audit."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).parent;BASE=HERE.parents[1];sys.path.insert(0,str(HERE))
from charged_reservoir import BoundarySource


def run():
    source=BoundarySource(3);ops=source.reservoir_instrument();independent=[];errors=[];integer_guard=True
    for meta,k in ops:
        initial=meta['reservoir_initial'];observed=meta['reservoir_final'];out=np.zeros_like(k)
        for col,(mask,f) in enumerate(source.states):
            source_bit=mask&1
            if observed!=source_bit:continue
            # The exact7mode fermionic SWAP has phase−1 on bothoccupied,
            # and parity of intervening sourceorbitals for actual exchange.
            if initial==source_bit:target=mask;sign=-1 if initial else 1
            else:target=mask^1;sign=(-1)**((mask>>1).bit_count())
            row=source.index[(target,f)];out[row,col]=sign/np.sqrt(2)
            before_B=3-mask.bit_count();after_B=3-target.bit_count()
            integer_guard &= initial-before_B==observed-after_B
        errors.append(float(np.max(abs(out-k))));independent.append((meta,out))
    h=source.hamiltonian;n=np.array([mask.bit_count() for mask,_ in source.states]);ids=np.flatnonzero(n==3)
    values,vectors=np.linalg.eigh(h[np.ix_(ids,ids)]);ground=float(values[0]);selected=abs(values-ground)<1e-10
    rho=np.zeros((source.n,source.n));rho[np.ix_(ids,ids)]=vectors[:,selected]@vectors[:,selected].T/sum(selected)
    receipt=HERE/'CHARGED_RESERVOIR.json';saved=json.loads(receipt.read_text());complete=np.zeros_like(rho);brancherror=0.;probabilities={str(q):0. for q in range(7)};work=0.;energyerror=0.
    for (meta,k),record in zip(independent,saved['records']):
        branch=k@rho@k.T;complete+=branch
        brancherror=max(brancherror,float(np.max(abs(branch-np.array(record['poststate_real'])))))
        probability=float(np.trace(branch));energy=float(np.trace(h@branch))
        for q in range(7):probabilities[str(q)]+=float(np.trace(branch[np.ix_(n==q,n==q)]))
        work+=energy-ground*probability+.23*(meta['reservoir_final']-meta['reservoir_initial'])*probability
        energyerror=max(energyerror,abs(sum(x['energy']*x['probability'] for x in record['complete_number_energy_eigenspace_records'])-energy))
    gram=sum((k.T@k for _,k in independent),start=np.zeros_like(h))
    # Opposing CARphase convention is physically distinct even if normalized.
    coherent=np.arange(1,source.n+1,dtype=complex);coherent+=1j*coherent[::-1];coherent/=np.linalg.norm(coherent);coherent=np.outer(coherent,coherent.conj())
    correct=sum((k@coherent@k.T for _,k in independent),start=np.zeros_like(coherent));unsigned=sum((abs(k)@coherent@abs(k).T for _,k in independent),start=np.zeros_like(coherent));phaseeffect=float(np.linalg.norm(correct-unsigned,ord='nuc')/2)
    checks={'allfourCARswapblocks_exact':max(errors)<1e-12,'complete_all192_source_domain':bool(np.allclose(gram,np.eye(source.n))),
       'reservoir_boundary_integer_Gauss':bool(integer_guard),'all_source_modular_Gauss':all(source.gauss(s)==(1,1,1,0,0,0) for s in source.states),
       'complete_actual_ground_reference_poststates':brancherror<1e-12,'allnumber_probabilities_not_postselected':abs(sum(probabilities.values())-1)<1e-12,
       'outside234_actual_probability_zero':sum(probabilities[str(q)] for q in (0,1,5,6))<1e-12,
       'complete_global_energy_readout':energyerror<1e-12,'same_initial_reference_work':abs(work-saved['mean_apparatus_work'])<1e-12,
       'unsigned_swap_phase_control_rejected':phaseeffect>.01}
    if not all(checks.values()):raise ValueError(checks)
    paths=[Path(__file__),HERE/'charged_reservoir.py',receipt,BASE/'next/C9/fermionic_shared_records.py']
    return {'schema':'oph.c10.charged_reservoir_independent_review.v1','checks':checks,'Kraus_max_error':max(errors),'poststate_max_error':brancherror,'energy_ledger_error':energyerror,'actual_allN_probabilities':probabilities,'mean_source_reservoir_work':work,'phase_omission_trace_distance':phaseeffect,'phase_control_domain':'CoherentallN mathematicallyadmittedinput, includingcrossfermionparity. No physicalisolatedparitysuperpositionclaim; fixedN3 signed/unsignedK differonlyoutcomeglobalphase.',
       'scope':'Exact finite supplied reservoir pulse, fairinitialoccupationlabel retained. Sourcecharge exchanged withreservoir; totalcharge conserved includingexplicitboundary. Work uses unchangedsourceH plus.23nr, with no independentboundaryHamiltonian supplied. No physicalelectron/mass/vacuum or autonomousbattery.',
       'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'full_criterion_met':False,'official_rows_changed':[]}

if __name__=='__main__':
 r=run();HERE.joinpath('CHARGED_RESERVOIR_INDEPENDENT_REVIEW.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
