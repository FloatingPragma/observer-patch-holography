"""Independent analytic spectrum and complete readout audit; no producer rerun."""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np
from scipy.linalg import eigh
from sympy import Matrix,Symbol,Poly,expand

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[2]
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.geometry import CARRIER_SEAMS,RING_PORTS


def require(condition,message):
    if not condition:raise ValueError(message)


def independent_incidence():
    b=np.zeros((24,63),int)
    for sub in (0,1):
        for edge,(a,c) in enumerate(CARRIER_SEAMS):b[12*sub+a,30*sub+edge]=1;b[12*sub+c,30*sub+edge]=-1
    for edge,port in enumerate(RING_PORTS):b[port,60+edge]=1;b[12+port,60+edge]=-1
    return b


def native_comparison(plus,minus,ep,em,time,p,seams,ports):
    psi=(np.exp(-1j*ep*time)*plus+np.exp(-1j*em*time)*minus)/math.sqrt(2)
    rho=np.outer(psi,psi.conj());output=(1-p)*rho;records={('no_commit',):1-p}
    for e,(a,b) in enumerate(seams):
        for law in ('C','E'):
            for coin in (0,1):
                groups={}
                for pos in range(len(psi)):
                    d=int(pos==a)-int(pos==b);current=(d-(2*coin-1)*(abs(d)%2))//2
                    value=current if law=='C' else d;target=b if current==1 else a if current==-1 else pos
                    groups.setdefault(value,[]).append((pos,target))
                for value,pairs in groups.items():
                    src,dst=map(np.array,zip(*pairs));branch=p*rho[np.ix_(src,src)]/(4*len(seams))
                    output[np.ix_(dst,dst)]+=branch;records[(e,law,coin,value)]=float(np.trace(branch).real)
    return output,records


def run():
    path=HERE/'ENERGY_GEV.json';receipt=json.loads(path.read_text())
    for relative,pin in receipt['source_sha256'].items():
        require(hashlib.sha256((BASE/relative).read_bytes()).hexdigest()==pin,'producer source pin stale: '+relative)
    b=independent_incidence();lap=b@b.T;x=Symbol('x')
    expected=x*(x*x-10*x+20)**3*(x-6)**7*(x-8)**2*(x*x-10*x+22)**2*(x**4-18*x**3+102*x*x-188*x+60)
    polynomial=Matrix(lap.tolist()).charpoly(x).as_expr()
    require(Poly(expand(polynomial-expected),x).is_zero,'canonical exact characteristic polynomial mismatch')
    roots=sorted(float(z.real) for z in np.roots([1,-18,102,-188,60]))
    eigenvalues=[0.]+roots+[5-math.sqrt(5)]*3+[5+math.sqrt(5)]*3+[6.]*7+[8.]*2+[5-math.sqrt(3)]*2+[5+math.sqrt(3)]*2
    mass=receipt['parameters']['mass_native'];positive=[math.sqrt(mass*mass+v) for v in eigenvalues if v>1e-12]
    expected_H=np.array(sorted([mass]+[-mass]*40+positive+[-e for e in positive]))
    reported=np.array(receipt['all_energy_levels_native'])
    spectrum_error=float(np.max(abs(expected_H-reported)))
    require(spectrum_error<1e-12 and len(reported)==87,'all-branch analytic H spectrum mismatch')
    h=np.block([[mass*np.eye(24),b],[b.T,-mass*np.eye(63)]])
    archived=np.asarray(receipt['energy_eigenvectors_columns'])
    require(archived.shape==(87,87),'missing complete declared energy basis')
    basis_error=max(float(np.max(abs(archived.T@archived-np.eye(87)))),float(np.max(abs(h@archived-archived*reported))))
    require(basis_error<1e-11,'archived spectral PVM not complete or not eigenbasis')
    independent_energy,independent_vectors=eigh(h,driver='evr')
    scale=receipt['calibration']['gev_per_native_unit'];anchor=receipt['calibration']['reference_gev'];reference=receipt['calibration']['reference_native']
    require(abs(anchor-.000511)<1e-15,'declared benchmark is not511keV inGeV')
    require(abs(reference-mass)<1e-12 and abs(scale*reference-anchor)<1e-15,'anchor ancestry mismatch')
    require(np.max(abs(np.array(receipt['all_energy_levels_gev'])-scale*reported))<1e-14,'all-branch GeV conversion mismatch')
    seams=[tuple(e) for e in CARRIER_SEAMS]+[(12+a,12+b) for a,b in CARRIER_SEAMS]+[(p,12+p) for p in RING_PORTS]
    apparatus_error=0.;interference_error=0.;native_energy_sum_error=0.;native_record_sum_error=0.;native_record_error=0.;aggregate_energy_error=0.;mean_energy_error=0.;gap_error=0.
    pair_summaries=[]
    for pair in receipt['pairs']:
        lam=pair['laplacian_eigenvalue'];ep=pair['positive_energy_native'];em=pair['negative_energy_native'];gap=pair['gap_native']
        plus=np.array(pair['positive_state_amplitudes']);minus=np.array(pair['negative_state_amplitudes'])
        require(np.max(abs(h@plus-ep*plus))<1e-12 and np.max(abs(h@minus-em*minus))<1e-12,'archived prepared eigenstate fails independent H')
        initial=(plus+minus)/math.sqrt(2);opposite=(plus-minus)/math.sqrt(2)
        complement=np.eye(87)-np.outer(initial,initial)-np.outer(opposite,opposite)
        require(np.linalg.eigvalsh(complement).min()>-1e-12,'apparatus complement not positive')
        for sample in pair['interference_samples']:
            t=sample['native_time'];actual=(np.exp(-1j*ep*t)*plus+np.exp(-1j*em*t)*minus)/math.sqrt(2)
            probabilities=np.array([abs(np.vdot(initial,actual))**2,abs(np.vdot(opposite,actual))**2,np.vdot(actual,complement@actual).real])
            interference_error=max(interference_error,float(np.max(abs(probabilities-sample['plus_minus_complement_probabilities']))))
            apparatus_error=max(apparatus_error,abs(float(sum(probabilities))-1))
        native,records=native_comparison(plus,minus,ep,em,.37,.29,seams,24)
        recorded=pair['native_step_complete_record_probabilities'];worker_records={}
        for r in recorded:
            key=('no_commit',) if r.get('status')=='no_commit' else (r['seam'],r['law'],r['coin'],r['record_value'])
            require(key not in worker_records,'duplicate native outcome record')
            worker_records[key]=r['probability']
        require(records.keys()==worker_records.keys(),'complete local native outcome keys differ')
        native_record_error=max(native_record_error,max(abs(records[k]-worker_records[k]) for k in records))
        native_record_sum_error=max(native_record_sum_error,abs(sum(worker_records.values())-1))
        native_probs=np.array(pair['native_step_energy_probabilities'])
        require(len(native_probs)==87 and native_probs.min()>-1e-12,'invalid complete native energy distribution')
        native_energy_sum_error=max(native_energy_sum_error,abs(float(sum(native_probs))-1))
        mean_energy_error=max(mean_energy_error,abs(float(reported@native_probs)-float(np.trace(h@native).real)))
        other=np.einsum('ij,ji->i',independent_vectors.conj().T@native,independent_vectors).real
        centers=[]
        for value in reported:
            if not centers or abs(value-centers[-1])>1e-8:centers.append(value)
        for center in centers:
            aggregate_energy_error=max(aggregate_energy_error,abs(float(native_probs[abs(reported-center)<1e-8].sum())-float(other[abs(independent_energy-center)<1e-8].sum())))
        gap_error=max(gap_error,abs(gap-2*math.sqrt(mass*mass+lam)),abs(pair['gap_gev']-scale*gap),abs(pair['half_gap_gev']-scale*gap/2))
        pair_summaries.append({'lambda':lam,'gap_native':gap,'half_gap_GeV':pair['half_gap_gev'],'native_record_outcomes':len(recorded)})
    maximum=max(spectrum_error,apparatus_error,interference_error,native_energy_sum_error,native_record_sum_error,native_record_error,aggregate_energy_error,mean_energy_error,gap_error)
    require(maximum<1e-10,'independent readout concordance failed')
    require(receipt['checks']['scalar_shift_preserves_gaps'],'producer scalar-shift check false')
    sources=[path,Path(__file__).resolve(),BASE/'candidate/oph_openwave/geometry.py']
    return {'schema':'oph.native.energy_gev.independent_review.v1','outcome':'PASS_WITH_SCOPE_NOTES',
            'independent_method':'Exact symbolic canonical24port characteristic polynomial; analytic ±sqrt(m²+lambda)/−m multiplicities; archived-state phases; independently constructed coarse C/E channel and SciPy evr degenerate-energy projectors.',
            'exact_L_characteristic_polynomial':str(expected),'native_H_spectrum_error':spectrum_error,
            'native_H_negative_m_multiplicity':int(np.count_nonzero(abs(reported+mass)<1e-8)),
            'calibration':{'benchmark_511_keV_in_GeV':.000511,'measured_half_gap_native':reference,'GeV_per_native_unit':scale,'supplied_empirical_anchor':True},
            'errors':{'interference_probabilities':interference_error,'apparatus_normalization':apparatus_error,
                      'complete_native_energy_normalization':native_energy_sum_error,'complete_native_record_normalization':native_record_sum_error,
                      'independent_native_record_probabilities':native_record_error,'independent_energy_projector_probabilities':aggregate_energy_error,
                      'native_energy_expectation':mean_energy_error,'archived_energy_basis':basis_error,'gap_and_GeV_conversion':gap_error},
            'selected_pair_summary':pair_summaries,
            'scope_notes':['The511keV benchmark is supplied calibration, not a predicted electron mass.',
                           'All negative/positive energies refer to a declared generator zero; gap/half-gap anchor is scalar-shift invariant.',
                           'Explicit coherent ± eigenvector superpositions and rankone probes make an interferometric gap observable; no claim that port-diagonal readout sees2m.',
                           'Per-vector energy probabilities within a degeneracy require a chosen fine basis; independent audit compares basis-invariant degenerate-energy projectors.',
                           'No physical vacuum, particle species, infinite-volume pole, relativistic mass or self-starting clock is identified.',
                           'Native C/E comparison is a separate declared Bernoulli step; no energy-conservation assertion through monitoring.'],
            'full_physical_criterion_met':False,'official_rows_promoted':0,'producer_rerun':False,
            'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}


if __name__=='__main__':
    result=run();(HERE/'INDEPENDENT_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'outcome':result['outcome'],'spectrum_error':result['native_H_spectrum_error'],'errors':result['errors']}))
