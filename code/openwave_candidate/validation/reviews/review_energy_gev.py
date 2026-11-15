"""Independent apparatus, calibration-input and one-pair native-outcome audit."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent;BASE=HERE.parents[1]
PRODUCER=BASE/'validation/observables/energy_gev/run_energy_gev.py'
spec=importlib.util.spec_from_file_location('energy_gev_review_source',PRODUCER)
e=importlib.util.module_from_spec(spec);spec.loader.exec_module(e)


def run():
    receipt=PRODUCER.with_name('ENERGY_GEV.json');data=json.loads(receipt.read_text());pair=data['pairs'][0]
    model=e.UniformFiberReduction(e.flat_ring(1,1,data['parameters']['mass_native']),6)
    energies,vectors=np.linalg.eigh(model.hamiltonian)
    plus=np.asarray(pair['positive_state_amplitudes'],complex);minus=np.asarray(pair['negative_state_amplitudes'],complex)
    initial=(plus+minus)/np.sqrt(2);opposite=(plus-minus)/np.sqrt(2)
    pplus=np.outer(initial,initial.conj());pminus=np.outer(opposite,opposite.conj());complement=np.eye(model.n)-pplus-pminus
    complement_min=float(np.linalg.eigvalsh(complement).min())
    basis_overlap=float(abs(np.vdot(plus,minus)))
    # One native step with complete explicit joint(record, spectral-basis-label).
    duration=.37;probability=.29
    unitary=(vectors*np.exp(-1j*duration*energies))@vectors.T
    evolved=unitary@initial
    branches=[np.sqrt(1-probability)*evolved]
    for seam in range(len(model.maps)):
        for _,src,dst in model.maps[seam]:
            branch=np.zeros(model.n,complex);branch[dst]=math.sqrt(probability/(4*len(model.maps)))*evolved[src]
            branches.append(branch)
    joint=np.array([abs(vectors.T@branch)**2 for branch in branches])
    record_marginal=joint.sum(axis=1);spectral_marginal=joint.sum(axis=0)
    observed,ledger=model.step(np.outer(initial,initial.conj()),duration,probability)
    independently_measured=np.diag(vectors.T@observed@vectors).real
    error=float(np.max(abs(spectral_marginal-independently_measured)))
    max_correlation=float(np.max(abs(joint-record_marginal[:,None]*spectral_marginal[None,:])))
    if error>1e-12 or abs(joint.sum()-1)>1e-12 or complement_min<-1e-12:raise ValueError('spectral/native apparatus control failed')
    scale_probes=[]
    for native,gev in ((1e-300,1e300),(1e300,1e-300),(True,.000511)):
        try:
            with np.errstate(all='ignore'):value=e.gev_scale(native,gev)
            valid=math.isfinite(value) and value>0 and type(native) is not bool
            scale_probes.append({'inputs':[native,gev],'status':'finite_positive_domain' if valid else 'unsupported_input_accepted','value':str(value)})
        except (ValueError,OverflowError,TypeError) as exc:scale_probes.append({'inputs':[native,gev],'status':'rejected','type':type(exc).__name__})
    return {'schema':'oph.energy_gev.adversarial_review.v1','source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (PRODUCER,receipt,Path(__file__))},
       'povm':{'basis_overlap':basis_overlap,'complement_min_eigenvalue':complement_min,'dimension':model.n,
               'scope':'Rank-one cross-branch superposition apparatus, not a port/current measurement.'},
       'one_pair_native_probe':{'joint_shape':list(joint.shape),'normalization':float(joint.sum()),
          'joint_vs_channel_spectral_marginal_error':error,'joint_vs_product_of_marginals_max_difference':max_correlation,
          'scope':'A complete sequential native-record plus spectral apparatus has a joint law; two separate marginals do not determine it.'},
       'calibration_input_probes':scale_probes,'producer_scope_supported':['GeV outputs explicitly calibrated','negative branches retained','vacuum/species/mass prediction not claimed'],
       'unresolved_resources':['Uniform modular-fiber preparation','relative positive/negative branch coherence','complete nonlocal spectral-basis or energy-projector apparatus','physical vacuum-to-particle and common clock/length dictionary'],
       'official_rows_promoted':0,'full_producer_run_repeated':False}

if __name__=='__main__':
    result=run();(HERE/'ENERGY_GEV_ADVERSARIAL_REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))
