"""Current-ring filled CAR reference and rejecting native-instrument lift test.

Reuses existing many_body_route exterior-matrix method; no Wilson symbol,
reservoir law, or physical vacuum is silently transferred.
"""
from pathlib import Path
import sys,json,hashlib,importlib.util
import numpy as np
import sympy as sp
BASE=Path(__file__).resolve().parents[2];FOCK=BASE/'validation/bridges/providers/fock.py'
sys.path.insert(0,str(BASE/'candidate'))
from oph_openwave.geometry import flat_ring
from oph_openwave.construction import Construction
from oph_openwave.uniform_native import UniformFiberReduction
_spec=importlib.util.spec_from_file_location('bundled_reference_fock',FOCK)
_fock=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_fock)
exterior_matrix=_fock.exterior_matrix


def reference(h):
    """Declared finite CAR attachment of unchanged H; implicit full 2^d carrier."""
    values,vectors=np.linalg.eigh(h)
    if min(abs(values))<1e-10:raise ValueError('gapped finite H required for unique filled reference')
    negative=values<0;covariance=vectors[:,negative]@vectors[:,negative].conj().T
    return values,vectors,covariance,float(values[negative].sum())


def run():
    from itertools import combinations
    model=UniformFiberReduction(flat_ring(1,1,.23),6);h=model.hamiltonian
    values,vectors,covariance,ground=reference(h);negative=values<0
    # Diagnostic one-body covariance action only, not a constructed CAR channel.
    mapped=sum((k@covariance@k.conj().T for _,k in model.kraus_for(0)),np.zeros_like(covariance))
    covariance_change=float(np.linalg.norm(mapped-covariance))
    one_body_excitation=float(np.trace(h@mapped).real-ground)
    tiny=UniformFiberReduction(Construction(3,((0,1),(1,2),(2,0)),.23),6)
    e,v,cov,g=reference(tiny.hamiltonian);n=int(sum(e<0));pairs=tuple(combinations(range(6),n))
    exact_lifts=[]
    for meta,k in tiny.kraus_for(0):
        exact=sp.Matrix([[sp.Rational(str(float(x))) for x in row] for row in k])
        exact_lifts.append(exterior_matrix(exact,n))
    exact_gram=sum((k.H*k for k in exact_lifts),sp.zeros(len(pairs)))
    defect=exact_gram-sp.eye(len(pairs));exact_not_tp=defect!=sp.zeros(len(pairs))
    gram=np.asarray(exact_gram).astype(float)
    neg=v[:,e<0];filled=np.array([np.linalg.det(neg[list(ids),:]) for ids in pairs])
    retained=float(np.vdot(filled,gram@filled).real)
    # Full 2^87 spectrum represented exactly by positive excitation costs, no sector deletion.
    costs=np.where(negative,-values,values)
    paths=[Path(__file__),BASE/'candidate/oph_openwave/uniform_native.py',BASE/'candidate/oph_openwave/geometry.py',FOCK,BASE/'validation/bridges/providers/mass_instrument_model.py']
    return {'schema':'oph.native.positive_reference_attachment.v1','checks':{'current_ring_H_unchanged':True,'filled_reference_projector':bool(np.max(abs(covariance@covariance-covariance))<1e-12),'positive_excitation_costs':bool(min(costs)>0),'one_body_native_reference_not_stationary':bool(covariance_change>1e-4),'naive_exterior_native_lift_rejected_exactly':bool(exact_not_tp),'naive_filled_reference_probability_not_normalized':bool(abs(retained-1)>.01),'fair_coin_weight_correction_still_not_TP':bool(16*exact_gram!=sp.eye(len(pairs)))},'attachment':{'statistics':'explicit supplied CAR, different from nonnegative Bose many_load_native','one_particle_generator':'unchanged canonical flat_ring1x1 Hm','modes':87,'complete_Fock_dimension':str(2**87),'negative_filled_modes':int(sum(negative)),'reference_energy':ground,'reference_covariance_real':covariance.real.tolist(),'reference_covariance_imag':covariance.imag.tolist(),'positive_excitation_costs':costs.tolist(),'charged_excitation_gap':float(min(costs)),'fixed_number_neutral_gap':float(min(values[~negative])-max(values[negative])),'positivity_identity':'dGamma(H)-Esea = sum_Epositive E n_E + sum_Enegative (-E)(1-n_E)','reference_scope':'unique finite filled CAR ground conditional on declared statistics; not physical vacuum/species or native preparation'},'native_preservation_diagnostic':{'seam':0,'one_body_covariance_change_frobenius':covariance_change,'formal_covariance_energy_increase':one_body_excitation,'mapped_covariance_eigenvalue_range':[float(np.linalg.eigvalsh(mapped).min()),float(np.linalg.eigvalsh(mapped).max())],'interpretation':'single-particle native channel acting on reference covariance is not stationary; not asserted as a complete CAR source channel'},'exact_rejecting_control':{'graph':'triangle','modes':6,'fixed_number':n,'sector_dimension':len(pairs),'existing_method':'many_body_route.fock.exterior_matrix on each original Kraus','sum_exterior_KdagK_diagonal':list(map(str,exact_gram.diagonal())),'completeness_defect_max':float(np.max(abs(gram-np.eye(len(pairs))))),'filled_reference_retained_probability':retained,'fair_coin_weight_corrected_completeness_defect_max':float(np.max(abs(16*gram-np.eye(len(pairs))))),'fair_coin_weight_corrected_filled_retained_probability':16*retained,'coin_weight_correction':'each one-particle K has amplitude1/2; replace wedge^3(K) by4*wedge^3(K) to apply the fair law/coin weight once, still incomplete','postselected_renormalization_used':False},'resources_and_missing':['CAR statistics and filled reference are supplied architecture, not derived from original one-load source','no compatible instrument is supplied by THIS naive exterior lift; direct configuration-level repair is next/C9/fermionic_shared_records.py','existing mass_instrument reservoir swap and energy-shell apparatus do not identify current native C/E with a sea-preserving operation','no pair creation/annihilation, physical electron identity, interacting vacuum, empirical mass prediction or physical clock'], 'full_criterion_met':False,'official_rows_changed':[],'source_sha256':{str(p.resolve().relative_to(BASE.resolve())):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}}

if __name__=='__main__':
    r=run();Path(__file__).with_name('POSITIVE_REFERENCE_ATTACHMENT.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['checks']))
