"""Exact neutral scalar pole and complete occupation/energy records."""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sympy as s
from model.defect_scalar_response.response import ground_pair, scalar_expectation
from model.defect_family_vertices.source import encoded_value
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent


@lru_cache(maxsize=1)
def build():
    data=ground_pair();h=data['block']['H'];n=data['n'];v=data['v'];field=data['field']
    assert n*n==n
    assert n*(v*v)==(v*v)*n
    assert h*h-h==v*v
    occupation=scalar_expectation(n)['value'];delta=data['delta'];lam=data['lam']
    assert occupation==(field.one-field.one/delta)/2
    weight=occupation*(1-occupation)
    assert weight==lam/(1+4*lam)
    first_moment=2*delta*weight
    # The literal capped scalar fields, rather than a continuum identity.
    phi=higgs_memory()[0]
    density=sum((x.H*x+x*x.H for x in phi),s.zeros(5))
    ns=s.diag(0,1,1,1,1)
    assert density==4*s.eye(5)-3*ns
    rows=[]
    for flag,p in [('occupied',occupation),('blank',1-occupation)]:
        for outcome,probability,energy in [('ground',p*p,field.zero),
             ('excited',weight,delta),('outside',field.zero,None)]:
            rows.append({'occupation_outcome':flag,'energy_outcome':outcome,
                'original_probability':encoded_value(probability),
                'source_energy_transfer':None if energy is None else encoded_value(energy)})
    assert occupation**2+(1-occupation)**2+2*weight==field.one
    # Connected density insertion has one positive line; its retarded
    # normalization is not the elementary charged scalar normalization.
    result={'schema':'oph.native_scalar_density_pole.v1',
        'source':{'CAR_modes':30,'scalar_states':5,
            'Hamiltonian':'full local unit native Yukawa action plus N_H',
            'reference':'local native-current singlet ground projector divided by two',
            'clock':'exp(-itH), hbar_model=1',
            'physical_source_dimension':'2^30*5','measured_mass_inputs':[]},
        'exact_identities':{'N_is_projection':True,'N_commutes_with_V_squared':True,
            'H_squared_minus_H_equals_V_squared':True,
            'pair_projector_is_exact_spectral_projector':True,
            'scalar_field_identity':'sum_k(phi_k^dagger phi_k+phi_k phi_k^dagger)=4I-3N_H'},
        'scalar_occupation':encoded_value(occupation),
        'connected_positive_pole_weight':encoded_value(weight),
        'positive_pole_source_energy':encoded_value(delta),
        'retarded_first_moment':encoded_value(first_moment),
        'retarded_response':'2*Delta*w/(z^2-Delta^2)',
        'normalized_coordinate':'(N_H-<N_H>)/sqrt(2*Delta*w)',
        'normalized_response':'1/(z^2-Delta^2)',
        'complete_occupation_energy_records':rows,
        'outside_pair_weight_exact':'0',
        'readout':{'instrument':'measure {N_H,I-N_H}, then complete H energy including outside-pair outcome',
            'postselection':False,'all_original_probabilities_sum_exactly_one':True,
            'unchanged_dual_reference_retained':True,
            'mean_energy_injected_by_nonselective_occupation_measurement':encoded_value(first_moment)},
        'scope':{'neutral_scalar_response_calculated':True,'propagation_calculated':False,
            'observed_Higgs_identification':False,'physical_rest_mass_GeV':None,
            'elementary_scalar_phi_poles_identical':False,
            'W_Z_photon_gluon_operators_added':False},
        'source_bindings':bindings()}
    return result


def bindings():
    names=['model/native_scalar_density_pole/pole.py','model/native_scalar_density_pole/PROOF.md',
           'model/defect_scalar_response/response.py','model/defect_scalar_response/results.json',
           'model/defect_family_vertices/source.py','model/defect_family_vertices/results.json',
           'model/coherent_defect_wave/chiral_higgs/model.py']
    out={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names}
    # The exact ground projectors and literal CAR operators are transitive
    # dependencies, checked against the source packages' declared closures.
    for package in ['defect_scalar_response','defect_family_vertices']:
        receipt=json.loads((ROOT/f'model/{package}/results.json').read_text())
        for p,digest in receipt['input_sha256'].items():
            candidate=ROOT/p
            if not candidate.is_file():raise ValueError('missing scalar ground input: '+p)
            actual=hashlib.sha256(candidate.read_bytes()).hexdigest()
            if actual!=digest:raise ValueError('stale scalar ground input: '+p)
            out[p]=actual
    return dict(sorted(out.items()))


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args()
    data=build();dest=HERE/'results.json'
    if a.check:
        if data!=json.loads(dest.read_text()):raise SystemExit('scalar density pole receipt differs')
    else:dest.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print('one exact neutral scalar pole; six complete occupation/energy records; no physical Higgs assignment')

if __name__=='__main__':main()
