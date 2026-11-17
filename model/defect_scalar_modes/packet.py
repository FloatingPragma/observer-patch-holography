"""Source-bound complete scalar modes and native filled-sea response."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy import linalg as la

from .modes import build_scalar_source, selected_response
from .backreaction import calculate_backreaction, small_car_control

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TOLERANCE=2e-8


def source_bindings():
    names={str(p.relative_to(ROOT)) for p in HERE.iterdir()
           if p.is_file() and p.suffix in ('.py','.md')}
    names.update(('model/compact_defect_pair/PROOF.md',
                  'model/defect_mode_overlap/PROOF.md','model/requirements.txt','requirements.txt'))
    pending=list(names)
    def include(module):
        if not module.startswith(('model.','tools.')):
            return
        stem=module.replace('.','/')
        pieces=stem.split('/')
        candidates=[stem+'.py',stem+'/__init__.py']
        candidates += ['/'.join(pieces[:n])+'/__init__.py' for n in range(1,len(pieces))]
        for candidate in candidates:
            if (ROOT/candidate).is_file() and candidate not in names:
                names.add(candidate);pending.append(candidate)
    while pending:
        name=pending.pop()
        if not name.endswith('.py'):
            continue
        module=name[:-3].replace('/','.')
        package=module.removesuffix('.__init__') if module.endswith('.__init__') else module.rpartition('.')[0]
        for node in ast.walk(ast.parse((ROOT/name).read_text())):
            if isinstance(node,ast.Import):
                for alias in node.names:
                    include(alias.name)
            elif isinstance(node,ast.ImportFrom):
                module='.'*node.level+(node.module or '')
                target=importlib.util.resolve_name(module,package) if node.level else module
                include(target)
                for alias in node.names:
                    include(target+'.'+alias.name)
    return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(names)}


def archive_arrays(data=None,back=None):
    data=build_scalar_source() if data is None else data
    back=calculate_backreaction() if back is None else back
    return {'points':np.asarray(data['points']),'K':data['K'],'Omega':data['Omega'],
            'modes':data['modes'],'frequencies':data['frequencies'],
            'scalar_field_to_canonical':data['scalar_field_to_canonical'],
            'pair_Gram':back['kernel'],'particle_Gram':back['particle_Gram'],
            'antiparticle_Gram':back['antiparticle_Gram'],
            'pair_Gram_spectrum':back['pair_Gram_spectrum'],
            'dominant_pair_profile':back['dominant_pair_profile']}


def calculate():
    pins=source_bindings()
    data=build_scalar_source();back=calculate_backreaction()
    public={key:value for key,value in back.items() if not isinstance(value,np.ndarray)}
    uncoupled=build_scalar_source(False)
    frequencies=data['frequencies']
    result={'schema':'oph.defect_scalar_modes.v1','source_bindings':pins,
        'source':{'radius':3,'mass':64,'half_separation':1,'sites':343,
            'native_internal_labels':15,'fermion_CAR_modes':10290,
            'scalar_total_occupation_cap':1,'scalar_dimension':1373,
            'stiffness':'K=I+L_graph+M^2 diag(|n_x|^2)',
            'scalar_frequency':'Omega=positive_sqrt(K)',
            'primitive_Yukawa_channels':{'up':1,'down':1,'charged_lepton':1},
            'site_phase':'(-1)^(x1+x2+x3)','clock':'exp(-it(H0+H_Y)), hbar=1'},
        'complete_scalar_stiffness_eigenvalues':data['stiffness_eigenvalues'].tolist(),
        'complete_scalar_frequencies':frequencies.tolist(),
        'complete_scalar_modes_count':len(frequencies),
        'lowest_mode':selected_response(data),
        'even_odd_stiffness_split':float(data['stiffness_eigenvalues'][1]-data['stiffness_eigenvalues'][0]),
        'spectral_identities_numerical_errors':data['errors'],
        'vacuum_response':{'matrix':'2 Omega (z^2-K)^(-1)',
            'equal_time_commutator':0,'first_moment':'2 Omega',
            'canonical_coordinate':'q=(2 Omega)^(-1/2) phi',
            'canonical_retarded_response':'(z^2-K)^(-1)',
            'canonical_contact':'sqrt(2 omega_j) times primitive mode overlap',
            'global_cap_canonical_CCR_claimed':False,
            'reference':'scalar blank under free H_H; distinct from the interacting ground'},
        'texture_decoupled_control':{'stiffness':'I+L_graph','lowest_frequency':1.,
            'constant_normalized_profile':float(1/np.sqrt(343)),
            'response':selected_response(uncoupled)},
        'native_filled_sea_first_image':public,'independent_full_CAR_control':small_car_control(),
        'counting':{'complete_spatial_band_rank':4,'native_band_CAR_coordinates':60,
            'maximum_distinct_singular_strengths_of_one_real_profile':2,
            'physical_generation_count_derived':False},
        'scope':{'all_scalar_modes_and_all_fifteen_native_labels_retained':True,
            'full_negative_spatial_sea_used':True,'global_native_gauge_singlet':True,
            'local_quantum_gauge_links_or_Gauss_law':False,
            'empirical_inputs':[],'inserted_family_coefficients':False,
            'source_selected_free_scalar_mode':True,'canonical_free_vacuum_response':True,
            'strict_full_interacting_variational_inequality':True,
            'numeric_two_vector_Ritz_value_outward_certified':False,
            'full_interacting_ground_reconstructed':False,'dressed_Higgs_pole_reconstructed':False,
            'relativistic_particle_mass_or_physical_Yukawa_hierarchy':False,
            'three_generations_derived':False}}
    if source_bindings()!=pins:
        raise ArithmeticError('source changed during scalar calculation')
    return result


def validate_archive(packet):
    expected=archive_arrays()
    if set(packet)!=set(expected):
        raise ArithmeticError('complete scalar archive inventory differs')
    for key,target in expected.items():
        actual=packet[key]
        if (not isinstance(actual,np.ndarray) or actual.shape!=target.shape
                or not np.isfinite(actual).all()):
            raise ArithmeticError('invalid scalar array: '+key)
        if not np.allclose(actual,target,atol=TOLERANCE,rtol=0):
            raise ArithmeticError('source scalar array differs: '+key)
    modes=packet['modes'];w=packet['frequencies']
    if (la.norm(modes.T@modes-np.eye(343))>TOLERANCE
            or la.norm(packet['K']@modes-modes*w*w)>TOLERANCE):
        raise ArithmeticError('complete scalar eigensystem failed')


def compare(saved,fresh,path=()):
    if isinstance(fresh,dict):
        return isinstance(saved,dict) and set(saved)==set(fresh) and all(compare(saved[k],v,path+(k,)) for k,v in fresh.items())
    if isinstance(fresh,list):
        return isinstance(saved,list) and len(saved)==len(fresh) and all(compare(a,b,path) for a,b in zip(saved,fresh,strict=True))
    if isinstance(fresh,float):
        probability=any('probability' in key for key in path)
        allowed=2e-14+TOLERANCE*abs(fresh) if probability else TOLERANCE*max(1,abs(fresh))
        return type(saved) in (float,int) and np.isfinite(saved) and abs(saved-fresh)<=allowed
    return type(saved) is type(fresh) and saved==fresh


def write():
    result=calculate()
    (HERE/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(HERE/'scalar_modes.npz',**archive_arrays())
    return result


def check():
    saved=json.loads((HERE/'results.json').read_text())
    fresh=calculate()
    if not compare(saved,fresh):
        raise ArithmeticError('complete source-bound scalar receipt differs')
    with np.load(HERE/'scalar_modes.npz',allow_pickle=False) as packet:
        validate_archive(packet)
    return fresh


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=write() if args.write else check()
    print('Complete scalar modes, native sea response and strict variational bound verified:',result['complete_scalar_modes_count'],'modes')
