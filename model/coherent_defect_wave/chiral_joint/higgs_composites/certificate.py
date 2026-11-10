"""Reproducible exact operator checks and complete periodic spectral diagnostics."""
from pathlib import Path
import hashlib
import json
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from .model import (tensor_certificate,compiler_certificate,noncanonical_certificate,
                    composite_dynamics_certificate,star_selection_certificate)
from .periodic_bridge import residue_certificate


ROOT=Path(__file__).resolve().parents[4]
CONSUMED_RECEIPT='model/coherent_defect_wave/chiral_joint/periodic_spectrum/results.json'


def validate_pins(pins):
    for name,expected in pins.items():
        path=ROOT/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('source custody mismatch: '+name)


def input_hashes():
    receipt=json.loads((ROOT/CONSUMED_RECEIPT).read_text())
    validate_pins(receipt['input_sha256'])
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(receipt['input_sha256'])
    seeds.update((CONSUMED_RECEIPT,'model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/field_limit/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/periodic_family/PROOF.md'))
    paths=set(import_closure(seeds))|seeds
    missing=sorted(name for name in paths if not (ROOT/name).is_file())
    if missing:raise ValueError('missing source inputs: '+', '.join(missing))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    pins=input_hashes()
    return {'schema':'oph.native_quantum_higgs_composites.v1',
            'input_sha256':pins,
            'consumed_periodic_spectrum_receipt_pins_validated':True,
            'exact_operator_identities':tensor_certificate(),'compiler_ancestry':compiler_certificate(),
            'noncanonical_CAR':noncanonical_certificate(),'local_star_selection':star_selection_certificate(),'quantum_composite_dynamics':composite_dynamics_certificate(),
            'same_periodic_H_spectral_bridge':residue_certificate(),
            'parameter_ancestry':{'operator_tensors':'actual native R15 representation',
                'local_Yukawa_coefficients':'independent symbolic real architectural parameters',
                'periodic_spectral_coefficients':'inherited y_e=omega=kappa=1, y_u=y_d=0; two quantum links',
                'empirical_mass_inputs':[],'empirical_mass_prediction':False},
            'scope':{'one_native_generation':True,'three_observed_generations_derived':False,
                     'full_3D_periodic_quantum_Higgs_EOM_operator_domain':True,
                     'periodic_residues_are_full3D_ground_residues':False,
                     'physical_interacting_continuum':False,'observed_particles_identified':False,
                     'GeV_masses':False}}


def write(path=None):
    path=Path(path) if path is not None else Path(__file__).with_name('results.json')
    result=certificate();path.write_text(json.dumps(result,indent=2)+'\n');return result


def replay_differences(expected,actual):
    """Exact custody/certified data; declared tolerance only for binary64 diagnostics."""
    expected=json.loads(json.dumps(expected))
    actual=json.loads(json.dumps(actual))
    differences=[]
    def compare(left,right,path=''):
        if type(left)!=type(right):differences.append(path+' (type)');return
        if isinstance(left,dict):
            if set(left)!=set(right):differences.append(path+' (keys)');return
            for key in left:compare(left[key],right[key],path+'/'+key)
        elif isinstance(left,list):
            if len(left)!=len(right):differences.append(path+' (length)');return
            for j,(a,b) in enumerate(zip(left,right,strict=True)):compare(a,b,path+'/'+str(j))
        elif isinstance(left,float):
            # Exact interval endpoints are strings and never enter this branch.
            import math
            if not math.isfinite(left) or not math.isfinite(right) or abs(left-right)>1e-10:
                differences.append(path+' (binary64 diagnostic)')
        elif left!=right:differences.append(path)
    compare(expected,actual)
    return differences
