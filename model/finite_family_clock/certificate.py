"""Independent calculation, immutable parent validation and source closure."""
import hashlib
import json
import math
from pathlib import Path

from model.common_source.chiral_candidate.action import SOURCE_FILES
from model.particle_common_contract.provenance import import_closure
from .family import calculate as family
from .memory import calculate as memory

ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
HIGGS_PARENT='model/coherent_defect_wave/chiral_higgs/results.json'
SEAM_PARENT='model/coherent_seam_repair/manifest.json'


def validate_pins(pins):
    for name,expected in pins.items():
        relative=Path(name)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('source custody requires repository-relative paths')
        path=ROOT/relative
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('source custody mismatch: '+name)


def input_hashes():
    higgs=json.loads((ROOT/HIGGS_PARENT).read_text())['input_sha256']
    seam=json.loads((ROOT/SEAM_PARENT).read_text())['files']
    seam={'model/coherent_seam_repair/'+name:value for name,value in seam.items()}
    validate_pins(higgs)
    validate_pins(seam)
    seeds={str(path.relative_to(ROOT)) for path in HERE.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(higgs)
    seeds.update(seam)
    seeds.update(SOURCE_FILES)
    seeds.update((HIGGS_PARENT,SEAM_PARENT,'model/requirements.txt'))
    names=import_closure(seeds)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(names)}


def certificate():
    pins=input_hashes()
    return {'schema':'oph.finite_family_clock.v1','input_sha256':pins,
            'immutable_parent_pins_validated':[HIGGS_PARENT,SEAM_PARENT],
            'family_instrument_and_collision_clock':family(),
            'same_H_quantum_memory':memory(),
            'scope':{'primary_common_source_transfer':False,
                     'observed_generations_derived':False,'physical_rest_mass_or_GeV':False,
                     'neutral_background_is_quantum_Higgs_state':False,
                     'empirical_parameters_fitted_here':False}}


def replay_differences(expected,actual):
    # Compare the serialized receipt types, including NumPy binary64
    # diagnostics whose JSON representation is an ordinary Python float.
    expected=json.loads(json.dumps(expected))
    actual=json.loads(json.dumps(actual))
    differences=[]
    def compare(left,right,path=''):
        if type(left)!=type(right):
            differences.append(path+' (type)');return
        if isinstance(left,dict):
            if set(left)!=set(right):
                differences.append(path+' (keys)');return
            for key in left:compare(left[key],right[key],path+'/'+key)
        elif isinstance(left,list):
            if len(left)!=len(right):
                differences.append(path+' (length)');return
            for i,(a,b) in enumerate(zip(left,right,strict=True)):compare(a,b,path+'/'+str(i))
        elif isinstance(left,float):
            # Strings carry exact symbolic and high-precision deterministic
            # values. Only binary64 diagnostics reach this tolerance branch.
            if not math.isfinite(left) or not math.isfinite(right) or abs(left-right)>1e-10:
                differences.append(path+' (binary64 diagnostic)')
        elif left!=right:
            differences.append(path)
    compare(expected,actual)
    return differences


def write(path=None):
    value=certificate()
    path=Path(path) if path is not None else HERE/'results.json'
    path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    return value
