"""Complete atomic-bound replay and fail-closed quartic-remainder custody."""
import hashlib
from pathlib import Path

from model.common_source.chiral_candidate import action
from model.particle_common_contract.provenance import import_closure
from ..atomic_reference_control.lattice import atomic_certificate,hopping_certificate,gap_data,hopping_norm
from .bounds import coefficient_certificate,quartic_remainder
from .symmetry import certificate as symmetry_certificate

ROOT=Path(__file__).resolve().parents[4]
PACKAGE='model/coherent_defect_wave/chiral_joint/atomic_even_remainder/'


def input_hashes():
    seeds={PACKAGE+name for name in ('__init__.py','__main__.py','bounds.py','symmetry.py',
        'certificate.py','test_remainder.py','PROOF.md','README.md')}
    seeds.update(action.SOURCE_FILES)
    seeds.update(('model/requirements.txt',
        'model/coherent_defect_wave/chiral_joint/atomic_reference_control/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/atomic_reference_control/results.json',
        'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
        'model/coherent_defect_wave/chiral_joint/local_reference/results.json'))
    for name in sorted(seeds):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing declared even remainder input: '+name)
    paths=import_closure(seeds)
    for name in sorted(paths):
        if not (ROOT/name).is_file():raise FileNotFoundError('missing even remainder import input: '+name)
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate():
    atomic=atomic_certificate();hopping=hopping_certificate();symmetry=symmetry_certificate()
    gap=gap_data()['gap'];tau=hopping_norm();endpoint=gap/(5*tau)
    cases=[]
    for coupling in (0,endpoint,-endpoint,endpoint/16):
        result=quartic_remainder(gap,tau,coupling)
        cases.append({key:str(value) for key,value in result.items()})
    return {'schema':'oph.full_Q_even_canonical_remainder.v1',
        'atomic_physical_reference':atomic,'full_carrier_hopping_control':hopping,
        'actual_source_parity':symmetry,'uniform_coefficient':coefficient_certificate(),
        'canonical_operator':'K_lambda=J_lambda^dagger (H_lambda-volume*E0) J_lambda',
        'complete_physical_complement':'Q=I_physical-P0',
        'second_order_operator':'K2=-P0 T Q [Q(H0-volume*E0)Q]^(-1) Q T P0',
        'full_Q_even_resolvent_identity':'Ravg-R0=D0^(-1/2) C^2(I-C^2)^(-1) D0^(-1/2), C=D0^(-1/2) Delta D0^(-1/2)',
        'controlled_remainder':'||K_lambda-lambda^2 K2|| <= v^4[1/(g(g^2-v^2))+12/(g-v)^3] <=25v^4/g^3',
        'norm_domain':'v=|lambda|tau<=g/5','admitted_source_hopping_interval':[str(-endpoint),str(endpoint)],
        'exact_bound_cases':cases,
        'stable_bound_API':'bounds.uniform_quartic_remainder(gap,tau,coupling)',
        'stable_gap_coupling_API':'bounds.choose_quartic_coupling(gap,tau,second_order_gap)',
        'gap_coupling_hypothesis':'the caller supplies a proved positive second-order spectral gap delta2',
        'gap_coupling_conclusion':'lambda^2 delta2-2 R4 >=lambda^2 delta2/2>0',
        'unique_interacting_ground_from_parity_alone':False,
        'unit_hopping_or_larger_spectral_lattice_transferred':False,
        'all_charged_intermediate_states_retained':True,'numerical_eigensolver_used':False,
        'input_sha256':input_hashes()}
