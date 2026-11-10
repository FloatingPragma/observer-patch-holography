"""One finite native Wilson–Maxwell source with complete interaction records."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from . import current,matter,bounds,field,gauss,exact_current as current_audit
from model.native_higgs_filtered_state import model as filtered
from model.native_block_dirac.resources import bounds as interval_bounds
from model.particle_common_contract.provenance import import_closure

PARENT='model/native_higgs_filtered_state/results.json'
PARENT_SHA='19527188c4ffa855df73330c5fc1d9efead08308e0c7719e6f67fd845e901029'
PEER_SHA='f646ac992fc589a042b1ddb635f6d7f288943d7bec15aeca846badda4608cb7e'
RESEARCH=HERE/'research_reference.json'
RESEARCH_SHA='2fcb312421565591e3a62308986b98d3803519320949185ef72e9c8f19c4495e'
FIELD_REFERENCE=HERE/'field_research_reference.json'


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def serial(value):return json.loads(json.dumps(value,default=str))


@lru_cache(None)
def admit():
    if digest(ROOT/PARENT)!=PARENT_SHA:raise ValueError('native shared-energy source differs')
    saved=json.loads((ROOT/PARENT).read_text());filtered.higgs.validate(saved['input_sha256'])
    if filtered.certificate()!=saved:raise ValueError('fresh same-source shared-energy result differs')
    return saved


@lru_cache(None)
def admit_field():
    path=FIELD_REFERENCE
    if digest(path)!=PEER_SHA:raise ValueError('frozen field/Gauss/current prerequisite differs')
    saved=json.loads(path.read_text());filtered.higgs.validate(saved['input_sha256'])
    return saved


def resources():
    source=admit();energy=source['typed_free_Dirac_corollary']['source_excitation_interval'];upper=F(energy[1])
    g=F(1,64);raw=32*(600+F(360,2**34))+99*upper+18*upper*g+F(923,4)*upper*g*g+128
    mp.iv.dps=90;log=mp.iv.ln
    capacity=32*(256*log(2)+8*log(5)+24*log(119)+log(2**31+1))+7*log(9)+log(65)+log(56)+log(1161)
    return {'source_selected_excitation_interval':energy,
        'selector':'positive elementary/Higgs cross residue maximum within the refined radius2v absolute source bin',
        'source_constituents':32,'global_exact_code_dimension':16**8,
        'full_native_constituent_kinematic_dimension':'2^240*5^8*119^24*(2^31+1)',
        'full_native_source_reference_trial_capacity_nats_formula':'32*(240ln2+8ln5+24ln119+ln(2^31+1)+16ln2)+7ln9+ln65+ln56+ln1161',
        'full_native_source_reference_trial_capacity_nats_interval':interval_bounds(capacity),
        'global_spectral_code_compression_is_declared_architecture':True,
        'full_raw_H_ext_norm_upper':str(raw),
        'Gauss_orthogonal_H_extension':'zero; hence identity autonomous unitary on those retained inputs',
        'native_physical_code_complement_H':'unchanged native sum plus finite field, full4|R| and charge-diagonal Coulomb; no kinetic/current coupler',
        'preparation_duration_source_units':'1','autonomous_wait':'1/(64m)',
        'readout_duration_source_units':'1','full_elapsed_source_time':'2+1/(64m)',
        'preparation_control':'i Log(S)-H_ext, explicitly cancelling autonomous H_ext during the declared pulse',
        'readout_control':'i Log(exp(-i H_ext) U_count)-H_ext; final common free source evolution retained',
        'Hermitian_log_branch':'phases in (-pi,pi]',
        'each_control_norm_upper':str(raw+F(22,7)),
        'net_source_work_absolute_upper':str(2*raw),
        'gross_two_pulse_switching_work_upper':str(4*(raw+F(22,7))),
        'arbitrary_pulse_error_target_probability_lower':'(sqrt(p0^32)*a_star-eta)_+^2',
        'arbitrary_pulse_error_connected_lower':'(sqrt(p0^32)*a_star-eta)_+^2-(tau*nu+eta)^4',
        'controller_coefficient_workspace_in_capacity':False,
        'efficient_native_controller_or_passive_repair_compilation_claimed':False}


def calculate():
    admit();field_data=serial(field.certificate());gauss_data=serial(gauss.calculate())
    independent_H,independent_currents=current_audit.literal_matter()[:2]
    actual_H,actual_currents=matter.cube_matter()[:2]
    assert independent_H==actual_H and tuple(independent_currents)==tuple(actual_currents)
    intervals=matter.filled_current();exact=serial(current_audit.calculate())
    peer=admit_field()
    assert field_data==peer['finite_field']
    assert gauss_data==peer['full_Gauss_embedding']
    assert exact==peer['independent_exact_current']
    ilow,ihigh=map(F,intervals['pair_one_photon_V1_vacuum_norm_squared_over_m_squared_interval'])
    elow,ehigh=map(F,exact['outward_interval'])
    assert ilow<elow<ehigh<ihigh
    return {'source':'actual native constituents plus explicitly declared open Wilson, finite cochain Maxwell and common Gauss dressing',
        'shared_field_prerequisite_sha256':PEER_SHA,
        'source_choices':{'g':'1/64','a':'1/m','clock':'one full H(g), including Coulomb and complete finite Peierls matrices',
            'reference_H':'4|R| on all65 retained static compensator states',
            'mass_m':'actual refined positive-cross source eigenenergy minus actual native Eg',
            'source_empirical_mass_or_coupling_calibration':False},
        'exact_current_interface':current.calculate(),'field_source':field_data,'full_Gauss_reference':gauss_data,
        'literal_trial_support':bounds.support_control(),'actual_preparation':bounds.preparation(),
        'full_charge_frame_convention':bounds.charge_frame_control(),
        'independent_current_coefficient':{'uniform_polynomial_method':intervals,
            'exact_CRT_method':exact,'strict_outward_interval_inclusion':True,
            'same_full_current_and_cochain_matrices_used':True},
        'complete_interaction_records':bounds.response(),
        'count_instrument':{'code_labels':'(N_a,N_b,N_ph) in {0,...,16}^2 x {0,...,3}',
            'native_complement_labels':'(complement,N_ph), N_ph=0,...,3',
            'physical_count_alphabet_size':1160,'full_pointer_dimension':1161,'target_label':[1,1,1],
            'projector_sum':'sum_physical Y P_label Y^dagger+(I-Y Y^dagger)=I on the complete kinematic carrier',
            'physical_detector':'Y P_label Y^dagger plus the distinct Gauss-orthogonal failure effect',
            'Gauss_failure_probability_on_declared_preparation':'0',
            'coupled_ground_vacuum_scattering_or_radiation_asserted':False},
        'clock_and_capacity':resources(),
        'boundaries':{'actual_finite_source_interaction_proved':True,
            'all_state_native_spatial_locality_transferred':False,
            'interacting_Lorentz_or_four_dimensional_QED_limit_proved':False,
            'SM_particle_or_weak_chirality_identification':False,
            'infinite_photon_carrier_used':False,
            'full_Fock_Hamiltonian_numerically_enumerated':False}}


def certificate():
    result=calculate();parent=admit();seeds=set(parent['input_sha256'])|{PARENT}
    seeds.update(admit_field()['input_sha256'])
    seeds.update(('model/particle_common_contract/provenance.py',
                  'model/particle_physical_causality/PROOF.md',
                  'model/particle_field_causality/PROOF.md',
                  'model/maxwell_local_selection/PROOF.md','model/maxwell_weak_field/PROOF.md'))
    seeds.update(str(p.relative_to(ROOT)) for p in HERE.iterdir() if p.suffix in ('.py','.md'))
    seeds.update((str(RESEARCH.relative_to(ROOT)),str(FIELD_REFERENCE.relative_to(ROOT))))
    result.update({'schema':'oph.native_wilson_maxwell.v1',
        'input_sha256':{name:digest(ROOT/name) for name in sorted(import_closure(seeds))},
        'owned_source_sha256':{p.name:digest(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')}})
    if digest(RESEARCH)!=RESEARCH_SHA:raise ValueError('immutable combined source calculation differs')
    ignored={'schema','input_sha256','owned_source_sha256','shared_field_source_sha256'}
    science=lambda value:{k:v for k,v in value.items() if k not in ignored}
    if science(result)!=science(json.loads(RESEARCH.read_text())):
        raise ValueError('fresh public coupled-source science differs from frozen research')
    result['independent_research_comparison']={
        'reference_sha256':RESEARCH_SHA,'shared_field_reference_sha256':PEER_SHA,
        'all_scientific_fields_equal':True,'ignored_only_provenance_fields':sorted(ignored),
        'both_independent_current_methods_and_all_source_record_bounds_replayed':True}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    args=parser.parse_args();path=HERE/'results.json'
    if args.check:
        saved=json.loads(path.read_text());filtered.higgs.validate(saved['input_sha256'])
        for name,value in saved['owned_source_sha256'].items():
            if digest(HERE/name)!=value:raise ValueError('declared source file differs: '+name)
    actual=certificate()
    if args.write:path.write_text(json.dumps(actual,indent=2,sort_keys=True)+'\n')
    elif actual!=saved:raise ValueError('complete native Wilson–Maxwell record differs')
    print('PASS: shared finite native/Dirac/Maxwell source, complete Gauss and nonzero connected records')


if __name__=='__main__':main()
