"""Complete charged-sector local bound, lattice reference and source custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from model.common_source.chiral_candidate import action
from ..higgs_composites.certificate import validate_pins
from . import local_bound,lattice,perturbation,ground_exchange,ground_records,resolved_exchange

ROOT=Path(__file__).resolve().parents[4]
CONSUMED=('model/coherent_defect_wave/chiral_joint/local_reference/results.json',
          'model/coherent_defect_wave/chiral_joint/peter_weyl_links/results.json')


def input_hashes():
    here=Path(__file__).parent
    seeds={str(path.relative_to(ROOT)) for path in here.iterdir() if path.suffix in ('.py','.md')}
    seeds.update(action.SOURCE_FILES);seeds.update(CONSUMED)
    for name in CONSUMED:
        source=json.loads((ROOT/name).read_text());validate_pins(source['input_sha256'])
        seeds.update(source['input_sha256'])
    seeds.update(('model/requirements.txt','model/coherent_defect_wave/chiral_joint/field_limit/PROOF.md'))
    paths=set(import_closure(seeds))|seeds
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def certificate_data():
    return {'schema':'oph.PW_atomic_reference_control.v1',
            'full_unconstrained_onsite_bound':local_bound.certificate(),
            'atomic_lattice_reference':lattice.atomic_certificate(),
            'larger_volume_atomic_control':lattice.atomic_certificate(side=4),
            'full_PW_hopping_controls':lattice.hopping_certificate(),
            'controlled_hopping_comparison':perturbation.certificate(),
            'actual_ground_two_hop_interaction':ground_exchange.second_moment_certificate(),
            'actual_denominator_resolved_spin_exchange':resolved_exchange.certificate(),
            'complete_unit_hopping_ground_quench_records':ground_records.certificate(),
            'parameter_ancestry':{'onsite_and_link_coefficients':'unchanged unit PW values',
                'controlled_comparison':'explicit lambda multiplying the complete actual kinetic operator',
                'quench_readout':'actual lambda=1, separate finite-time bound around H0 reference',
                'empirical_inputs':[],'mass_or_clock_calibration':False},
            'scope':{'all_unconstrained_local_Gauss_charge_sectors_bounded':True,
                'unit_link_penalty_atomic_ground_isolated':True,'all_source_states_retained':True,
                'small_hopping_stationary_low_band_control':True,
                'unit_hopping_stationary_interacting_vacuum':False,'thermodynamic_hopping_stability':False,
                'physical_particle_mass_or_relativistic_phase':False,'GeV_units':False}}


def certificate():
    pins=input_hashes()
    return {**certificate_data(),'input_sha256':pins,'consumed_receipt_pins_validated':True}


def write():
    result=certificate()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
