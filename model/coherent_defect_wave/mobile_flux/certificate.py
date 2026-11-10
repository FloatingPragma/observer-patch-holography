"""Exact finite records and source custody for mobile paired flux."""
import hashlib
from pathlib import Path
import sympy as s

from model.particle_common_contract.provenance import import_closure
from .model import (observer,PairState,positive,insert_pair,exchange,charge,
                    charge_measurement,constraint_branches,onsite,MODES,mode,STAR)

ROOT=Path(__file__).resolve().parents[3]


def hashes():
    here=Path(__file__).parent
    paths={str(p.relative_to(ROOT)) for p in here.iterdir()
           if p.suffix in ('.py','.md') and not p.name.startswith('test_')}
    paths.update(('model/requirements.txt','model/coherent_defect_wave/species/PROOF.md',
                  'model/coherent_defect_wave/spatial/PROOF.md',
                  'model/coherent_defect_wave/spatial/ROUTE.md',
                  'evidence/native_geometric_source_20260925/spec.json',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
                  'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py'))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in sorted(import_closure(paths)) if (ROOT/name).is_file()}


def experiment():
    law=observer().base.law
    words=[{'word':g,'regional_charge':charge(g),'sphere_defect_count':law.syndrome(g).bit_count()}
           for g in observer().base.link_words]
    assert all(row['sphere_defect_count']==2*row['regional_charge'] for row in words)
    raw=PairState({((0,0),()):1});vac=positive(raw);pair=insert_pair(vac,0)
    moved=exchange(pair)
    charge_records={''.join(map(str,q)):str(v.norm_squared()) for q,v in charge_measurement(moved)}
    assert charge_records=={'00':'0','01':'1/2','10':'1/2','11':'0'}
    final=charge_measurement(exchange(moved))
    assert dict((q,z.norm_squared()) for q,z in final)[0,1]==1
    flags={str(q):str(z.norm_squared()) for q,z in constraint_branches(raw)}
    assert all(value=='1/4' for value in flags.values())
    # Full native wave, on the receiver, accesses the very link word moved by E.
    probe=PairState({((1,0),(MODES+mode(p,3,0,a),)):1 for p in (0,1) for a in range(60)},120)
    result=onsite(exchange(positive(probe)),1)
    assert result.norm_squared()==1
    records={}
    for receiver_charge in (0,1):
        for yes in (False,True):
            selected=PairState({key:v for key,v in result.amplitudes
                               if charge(key[0][1])==receiver_charge
                               and ((MODES+mode(0,3,0,0) in key[1])==yes)},result.norm_denominator)
            records[f'Q_receiver={receiver_charge},pointer={yes}']=str(selected.norm_squared())
    assert list(records.values())==['119/240','1/240','1/2','0']
    return {'program':'coherent-mobile-paired-flux-v1','link_words':words,
            'local_one_wave_dimension':observer().base.dimension,
            'old_traces_or_reference_constants_transferred':False,
            'two_cell_Gauss_preparation_flags':flags,
            'successful_records_require_preparation_probability_factor':'1/4',
            'preparation_failure_records_retained':True,
            'charge_position_records_after_one_exchange':charge_records,
            'two_exchanges_move_pair_with_probability':'1',
            'receiver_native_wave_joint_records':records,
            'receiver_output_amplitude_support':len(result.amplitudes),
            'transported_object':'inside flux plus exterior compensator and whole carrier boundary data',
            'individual_old_carrier_boundary_kept_fixed':False,
            'separate_optional_H_flux_band':'2+(pi/2)*sum_i(1-cos(k_i))',
            'H_flux_GeV_or_relativistic_mass_claim':False}


def certificate():
    from .routing import certificate as route_certificate
    from .lane_collision import certificate as collision_certificate
    return {'schema':'oph.coherent_defect_wave.mobile_flux.v1',
            'paired_flux':experiment(),'protected_object_route':route_certificate(),
            'complete_lane_collision':collision_certificate(),
            'input_sha256':hashes()}
