"""Exact charged response, pair resolvent and transitive source custody."""
import hashlib
import json
from pathlib import Path
from model.particle_common_contract.provenance import import_closure
from ..higgs_composites.certificate import validate_pins
from ..atomic_reference_control.resolved_exchange import certificate as exchange_control
from .model import separate_spin
from .resolvent import field_data,encoded,response,pair_integral,exchange

ROOT=Path(__file__).resolve().parents[4]
CONSUMED=('model/coherent_defect_wave/chiral_joint/atomic_reference_control/results.json',
          'model/coherent_defect_wave/chiral_joint/ground_composite_response/results.json',
          'model/coherent_defect_wave/chiral_joint/local_reference/results.json')


def input_hashes():
    seeds={str(p.relative_to(ROOT)) for p in Path(__file__).parent.iterdir() if p.suffix in ('.py','.md')}
    seeds.update(CONSUMED)
    for name in CONSUMED:
        source=json.loads((ROOT/name).read_text());validate_pins(source['input_sha256'])
        seeds.update(source['input_sha256'])
    seeds.add('model/requirements.txt');paths=seeds|set(import_closure(seeds))
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}


def matrix(value):return [[str(v) for v in row] for row in value.tolist()]


def polynomial(value):return [encoded(c) for c in value.rep.to_list()]


def certificate_data():
    e,f,field=field_data();rows=[]
    for family in ('ec','L'):
        data=response(family);result,proof=pair_integral(family)
        blocks=[]
        for spin,block in enumerate(data['blocks']):
            blocks.append({'native_removed_mode':block['native_removed_mode'],
                'target_twice_spin_weight':2*spin,
                'full_weight_dimension':block['full_weight_dimension'],
                'complete_highest_weak_irrep_dimension':block['complete_charged_irrep_dimension'],
                'Hamiltonian':matrix(block['H']),'metric':matrix(block['gram']),
                'complete_elementary_transfer':matrix(block['map']),
                'actual_metric_adjoint':matrix(block['adjoint'])})
        rows.append({'family':family,'complete_blocks':blocks,
            'complete_target_annihilator':str(data['q'].as_expr()),
            'complete_target_annihilator_degree':data['q'].degree(),
            'ground_response_reduced_pole_degree':data['reduced_denominator'].degree(),
            'signed_spin_response_moments':[encoded(v) for v in data['moments']],
            'rational_response_numerator':polynomial(data['numerator']),
            'rational_response_denominator':polynomial(data['denominator']),
            'reflected_denominator':polynomial(proof['reflected_denominator']),
            'exact_modular_inverse':polynomial(proof['inverse']),
            'denominator_resolved_response':polynomial(proof['evaluated_resolvent']),
            'pair_integral_coefficients_descending_in_E0':encoded(result),
            'exact_modular_pair_resolvent_identity':True})
    controls=exchange_control();actual=exchange()
    import sympy as s
    lower,upper=map(s.Rational,actual['exact_J_interval'])
    oldlo,oldhi=map(s.Rational,controls['certified_J_edge_interval'])
    assert oldlo<lower<upper<oldhi
    return {'schema':'oph.PW_charged_ground_exchange.v1',
        'actual_full_source_separate_spin_identities':separate_spin(),
        'algebraic_ground_field':{'minimal_polynomial':str(f.as_expr()),'root':'smallest real root',
            'coefficient_encoding':'descending powers of E0, exact rational strings'},
        'complete_charged_local_response':rows,'exact_weighted_exchange':actual,
        'strict_refinement_of_full_source_exchange_interval':True,
        'complete_effective_edge_operator':'-A_edge I+J_edge sum_i R_axis,ii S_x,i S_y,i',
        'scalar_part':'unchanged complete five-family integral and bound from atomic_reference_control',
        'all_other_edges_and_charged_intermediates':'retained; cross-edge second-order terms vanish by exact link-sector orthogonality',
        'full_effective_H_remainder':'inherited canonical graph bound 6 (abs(lambda) tau)^3/gap^2 on its stated range',
        'source_Hamiltonian_parameters_changed':False,'empirical_inputs':[],
        'scope':{'exact_all_time_quark_spin_response_cancellation':True,
            'exact_actual_ground_pair_energy_denominator':True,
            'full_unconstrained_local_source_retained':True,
            'quark_contribution_to_scalar_energy_discarded':False,
            'unit_hopping_ground_theorem_transferred':False,
            'particle_mass_or_relativistic_phase':False,'GeV_units':False}}


def certificate():
    pins=input_hashes()
    return {**certificate_data(),'input_sha256':pins,'consumed_receipt_pins_validated':True}


def write():
    result=certificate()
    Path(__file__).with_name('results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
