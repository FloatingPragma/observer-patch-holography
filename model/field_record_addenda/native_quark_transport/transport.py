"""Exact quark-spin-singlet first-hop calculation on the existing native cube."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
import sympy as s

from model.coherent_defect_wave.chiral_joint.interacting_charged_transport.algebra import routes, spin_effects
from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
from model.coherent_defect_wave.chiral_joint.atomic_momentum_response.model import reference_data
from model.coherent_defect_wave.chiral_joint.collective_spin_response.transfer import data as neutral_data
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import sqrt_interval, outward
from model.coherent_defect_wave.chiral_joint.peter_weyl_links.clebsch import link
from model.coherent_defect_wave.physical_attachment.chiral import representation
from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.model import certificate as preparation

CHANNELS = ('uc', 'dc', 'Q')
TIME = s.Rational(1, 2**22)
PARENTS = (
    'model/native_quark_response/results.json',
    'model/coherent_defect_wave/chiral_joint/interacting_charged_transport/results.json',
    'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/results.json',
    'model/coherent_defect_wave/chiral_joint/local_reference/results.json',
)


def pin(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_algebra():
    ref, _ = reference_data()
    local = json.loads((ROOT/PARENTS[-1]).read_text())
    irreps = [(row['transformed_quark_number'], module)
              for row in local['complete_quark_character_bound']['complete_occupation_characters']
              for module in row['modules'] if module['Casimir'] == '54']
    assert len(irreps) == 1
    number, module = irreps[0]
    assert number == 12 and module['dimension'] == 140 and module['multiplicity'] == 1
    assert local['exact_particle_hole_algebras']['quark_collective_color_and_spin_commutators'] == 11
    b, effects = spin_effects()
    assert all(effect.subs(b, 0) == s.eye(4)/128 for effect in effects)
    hyper = representation()[2][-1]
    assert all(hyper[a,a]+hyper[c,c] != 0 for a in range(15) for c in range(15))
    result = []
    for family in CHANNELS:
        indices = FAMILIES[family]; d = len(indices)
        assert ref['native_irrep_ground_occupations'][family]['exact_ground_mean_occupation'] == str(d)
        columns = [v for a in indices for c in indices
                   for v in (link(a,c)[:,0], link(a,c).H[:,0])]
        matrix = s.ImmutableSparseMatrix.hstack(*columns)
        assert matrix.H*matrix == s.eye(2*d*d)/d
        y6 = int(6*hyper[indices[0], indices[0]])
        result.append(dict(family=family, internal_dimension=d,
            source_hypercharge_change=str(-hyper[indices[0],indices[0]]),
            reference_six_Y=y6, arrival_six_Y=-y6,
            source_raw_B3_change=1 if family in ('uc','dc') else -1,
            source_raw_L_change=0, link_gram_dimension=2*d*d,
            link_gram='I/d', complete_one_body_compression='delta_spin delta_internal P_ground/2',
            complete_first_hop_effect='P0/128', moment='1/128'))
    return dict(unique_quark_module=module, transformed_number=number,
                three_axis_spin_effects='I_4/128', channels=result,
                actual_routes=routes())


def bounds(time=TIME):
    if not isinstance(time, s.Rational) or time < 0:
        raise ValueError('nonnegative exact rational source time required')
    n = neutral_data(); prep = preparation()['finite_energy_agreement_filter']
    lam, tau = n['coupling'], n['tau']; B0=s.Integer(576)
    assert lam == s.Rational(1, 2**34) and tau == 360
    eps = 2*s.Rational(prep['conditional_one_source_trace_distance_upper'])
    a_lower = 1-(n['errors']['invariant_graph_norm']+eps)**2
    q_lower = s.Rational(prep['success_probability_lower'])
    rt_lo,rt_hi = sqrt_interval((s.Rational(1,128),)*2,192)
    half_lo,half_hi = sqrt_interval((s.Rational(1,2),)*2,192)
    error = tau*time*(2*B0+lam*tau)/2
    margin = rt_lo-error*half_hi
    low = q_lower*a_lower*(lam*time)**2*max(0,margin)**2
    high = (lam*time)**2*(rt_hi+error*half_hi)**2
    assert 0 < a_lower <= 1 and 0 < q_lower <= 1
    nonarrival = q_lower*a_lower/2-6*high
    assert nonarrival > 0
    return dict(source_time=str(time), coupling=str(lam), full_hopping_norm=str(tau),
        centered_atomic_norm=str(B0), physical_preparation_success_lower=str(q_lower),
        localization_success_lower=str(a_lower), scaled_time_remainder=str(error),
        strict_amplitude_margin=str(margin),
        original_single_route_interval=outward((low,high),192),
        original_six_route_interval=outward((6*low,6*high),192),
        complete_nonarrival_lower=str(nonarrival),
        conditional_pump_probability='a/2',
        full_original_law=['1-q','q(1-a)','qa/2','q p_e for all six routes','q(a/2-sum_e p_e)'],
        coarse_records=10, all_environment_and_failure_records_retained=True,
        ground_trial_vector_subtraction_needed_for_first_hop=False,
        all_source_Gauss_constraints_with_fixed_reference=True,
        exact_idle_arrival='0', source_clock='exp(-it (H0+2^-34 T))',
        reference_clock='exp(-4it N_ref)')


def certificate():
    parents={}; inputs={}
    for name in PARENTS:
        value=json.loads((ROOT/name).read_text()); parents[name]=pin(ROOT/name)
        for path,expected in value['input_sha256'].items():
            assert pin(ROOT/path)==expected, path
            if path in inputs:assert inputs[path]==expected
            inputs[path]=expected
    from ..provenance import import_closure
    seeds=set(inputs)
    seeds.update(('model/coherent_defect_wave/chiral_joint/interacting_charged_transport/algebra.py',
                  'model/coherent_defect_wave/chiral_joint/interacting_ground_preparation/model.py'))
    for path in import_closure(seeds): inputs[path]=pin(ROOT/path)
    return dict(schema='oph.research.native_quark_transport.v1', algebra=source_algebra(),
        finite_records=bounds(), parent_sha256=parents, input_sha256=dict(sorted(inputs.items())),
        scope=dict(same_full_H_compensated_quark_type_transport=True,
                   energy_resolved_propagation=False, relativistic_dispersion=False,
                   measured_quark_identity_or_mass=False, GeV_calibration=False),
        owned_source_sha256={p.name:pin(p) for p in sorted(HERE.iterdir()) if p.suffix in ('.py','.md')})


if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true'); args=parser.parse_args()
    value=certificate(); path=HERE/'results.json'
    if args.write:path.write_text(json.dumps(value,indent=2,sort_keys=True)+'\n')
    else:assert json.loads(path.read_text())==value
    print('PASS: native quark transport; exact moment 1/128; every complete record retained')
    print(json.dumps(value['finite_records']['original_single_route_interval']))
