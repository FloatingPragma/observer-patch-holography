"""Exact scalar-cost derivatives of full native CAR energy-transfer moments."""
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import sympy as s

from model.coherent_defect_wave.chiral_joint.ground_composite_response.model import source
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import (
    encoded, field_data, polynomial,
)
from model.coherent_defect_wave.chiral_joint.charged_ground_exchange import model as lepton
from model.coherent_defect_wave.chiral_joint.family_cell.quotient import evaluate, outward
from model.native_quark_response import local as quark

ROOT = Path(__file__).resolve().parents[2]
CENTRE_CLASS = {'uc': 2, 'dc': 2, 'Q': 1, 'ec': 0, 'L': 3}


def scalar_projector(block):
    """Restrict the literal five-state scalar occupation, checking its full image."""
    raw = s.diag(*(int(key[2] != 0) for key in block['keys']))
    image = raw*block['basis']
    left = block.get('left', s.eye(block['basis'].cols))
    reduced = s.ImmutableSparseMatrix(left*image.extract(
        block['coordinates'], range(block['basis'].cols)))
    assert block['basis']*reduced == image
    assert reduced**2 == reduced
    assert block['gram']*reduced == reduced.H*block['gram']
    return reduced


@lru_cache(maxsize=1)
def ground_pair():
    """The complete rank-two scalar pair containing the actual local ground."""
    e, minimal, field = field_data()
    block = source(1)
    n = scalar_projector(block)
    v = block['H']-n
    u = 2*n-s.eye(n.rows)
    assert u*v+v*u == s.zeros(n.rows)
    squared = v*v
    x = s.Symbol('x')
    q = s.Poly(squared.charpoly(x).as_expr(), x).sqf_part()
    lam = field.unit*(field.unit-field.one)
    qf = s.Poly(q.as_expr(), x, domain=field)
    linear = s.Poly.from_dict({(1,): field.one, (0,): -lam}, x, domain=field)
    quotient = qf.exquo(linear)
    derivative = field.zero
    for c in q.diff().all_coeffs():
        derivative = derivative*lam+field.convert(c)
    projector = quotient.mul_ground(field.one/derivative)
    assert (projector*projector-projector).rem(qf).is_zero
    assert (linear*projector).rem(qf).is_zero
    powers = [s.eye(n.rows)]
    for _ in range(projector.degree()):
        powers.append(squared*powers[-1])
    coefficients = tuple(projector.rep.nth(i) for i in range(projector.degree()+1))

    def pair_trace(operator):
        return sum((c*field.convert(s.trace(power*operator))
                    for c, power in zip(coefficients, powers)), field.zero)

    delta = field.one-2*field.unit
    assert pair_trace(s.eye(n.rows)) == field.convert(2)
    assert pair_trace(u) == field.zero
    assert pair_trace(v) == field.zero
    assert delta*delta == field.one+4*lam
    root_intervals = minimal.intervals(eps=s.Rational(1, 10**90))
    ground_interval = tuple(map(s.Rational, root_intervals[0][0]))
    assert ground_interval[1] < 0
    return dict(block=block, n=n, v=v, u=u, lam=lam, delta=delta,
                field=field, ground_interval=ground_interval, trace=pair_trace,
                squared_annihilator=q, projector_polynomial=projector)


def scalar_expectation(operator):
    """Exact value and ground-following derivative at scalar cost mu=1."""
    data = ground_pair()
    trace, delta, u, v = data['trace'], data['delta'], data['u'], data['v']
    a = trace(operator)
    b = trace(u*operator)
    c = trace(2*v*operator)
    value = a/2-(b+c)/(2*delta)
    derivative = -b/(2*delta)+(b+c)/(2*delta**3)
    return dict(value=value, derivative=derivative, pair_trace=a,
                scalar_polarization_trace=b, transition_trace=c)


def response_operator(block, target):
    data = ground_pair()
    ns = data['n']
    nt = scalar_projector(target)
    a, ad = block['map'], block['adjoint']
    assert nt*a == a*ns
    hs, ht = data['block']['H'], block['H']
    assert ht == target['H']
    vt = ht-nt
    ut = 2*nt-s.eye(nt.rows)
    assert ut*vt+vt*ut == s.zeros(nt.rows)
    commutator = ht*a-a*hs
    # The entire commutator is independent of the scalar cost.
    assert commutator == vt*a-a*data['v']
    commutator_adjoint = s.ImmutableSparseMatrix(
        data['block']['gram'].inv()*commutator.H*block['gram'])
    moment = s.ImmutableSparseMatrix(commutator_adjoint*commutator)
    assert moment*ns == ns*moment
    effect = s.ImmutableSparseMatrix(ad*a)
    assert effect*ns == ns*effect
    return effect, moment


def encode_statistic(result):
    data = ground_pair()
    result = dict(result)
    result['value_interval'] = outward(evaluate(polynomial(result['value']), data['ground_interval']), 45)
    result['derivative_interval'] = outward(evaluate(polynomial(result['derivative']), data['ground_interval']), 45)
    for key in ('value', 'derivative', 'pair_trace', 'scalar_polarization_trace', 'transition_trace'):
        result[key+'_coefficients_descending_in_ground_energy'] = encoded(result.pop(key))
    return result


@lru_cache(maxsize=None)
def channel(family):
    if family in quark.CHANNELS:
        blocks = [(quark.transfer(family), quark.target(family))]
    elif family in ('ec', 'L'):
        # Spin averaging is equivalent to the physical mixed ground doublet.
        blocks = [(lepton.transfer(family, spin), lepton.target(family, 2*spin))
                  for spin in (0, 1)]
    else:
        raise ValueError('native channel required')
    source_dimension = ground_pair()['block']['H'].rows
    effect, moment = s.zeros(source_dimension), s.zeros(source_dimension)
    for block, target in blocks:
        weight, second = response_operator(block, target)
        effect += weight/s.Integer(len(blocks))
        moment += second/s.Integer(len(blocks))
    success = scalar_expectation(effect)
    result = scalar_expectation(moment)
    field = ground_pair()['field']
    assert success['value'] == field.convert(s.Rational(1, 2))
    assert success['derivative'] == field.zero
    assert result['transition_trace'] == field.zero
    assert result['derivative'] == -2*ground_pair()['lam']*result['scalar_polarization_trace']/ground_pair()['delta']**3
    if family in quark.CHANNELS:
        retained = json.loads((ROOT/'model/native_quark_response/results.json').read_text())
        coefficients = retained['exact_local_response'][family]['exact_ground_moments'][:3]
        factor = field.one
    else:
        retained = json.loads((ROOT/'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json').read_text())
        coefficients = retained['exact_local_channels'][family]['positive_total_raw_moments_coefficients_in_E0'][:3]
        factor = field.convert(s.Rational(1, 2))
    moments = []
    for row in coefficients:
        value = field.zero
        for coefficient in row:
            value = value*field.unit+field.convert(s.Rational(coefficient))
        moments.append(value)
    assert result['value'] == factor*(moments[2]-2*field.unit*moments[1]+field.unit**2*moments[0])
    return dict(family=family, centre_class_of_removed_field=CENTRE_CLASS[family],
                centre_class_of_created_hole=(-CENTRE_CLASS[family]) % 6,
                complete_target_dimensions=[block['H'].rows for block, _ in blocks],
                native_removed_modes=[block['native_removed_mode'] for block, _ in blocks],
                spin_average_weights=[str(s.Rational(1, len(blocks)))]*len(blocks),
                success_probability='1/2', failure_probability='1/2',
                success_probability_derivative='0', outside_complete_target_probability='0',
                complete_unconditioned_second_energy_transfer_moment=encode_statistic(result))


def calculate():
    data = ground_pair()
    from model.coherent_defect_wave.physical_attachment.chiral import representation
    from model.coherent_defect_wave.chiral_joint.interacting_charged_response.model import FAMILIES
    hypercharge = representation()[2][-1]
    for family, expected in CENTRE_CLASS.items():
        assert {int(6*hypercharge[i, i]) % 6 for i in FAMILIES[family]} == {expected}
    rows = {family: channel(family) for family in CENTRE_CLASS}
    assert rows['uc']['complete_unconditioned_second_energy_transfer_moment'] == rows['dc']['complete_unconditioned_second_energy_transfer_moment']
    return dict(schema='oph.defect_scalar_cost_response.v1',
                construction=dict(local_car_modes=30, scalar_states=5,
                    original_scalar_cost='1', all_native_transition_coefficients='1',
                    intervention='H(mu)=V+mu N_H, N_H=diag(0,1,1,1,1)',
                    ground_family='actual physical local ground doublet, equal mixture',
                    source_clock='declared local Hamiltonian clock',
                    ground_selected_energy_polynomial=str(field_data()[1].as_expr()),
                    selected_energy_root='smallest real root',
                    selected_energy_interval=outward(data['ground_interval'], 60),
                    squared_transfer_eigenvalue='lambda=e(e-1)',
                    scalar_pair_dimension=2,
                    squared_transfer_complete_annihilator=str(data['squared_annihilator'].as_expr())),
                statistic=dict(definition='M2_r(mu)=sum_j w_rj(mu) [epsilon_rj(mu)-e(mu)]^2',
                    exact_sum_rule='tr rho_ground(mu) [H(mu),c_r]^dagger [H(mu),c_r]',
                    failure_record_value='0', failure_probability='1/2',
                    success_weight_is_not_divided_out=True,
                    every_target_energy_included_in_polynomial_sum_rule=True,
                    state_derivative_included=True,
                    scalar_cost_derivative_point='mu=1',
                    general_ground_formula='P(mu)=Pi/2-Pi [2V+mu(2N_H-I)]/(2 sqrt(mu^2+4 lambda))',
                    general_moment_formula='M2(mu)=a/2-mu b/(2 sqrt(mu^2+4 lambda))',
                    general_derivative_formula='dM2/dmu=-2 lambda b/(mu^2+4 lambda)^(3/2)'),
                channels=rows,
                scope=dict(exact_local_source_calculation=True,
                    target_particle_data_inputs=[],
                    fitted_Yukawa_parameters=[],
                    physical_mass_predictions=[],
                    scalar_cost_is_neutral_Higgs_field_amplitude=False,
                    centre_class_is_individual_flux_creation_map=False,
                    source_clock_converted_to_physical_time=False,
                    varied_ground_preparation_implemented=False,
                    perturbed_eight_cell_ground_or_error_bounds_claimed=False))


def certificate():
    from model.particle_common_contract.provenance import import_closure
    value = calculate()
    seeds = {'model/defect_scalar_response/response.py',
             'model/defect_scalar_response/PROOF.md',
             'model/native_quark_response/results.json',
             'model/native_quark_response/PROOF.md',
             'model/coherent_defect_wave/chiral_joint/interacting_charged_response/results.json',
             'model/coherent_defect_wave/chiral_joint/local_reference/PROOF.md',
             'model/coherent_defect_wave/chiral_joint/charged_ground_exchange/PROOF.md',
             'model/defect_particle_classes/model.py'}
    seeds.update(path.relative_to(ROOT).as_posix() for path in Path(__file__).parent.iterdir()
                 if path.suffix in ('.py', '.md'))
    paths = import_closure(seeds)
    value['input_sha256'] = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                             for name in sorted(paths)}
    return value
