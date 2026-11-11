"""Full90-CAR neutral-composite propagation with exact finite-time bounds."""
from functools import lru_cache
import sympy as s

from model.defect_mass_join.source import Source, car
from model.defect_mass_join.encounter import exact_apply, exact_norm_squared, expression


def packet(source, site, spin, kind):
    poly = car.neutral_composite(spin)
    if kind == 'addition':
        poly = car.adjoint(poly)
        local_mask = base = 0
    elif kind == 'removal':
        local_mask, base = (1 << 30)-1, (1 << source.modes)-1
    else:
        raise ValueError('addition or removal required')
    outer = base & ~(((1 << 30)-1) << (30*site))
    return {(outer | (mask << (30*site)), source.scalar_label(site, scalar)): (int(value), 0, 0, 0)
            for (mask, scalar), value in car.column(poly, local_mask, 0).items()}


def inner(left, right, power):
    return s.simplify(sum(s.conjugate(expression(value))*expression(right.get(key, (0, 0, 0, 0)))
                          for key, value in left.items())/s.Integer(6)**power/2)


@lru_cache(maxsize=1)
def calculate():
    source = Source(((-1, 0, 0), (0, 0, 0), (1, 0, 0)), texture=False)
    r = s.sqrt(2)
    omega = s.Matrix([[4+3*r, -2, 4-3*r], [-2, 10, -2], [4-3*r, -2, 4+3*r]])/6
    assert s.simplify(omega*omega-source.spatial['K']) == s.zeros(3)
    time = s.Rational(1, 10000)
    rows = []
    for kind in ('addition', 'removal'):
        initial = packet(source, 0, 0, kind)
        target = packet(source, 1, 1, kind)
        powers = [initial]
        for _ in range(3):
            powers.append(exact_apply(source, omega, powers[-1]))
        norms = [exact_norm_squared(state, n)/2 for n, state in enumerate(powers)]
        amplitudes = [inner(target, state, n) for n, state in enumerate(powers)]
        assert norms[0] == 1 and amplitudes[0] == amplitudes[1] == 0
        assert s.simplify(amplitudes[2]*s.conjugate(amplitudes[2])) == s.Rational(1, 9)
        leading = time**2/6
        remainder = time**3*s.sqrt(norms[3])/6
        assert s.simplify(leading**2-remainder**2).is_positive
        bounds = [(leading-remainder)**2, (leading+remainder)**2]
        rows.append({'operation': kind,
            'unperturbed_preparation': 'empty90mode fermion occupation and scalar blank' if kind == 'addition'
                                        else 'full90mode fermion occupation and scalar blank',
            'initial_packet': f'normalized J_left,spin0^{"dagger" if kind == "addition" else "identity"} applied to that preparation',
            'coherent_prep': 'X=J_kind tensor |odd compensating reference><blank|+h.c.; anglepi/(2sqrt2) prepares the normalized packet with probability1',
            'retained_reference': 'neutral gauge representation, opposite lepton number and odd parity on the prepared branch; intrinsic Hamiltonian0',
            'initial_state_is_interacting_ground': False,
            'exact_power_supports': [len(state) for state in powers],
            'normalized_exact_power_norm_squared': list(map(str, norms)),
            'normalized_target_matrix_elements_of_H_powers': list(map(str, amplitudes)),
            'time': str(time), 'Taylor_degree': 2,
            'exact_state_remainder_bound': str(remainder),
            'exact_receiver_probability_bounds': list(map(str, bounds)),
            'receiver_probability_bounds_numeric': list(map(float, bounds)),
            'receiver_projector': 'normalized J_middle,spin1 packet projector, with every other source outcome in its complement',
            'complete_no_receiver_probability': '[1-upper_yes,1-lower_yes]',
            'all_quark_occupations_and_scalar_states_retained': True,
            'fermion_number_cutoff': None,
            'same_source_343_ground_or_neutrino_rest_mass_inferred': False})
    return {'source': source.descriptor(), 'complete_addition_and_removal_transport_controls': rows,
        'original_reference_same_for_addition_and_removal': False,
        'actual_ground_both_branches': 'the separate one-cell ground-composite measure retains both J and J^dagger on the same physical ground',
        'lightcone_or_Lorentz_limit_claimed': False}
