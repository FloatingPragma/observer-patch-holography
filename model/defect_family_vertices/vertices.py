"""All native contact components and complete charged-scalar pair pointers."""
from functools import lru_cache
import sympy as s

from model.coherent_defect_wave.chiral_joint.higgs_composites import model as car
from model.coherent_defect_wave.chiral_higgs.model import higgs_memory
from model.coherent_defect_wave.physical_attachment.chiral import representation
from .source import (encoded_value, expectation, ground_pair, removal_moments,
                     scalar_population, scalar_readout)

CHANNELS = {'up': ('uc', 3), 'down': ('dc', 3), 'charged_lepton': ('ec', 1)}


@lru_cache(maxsize=1)
def action():
    symbolic, couplings = car.quantum_yukawa()
    substitutions = {symbol: 1 for symbol in couplings.values()}
    parts = {name: car.clean({key: value.subs(substitutions) for key, value in part.items()})
             for name, part in symbolic.items()}
    interaction = car.add(*parts.values())
    h = car.add(interaction, car.monomial(matrix=higgs_memory()[2]))
    assert h == car.adjoint(h)
    for current in range(12):
        assert not car.commutator(car.gauge_current(current), h)
    return h, parts, interaction


def scalar_matrix(operator, dimension=5):
    assert all(not word for word, _, _ in operator)
    return s.ImmutableSparseMatrix(dimension, dimension, {(row, col): value
                                         for (_, row, col), value in operator.items()})


@lru_cache(maxsize=1)
def scalar_double_commutator():
    h, _, interaction = action()
    fields = [car.monomial(matrix=p) for p in higgs_memory()[0]]
    double = car.add(*(car.commutator(phi, car.commutator(h, car.adjoint(phi)))
                       for phi in fields))
    symmetric = car.scale(car.add(double, car.adjoint(double)), s.Rational(1, 2))
    expected = car.add(car.scale(car.monomial(), 4),
        car.scale(car.monomial(matrix=higgs_memory()[2]), -5),
        car.scale(interaction, -s.Rational(3, 2)))
    assert symmetric == expected
    commutator_sum = car.add(*(car.commutator(phi, car.adjoint(phi)) for phi in fields))
    assert car.add(double, car.scale(car.adjoint(double), -1)) == car.commutator(h, commutator_sum)
    return {'full_CAR_scalar_symmetric_identity': 'sum_i [phi_i,[H,phi_i^dagger]] + h.c. = 8I-10N_H-3V',
            'antihermitian_part_is_ground_stationarity_commutator': True,
            'complete_scalar_normalization': encoded_value(scalar_population()['scalar_normalization']),
            'normalization_copied_from_lepton_only_source': False}


@lru_cache(maxsize=1)
def contacts():
    h, _, _ = action()
    phi = higgs_memory()[0]
    scalar_basis = (*phi, *(p.H for p in phi))
    population = scalar_population()
    probabilities = population['populations']
    covariance = population['scalar_covariance']
    field = ground_pair()['field']
    annihilators = [car.monomial(((i, False),)) for i in range(30)]
    rows, matrices = [], {}
    for i, ci in enumerate(annihilators):
        assert car.anticommutator(ci, car.adjoint(ci)) == car.monomial()
        derivative = car.commutator(ci, h)
        for j, cj in enumerate(annihilators):
            insertion = scalar_matrix(car.anticommutator(derivative, cj))
            coefficients = (insertion[0, 1], insertion[0, 2], insertion[1, 0], insertion[2, 0])
            assert insertion == sum((a*b for a, b in zip(coefficients, scalar_basis)), s.zeros(5))
            mixed = []
            for basis in scalar_basis:
                product = insertion*basis.H
                mixed.append(sum((probabilities[level]*field.convert(product[level, level])
                                  for level in range(5)), field.zero))
            assert all(value/covariance == field.convert(coefficient)
                       for value, coefficient in zip(mixed, coefficients))
            matrices[i, j] = insertion
            if any(coefficients):
                channel = next(name for name, indices in car.SINGLETS.items()
                               if i % 15 in indices or j % 15 in indices)
                rows.append({'left_mode': i, 'right_mode': j, 'channel': channel,
                             'scalar_basis_coefficients': list(map(str, coefficients)),
                             'mixed_three_point': [encoded_value(value) for value in mixed]})
    assert all(matrices[i, j] == -matrices[j, i] for i in range(30) for j in range(30))
    return {'all_fermion_pairs_checked': 900, 'full_CAR_modes': 30,
            'scalar_basis': ['phi_0', 'phi_1', 'phi_0^dagger', 'phi_1^dagger'],
            'fermion_complete_addition_removal_normalization': 'I30',
            'scalar_Gram': 'C I4', 'C': encoded_value(covariance),
            'scalar_spectral_normalization': encoded_value(population['scalar_normalization']),
            'canonical_contact_modulus': 'sqrt(Z_H) for each nonzero unit primitive component',
            'nonzero_components': rows, 'zero_components': 900-len(rows),
            'color_components_per_channel': {'up': 3, 'down': 3, 'charged_lepton': 1},
            'full_Fock_operator_reconstruction': True,
            'flavor_dependent_primitive_parameters_supplied': False,
            'generation_copies_inferred_from_color': False}


def lift_reference(operator, reference):
    return {(word, 3*row+a, 3*col+b): value*entry
            for (word, row, col), value in operator.items()
            for (a, b), entry in s.ImmutableSparseMatrix(reference).todok().items()}


def vacuum_reference_block(operator):
    return car.clean({(word, row//3, col//3): value
                      for (word, row, col), value in operator.items()
                      if row % 3 == col % 3 == 0})


@lru_cache(maxsize=None)
def channel_scalar_energies(family, colors):
    """Split each literal channel energy using its complete native probe map."""
    from model.native_quark_response import local as quark
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange import model as lepton
    from model.coherent_defect_wave.chiral_joint.local_reference import physical
    source = ground_pair()['block']
    raw_metric = physical.metric_on(source['keys'])
    selectors = []
    for levels in ((1, 2), (3, 4)):
        raw = s.diag(*(int(key[2] in levels) for key in source['keys']))
        compressed = source['gram'].inv()*source['basis'].H*raw_metric*raw*source['basis']
        assert source['basis']*compressed == raw*source['basis']
        selectors.append(compressed)
    blocks = ([quark.transfer(family)] if family in quark.CHANNELS else
              [lepton.transfer(family, spin) for spin in (0, 1)])
    symmetric = s.zeros(source['H'].rows)
    for block in blocks:
        k = block['adjoint']*block['H']*block['map']-block['adjoint']*block['map']*source['H']
        adjoint = source['gram'].inv()*k.H*source['gram']
        symmetric += (k+adjoint)/len(blocks)
    total = -2*colors*expectation(symmetric)
    result = tuple(-2*colors*expectation(p*symmetric+symmetric*p) for p in selectors)
    assert sum(result) == total == -4*colors*removal_moments()[family]
    return result


@lru_cache(maxsize=1)
def pulse_responses():
    h, parts, interaction = action()
    phi, scalar_currents, _, _ = higgs_memory()
    raising = [s.ImmutableSparseMatrix(3, 3, {(i+1, 0): 1}) for i in range(2)]
    x = car.monomial(matrix=sum((s.kronecker_product(p.H, r)+s.kronecker_product(p, r.H)
                                for p, r in zip(phi, raising)), s.zeros(15)))
    # Exact blank-reference blocks of exp(-i theta X). The independent
    # trigonometric identities enter through c^2+s^2=C^2+S^2=1.
    c, sn, C, S = s.symbols('c sn C S', real=True)
    blocks = [s.diag(C, 1, 1, c, c)]
    for component in range(2):
        block = s.zeros(5)
        block[component+1, 0] = -s.I*S/s.sqrt(2)
        block[0, component+3] = -s.I*sn
        blocks.append(block)
    gram = sum((block.H*block for block in blocks), s.zeros(5))
    assert gram == s.diag(C*C+S*S, 1, 1, c*c+sn*sn, c*c+sn*sn)
    blank_evolution = s.Matrix(15, 5, lambda row, col: blocks[row % 3][row//3, col])
    derivative = (blank_evolution.diff(c)*(-sn)+blank_evolution.diff(sn)*c
                  +blank_evolution.diff(C)*(-s.sqrt(2)*S)
                  +blank_evolution.diff(S)*(s.sqrt(2)*C))
    assert derivative == -s.I*scalar_matrix(x, 15)*blank_evolution
    initial = s.Matrix(15, 5, lambda row, col: int(row == 3*col))
    assert blank_evolution.subs({c: 1, C: 1, sn: 0, S: 0}) == initial
    moments = removal_moments()
    field = ground_pair()['field']
    rows, total_energy = {}, field.zero
    for channel, (family, colors) in CHANNELS.items():
        part = parts[channel]
        pairs = []
        for component in range(2):
            coefficients = {word: value for (word, row, col), value in part.items()
                            if (row, col) == (component+1, 0)}
            assert sum(abs(v) for v in coefficients.values()) == 2*colors
            pairs.append(car.add(*(car.scale(car.monomial(word), value)
                                   for word, value in coefficients.items())))
        positive = car.add(*(car.multiply(pair, car.monomial(matrix=p.H))
                             for pair, p in zip(pairs, phi)))
        assert car.add(positive, car.adjoint(positive)) == part
        b = car.add(*(lift_reference(pair, r.H) for pair, r in zip(pairs, raising)))
        divisor = 8*colors
        pointer = car.scale(car.add(b, car.scale(car.adjoint(b), -1)), s.I/divisor)
        pulled_b = car.add(*(car.multiply(car.monomial(matrix=blocks[0].H),
            car.multiply(pair, car.monomial(matrix=blocks[k+1]))) for k, pair in enumerate(pairs)))
        pulled_pointer = car.scale(car.add(pulled_b, car.scale(car.adjoint(pulled_b), -1)), s.I/divisor)
        particle_selector = car.monomial(matrix=s.diag(0, 1, 1, 0, 0))
        antiparticle_selector = car.monomial(matrix=s.diag(0, 0, 0, 1, 1))
        particle_part = car.anticommutator(particle_selector, part)
        antiparticle_part = car.anticommutator(antiparticle_selector, part)
        assert car.add(particle_part, antiparticle_part) == part
        assert pulled_pointer == car.add(car.scale(particle_part, S/(s.sqrt(2)*divisor)),
                                          car.scale(antiparticle_part, C*sn/divisor))
        derivative = car.scale(car.commutator(x, pointer), s.I)
        assert vacuum_reference_block(derivative) == car.scale(part, s.Rational(1, divisor))
        # Full CAR sum rule derives channel energy from complete removal moments.
        number_moment = {}
        for index in car.SINGLETS[channel]:
            for spin in (0, 1):
                c = car.monomial(((index+15*spin, False),))
                number_moment = car.add(number_moment,
                    car.multiply(car.adjoint(c), car.commutator(h, c)))
        assert part == car.scale(car.add(number_moment, car.adjoint(number_moment)), -1)
        channel_energy = -4*colors*moments[family]
        particle_energy, antiparticle_energy = channel_scalar_energies(family, colors)
        total_energy += channel_energy
        # For theta=pi/(2 sqrt(2)), the full blank block gives
        # <A>=<V_particle>/(8 d sqrt(2)) exactly, with both pointer outcomes.
        from .source import evaluate, polynomial, sqrt_interval, divide, outward
        finite_coefficient = particle_energy/divisor
        finite_interval = divide(evaluate(polynomial(finite_coefficient), ground_pair()['ground_interval']),
                                 sqrt_interval((s.Integer(2), s.Integer(2))))
        for current in range(12):
            t = scalar_currents[current][1:3, 1:3]
            qr = s.diag(s.zeros(1), -t.T)
            total_current = car.add(lift_reference(car.gauge_current(current), s.eye(3)),
                                   car.monomial(matrix=s.kronecker_product(s.eye(5), qr)))
            assert not car.commutator(total_current, x)
            assert not car.commutator(total_current, pointer)
        rows[channel] = {'all_colors_included': colors,
                         'pair_norm_upper_each_component': 2*colors,
                         'pointer_divisor': divisor, 'pointer_norm_upper': '1',
                         'pointer_effects': '(I+A)/2 and (I-A)/2',
                         'complete_singlet_removal_first_moment': encoded_value(moments[family]),
                         'full_channel_interaction_energy': encoded_value(channel_energy),
                         'scalar_particle_channel_energy': encoded_value(particle_energy),
                         'scalar_antiparticle_channel_energy': encoded_value(antiparticle_energy),
                         'pair_pointer_pulse_susceptibility': encoded_value(channel_energy/divisor),
                         'color_average_normalization': True,
                         'pointer_probabilities_at_zero_pulse': ['1/2', '1/2'],
                         'probability_derivatives_at_zero_pulse': [encoded_value(channel_energy/(2*divisor)),
                                                                    encoded_value(-channel_energy/(2*divisor))],
                         'exact_finite_pulse': {
                             'angle': 'pi/(2 sqrt(2))',
                             'pointer_mean': 'particle_channel_energy/(8 d sqrt(2))',
                             'coefficient_of_inverse_sqrt_two': encoded_value(finite_coefficient),
                             'pointer_mean_interval': outward(finite_interval, 60),
                             'all_pointer_probability_intervals': [
                                 outward(((1+finite_interval[0])/2, (1+finite_interval[1])/2), 60),
                                 outward(((1-finite_interval[1])/2, (1-finite_interval[0])/2), 60)]}}
    assert total_energy == scalar_population()['interaction_energy']
    assert rows['up']['pair_pointer_pulse_susceptibility'] == rows['down']['pair_pointer_pulse_susceptibility']
    return {'intervention': 'X=sum_k phi_k^dagger tensor |k dual><0| + h.c.',
            'reference_dimension': 3, 'reference_is_retained': True,
            'finite_pulse_blank_blocks': {
                'K_0': 'diag(cos(sqrt(2) theta),1,1,cos(theta),cos(theta))',
                'K_k': '-i sin(sqrt(2) theta)/sqrt(2) |p_k><0| - i sin(theta) |0><a_k|',
                'literal_Kraus_completeness': True,
                'pointer_mean_all_angles': '[sin(sqrt(2) theta)/sqrt(2) V_particle + cos(sqrt(2) theta) sin(theta) V_antiparticle]/(8 d)'},
            'global_gauge_covariance_checks': 72,
            'scalar_cost_varied': False, 'primitive_Yukawa_coefficients_varied': False,
            'source_state': 'actual full-thirty-mode physical ground doublet / 2',
            'all_pointer_outcomes_retained': True, 'channels': rows}
