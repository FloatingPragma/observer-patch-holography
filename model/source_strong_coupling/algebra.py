"""Exact coefficients of the configured native quantum-link source.

The auxiliary matter representation and link Hamiltonian are inherited
architecture. All source labels, Gauss generators and electric sectors are
retained. No QCD coefficient is inserted or inferred from dimension counts.
"""
from functools import lru_cache
import inspect

import sympy as s

from ..common_source import native_matter_attachment as matter
from ..common_source import native_mass_dressing as dressing
from ..native_quantum_gauss import algebra as native
from ..native_gauge_strings import strings


def _zero(matrix):
    return all(s.simplify(value) == 0 for value in matrix)


@lru_cache(maxsize=1)
def color_generators():
    """Hermitian color basis with the fundamental trace convention 1/2."""
    raw = matter.canonical_currents()
    normalized = tuple([s.I*g/2 for g in raw[:7]]
                       + [s.I*(raw[6]+2*raw[7])/(2*s.sqrt(3))])
    gram = s.Matrix(8, 8, lambda a, b: s.simplify(s.trace(normalized[a]*normalized[b])))
    if gram != s.eye(8)/2:
        raise ArithmeticError('fundamental color trace convention failed')
    # Each i*T is admitted through the actual native port control map.
    coefficients = tuple(matter.native_current.admit_generator(s.I*g) for g in normalized)
    selected = tuple(i for i in range(32) if i.bit_count() in (2, 4))
    lifted = tuple(matter.exterior_derivative(-s.I*matter.spin_lift(s.I*g))
                   .extract(selected, selected) for g in normalized)
    casimir = sum((g*g for g in lifted), s.zeros(15)).applyfunc(s.simplify)
    r = matter.representation()
    expected = r['operators']['color_C2'].extract(selected, selected)
    if casimir != expected or any(not _zero(casimir*g-g*casimir) for g in native.native_currents()):
        raise ArithmeticError('same native normalized color Casimir failed')
    matter_gram = s.Matrix(8, 8, lambda a, b: s.simplify(s.trace(lifted[a]*lifted[b])))
    if matter_gram != 2*s.eye(8):
        raise ArithmeticError('full fifteen-state trace index changed')
    source_frame=dressing.current_data()['matter_frame'][:8]
    if not _zero(sum((g.H*g for g in source_frame),s.zeros(15))-2*casimir):
        raise ArithmeticError('actual native trace-frame color contraction changed')
    return normalized, lifted, coefficients, casimir


@lru_cache(maxsize=1)
def source_certificate():
    """Replay every native multiplet and both endpoints in the ambient source."""
    normalized, lifted, coefficients, c3 = color_generators()
    r = matter.representation(); selected = list(r['selected'])
    c2 = r['operators']['weak_C2'].extract(selected, selected)
    q = r['operators']['q'].extract(selected, selected)
    currents = native.native_currents(); blocks = native.native_blocks()
    for action in (native.hamiltonian_action,strings.hamiltonian_action):
        parameters = inspect.signature(action).parameters
        if parameters['hopping'].default != 1 or parameters['link_energy'].default != 1:
            raise ArithmeticError('the committed native action defaults changed')
    block_indices = tuple(indices for _, indices in blocks)
    rows = []; gauss = hamiltonian = 0
    for charge, indices in blocks:
        a, b = native.reference_states(15, indices)
        d = len(indices)
        for j, energy in ((s.Rational(1), s.Rational(1)),
                          (s.Rational(2), s.Rational(1)),
                          (s.Rational(1), s.Rational(2))):
            if native.hamiltonian_action(a, 15, j, energy, blocks=block_indices) != native.scale(b, j):
                raise ArithmeticError('full source action on vacuum-link singlet failed')
            expected = native.sum_states(native.scale(a, d*j), native.scale(b, energy))
            if native.hamiltonian_action(b, 15, j, energy, blocks=block_indices) != expected:
                raise ArithmeticError('full source action on occupied-link singlet failed')
            hamiltonian += 2
        for generator in currents:
            for vertex in (0, 1):
                for state in (a, b):
                    if native.gauss_action(generator, state, vertex):
                        raise ArithmeticError('a complete native singlet violates Gauss')
                    gauss += 1
        color, weak, charge6 = (matrix[indices[0], indices[0]] for matrix in (c3, c2, q))
        if any(matrix[i, i] != value for matrix, value in ((c3, color),(c2,weak),(q,charge6))
               for i in indices):
            raise ArithmeticError('a source block lacks a scalar sector readout')
        rows.append({'central_label': charge, 'dimension': d, 'color_C2': str(color),
                     'weak_C2': str(weak), 'Y_squared': str(charge6**2/36),
                     'configured_electric_energy': '1',
                     'squared_norm_A': str(d), 'squared_norm_B': str(d*d)})
    # These are the actual five nonvacuum sectors, with vacuum energy zero.
    # The same electric law cannot be a three-factor quadratic Casimir law.
    table = s.Matrix([[s.Rational(row[key]) for key in ('color_C2','weak_C2','Y_squared')]
                      for row in rows])
    witness = next(v for v in table.T.nullspace() if sum(v) != 0)
    witness = witness * s.ilcm(*(x.q for x in witness))
    if (witness.T*table != s.zeros(1,3) or sum(witness) == 0
            or table.rank() != 3 or table.row_join(s.ones(5,1)).rank() != 4):
        raise ArithmeticError('native electric-law obstruction failed')
    return {'native_current_count': len(currents), 'matter_dimension': 15,
            'complete_link_dimension': 1+sum(len(i)**2 for _, i in blocks),
            'fundamental_generator_trace': 'Tr(T_a T_b)=delta_ab/2',
            'matter_generator_trace': 'Tr_15(T_a T_b)=2 delta_ab',
            'native_trace_frame_color_casimir': '2 C3',
            'normalization_is_a_physical_coupling_selection': False,
            'port_admission_count': len(coefficients),
            'full_source_Gauss_checks': gauss, 'full_source_action_checks': hamiltonian,
            'blocks': rows, 'quadratic_Casimir_table_rank': table.rank(),
            'electric_augmented_rank': table.row_join(s.ones(5,1)).rank(),
            'electric_incompatibility_left_witness': list(map(str,witness)),
            'witness_dot_actual_electric_spectrum': str(sum(witness)),
            'same_source_reduced_hamiltonian': '[[0,J sqrt(d)],[J sqrt(d),epsilon]]',
            'free_coefficient_controls': ['(J,epsilon)=(1,1)', '(2,1)', '(1,2)'],
            'committed_source_defaults': {'hopping_J':'1','electric_epsilon':'1',
                                          'epsilon_over_abs_J':'1'},
            'fixed_source_defaults_are_ambiguous': False,
            'covariance_alone_selects_the_defaults': False,
            'physical_QCD_color_action_selected': False,
            'alpha_s_MSbar_prediction': None}
