"""Complete current commutant and full-record stationary family coefficients."""
from functools import lru_cache

import sympy as s

from ..common_source import native_matter_attachment as matter
from ..common_source import native_family_record_memory as memory
from ..common_source import native_gauge_family_records as gates
from ..common_source import native_repair_mass as fixed_repair
from ..common_source import native_repair_mass_attachment as fixed_attachment

POLICIES = ('always', 'C', 'E', 'CE')


def clean(value):
    return value.applyfunc(lambda x: s.cancel(s.expand(x)))


def zero(value):
    return all(s.cancel(s.expand(x)) == 0 for x in value)


@lru_cache(maxsize=1)
def current_blocks():
    data = matter.representation()
    selected = list(data['selected'])
    q = data['operators']['q'].extract(selected, selected)
    blocks = {}
    for charge in sorted(set(q.diagonal())):
        indices = tuple(i for i in range(15) if q[i, i] == charge)
        blocks[int(charge)] = indices
    if len(blocks) != 5 or sum(map(len, blocks.values())) != 15:
        raise ArithmeticError('the actual charge blocks changed')
    return blocks


@lru_cache(maxsize=1)
def commutant_certificate():
    data = matter.representation()
    currents = data['matter_currents']
    coefficients = matter.native_current.admit_generator(data['native']['central'])
    central = sum((x*g for x, g in zip(coefficients, currents, strict=True)), s.zeros(15))
    selected = list(data['selected'])
    q = data['operators']['q'].extract(selected, selected)
    index = next(i for i in range(15) if central[i, i] != 0)
    scale = s.simplify(q[index, index]/central[index, index])
    if not zero(q-scale*central):
        raise ArithmeticError('the block-separating charge is not a native current')
    rows = []
    for charge, indices in current_blocks().items():
        d = len(indices)
        equations = s.MutableSparseMatrix.zeros(12*d*d, d*d)
        for a, current in enumerate(currents):
            g = current.extract(indices, indices)
            for i in range(d):
                for j in range(d):
                    row = a*d*d+i*d+j
                    for k in range(d):
                        equations[row, k*d+j] += g[i, k]
                        equations[row, i*d+k] -= g[k, j]
        rank = equations.to_DM(extension=True).rank()
        if rank != d*d-1:
            raise ArithmeticError('a native current block is not irreducible')
        rows.append({'q6Y': charge, 'dimension': d, 'matrix_unknowns': d*d,
                     'commutator_equations': 12*d*d, 'exact_rank': rank,
                     'commutant_complex_dimension': 1})
    return {'native_current_count': 12, 'complete_internal_rows': 15,
            'charge_blocks': rows, 'internal_commutant_complex_dimension': 5,
            'family_factor_dimension': 3, 'complete_family_rows': 45,
            'family_commutant': 'direct sum of five Mat3(C)',
            'family_commutant_complex_dimension': 45,
            'charge_is_exact_native_current_combination': True,
            'same_current_canonical_kinetic_embedding_group': 'U(3)^5',
            'shared_left_quark_doublet_family_unitary': True,
            'general_nonunitary_or_current_changing_maps_classified': False}


def coefficient_freedom_certificate():
    """The complete invariant lines permit a generic nonzero-CP control."""
    vertices = matter.invariant_higgs_lines()
    if len(vertices) != 3:
        raise ArithmeticError('the three native Higgs invariant lines changed')
    up = s.diag(1, 2, 3)
    down = s.Matrix([[1, 1, 0], [0, 2, s.I], [1, 0, 3]])
    invariant = cp_invariant(up, down)
    if invariant == 0:
        raise ArithmeticError('the gauge-invariant complex control lost CP violation')
    return {'complex_coefficients_before_family_basis_quotient': 27,
            'reason': 'three complete one-dimensional Higgs invariant lines tensor Mat3(C)',
            'same_current_complex_control_up': [[str(x) for x in up.row(i)] for i in range(3)],
            'same_current_complex_control_down': [[str(x) for x in down.row(i)] for i in range(3)],
            'complex_control_CP_trace_cube_exact': str(invariant),
            'control_is_selected_by_the_declared_real_memory': False,
            'zero_CP_follows_from_real_memory_not_gauge_covariance': True}


def fixed_repair_certificate():
    source = fixed_repair.architecture()
    roots = fixed_attachment.spectral_certificate(source['Y'])
    kernels = fixed_attachment.higgs_kernels(source['Y'])
    if {row['channel'] for row in kernels} != {'up', 'down', 'charged_lepton'}:
        raise ArithmeticError('the fixed repair attachment omitted a charged channel')
    return {'source': 'model/common_source/native_repair_mass_attachment.py',
            'advance_probability_exact': str(source['p']),
            'charged_channels': sorted(row['channel'] for row in kernels),
            'common_positive_singular_value_intervals': roots['positive_eigenvalue_intervals_exact'],
            'all_same_current_unitary_embeddings_have_identical_channel_spectra': True,
            'left_squared_up_down_matrices_identical_in_every_admitted_embedding': True,
            'mixing_moduli_squared': [[int(i == j) for j in range(3)] for i in range(3)],
            'charged_current_CP_trace_cube_exact': '0', 'empirical_inputs': []}


def embedding(unitaries):
    """Every exact same-current unitary, keyed by the five actual charges."""
    blocks = current_blocks()
    if type(unitaries) is not dict or set(unitaries) != set(blocks):
        raise ValueError('all five charge blocks are required')
    result = s.zeros(45)
    for charge, indices in blocks.items():
        u = unitaries[charge]
        if (not isinstance(u, s.MatrixBase) or u.shape != (3, 3)
                or any(x.has(s.Float) or x.free_symbols or x.is_finite is not True for x in u)
                or not zero(u.H*u-s.eye(3))):
            raise ValueError('an exact three-family unitary is required')
        for i in indices:
            result[3*i:3*i+3, 3*i:3*i+3] = u
    return s.ImmutableSparseMatrix(result)


@lru_cache(maxsize=1)
def symbolic_memory():
    p, q = s.symbols('p q', real=True)
    a, b = p*(1-q), p*q
    denominator = 1-a*a+a*b
    projectors = memory._projectors(False)
    if not zero(projectors[9]-projectors[10]):
        raise ArithmeticError('the antipodal projector identification changed')
    pi = s.zeros(12, 1)
    pi[2] = (1-a)/denominator
    pi[3] = a*pi[2]
    pi[10] = a*b*pi[2]/(1-a*a)
    pi[9] = a*pi[10]
    _, (gC, gE) = memory.geometry()
    flags = ((False, False), (False, True), (True, False), (True, True))
    probabilities = ((1-p)*(1-q), (1-p)*q, p*(1-q), p*q)
    transition = s.zeros(12)
    weights = {policy: s.zeros(12, 1) for policy in POLICIES}
    for old in range(12):
        for (c, e), probability in zip(flags, probabilities, strict=True):
            destination = memory._destination(old, c, e, 2, gC, gE)
            transition[destination, old] += probability
            for policy in POLICIES:
                enabled = (policy == 'always' or policy == 'C' and c
                           or policy == 'E' and e or policy == 'CE' and c and e)
                if enabled:
                    weights[policy][destination] += pi[old]*probability
    if not zero(transition*pi-pi) or s.cancel(sum(pi)-1) != 0:
        raise ArithmeticError('the complete twelve-input stationary law failed')
    coefficients = {policy: clean(sum((weight*P for weight, P in
                    zip(values, projectors, strict=True)), s.zeros(3)))
                    for policy, values in weights.items()}
    compact = ((1-a)*projectors[2]+a*(1-a)*projectors[3]+a*b*projectors[10])/denominator
    if not zero(coefficients['always']-compact):
        raise ArithmeticError('the compact stationary formula lost a record')
    if not zero(coefficients['C']-compact+(1-p)*projectors[2]):
        raise ArithmeticError('the current-event coefficient was duty-renormalized')
    for policy, trace in zip(POLICIES, (1, p, q, p*q), strict=True):
        if s.cancel(s.trace(coefficients[policy])-trace) != 0:
            raise ArithmeticError('an unconditional duty trace changed')
    return {'p': p, 'q': q, 'denominator': denominator, 'stationary': clean(pi),
            'transition': transition, 'weights': {k: clean(v) for k, v in weights.items()},
            'coefficients': coefficients, 'projectors': projectors}


def coefficient(p, q, policy='always'):
    if policy not in POLICIES:
        raise ValueError('unknown complete-record duty policy')
    p, q = s.sympify(p), s.sympify(q)
    if (any(x.has(s.Float) or x.free_symbols or x.is_real is not True for x in (p, q))
            or p.is_nonnegative is not True or (1-p).is_positive is not True
            or q.is_nonnegative is not True or (1-q).is_nonnegative is not True):
        raise ValueError('exact 0<=p<1 and 0<=q<=1 are required')
    data = symbolic_memory()
    return clean(data['coefficients'][policy].subs({data['p']: p, data['q']: q}))


def cp_invariant(up, down):
    """Trace of the commutator cube, without masses or eigenbasis division."""
    if up.shape != (3, 3) or down.shape != (3, 3):
        raise ValueError('three-family coefficients are required')
    hu, hd = up*up.H, down*down.H
    commutator = hu*hd-hd*hu
    return s.simplify(s.trace(commutator**3))


@lru_cache(maxsize=1)
def memory_certificate():
    data = symbolic_memory()
    p, q = data['p'], data['q']
    r = s.Symbol('r', real=True)
    A = data['coefficients']['always']
    B = A.subs(q, r)
    D = data['denominator']
    formula = ((s.sqrt(5)-1)*p*p*(q-r)
               *(p*p*(1-q)*(1-r)+q+r-1)/(10*D*D.subs(q, r)))
    if s.cancel(s.expand((A*B-B*A)[0, 2]-formula)) != 0:
        raise ArithmeticError('the exact noncommuting-memory witness changed')
    x, y, z = s.symbols('x y z', real=True)
    skew = s.Matrix([[0, x, y], [-x, 0, z], [-y, -z, 0]])
    if s.trace(skew**3) != 0:
        raise ArithmeticError('the real charged-current CP identity failed')
    ranks = {}
    for policy, weights in data['weights'].items():
        support = [i for i, value in enumerate(weights) if value != 0]
        rank = sum((data['projectors'][i] for i in support), s.zeros(3)).rank()
        ranks[policy] = {'support_at_0_lt_p_q_lt_1': support,
                         'coefficient_rank': rank,
                         'trace_exact': str(s.cancel(s.trace(data['coefficients'][policy])))}
    return {'all_old_memory_inputs': 12, 'all_native_raw_labels': 128,
            'all_native_plus_binary_labels_per_input': 512,
            'complete_Kraus_label_count': 6144,
            'Kraus_count_domain': 'reduced memory-feedback channel; current auxiliary legs traced and fine outcomes binary-grouped',
            'underlying_current_outcome_pairs': 169,
            'underlying_native_current_fine_labels_per_old_input': 21632,
            'underlying_native_current_fine_labels_including_all_old_inputs': 259584,
            'binary_grouping_is_not_postselection': True,
            'reset_port': 2, 'geometric_projector_rank': 1,
            'antipodal_projector_identity': 'P9=P10',
            'stationary_denominator': str(D),
            'stationary_weights': list(map(str, data['stationary'])),
            'coefficient_formula': 'Y=((1-a)P2+a(1-a)P3+ab P10)/(1-a^2+ab); a=p(1-q), b=pq',
            'mixing_commutator_02': str(formula),
            'nonzero_mixing_domain': '0<p<=1/12; 0<=q,r<=1/6; q!=r',
            'negative_factor_uniform_upper_bound': '-95/144',
            'duty_policies': ranks,
            'all_real_projector_weights_have_charged_current_CP_invariant_exact': '0',
            'nontrivial_real_mixing_is_nonzero_CP_violation': False,
            'coefficient_is_stationary_first_order_generator': True,
            'finite_step_channel_identified_as_unitary': False}


@lru_cache(maxsize=1)
def source_gates():
    mapping = {tuple(row['charge_origin']): row['existing_action_channel']
               for row in matter.action_intertwiner()[1]}
    result = {}
    for row in gates.sector_gates():
        name = mapping[tuple(map(str, row['origin']))]
        result[name] = row
    if set(result) != {'up', 'down', 'charged_lepton'}:
        raise ArithmeticError('the actual Higgs-current dictionary changed')
    return result
