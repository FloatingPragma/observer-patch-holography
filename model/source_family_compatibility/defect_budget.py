"""Exact full-source defect budgets for complete globally sampled repairs."""
from fractions import Fraction as F
from functools import lru_cache

import sympy as s

from model.defect_field_bridge.full_domain import finite_branches, RAW_SLOTS
from model.strong_coupling_refinement.algebra import operators, validate_graph


def rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational required')
    value = F(value)
    if value < 0 or (positive and not value):
        raise ValueError('positive rational required' if positive else 'negative input')
    return value


def defect_counts(sites):
    if type(sites) is not int or not 2 <= sites <= 4:
        raise ValueError('explicit whole-Fock oracle uses two through four sites')
    return tuple(sum(((mask >> (2*x)) & 3).bit_count() != 1
                     for x in range(sites)) for mask in range(1 << (2*sites)))


@lru_cache(maxsize=None)
def exact_source(sites, edges, coupling, contact):
    """Reassemble every actual number block; do not compress the repair domain."""
    validate_graph(sites, edges)
    counts = defect_counts(sites)
    c, u = rational(coupling, positive=True), rational(contact, positive=True)
    dim = len(counts)
    hopping = s.zeros(dim)
    for number in range(2*sites+1):
        masks, twice, doubles = operators(sites, edges, number)
        coo = twice.tocoo()
        for row, col, value in zip(coo.row, coo.col, coo.data):
            hopping[masks[row], masks[col]] = s.Rational(c.numerator, 2*c.denominator)*int(value)
        for mask, double in zip(masks, doubles):
            if counts[mask] != 2*int(double)-number+sites:
                raise ArithmeticError('whole-source defect/occupation identity failed')
    defect = s.diag(*counts)
    h = hopping+s.Rational(u.numerator, 2*u.denominator)*defect
    if hopping != hopping.H or h != h.H:
        raise ArithmeticError('actual signed source is not Hermitian')
    return s.ImmutableMatrix(h), s.ImmutableMatrix(hopping), s.ImmutableMatrix(defect)


def seam_gram(sites, seam, theta):
    """Compute sum over all actual delta-Kraus columns using rational weights."""
    theta = rational(theta)
    if theta > 1:
        raise ValueError('law weight must be at most one')
    dim = len(defect_counts(sites)); gram = s.zeros(dim)
    for law, weight in (('C', 1-theta), ('E', theta)):
        rows = finite_branches(sites=sites, flavors=2, law=law, seam=seam)
        if len(rows) != 2*RAW_SLOTS:
            raise ArithmeticError('padded raw records omitted')
        for (_, raw), entries in rows.items():
            delta = s.zeros(dim)
            for initial, final, sign in entries:
                delta[final, initial] = sign
            if raw == 0:
                delta -= s.eye(dim)
            gram += s.Rational(weight.numerator, 2*weight.denominator)*delta.H*delta
    expected = []
    q = []
    for mask in range(dim):
        n = [((mask >> (2*x)) & 3).bit_count() for x in seam]
        expected.append(2 if 2 in n else s.Rational(1+theta) if sum(n) == 1 else 0)
        q.append(int(n != [1, 1]))
    if gram != s.diag(*expected):
        raise ArithmeticError('whole-source repair Gram formula differs')
    if any(s.Rational(2*q[j])-gram[j, j] < 0 for j in range(dim)):
        raise ArithmeticError('defective endpoint bound failed')
    return s.ImmutableMatrix(gram), tuple(q)


def energy_budget(vertices, edges, degree, coupling, contact, energy):
    """Rational expectation bound, valid for an H or H-E0 budget C>=0."""
    if any(type(x) is not int or x < 1 for x in (vertices, edges, degree)):
        raise ValueError('positive integer graph counts required')
    if vertices < 2 or degree >= vertices or edges > vertices*degree//2:
        raise ValueError('graph counts incompatible with a simple maximum-degree bound')
    c, u, budget = rational(coupling, positive=True), rational(contact, positive=True), rational(energy)
    crude = 2*(budget+c*edges)/u
    a = u/2-c*degree
    local = 2*budget/a+4*c*c*edges*degree/(a*a) if a > 0 else None
    bound = min(F(vertices), crude, *([local] if local is not None else []))
    return {'defect_expectation_upper': bound, 'global_hopping_bound': crude,
            'local_form_bound': local, 'local_coefficient_a': a,
            'ground_relative_energy_budget_allowed': True,
            'ground_sector_or_projection_assumed': False}


def uniform_history(vertices, edges, degree, coupling, contact, energy, attempts):
    if type(attempts) is not int or attempts < 0:
        raise ValueError('nonnegative integer attempt count required')
    budget = energy_budget(vertices, edges, degree, coupling, contact, energy)
    b = budget['defect_expectation_upper']
    # The trace bound is 2m sqrt(2 d b/E), capped at two.
    squared = min(F(4), 8*attempts*attempts*degree*b/edges)
    return {**budget, 'attempts': attempts, 'complete_trace_bound_squared': squared,
            'bound_controls_entire_source_reference_archive': True,
            'ideal_energy_preservation_or_work_budget_required': True}


def finite_controls():
    rows = []
    for sites, edges in ((2, ((0, 1, 1),)), (3, ((0, 1, 1), (1, 2, -1)))):
        h, hopping, defect = exact_source(sites, edges, F(3, 2), F(4))
        aggregate = s.zeros(h.rows)
        for x, y, sign in edges:
            local_h, local_t, _ = exact_source(sites, ((x, y, sign),), F(3, 2), F(4))
            gram, q = seam_gram(sites, (x, y), F(1, 2))
            p = s.eye(h.rows)-s.diag(*q)
            if p*local_t*p != s.zeros(h.rows):
                raise ArithmeticError('hopping compression on nondefective endpoints is nonzero')
            aggregate += gram/len(edges)
        degree = max(sum(x in e[:2] for e in edges) for x in range(sites))
        if any(2*degree*defect[j, j]/len(edges)-aggregate[j, j] < 0 for j in range(h.rows)):
            raise ArithmeticError('global sparse-defect bound failed')
        rows.append({'sites': sites, 'whole_Fock_dimension': h.rows,
                     'complete_record_slots_per_attempt': len(edges)*2*2*RAW_SLOTS,
                     'all_matrix_units_and_reference_coefficients_covered': True,
                     'signed_local_hopping_compression_zero': True})
    h, _, defect = exact_source(2, ((0, 1, 1),), F(3, 2), F(4))
    masks = [m for m in range(16) if m.bit_count() == 2]
    block = h.extract(masks, masks); identity = s.eye(6)
    eigenvalues = (s.Rational(-1, 2), s.Integer(0), s.Integer(4), s.Rational(9, 2))
    phases = (s.I, 1, 1, -s.I)
    evolution = s.zeros(6)
    for eigenvalue, phase in zip(eigenvalues, phases):
        projection = identity
        for other in eigenvalues:
            if other != eigenvalue:
                projection = projection*(block-other*identity)/(eigenvalue-other)
        evolution += phase*projection
    if ((evolution.H*evolution-identity).applyfunc(s.expand) != s.zeros(6)
            or (evolution*block-block*evolution).applyfunc(s.expand) != s.zeros(6)):
        raise ArithmeticError('actual Hubbard time-pi control failed')
    initial = s.zeros(6, 1); initial[masks.index(9)] = 1
    evolved = evolution*initial
    state = s.zeros(16, 1)
    for j, mask in enumerate(masks): state[mask] = evolved[j]
    b = s.expand((state.H*defect*state)[0])
    if b != s.Rational(9, 25) or s.expand((state.H*h*state)[0]) != 0:
        raise ArithmeticError('interacting leakage or conserved energy changed')
    grams = [s.expand((state.H*seam_gram(2, (0, 1), theta)[0]*state)[0])
             for theta in (F(0), F(1, 2), F(1))]
    if grams != [s.Rational(9, 25)]*3:
        raise ArithmeticError('complete repair did not retain collision disturbance')
    return {'whole_source_controls': rows,
            'actual_interacting_time_control': {
                'coupling': '3/2', 'contact': '4', 'duration': 'pi',
                'initial_mask': 9, 'whole_Fock_dimension': 16,
                'initial_and_final_H_expectation': '0',
                'doublon_hole_probability': '9/50', 'defect_expectation': str(b),
                'complete_dilation_difference_squared': [str(g) for g in grams],
                'single_occupation_invariance_assumed': False}}


def certificate():
    def encode(value):
        if isinstance(value, F): return str(value)
        if isinstance(value, dict): return {k: encode(v) for k, v in value.items()}
        if isinstance(value, list): return [encode(v) for v in value]
        return value
    return encode({'finite_source': finite_controls(),
                   'uniform_energy_history_examples': [
                       uniform_history(8, 12, 3, 1, u, 0, 4) for u in (64, 256, 1024)],
                   'full_Hubbard_propagation_in_same_repair_domain': True,
                   'capacity_action_or_physical_clock_selected': False,
                   'interacting_relativistic_particle_limit_proved': False})
