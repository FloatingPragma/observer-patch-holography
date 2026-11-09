"""Complete native records identify the Hubbard hopping/contact ratio."""
from fractions import Fraction as F
from math import ceil

import sympy as s

from model.defect_field_bridge.full_domain import finite_branches, RAW_SLOTS
from model.full_fock_field.field import FieldGeometry
from model.field_symmetry import symmetry
from model.native_program import interpreter as native
from model.source_family_compatibility.defect_budget import exact_source
from model.source_strong_coupling.readout import Interval, probability, infer_coefficients


def _rational(value, *, positive=False):
    if type(value) not in (int, F):
        raise ValueError('exact rational input required')
    value = F(value.numerator, value.denominator)
    if value < 0 or (positive and value == 0):
        raise ValueError('positive input required' if positive else 'negative input')
    return value


def raw_matrices(sites, seam, law):
    """Unweighted integer branches; each coin has probability one half."""
    if type(sites) is not int or not 2 <= sites <= 3:
        raise ValueError('bounded exact oracle uses two or three sites')
    rows = finite_branches(sites=sites, flavors=2, seam=seam, law=law)
    result = {}
    for key, columns in rows.items():
        matrix = s.zeros(4**sites)
        for initial, final, sign in columns:
            matrix[final, initial] = sign
        result[key] = s.ImmutableMatrix(matrix)
    if len(result) != 2*RAW_SLOTS:
        raise ArithmeticError('padded source records omitted')
    return result


def branch_invariants():
    """Full-domain conservation, including nonadjacent CAR signs."""
    checks = []
    for sites, seam in ((2, (0, 1)), (3, (0, 2))):
        dimension = 4**sites
        for law in ('C', 'E'):
            rows = finite_branches(sites=sites, flavors=2, seam=seam, law=law)
            counts = [0]*dimension
            for columns in rows.values():
                for initial, final, sign in columns:
                    counts[initial] += 1
                    if (initial.bit_count() != final.bit_count() or
                            any((((initial >> (2*x)) & 3) == 3) !=
                                (((final >> (2*x)) & 3) == 3) for x in range(sites))):
                        raise ArithmeticError('source branch changed a doublon label or total number')
            if counts != [2]*dimension:
                raise ArithmeticError('source instrument is incomplete')
            checks.append({'sites': sites, 'seam': list(seam), 'law': law,
                           'full_Fock_dimension': dimension, 'raw_slots': len(rows),
                           'each_doublon_projector_and_total_number_commute': True})
    return checks


def dimer_algebra():
    """Rebuild every source sector and its actual two-dimensional response."""
    h, hopping, defect = exact_source(2, ((0, 1, 1),), F(1), F(1))
    c, u = s.symbols('c u', positive=True)
    h = c*hopping+u*defect/2
    initial = s.eye(16)[:, 9]  # site0 up, site1 down
    bright = s.sqrt(2)*hopping*initial
    partner = ((h-u*s.eye(16))*bright/c).applyfunc(s.simplify)
    dark = initial-partner*(partner.H*initial)[0]
    paired = s.diag(*(int(mask in (3, 12)) for mask in range(16)))
    protected = s.diag(*(int(mask in (5, 6, 9, 10)) for mask in range(16)))
    if (bright.H*bright)[0] != 1 or (partner.H*partner)[0] != 1:
        raise ArithmeticError('source cyclic vectors are not normalized')
    if (h*partner-c*bright).applyfunc(s.expand) != s.zeros(16, 1):
        raise ArithmeticError('source hopping coefficient differs')
    if (h*bright-c*partner-u*bright).applyfunc(s.expand) != s.zeros(16, 1):
        raise ArithmeticError('source contact coefficient differs')
    if h*dark != s.zeros(16, 1) or (dark.H*dark)[0] != s.Rational(1, 2):
        raise ArithmeticError('original dark probability differs')
    if paired*bright != bright or paired*partner != s.zeros(16, 1):
        raise ArithmeticError('collision event differs from source pair projection')
    record_rows = []
    for law in ('C', 'E'):
        for (coin, raw), k in raw_matrices(2, (0, 1), law).items():
            expected = protected if raw == 0 else paired if raw == 8 else s.zeros(16)
            if k*(paired+protected) != expected:
                raise ArithmeticError('fine repair backaction differs on total number two')
            record_rows.append({'law': law, 'coin': coin, 'raw': raw,
                                'probability': '(1-p)/4' if raw == 0 else 'p/4' if raw == 8 else '0'})
    return {'full_Fock_dimension': 16, 'initial_mask': 9,
            'coupled_block': [['0', 'c'], ['c', 'U']],
            'original_coupled_weight': '1/2', 'original_dark_weight': '1/2',
            'complete_records': record_rows, 'complete_record_count': 64,
            'postselection_used': False}


def original_weyl_attachment():
    """The same number-two detector is invariant under all three source frames."""
    selected = s.diag(*(int(mask.bit_count() == 2) for mask in range(16)))
    original_masks = []
    def exact(value):
        return s.Matrix(value.real.toarray())+s.I*s.Matrix(value.imag.toarray())
    for axis in range(3):
        sides = tuple(2 if j == axis else 1 for j in range(3))
        source = symmetry.operators(FieldGeometry(sides))
        frame, original = exact(source.frame), exact(source.K)/source.scale
        positive, _, _ = exact_source(2, ((0, 1, 1),), F(1, 8), F(1))
        if frame*original*frame.H != positive-s.eye(16):
            raise ArithmeticError('actual original Weyl and signed source clocks differ')
        _, unit_hopping, unit_defect = exact_source(2, ((0, 1, 1),), F(1), F(1))
        for coupling, contact, target in ((1., 0., unit_hopping),
                                          (0., 1., unit_defect/2-s.eye(16))):
            piece = symmetry.operators(FieldGeometry(sides, coupling, contact))
            if frame*exact(piece.K)*frame.H/piece.scale != target:
                raise ArithmeticError('independent source hopping/contact frame identity failed')
        for law in ('C', 'E'):
            for matrix in raw_matrices(2, (0, 1), law).values():
                if (matrix*frame-frame*matrix)*selected != s.zeros(16):
                    raise ArithmeticError('number-two bare repair does not intertwine this frame')
        initial = frame.H*s.eye(16)[:, 9]
        support = [i for i in range(16) if initial[i]]
        if len(support) != 1 or abs(initial[support[0]]) != 1:
            raise ArithmeticError('source frame requires an uncounted coherent initial superposition')
        original_masks.append(support[0])
    return {'original_axes': [0, 1, 2], 'original_initial_masks': original_masks,
            'all_16_source_rows_and_all_64_fine_record_maps_checked': True,
            'independent_hopping_and_contact_coefficients_checked': True,
            'bare_repair_covariance_scope': 'total occupation two on this dimer',
            'full_domain_bare_repair_covariance_asserted': False}


def collision_probability(coupling, contact, time, terms=8):
    c, u, t = (_rational(x, positive=True) for x in (coupling, contact, time))
    bright = probability(1, c, u, t, terms)
    return Interval(bright.low/2, bright.high/2)


def identify_action(time, p1, p2, coupling_upper, contact_upper):
    """Deterministic complete-record intervals, with no zero-time derivative."""
    t, cmax, umax = (_rational(x, positive=True)
                     for x in (time, coupling_upper, contact_upper))
    doubled = []
    for interval in (p1, p2):
        if not isinstance(interval, Interval) or interval.low < 0 or interval.high > 1:
            raise ValueError('complete collision probability interval in [0,1] required')
        lo, hi = max(F(), interval.low), min(F(1, 2), interval.high)
        if lo > hi:
            raise ValueError('collision data contradict the declared dimer preparation')
        doubled.append(Interval(2*lo, 2*hi))
    result = infer_coefficients(1, t, *doubled, cmax, umax)
    return {'coupling_squared': result['hopping_squared'],
            'contact_squared': result['electric_squared'],
            'contact_over_coupling_squared': result['electric_over_hopping_squared'],
            'frequency_squared': result['omega_squared']}


def sampling_budget(probability_error, failure_probability):
    """Rational sufficient Hoeffding count for two independent time settings.

    log(4/delta) <= ceil(log2(4/delta)); no floating logarithm is used.
    """
    epsilon = _rational(probability_error, positive=True)
    delta = _rational(failure_probability, positive=True)
    if epsilon >= 1 or delta >= 1:
        raise ValueError('probability error and confidence failure must be below one')
    target = 4/delta
    exponent = max(0, target.numerator.bit_length()-target.denominator.bit_length())
    while F(2)**exponent < target:
        exponent += 1
    count = ceil(F(exponent, 2)/epsilon**2)
    return {'trials_per_time': count, 'total_trials': 2*count,
            'simultaneous_failure_probability_upper': str(delta),
            'original_collision_indicator': True,
            'independent_identically_prepared_trials_required': True}


def fifo_archive_equivalence():
    """Actual native executions differ only in retained/exported partition."""
    basis = native.local.basis(True)
    x = basis.index(tuple(int(k in (0, 4)) for k in range(12)))
    y = basis.index(tuple(int(k in (2, 4)) for k in range(12)))
    density = {(a, r, b, q): F(1, 2) for a, r in ((x, 0), (y, 1))
               for b, q in ((x, 0), (y, 1))}
    histories = []
    for capacity in (1, 2):
        p = native.program(biased=False, capacity=capacity, commands=[{
            'op': 'repair', 'selector': {'kind': 'fixed', 'seam': [0, 1]},
            'law': 'mixture', 'repeat': 2}])
        h = native.execute(p, density)
        if h.status != 'COMPLETE':
            raise ArithmeticError('archive comparison did not execute completely')
        histories.append(h)
    for left, right in zip(histories[0].branches, histories[1].branches, strict=True):
        if left.archive+left.memory != right.archive+right.memory or left.density != right.density:
            raise ArithmeticError('FIFO relocation changed a complete source/reference history')
    return {'capacities': [1, 2], 'attempts': 2,
            'allocated_history_blocks': len(histories[0].branches),
            'complete_source_reference_record_equality': True,
            'policy_does_not_depend_on_memory_capacity': True,
            'overflow_policy': 'export_oldest', 'full_record_order_retained': True}


def certificate():
    geometry = FieldGeometry((2, 1, 1))
    if geometry.speed != 1/8 or geometry.interaction != 1:
        raise ArithmeticError('configured finite source defaults differ')
    rows = []
    for denominator in (32, 64, 128):
        t = F(1, denominator)
        p1 = collision_probability(F(1, 8), F(1), t)
        p2 = collision_probability(F(1, 8), F(1), 2*t)
        inferred = identify_action(t, p1, p2, F(1, 4), F(2))
        if not (inferred['coupling_squared'].contains(F(1, 64))
                and inferred['contact_squared'].contains(F(1))
                and inferred['contact_over_coupling_squared'].contains(F(64))):
            raise ArithmeticError('complete source records exclude the configured action')
        rows.append({'time': str(t), 'probabilities': [p1.encoded(), p2.encoded()],
                     'inferred': {key: value.encoded() for key, value in inferred.items()}})
    return {'source': dimer_algebra(), 'original_Weyl': original_weyl_attachment(),
            'branch_invariants': branch_invariants(),
            'FIFO_archive': fifo_archive_equivalence(), 'default_action_readout': rows,
            'sampling_example': sampling_budget(F(1, 2**20), F(1, 100)),
            'configured_c': '1/8', 'configured_U': '1',
            'source_clock_and_known_initial_preparation_required': True,
            'thermodynamic_ground_selected_by_this_readout': False,
            'identified_quantity': 'declared Hubbard contact U over hopping parameter c',
            'physical_strong_coupling_or_QCD_alpha_selected': False,
            'coherent_clock_pulse_tape_executed_here': False,
            'capacity_selects_action_or_physical_clock': False,
            'physical_particle_or_mass_identified': False}
