"""Exact local proposal guards; every accepted and rejected record is retained."""
from fractions import Fraction as F
from itertools import product

from model.compact_defect_pair.topology import direction, separation
from model.lattice_texture_selection.shells import admissible_directions, degree, dot, shell


def rational_vector(value):
    if type(value) not in (tuple, list) or len(value) != 3:
        raise ValueError('three exact rational components required')
    if any(type(x) not in (int, F) for x in value):
        raise ValueError('three exact rational components required')
    return tuple(F(x) for x in value)


def initial(half_separation=1):
    a = separation(half_separation)
    return {p: rational_vector(direction(p, a))
            for p in product(range(-a-1, a+2), repeat=3)}


def surfaces(half_separation=1):
    a = separation(half_separation)
    return tuple(shell(1, (0, 0, z))[0] for z in (-a, a))


def parameters(eta):
    if type(eta) not in (int, F) or not 0 < eta < F(1, 4):
        raise ValueError('an exact relative budget in (0,1/4) is required')
    return F(eta)


def validate(state, half_separation=1, eta=F(1, 512)):
    eta = parameters(eta)
    baseline = initial(half_separation)
    if type(state) is not dict or set(state) != set(baseline):
        raise ValueError('complete compact support required')
    values = {p: rational_vector(v) for p, v in state.items()}
    for p, f in baseline.items():
        v = values[p]
        if not any(f):
            if any(v):
                raise ValueError('both mass-zero cores are fixed')
        elif not any(v) or dot(tuple(v[i]-f[i] for i in range(3)),
                              tuple(v[i]-f[i] for i in range(3))) > eta*eta*dot(f, f):
            raise ValueError('noncore direction violates reference budget')
    for triangles in surfaces(half_separation):
        admissible_directions(triangles, {p: values[p] for t in triangles for p in t})
    return values


def charges(state, half_separation=1):
    return [degree(t, {p: state[p] for triangle in t for p in triangle})['degree']
            for t in surfaces(half_separation)]


def update(state, point, proposal, half_separation=1, eta=F(1, 512)):
    """A single local record; rejection is an identity transition."""
    old = validate(state, half_separation, eta)
    candidate = rational_vector(proposal)
    baseline = initial(half_separation)
    reason = 'accepted'
    incident = [triangle for surface in surfaces(half_separation)
                for triangle in surface if point in triangle]
    if point not in old:
        reason = 'fixed_exterior'
    elif not any(baseline[point]):
        reason = 'fixed_core'
    elif not any(candidate):
        reason = 'zero_direction'
    else:
        f = baseline[point]
        difference = tuple(candidate[i]-f[i] for i in range(3))
        if dot(difference, difference) > eta*eta*dot(f, f):
            reason = 'reference_budget'
        elif dot(candidate, old[point]) < 0 or any(
                dot(candidate, old[q]) < 0 for triangle in incident
                for q in triangle if q != point):
            reason = 'shell_homotopy'
    result = dict(old)
    if reason == 'accepted':
        result[point] = candidate
    return result, {'point': list(point), 'proposal': [str(x) for x in candidate],
                    'accepted': reason == 'accepted', 'reason': reason,
                    'incident_triangles_checked': len(incident)}


def complete_protocol(half_separation=1, eta=F(1, 512)):
    """Four externally chosen proposals form a complete finite instrument."""
    a = separation(half_separation)
    eta = parameters(eta)
    state = initial(a)
    point = (0, 0, 0)
    f = state[point]
    proposals = [(point, (eta, 0, f[2])), (point, (-eta, 0, f[2])),
                 (point, (0, 0, 0)), ((0, 0, a), (0, 0, 1))]
    branches = []
    for label, (site, proposal) in enumerate(proposals):
        target, record = update(state, site, proposal, a, eta)
        branches.append((target, dict(label=label, probability='1/4',
                                      local_shell_degrees=charges(target, a), **record)))
    if sum(F(row['probability']) for _, row in branches) != 1:
        raise ArithmeticError('incomplete instrument')
    if any(row['local_shell_degrees'] != [-1, 1] for _, row in branches):
        raise ArithmeticError('an admitted record changed shell degree')
    return state, branches


def continuing_protocol():
    state, initial_branches = complete_protocol()
    states = [state, initial_branches[0][0], initial_branches[1][0]]
    rows, transition = [], [[F(0) for _ in states] for _ in states]
    for before, current in enumerate(states):
        for _, proposal in initial_branches:
            target, record = update(current, tuple(proposal['point']),
                                    tuple(F(x) for x in proposal['proposal']))
            after = next(i for i, value in enumerate(states) if value == target)
            transition[before][after] += F(1, 4)
            rows.append(dict(record, source_state=before, target_state=after,
                             proposal_label=proposal['label'], probability='1/4',
                             local_shell_degrees=charges(target)))
    if any(sum(row) != 1 for row in transition):
        raise ArithmeticError('every retained source state needs a complete law')
    if any(row['local_shell_degrees'] != [-1, 1] for row in rows):
        raise ArithmeticError('a continuing record changed shell degree')
    return states, rows, transition


def exact_control():
    state, branches = complete_protocol()
    _, all_rows, transition = continuing_protocol()
    eta = F(1, 512)
    _, bad_guard = update(state, (-1, -1, -1), (-1, -1, -F(1003, 3000)),
                          eta=F(6, 25))
    if bad_guard['reason'] != 'shell_homotopy':
        raise ArithmeticError('a shell-only negative control is required')
    erased = {p: ((F(0), F(0), F(1)) if any(v) else v) for p, v in state.items()}
    erased_charges = charges(erased)
    if erased_charges != [0, 0]:
        raise ArithmeticError('unguarded texture-erasure control failed')
    return {'relative_budget': str(eta),
            'all_instrument_records': [row for _, row in branches],
            'all_continuing_instrument_records': all_rows,
            'complete_texture_transition_matrix': [[str(x) for x in row] for row in transition],
            'stationary_texture_law': ['0', '1/2', '1/2'],
            'accepted_records': sum(row['accepted'] for _, row in branches),
            'complete_probability': '1', 'shared_face_update_checked_on_both_shells': True,
            'shell_guard_control_with_relative_budget_6_over_25': bad_guard,
            'unguarded_erasure_degrees': erased_charges,
            'erasure_retains_four_zero_onsite_coordinates': True,
            'degree_is_fermion_number': False}
