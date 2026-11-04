"""Exact local Gauss identities and complete records on a finite link.

The link representation and coherent hopping are declared simulator design.
They use the existing native matter representation without identifying its
finite transport frequency with a physical mass or continuum gauge field.
"""
from functools import lru_cache
from itertools import combinations

import sympy as s

from model.common_source import native_matter_attachment as matter


def clean(value):
    return s.expand(s.simplify(value))


def exact_real(value, *, nonnegative=False):
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols or value.is_real is not True:
        raise ValueError('an exact finite real scalar is required')
    if value.is_finite is not True or (nonnegative and value.is_nonnegative is not True):
        raise ValueError('a finite nonnegative scalar is required')
    return value


@lru_cache(maxsize=1)
def native_currents():
    """The actual port-admitted basis, lifted to the fifteen selected masks."""
    selected = tuple(i for i in range(32) if i.bit_count() in (2, 4))
    currents = tuple(matter.exterior_derivative(matter.spin_lift(g)).extract(selected, selected)
                     for g in matter.canonical_currents())
    if any(g + g.H != s.zeros(15) for g in currents):
        raise ArithmeticError('anti-Hermitian native currents required')
    return currents


@lru_cache(maxsize=1)
def native_blocks():
    """Use the actual central current to retain every native multiplet."""
    currents = native_currents()
    center = currents[-1]/s.I
    if center != s.diag(*center.diagonal()):
        raise ArithmeticError('the source center is not diagonal in its mask basis')
    blocks = {}
    for index, charge in enumerate(center.diagonal()):
        blocks.setdefault(charge, []).append(index)
    result = tuple((str(charge), tuple(indices)) for charge, indices in
                   sorted(blocks.items(), key=lambda item: item[0]))
    for _, indices in result:
        outside = set(range(15))-set(indices)
        if any(g[i, j] != 0 or g[j, i] != 0 for g in currents for i in indices for j in outside):
            raise ArithmeticError('the actual source currents mix the declared blocks')
    return result


def add(out, key, value):
    out[key] = out.get(key, 0) + value


def normalized(state):
    return {key: value for key, raw in state.items() if (value := clean(raw)) != 0}


def reference_states(d, indices=None):
    """Unnormalized singlets A and B, of squared norms d and d**2.

    Keys are (matter endpoint, matter label, reference label, link a, link b).
    The reference transforms in the dual representation at endpoint one.
    Link (-1,-1) is its invariant vacuum; (a,b) transforms in V* at zero
    and V at one. The nontrivial static reference is an explicit input.
    """
    if type(d) is not int or d < 1:
        raise ValueError('a positive integer representation dimension is required')
    indices = tuple(range(d)) if indices is None else tuple(indices)
    if not indices or len(set(indices)) != len(indices) or any(type(i) is not int or not 0 <= i < d for i in indices):
        raise ValueError('a nonempty set of distinct representation indices is required')
    a = {(1, i, i, -1, -1): s.Integer(1) for i in indices}
    b = {(0, i, j, i, j): s.Integer(1) for i in indices for j in indices}
    return a, b


def gauss_action(generator, state, vertex, *, link_sign=1, include_reference=True):
    """Literal generator action on sparse complete one-particle/link states."""
    if vertex not in (0, 1) or generator.rows != generator.cols:
        raise ValueError('two endpoints and a square generator are required')
    out = {}
    entries = {(i, j): x for (i, j), x in generator.todok().items() if x != 0}
    for (site, i, r, a, b), amplitude in state.items():
        for (row, col), x in entries.items():
            if site == vertex and col == i:
                add(out, (site, row, r, a, b), amplitude*x)
            if vertex == 1 and include_reference and row == r:
                add(out, (site, i, col, a, b), -amplitude*x)
            if a >= 0:
                if vertex == 0 and row == a:
                    add(out, (site, i, r, col, b), -link_sign*amplitude*x)
                if vertex == 1 and col == b:
                    add(out, (site, i, r, a, row), link_sign*amplitude*x)
    return normalized(out)


def hamiltonian_action(state, d, hopping=1, link_energy=1, *, omit_link=False, blocks=None):
    """One-particle restriction of J sum c0† U c1 + h.c. + epsilon Pe."""
    hopping, link_energy = exact_real(hopping), exact_real(link_energy, nonnegative=True)
    out = {}
    blocks = (tuple(range(d)),) if blocks is None else blocks
    partners = {i: block for block in blocks for i in block}
    if set(partners) != set(range(d)) or sum(map(len, blocks)) != d:
        raise ValueError('link blocks must partition the representation')
    for (site, i, r, a, b), amplitude in state.items():
        if a >= 0:
            add(out, (site, i, r, a, b), link_energy*amplitude)
        if omit_link:
            add(out, (1-site, i, r, a, b), hopping*amplitude)
        elif site == 1 and a == -1:
            for j in partners[i]:
                add(out, (0, j, r, j, i), hopping*amplitude)
        elif site == 0 and i == a:
            add(out, (1, b, r, -1, -1), hopping*amplitude)
    return normalized(out)


def scale(state, value):
    return normalized({key: value*x for key, x in state.items()})


def sum_states(*states):
    out = {}
    for state in states:
        for key, value in state.items():
            add(out, key, value)
    return normalized(out)


def exact_certificate():
    currents = native_currents()
    d = currents[0].rows
    blocks = tuple(indices for _, indices in native_blocks())
    for indices in blocks:
        a, b = reference_states(d, indices)
        if hamiltonian_action(a, d, blocks=blocks) != b:
            raise ArithmeticError('the literal hopping lost A_r -> B_r')
        if hamiltonian_action(b, d, blocks=blocks) != sum_states(scale(a, len(indices)), b):
            raise ArithmeticError('the literal hopping lost B_r -> dim(r) A_r + B_r')
        for current in currents:
            for vertex in (0, 1):
                if gauss_action(current, a, vertex) or gauss_action(current, b, vertex):
                    raise ArithmeticError('the actual native Gauss generator fails')
    rank = s.Matrix.hstack(*(g.reshape(d*d, 1) for g in currents)).rank()
    noncommuting = sum(g*h != h*g for g, h in combinations(currents, 2))
    if rank != 12 or noncommuting == 0:
        raise ArithmeticError('the full native non-Abelian current was lost')
    a = sum_states(*(reference_states(d, indices)[0] for indices in blocks))
    b = sum_states(*(reference_states(d, indices)[1] for indices in blocks))
    omitted = hamiltonian_action(a, d, omit_link=True, blocks=blocks)
    broken_no_link = sum(bool(gauss_action(g, omitted, v)) for g in currents for v in (0, 1))
    broken_sign = sum(bool(gauss_action(g, b, v, link_sign=-1))
                      for g in currents for v in (0, 1))
    broken_reference = sum(bool(gauss_action(g, a, 1, include_reference=False)) for g in currents)
    if min(broken_no_link, broken_sign, broken_reference) == 0:
        raise ArithmeticError('a required negative control failed to fail')
    link_dimension = 1+sum(len(indices)**2 for indices in blocks)
    return {'matter_dimension': d, 'native_current_rank': rank,
            'noncommuting_basis_pairs': noncommuting,
            'link_dimension': link_dimension, 'one_particle_with_reference_dimension': 2*d*d*link_dimension,
            'literal_Gauss_checks': len(currents)*2*2*len(blocks),
            'native_blocks': [{'center_eigenvalue': charge, 'indices': list(indices), 'dimension': len(indices)}
                              for charge, indices in native_blocks()],
            'H_on_A_r': 'B_r', 'H_on_B_r': 'dim(r) A_r + B_r',
            'failures_when_link_omitted': broken_no_link,
            'failures_when_link_orientation_reversed': broken_sign,
            'failures_when_static_reference_omitted': broken_reference}


def complete_history(d=6, *, hopping=1, link_energy=1, duration=s.Rational(1, 8), rounds=4):
    """All binary site records, with exact detector energy accounting.

    Coefficients refer to A/sqrt(d), B/sqrt(d); their Gram matrix is
    diag(1,d). A site measurement acts on the same gauge-singlet states.
    Cayley steps are the declared finite evolution, not a continuum estimate.
    """
    reference_states(d)
    if type(rounds) is not int or not 1 <= rounds <= 12:
        raise ValueError('the retained complete-history domain is one to twelve rounds')
    j = exact_real(hopping)
    e = exact_real(link_energy, nonnegative=True)
    t = exact_real(duration, nonnegative=True)
    h = s.Matrix([[0, d*j], [j, e]])
    gram = s.diag(1, d)
    u = ((s.eye(2)+s.I*t*h/2).inv()*(s.eye(2)-s.I*t*h/2)).applyfunc(clean)
    if (u.H*gram*u-gram).applyfunc(clean) != s.zeros(2):
        raise ArithmeticError('the full Cayley step is not unitary')
    if (h.H*gram-gram*h) != s.zeros(2):
        raise ArithmeticError('the restricted generator is not Hermitian')
    norm = lambda v: clean((v.H*gram*v)[0])
    energy = lambda v: clean((v.H*gram*h*v)[0])
    states = {'': s.Matrix([1, 0])}
    record_rounds = []
    initial_energy = sum((energy(v) for v in states.values()), s.Integer(0))
    work_total = s.Integer(0)
    for index in range(rounds):
        before = clean(sum((energy(v) for v in states.values()), s.Integer(0)))
        evolved = {word: (u*v).applyfunc(clean) for word, v in states.items()}
        after_unitary = clean(sum((energy(v) for v in evolved.values()), s.Integer(0)))
        if after_unitary != before:
            raise ArithmeticError('Cayley evolution changed same-Hamiltonian energy')
        states = {word+label: projection*v for word, v in evolved.items()
                  for label, projection in [('R', s.diag(1, 0)), ('L', s.diag(0, 1))]}
        if clean(sum((norm(v) for v in states.values()), s.Integer(0))) != 1:
            raise ArithmeticError('a complete record branch was dropped')
        after = clean(sum((energy(v) for v in states.values()), s.Integer(0)))
        work = clean(after-before)
        work_total = clean(work_total+work)
        record_rounds.append({'round': index+1, 'energy_before': str(before),
                              'energy_after_unitary': str(after_unitary),
                              'energy_after_readout': str(after), 'detector_work': str(work)})
    if clean(after-initial_energy-work_total) != 0:
        raise ArithmeticError('detector energy accounting failed')
    return {'dimension': d, 'hopping': str(j), 'link_energy': str(e), 'Cayley_duration': str(t),
            'rounds': rounds, 'outcomes': [{'word': w, 'probability': str(norm(v))}
                                          for w, v in sorted(states.items())],
            'normalization': '1', 'energy_accounting': record_rounds,
            'total_detector_work': str(work_total),
            'normalized_H': [['0', f'{j}*sqrt({d})'], [f'{j}*sqrt({d})', str(e)]],
            'squared_continuous_gap': str(clean(e*e+4*d*j*j)),
            'physical_particle_mass': False}
