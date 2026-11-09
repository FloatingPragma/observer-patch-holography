"""Exact local sixth/eighth moments of the specified centered cubic source.

The sparse words determine observables. Their selected-ground expectations
are inputs to the moment bound, not values supplied by this calculation.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
import json
from math import gcd

import numpy as np

from . import composite, momentum_response as momentum


def _copy(value):
    if type(value) is F:
        return F(value.numerator, value.denominator)
    if type(value) is dict:
        return {k: _copy(v) for k, v in value.items()}
    if type(value) is tuple:
        return tuple(_copy(v) for v in value)
    if type(value) is list:
        return [_copy(v) for v in value]
    return value


def _sites(radius):
    if type(radius) is not int or not 1 <= radius <= 6:
        raise ValueError('integer local source radius in 1..6 required')
    return tuple(sorted(x for x in product(range(-radius, radius+1), repeat=3)
                        if sum(map(abs, x)) <= radius))


def _source(sites):
    """Numerator4W at c=1,U=8, in actual signed Majorana coordinates."""
    sites = tuple(sorted(sites))
    if len(set(sites)) != len(sites):
        raise ValueError('distinct cubic source sites required')
    index = {x: j for j, x in enumerate(sites)}
    adjacency = [[] for _ in range(4*len(sites))]
    terms = {}

    def include(mask, re, im):
        terms[mask] = (re, im)
        word = mask
        while word:
            bit = word & -word
            adjacency[bit.bit_length()-1].append(mask)
            word ^= bit

    for x, j in index.items():
        include(15 << (4*j), -8, 0)
        for axis in range(3):
            y = list(x)
            y[axis] += 1
            y = tuple(y)
            if y not in index:
                continue
            k, phase = index[y], int(momentum.signed_link(x, axis))
            for spin in range(2):
                include((1 << (4*j+2*spin)) | (1 << (4*k+2*spin+1)), 0, phase)
                include((1 << (4*j+2*spin+1)) | (1 << (4*k+2*spin)), 0, -phase)
    return sites, index, adjacency, terms


def _derivative(operator, adjacency, terms):
    """Exact numerator of [O,4W]; only overlapping even source words act."""
    out = {}
    for mask, (ar, ai) in operator.items():
        word, seen = mask, set()
        while word:
            bit = word & -word
            word ^= bit
            for term in adjacency[bit.bit_length()-1]:
                if term in seen:
                    continue
                seen.add(term)
                if (mask & term).bit_count() % 2 == 0:
                    continue
                br, bi = terms[term]
                factor = 2*momentum._product_sign(mask, term)
                key = mask ^ term
                re, im = out.get(key, (0, 0))
                re += factor*(ar*br-ai*bi)
                im += factor*(ar*bi+ai*br)
                if re or im:
                    out[key] = (re, im)
                else:
                    out.pop(key, None)
    return out


def _b_contraction(operator, sites):
    """Hermitian nodal contraction with b*; its common denominator is external."""
    eligible = {j: (-1 if (sum(x)//2) % 2 else 1)
                for j, x in enumerate(sites) if all(v % 2 == 0 for v in x)}
    out = {}
    for mask, (ar, ai) in operator.items():
        word, seen = mask, set()
        while word:
            bit = word & -word
            word ^= bit
            site = (bit.bit_length()-1)//4
            if site in seen or site not in eligible:
                continue
            seen.add(site)
            for local in (13, 14):
                target = local << (4*site)
                if (mask & target).bit_count() % 2 == 0:
                    continue
                factor = eligible[site]*momentum._product_sign(mask, target)
                re, im = (-ai, ar) if local == 13 else (ar, ai)
                key = mask ^ target
                degree = key.bit_count()
                phase = (degree*(degree-1)//2) % 4
                coefficient = factor*(re, im, -re, -im)[phase]
                if coefficient:
                    value = out.get(key, 0)+coefficient
                    if value:
                        out[key] = value
                    else:
                        out.pop(key, None)
    return out


def transform_word(word, perm=(0, 1, 2), signs=(1, 1, 1), shift=(0, 0, 0)):
    """Actual cubic lift followed by a magnetic translation of an even word.

    Returns sorted coordinates and the complete gauge/fermionic sign. It
    implements an algebra automorphism, before any invariant-state quotient.
    """
    word = tuple(tuple(x) for x in word)
    if (len(word) % 2 or len(set(word)) != len(word)
            or any(len(x) != 4 or any(type(v) is not int for v in x)
                   or not 0 <= x[3] < 4 for x in word)
            or any(type(v) is not int for v in perm) or tuple(sorted(perm)) != (0, 1, 2)
            or len(signs) != 3 or any(type(v) is not int or v not in (-1, 1) for v in signs)
            or len(shift) != 3 or any(type(v) is not int for v in shift)):
        raise ValueError('distinct even Majorana word and exact cubic transform required')
    moved, exponent = [], 0
    for x, y, z, local in word:
        point = (x, y, z)
        target = momentum.rotate(point, perm, signs)
        exponent += sum(point[i]*point[j] for i in range(3) for j in range(i+1, 3)
                        if perm[i] > perm[j])
        exponent += shift[1]*target[0]+shift[2]*(target[0]+target[1])
        moved.append(tuple(target[j]+shift[j] for j in range(3))+(local,))
    exponent += sum(moved[i] > moved[j] for i in range(len(moved)) for j in range(i+1, len(moved)))
    return tuple(sorted(moved)), (-1 if exponent % 2 else 1)


def _anchor(word):
    if not word:
        return (), 1
    return transform_word(word, shift=tuple(-v for v in word[0][:3]))


def _translated(operator, sites):
    out = {}
    for mask, coefficient in operator.items():
        word, bits = [], mask
        while bits:
            bit = bits & -bits
            bits ^= bit
            j = bit.bit_length()-1
            word.append(sites[j//4]+(j % 4,))
        key, factor = _anchor(tuple(word))
        value = out.get(key, 0)+factor*coefficient
        if value:
            out[key] = value
        else:
            out.pop(key, None)
    return out


_ROTATIONS = tuple((p, s) for p in permutations(range(3))
                   for s in product((-1, 1), repeat=3))


def _invariant_reduction(operator):
    """Equal expectation on the specified state; not an operator identity."""
    cache, out, zero = {}, {}, 0
    for word, coefficient in operator.items():
        if word not in cache:
            orbit = {}
            for perm, signs in _ROTATIONS:
                moved, factor = transform_word(word, perm, signs)
                moved, translated = _anchor(moved)
                orbit.setdefault(moved, set()).add(factor*translated)
            if any(len(values) != 1 for values in orbit.values()):
                cache.update((key, ((), 0)) for key in orbit)
            else:
                canonical = min(orbit)
                sign = next(iter(orbit[canonical]))
                cache.update((key, (canonical, sign*next(iter(values))))
                             for key, values in orbit.items())
        key, factor = cache[word]
        if factor:
            value = out.get(key, 0)+factor*coefficient
            if value:
                out[key] = value
            else:
                out.pop(key, None)
        else:
            zero += 1
    return out, zero


def _digest(operator, denominator):
    digest = sha256()
    digest.update((str(denominator)+'\n').encode('ascii'))
    for word, coefficient in sorted(operator.items()):
        digest.update((json.dumps([word, coefficient], separators=(',', ':'))+'\n').encode('ascii'))
    return digest.hexdigest()


@lru_cache(maxsize=1)
def _readouts():
    sites, index, adjacency, terms = _source(_sites(6))
    origin = 4*index[(0, 0, 0)]
    operator = {13 << origin: (0, 1), 14 << origin: (-1, 0)}  # numerator2b
    rows, counts = {}, []
    expected = {4: (211, 181, 19, 328, 1264),
                6: (11857, 9781, 493, 8326, 277408),
                8: (997010, 843537, 27965, 278062, 133578340)}
    for n in range(1, 7):
        operator = _derivative(operator, adjacency, terms)
        counts.append(len(operator))
        if n % 2:
            continue
        order = n+2
        raw = _b_contraction(operator, sites)
        raw_count = len(raw)
        translated = _translated(raw, sites)
        del raw
        translated_count = len(translated)
        reduced, zero = _invariant_reduction(translated)
        del translated
        denominator = 2*4**n//16
        common = denominator
        for coefficient in reduced.values():
            common = gcd(common, abs(coefficient))
        reduced = {word: coefficient//common for word, coefficient in reduced.items()}
        denominator //= common
        bound = F(sum(abs(c) for c in reduced.values()), denominator)
        constant = F(reduced.get((), 0), denominator)
        observed = (raw_count, translated_count, len(reduced), constant, bound)
        if observed != expected[order] or zero:
            raise ArithmeticError('actual local nodal moment reduction changed')
        rows[order] = (reduced, denominator, {
            'moment_order': order, 'c': F(1), 'U': F(8),
            'local_derivative_order_of_b': n,
            'local_derivative_word_counts': tuple(counts),
            'unsymmetrized_Hermitian_words': raw_count,
            'magnetic_translation_words': translated_count,
            'invariant_Hermitian_words': len(reduced),
            'maximum_occupied_sites_per_word': max(len({w[:3] for w in word}) for word in reduced),
            'source_observable_norm_upper': bound,
            'tracial_control_constant': constant,
            'selected_ground_moment': None,
            'common_denominator': denominator,
            'sparse_observable_sha256': _digest(reduced, denominator),
            'quotient_is_invariant_state_expectation_equality': True,
            'quotient_is_operator_equality': False})
    return rows


def nodal_moment(order):
    if type(order) is not int or order not in (4, 6, 8):
        raise ValueError('source moment order4,6,or8 required')
    return _copy(_readouts()[order][2])


def score_plan(order):
    nodal_moment(order)
    operator, denominator, _ = _readouts()[order]
    total = sum(abs(c) for c in operator.values())
    return {'moment_order': order, 'score_absolute_bound': F(total, denominator),
            'settings': tuple({'word': word, 'probability': F(abs(c), total),
                               'plus_score': F(total if c > 0 else -total, denominator)}
                              for word, c in sorted(operator.items())),
            'both_pointer_outcomes_per_setting_retained': True,
            'identity_word_negative_outcome_is_zero': True,
            'large_native_readout_executed': False,
            'state_or_service_error_implied_zero': False}


def christoffel(s, t, u, alpha=4):
    """Exact source-moment bound, including the singular memory cases."""
    if any(type(x) not in (int, F) for x in (s, t, u, alpha)):
        raise ValueError('exact rational source moments required')
    s, t, u, alpha = (F(x.numerator, x.denominator) if type(x) is F else F(x)
                       for x in (s, t, u, alpha))
    if (min(s, t, u) < 0 or alpha <= 0 or s*u < t*t
            or (not s and (t or u)) or (not t and u)):
        raise ValueError('positive source Gram and Cauchy constraints required')
    a, h = alpha*alpha, alpha*alpha+s
    if not u:
        if t:
            raise ValueError('zero fourth memory moment requires zero second moment')
        value, coefficients, case = s/h, (F(1), -1/h, F(0)), 'memory_in_kernel'
    else:
        denominator = h*u-t*t
        if denominator <= 0:
            raise ArithmeticError('source moment Gram denominator is not positive')
        value = (s*u-t*t)/denominator
        coefficients = (F(1), -(h*t+u)/denominator, t/denominator)
        case = 'positive_fourth_memory_moment'
    return {'original_moments_0_2_4_6_8': (F(1), a, a*h, a*(h*h+t),
                                           a*(h**3+2*h*t+u)),
            'polynomial_coefficients_1_E2_E4': coefficients,
            'original_zero_atom_upper': value, 'degree_two_upper': s/h,
            'case': case, 'selected_source_values_certified_by_helper': False,
            'actual_default_source_atom_positive': False}


def onsite_force_control():
    """Exact seven-site force identity on the full CAR algebra."""
    force = momentum.scale(momentum.add(*(momentum.bond_derivative(0, j)
                                          for j in range(1, 7))), F(1, 2))
    onsite = {15 << (4*j): (F(-2), F(0)) for j in range(7)}
    first = momentum.commutator(force, onsite)
    second = momentum.commutator(first, onsite)
    if second != momentum.scale(force, 16):
        raise ArithmeticError('source onsite force identity failed')
    difference = momentum.add(momentum.anticommutator(first, momentum.dagger(first)),
                              momentum.scale(momentum.anticommutator(force, momentum.dagger(force)), -16))
    if len(difference) != 300:
        raise ArithmeticError('full-CAR onsite norm distinction changed')
    return {'c': F(1), 'U': F(8), 'force_words': len(force),
            'onsite_derivative_squared_multiplier': F(16),
            'nonzero_anticommutator_norm_difference_words': len(difference),
            'operator_identity_before_state_selection': True,
            'onsite_derivation_selfadjoint_in_interacting_ground_assumed': False,
            'ground_norm_equality_or_charge_gap_inferred': False}


def _trace_product(a, b):
    return sum(int(x) for x in (a*b.T).flat)


@lru_cache(maxsize=1)
def _finite_control():
    """Independent full64-state actual-source moment contraction."""
    data = composite.source((3, 1, 1), F(2), F(3, 7))
    rho, denominator, rank = composite.ground(data, F(-45, 14))
    rho, h = rho.astype(object), data['H'].astype(object)
    scale, alpha = data['scale'], F(3, 14)
    # h(1,0,-1)=0. The factor2 below is the squared smearing norm.
    a = (data['a'][0]-data['a'][4]).astype(object)
    b = (data['b'][0]-data['b'][4]).astype(object)
    delta_a, delta_b = a.copy(), b.copy()
    b_moments = {0: F(_trace_product(rho, b@b.T+b.T@b), 2*denominator)}
    moments = {0: F(_trace_product(rho, a@a.T+a.T@a), 2*denominator)}
    for n in range(1, 9):
        delta_a = delta_a@h-h@delta_a
        if n <= 6:
            delta_b = delta_b@h-h@delta_b
        if n % 2 == 0:
            moments[n] = F(_trace_product(rho, delta_a@a.T+a.T@delta_a),
                           2*denominator*scale**n)
            if n <= 6:
                b_moments[n] = F(_trace_product(rho, delta_b@b.T+b.T@delta_b),
                                2*denominator*scale**n)
            if moments[n] != alpha*alpha*b_moments[n-2]:
                raise ArithmeticError('actual finite nodal source identity failed')
    s = moments[4]/alpha**2-alpha**2
    h0 = alpha**2+s
    t = moments[6]/alpha**2-h0*h0
    u = moments[8]/alpha**2-h0**3-2*h0*t
    result = christoffel(s, t, u, alpha)
    if result['original_moments_0_2_4_6_8'] != tuple(moments[j] for j in (0, 2, 4, 6, 8)):
        raise ArithmeticError('actual full-ground moment reconstruction differs')
    onsite = np.diag(np.diag(h))
    kinetic = h-onsite
    center_b = data['b'][2].astype(object)
    force = center_b@kinetic-kinetic@center_b
    onsite_force = force@onsite-onsite@force
    if np.any(onsite_force@onsite-onsite@onsite_force != int(alpha*scale)**2*force):
        raise ArithmeticError('finite source onsite-square identity differs')
    force_norm = F(_trace_product(rho, force@force.T+force.T@force), denominator*scale**2)
    onsite_norm = F(_trace_product(rho, onsite_force@onsite_force.T+onsite_force.T@onsite_force),
                   denominator*scale**4)
    if onsite_norm-alpha**2*force_norm != F(-216, 805):
        raise ArithmeticError('actual-ground onsite norm countercontrol changed')
    return {'whole_Fock_dimension': 64, 'full_ground_rank': len(h)-rank,
            'c': F(2), 'U': F(3, 7), 'ground_energy': F(-45, 14),
            'smearing_numerator': (1, 0, -1), 'smearing_squared_norm': 2,
            'original_signed_moments': tuple(moments[j] for j in (0, 2, 4, 6, 8)),
            'memory_s_t_u': (s, t, u), 'Christoffel_bound': result,
            'center_force_ground_norm_squared': force_norm,
            'onsite_center_force_ground_norm_squared': onsite_norm,
            'onsite_norm_difference_from_alpha_squared': F(-216, 805),
            'finite_graph_zero_hopping_mode_is_infinite_lattice_momentum': False,
            'finite_ground_replaces_default_infinite_ground': False}


def finite_control():
    return _copy(_finite_control())


def memory_budget(eta=F(1, 8)):
    """Actual nodal memory floor from the generated source moment bounds."""
    from . import nodal_memory
    if type(eta) not in (int, F) or eta <= 0:
        raise ValueError('exact positive imaginary source frequency required')
    eta = F(eta.numerator, eta.denominator) if type(eta) is F else F(eta)
    six = nodal_moment(6)['source_observable_norm_upper']
    eight = nodal_moment(8)['source_observable_norm_upper']
    t_upper, u_upper = six/16-16**2, eight/16-16**3
    lower = nodal_memory.default_budget(F(1), eta)['node_Gram_strict_lower']
    if (t_upper, u_upper, lower) != (F(17082), F(33378201, 4), F(36, 125)):
        raise ArithmeticError('source moment and actual nodal lower bounds differ')
    return {'c': F(1), 'U': F(8), 'eta': eta,
            'source_s_strict_lower': lower, 'source_t_upper': t_upper,
            'source_u_upper': u_upper,
            'M_eta_lower': eta*eta*lower*lower/(eta*eta*lower+t_upper),
            'lower_bound_is_order_eta_squared': True,
            'zero_frequency_atom_inferred': False}


def certificate():
    return {'same_source_normalization': {'c': F(1), 'U': F(8), 'alpha': F(4)},
            'actual_local_nodal_moment_readouts': tuple(nodal_moment(n) for n in (4, 6, 8)),
            'onsite_force_identity': onsite_force_control(),
            'independent_full_source_control': finite_control(),
            'actual_source_memory_bound': memory_budget(),
            'source_selected_ground_moments_evaluated': False,
            'actual_zero_atom_decided': False,
            'physical_species_identified': [], 'physical_masses': []}
