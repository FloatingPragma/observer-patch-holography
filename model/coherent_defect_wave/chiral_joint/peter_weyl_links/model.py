"""The complete periodic source with literal compressed group coordinates."""
from collections import defaultdict
from functools import lru_cache

import sympy as s

from ...chiral_source.model import PAULI
from ..model import clean, sequence
from ..periodic_family.model import Configuration, kinetic_terms, local_yukawa
from ..oriented_links.model import (
    labels, label_index, validate_key, empty_key, gauss_column,
    translation_column, endpoint_current, orientation,
)
from .clebsch import irreps, native_charges, coordinate


@lru_cache(maxsize=2)
def single_link_columns(adjoint=False):
    out = defaultdict(list)
    for charge in native_charges():
        rep = irreps()[charge]
        for a, ia in enumerate(rep['native_indices']):
            for b, ib in enumerate(rep['native_indices']):
                matrix = coordinate(charge, a, b)
                for (row, col), value in matrix.todok().items():
                    if adjoint:
                        row, col, value = col, row, s.conjugate(value)
                    out[col].append((ia, ib, row, value))
    return {key: tuple(values) for key, values in out.items()}


@lru_cache(maxsize=32768)
def path_transitions(values, adjoint=False):
    """Every coordinate of a path, with its complete occupied-link columns.

    For the adjoint branch a,b remain the original coordinate indices; the
    matter word reverses endpoints. Distinct edges commute with one another.
    """
    if not values or any(type(v) is not int or v not in range(119) for v in values):
        raise ValueError('a nonempty word of valid link states is required')
    columns = single_link_columns(adjoint)
    chains = [(a, b, (target,), value) for a, b, target, value in columns.get(values[0], ())]
    for label in values[1:]:
        chains = [(a, last, targets + (target,), s.expand(value * factor))
                  for a, middle, targets, value in chains
                  for first, last, target, factor in columns.get(label, ())
                  if middle == first]
    return tuple(chains)


def column(key, config=Configuration()):
    """Full 30-CAR/site, scalar-five, link-119 action, without an extra hop."""
    mask, links, scalars = validate_key(key, config)
    out = defaultdict(lambda: s.Integer(0))
    out[key] += (config.link_penalty * sum(v != 0 for v in links)
                 + config.scalar_penalty * sum(v != 0 for v in scalars))
    couplings = config.y_up, config.y_down, config.y_lepton
    for cell, scalar in enumerate(scalars):
        offset = 30 * cell
        local_mask = (mask >> offset) & ((1 << 30) - 1)
        for (target, h), coefficient in local_yukawa(local_mask, scalar, couplings).items():
            target_scalars = list(scalars)
            target_scalars[cell] = h
            target_mask = mask ^ ((local_mask ^ target) << offset)
            out[target_mask, links, tuple(target_scalars)] += coefficient
    for target, value in hopping_column(key, kinetic_terms(config)).items():
        out[target] += value
    return clean(out)


def hopping_column(key, terms):
    """Exact specified kinetic summands, used for supported matrix elements."""
    mask, links, scalars = key
    out = defaultdict(lambda: s.Integer(0))
    for term in terms:
        values = tuple(links[e] for e in term.path)
        for adjoint in (False, True):
            paths = path_transitions(values, adjoint)
            for (u, v), spin_weight in PAULI[term.axis].todok().items():
                weight = term.coefficient * spin_weight
                if adjoint:
                    weight = s.conjugate(weight)
                for a, b, targets, factor in paths:
                    p = 30 * term.left_cell + 15 * u + a
                    q = 30 * term.right_cell + 15 * v + b
                    if adjoint:
                        p, q = q, p
                    step = sequence(mask, ((p, True), (q, False)))
                    if step is None:
                        continue
                    target_links = list(links)
                    for edge, target in zip(term.path, targets, strict=True):
                        target_links[edge] = target
                    out[step[0], tuple(target_links), scalars] += weight * factor * step[1]
    return clean(out)


def descriptor(config=Configuration()):
    volume = config.side ** 3
    return {
        'schema': 'oph.chiral_peter_weyl_links.v1',
        'group': 'connected compact S(U(3) x U(2)), in the actual native exterior representation',
        'retained_irreducible_charges_six_Y': list(irreps()),
        'retained_matrix_coefficient_dimension': 119,
        'full_dimension': f'2^{30 * volume} * 119^{3 * volume} * 5^{volume}',
        'side': config.side, 'spacing': str(config.spacing),
        'source_coefficients': {name: str(getattr(config, name)) for name in
            ('y_up', 'y_down', 'y_lepton', 'scalar_penalty', 'link_penalty')},
        'coordinate': 'P multiplication(D_r,ab) P on normalized Haar matrix coefficients',
        'multiplication_coefficient': 'sqrt(dim(source)/dim(target)) C_left conjugate(C_right)',
        'kinetic_source': 'complete periodic_family.kinetic_terms with these coordinates',
        'onsite_source': 'complete periodic_family.local_yukawa and scalar/link excitation penalties',
        'clock': 'exp(-it H_PW) on the complete finite carrier',
        'coupling_ancestry': 'declared unit coefficients and spacing pi, independent of observed masses',
        'additional_neutral_hopping_term': False,
        'scalar_or_link_expectation_substitution': False,
        'canonical_CAR_assumed_for_operator_valued_transport': False,
        'matter_charge_conjugation_under_rotation': False,
        'ground_or_spectrum_transferred_from_other_link_architecture': False,
        'native_idle_link_implementation_proved': False,
        'relativistic_field_limit_or_measured_particle_identification': False,
    }
