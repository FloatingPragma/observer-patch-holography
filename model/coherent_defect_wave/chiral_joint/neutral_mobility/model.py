"""Exact neutral composite mobility without compressing the quantum fields."""
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache

import sympy as s

from ...chiral_source.model import PAULI
from ..model import clean, sequence
from ..periodic_family.model import (Configuration as BaseConfiguration,
    cells, cell_index, shift_cell, edges, apply_state, norm2)
from ..oriented_links.model import (
    column as base_column, validate_key, empty_key, gauss_column,
    descriptor as base_descriptor,
)
from ..higgs_composites.model import neutral_composite, adjoint


@dataclass(frozen=True)
class Configuration(BaseConfiguration):
    neutral_hopping: object = s.Integer(1)

    def __post_init__(self):
        super().__post_init__()
        value = s.sympify(self.neutral_hopping)
        if value.has(s.Float) or value.is_real is not True or value.is_finite is not True:
            raise ValueError('finite exact real neutral hopping coefficient required')
        object.__setattr__(self, 'neutral_hopping', value)


@lru_cache(maxsize=64)
def neutral_words(config=Configuration()):
    """Every normal-ordered CAR/scalar matrix word, including its adjoint.

    H_N=(i lambda/2) sum_(x,i) q(x,i) D_x^dagger sigma_i D_(x+e_i)
    + adjoint, where D=J/sqrt(2). Scalar memories are quantum operators.
    """
    out = defaultdict(lambda: s.Integer(0))
    for xyz, axis in edges(config.side):
        left = cell_index(config.side, xyz)
        right = cell_index(config.side, shift_cell(config.side, xyz, axis))
        seam = -1 if xyz[axis] == config.side - 1 else 1
        for (u, v), spin_weight in PAULI[axis].todok().items():
            weight = s.I * config.neutral_hopping * seam * spin_weight / 4
            for (word_x, row_x, col_x), value_x in adjoint(neutral_composite(u)).items():
                for (word_y, row_y, col_y), value_y in neutral_composite(v).items():
                    word = (tuple((30 * left + i, create) for i, create in word_x)
                            + tuple((30 * right + i, create) for i, create in word_y))
                    scalars = tuple(sorted(((left, row_x, col_x), (right, row_y, col_y))))
                    value = weight * value_x * value_y
                    out[word, scalars] += value
                    dagger_word = tuple((i, not create) for i, create in reversed(word))
                    dagger_scalars = tuple((cell, col, row) for cell, row, col in scalars)
                    out[dagger_word, dagger_scalars] += s.conjugate(value)
    return tuple((word, scalars, value) for (word, scalars), value in clean(out).items())


def hopping_column(key, config=Configuration()):
    mask, links, scalars = validate_key(key, config)
    out = defaultdict(lambda: s.Integer(0))
    for word, scalar_changes, coefficient in neutral_words(config):
        if any(scalars[cell] != col for cell, _, col in scalar_changes):
            continue
        step = sequence(mask, word)
        if step is None:
            continue
        target_scalars = list(scalars)
        for cell, row, _ in scalar_changes:
            target_scalars[cell] = row
        out[step[0], links, tuple(target_scalars)] += coefficient * step[1]
    return clean(out)


def column(key, config=Configuration()):
    out = defaultdict(lambda: s.Integer(0), base_column(key, config))
    for target, value in hopping_column(key, config).items():
        out[target] += value
    return clean(out)


def composite_column(key, cell, spin, *, create=False, config=Configuration()):
    """The full-domain D or D^dagger action, with the global CAR signs."""
    mask, links, scalars = validate_key(key, config)
    if cell not in range(config.side ** 3) or spin not in (0, 1):
        raise ValueError('one source cell and one Weyl component required')
    poly = adjoint(neutral_composite(spin)) if create else neutral_composite(spin)
    out = defaultdict(lambda: s.Integer(0))
    for (word, row, col), coefficient in poly.items():
        if scalars[cell] != col:
            continue
        global_word = tuple((30 * cell + i, creation) for i, creation in word)
        step = sequence(mask, global_word)
        if step is None:
            continue
        target_scalars = list(scalars)
        target_scalars[cell] = row
        out[step[0], links, tuple(target_scalars)] += coefficient * step[1] / s.sqrt(2)
    return clean(out)


def preparation(cell=0, spin=0, config=Configuration()):
    result = composite_column(empty_key(config), cell, spin, create=True, config=config)
    assert norm2(result) == 1
    return result


def star_signature(key, config=Configuration()):
    mask, links, _ = validate_key(key, config)
    parity = [((mask >> (30 * cell)) & ((1 << 30) - 1)).bit_count() % 2
              for cell in range(config.side ** 3)]
    for edge, (xyz, axis) in enumerate(edges(config.side)):
        if links[edge] != 0:
            left = cell_index(config.side, xyz)
            right = cell_index(config.side, shift_cell(config.side, xyz, axis))
            parity[left] ^= 1
            parity[right] ^= 1
    return tuple(-1 if value else 1 for value in parity)


def target_signature(cell, config=Configuration()):
    if cell not in range(config.side ** 3):
        raise ValueError('target cell outside source')
    return tuple(-1 if site == cell else 1 for site in range(config.side ** 3))


def star_arrival_projector(state, cell, config=Configuration()):
    target = target_signature(cell, config)
    return {key: value for key, value in state.items() if star_signature(key, config) == target}


def arrival_projector(state, cell, config=Configuration()):
    """Local neutral spin-doublet projector times the target's negative star.

    Every state outside this explicitly specified channel is retained in the
    complementary outcome. All other sites and link occupations are kept.
    """
    if cell not in range(config.side ** 3):
        raise ValueError('arrival cell outside source')
    out = defaultdict(lambda: s.Integer(0))
    for key, coefficient in state.items():
        mask, links, scalars = validate_key(key, config)
        if star_signature(key, config)[cell] != -1:
            continue
        local_mask = (mask >> (30 * cell)) & ((1 << 30) - 1)
        for spin in range(2):
            components = ((1 << (15 * spin + 11), 1, 1),
                          (1 << (15 * spin + 10), 2, -1))
            overlap = next((sign for local, scalar, sign in components
                            if local == local_mask and scalar == scalars[cell]), 0)
            if not overlap:
                continue
            for target, scalar, sign in components:
                target_mask = mask ^ ((local_mask ^ target) << (30 * cell))
                target_scalars = list(scalars)
                target_scalars[cell] = scalar
                out[target_mask, links, tuple(target_scalars)] += s.sympify(coefficient) * s.Rational(overlap * sign, 2)
    return clean(out)


def complete_readout(state, cell, config=Configuration()):
    arrival = arrival_projector(state, cell, config)
    no_arrival = defaultdict(lambda: s.Integer(0), state)
    for key, value in arrival.items():
        no_arrival[key] -= value
    no_arrival = clean(no_arrival)
    return arrival, no_arrival


def descriptor(config=Configuration()):
    return {
        'schema': 'oph.chiral_neutral_mobility.v1',
        'base_architecture': base_descriptor(config),
        'composite': '(phi_0 c_11 - phi_1 c_10)/sqrt(2), independently for each spin',
        'neutral_hopping': str(config.neutral_hopping),
        'additional_operator': '(i lambda/2) sum_xi AP_seam D_x^dagger sigma_i D_(x+e_i) + adjoint',
        'coupling_ancestry': 'declared source architecture, independent of mass data',
        'clock': 'exp(-it (H_oriented + H_neutral)) on the full carrier',
        'additional_term_spatial_range_cells': 1,
        'all_base_interactions_retained': True,
        'canonical_composite_CAR_assumed': False,
        'scalar_replaced_by_expectation': False,
        'old_source_ground_or_spectrum_transferred': False,
        'source_operator_identified_with_measured_particle': False,
        'relativistic_limit_or_mass_claim': False,
    }
