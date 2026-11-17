"""Exact translation characters and physical spin-reference purification."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
from math import gcd, lcm
import json
from pathlib import Path

import sympy as s

from ..periodic_family.model import Configuration, cells, cell_index
from ..higgs_composites.certificate import validate_pins
from ..model import clean

ROOT = Path(__file__).resolve().parents[4]
REFERENCE_RECEIPT = 'model/coherent_defect_wave/chiral_joint/local_reference/results.json'


def reference_data():
    receipt = json.loads((ROOT / REFERENCE_RECEIPT).read_text())
    validate_pins(receipt['input_sha256'])
    source = receipt['source']; reference = receipt['selected_stationary_reference']
    spectrum = receipt['global_physical_ground_and_gap']
    assert source['full_CAR_modes'] == 30 and source['full_scalar_states'] == 5
    assert source['coefficients'] == {'up': '1', 'down': '1', 'charged_lepton': '1', 'scalar_penalty': '1'}
    assert reference['rank'] == spectrum['global_physical_ground_rank'] == 2
    assert spectrum['global_physical_ground_twice_spin'] == 1
    assert reference['fermion_parity'] == 'odd' and reference['raw_L'] == '1'
    assert reference['physical_Gauss_invariant'] and reference['rotation_invariant']
    return reference, spectrum


def validate_side(side):
    Configuration(side=side)
    return side


def translation_order(side, displacement):
    validate_side(side)
    if len(displacement) != 3 or any(type(a) is not int for a in displacement):
        raise ValueError('three integer displacement coordinates required')
    return lcm(*(side // gcd(side, a) for a in displacement))


def source_character(side, displacement):
    """Normalized trace on the full product of physical local ground doublets."""
    order = translation_order(side, displacement)
    volume = side ** 3
    return s.Rational(2) ** (volume // order - volume)


def momentum_multiplicity(side, momentum):
    """Integer Fourier multiplicity, evaluated by divisor/Mobius inversion."""
    validate_side(side)
    if len(momentum) != 3 or any(type(k) is not int or not 0 <= k < side for k in momentum):
        raise ValueError('three momentum indices in 0..N-1 required')
    total = s.Integer(0)
    for order in s.divisors(side):
        exact_order_character_sum = sum(s.mobius(order // divisor) * divisor ** 3
            for divisor in s.divisors(order) if all(k % divisor == 0 for k in momentum))
        total += 2 ** (side ** 3 // order) * exact_order_character_sum
    result = total / side ** 3
    assert result.is_Integer and result >= 0
    return int(result)


def momentum_probability(side, momentum):
    return s.Rational(momentum_multiplicity(side, momentum), 2 ** (side ** 3))


@lru_cache(maxsize=256)
def site_permutation(side, displacement):
    translation_order(side, displacement)
    return tuple(cell_index(side, tuple((xyz[j] + displacement[j]) % side for j in range(3)))
                 for xyz in cells(side))


def shift_spins(mask, side, displacement):
    validate_side(side)
    if type(mask) is not int or not 0 <= mask < (1 << side ** 3):
        raise ValueError('complete atomic spin-product basis coordinate required')
    permutation = site_permutation(side, tuple(displacement))
    return sum(1 << permutation[j] for j in range(side ** 3) if mask >> j & 1)


def purification(side=2):
    if side != 2:
        raise ValueError('explicit all-spin purification execution is specified at N=2')
    volume = side ** 3; complement = (1 << volume) - 1
    return {(mask, mask ^ complement): s.Integer((-1) ** mask.bit_count()) / 2 ** (volume // 2)
            for mask in range(1 << volume)}


def translate_state(state, displacement, *, side=2, which='source', reference_character=(0, 0, 0)):
    if which not in ('source', 'reference', 'joint'):
        raise ValueError('source, reference or joint translation required')
    phase = s.Integer(1)
    if which != 'source':
        phase = s.expand_complex(s.exp(2 * s.pi * s.I * sum(a * b for a, b in
            zip(reference_character, displacement, strict=True)) / side))
    result = defaultdict(lambda: s.Integer(0))
    for (source, reference), coefficient in state.items():
        if which != 'reference':
            source = shift_spins(source, side, displacement)
        if which != 'source':
            reference = shift_spins(reference, side, displacement)
        result[source, reference] += phase * coefficient
    return clean(result)


def phase_projector(state, momentum, *, side=2, which='source', reference_character=(0, 0, 0)):
    """Restricted complete finite translation PVM, with all source/ref states kept."""
    if len(momentum) != 3 or any(type(k) is not int or not 0 <= k < side for k in momentum):
        raise ValueError('three momentum indices required')
    result = defaultdict(lambda: s.Integer(0))
    for displacement in product(range(side), repeat=3):
        phase = s.expand_complex(s.exp(-2 * s.pi * s.I * sum(k * a for k, a in
            zip(momentum, displacement, strict=True)) / side)) / side ** 3
        translated = translate_state(state, displacement, side=side, which=which,
                                     reference_character=reference_character)
        for key, value in translated.items():
            result[key] += phase * value
    return clean(result)


def inner(left, right):
    return s.expand(sum(s.conjugate(value) * right.get(key, 0) for key, value in left.items()))


def norm2(state):
    return inner(state, state)


def full_phase_probability(side, phase):
    """All (2N)^3 source phase outcomes, including the parity-forbidden zeros."""
    if len(phase) != 3 or any(type(k) is not int or not 0 <= k < 2 * side for k in phase):
        raise ValueError('complete source phase labels lie in 0..2N-1')
    if any(k % 2 for k in phase):
        return s.Integer(0)
    return momentum_probability(side, tuple(k // 2 for k in phase))
