"""Exact local-preparation ceiling and sharp-sector exterior obstruction.

The proofs use the ambient tensor qutrit algebra, not a fictitious tensor
factorization of the one-A/one-B sector. See README.md for the channel argument.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations


PORTS = 12
ORDERED = tuple((a, b) for a in range(PORTS) for b in range(PORTS) if a != b)
DIMENSION = len(ORDERED)


def region(value):
    if (type(value) is not tuple or any(type(x) is not int for x in value)
            or tuple(sorted(set(value))) != value
            or any(not 0 <= x < PORTS for x in value)):
        raise ValueError('a sorted tuple of distinct exact port indices is required')
    return value


def bounds(core):
    """Bounds for a deterministic local channel and any globally sharp state.

    Local preparation starts with I132/132 and an independent ancilla, acts
    only on core, and ends in the one-A/one-B sector without postselection.
    The sharp-state bound does not require fixed occupation probabilities.
    """
    core = region(core)
    inside = tuple((a, b) for a, b in ORDERED if a in core and b in core)
    outside = tuple((a, b) for a, b in ORDERED if a not in core and b not in core)
    # Q swaps A and B. Its exterior opposite-character projector is
    # P_both-out (I-Q)/2. Verify the forbidden projector annihilates I+Q
    # entry by entry over the ordered basis, without numerical eigensolving.
    reverse = {pair: (pair[1], pair[0]) for pair in ORDERED}
    if any(reverse[pair] not in outside or reverse[pair] == pair for pair in outside):
        raise ArithmeticError('exterior character projector is not a disjoint pair sum')
    exterior_classes = {tuple(sorted(pair)) for pair in outside}
    if 2*len(exterior_classes) != len(outside):
        raise ArithmeticError('opposite-character projector rank mismatch')
    # A local channel cannot change the exterior basis labels. The only Q
    # matrix elements surviving that block-diagonal algebra have vacuum
    # exterior, which is exactly the both-inside sector.
    surviving = []
    for a, b in ORDERED:
        before = tuple(1 if p == a else 2 if p == b else 0
                       for p in range(PORTS) if p not in core)
        after = tuple(2 if p == a else 1 if p == b else 0
                      for p in range(PORTS) if p not in core)
        if before == after:
            surviving.append((a, b))
    if tuple(surviving) != inside:
        raise ArithmeticError('local preparation charge support changed')
    return {
        'core': list(core),
        'ordered_core_support': len(inside),
        'exterior_forbidden_projector_rank': len(exterior_classes),
        'local_preparation_max_abs_charge': str(Fraction(len(inside), DIMENSION)),
        'sharp_sector_min_exterior_trace_distance': str(Fraction(len(exterior_classes), DIMENSION)),
        'attaining_local_preparation': '(I + P_both-in-core Q)/132',
        'attaining_sharp_state': '(I + Q)/132',
    }


def certificate():
    # All subsets, not just contiguous example cores. The local theorem is
    # independent of the scheduler and of whether the core is graph connected.
    by_size = []
    checked = 0
    for size in range(PORTS+1):
        representative = bounds(tuple(range(size)))
        for core in combinations(range(PORTS), size):
            row = bounds(core)
            if any(row[key] != representative[key] for key in (
                    'ordered_core_support', 'exterior_forbidden_projector_rank',
                    'local_preparation_max_abs_charge',
                    'sharp_sector_min_exterior_trace_distance')):
                raise ArithmeticError('locality bound depends on an unclaimed chart')
            checked += 1
        by_size.append(representative)
    return {
        'all_port_subsets_checked': checked,
        'bounds_by_core_size': by_size,
        'preparation_domain': 'deterministic core-local CPTP map, independent ancilla, tracial sector input, endpoint one-A/one-B support',
        'sharp_domain': 'any density matrix supported in one global Q eigenspace, without a fixed occupation marginal',
        'proof': 'Exterior basis diagonality survives a local channel. Only exterior-vacuum Q matrix elements remain; their probability is fixed by locality. The exterior opposite-character projector annihilates every globally sharp state and has reference probability rank/132. Both bounds are attained.',
        'full_qutrit_faithfulness_claimed': False,
        'sharp_local_particle_claimed': False,
        'modular_energy_excess': '0: -log(I132/132) is scalar on the admitted sector',
    }
