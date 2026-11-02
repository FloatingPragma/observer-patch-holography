"""Actual full-domain repair columns and exact twist/translation arithmetic."""
from itertools import combinations

from model.defect_field_bridge.full_domain import branch_action, RAW_SLOTS
from model.strong_coupling_refinement.algebra import geometry


def spin_mask(word, sites):
    return sum(1 << (2*x+((word >> x) & 1)) for x in range(sites))


def protected_columns():
    sites, edges = geometry((2, 2, 2))
    checked = 0
    for a, b, _ in edges:
        for law in ('C', 'E'):
            for coin in (0, 1):
                for word in range(1 << sites):
                    mask = spin_mask(word, sites)
                    if branch_action(mask, law, coin, sites=sites, flavors=2, seam=(a, b)) != (0, mask, 1, 0):
                        raise ArithmeticError('field repair changed the singly occupied spin carrier')
                    checked += 1
    # Independently verify completeness and injectivity on every sector of a smaller full carrier.
    full_count = 0
    for law in ('C', 'E'):
        for seam in ((0, 1), (0, 2), (1, 2)):
            buckets = {(coin, raw): [] for coin in (0, 1) for raw in range(RAW_SLOTS)}
            for coin in (0, 1):
                for mask in range(64):
                    raw, output, sign, _ = branch_action(mask, law, coin, sites=3, flavors=2, seam=seam)
                    if sign not in (-1, 1) or output.bit_count() != mask.bit_count():
                        raise ArithmeticError('full-sector repair is not number preserving')
                    buckets[coin, raw].append(output)
                    full_count += 1
            if any(len(row) != len(set(row)) for row in buckets.values()):
                raise ArithmeticError('full-sector recorded repair is not injective')
    return {'protected_cube_columns': checked, 'small_full_Fock_columns': full_count,
            'all_32_coin_raw_slots_retained': True, 'ideal_raw_record': 0,
            'identity_in_every_internal_spin_state': True,
            'bare_one_bit_occupation_law_reinterpreted': False,
            'large_carrier_native_pulse_execution_claimed': False}


def twist_translation(length, area):
    """Check the twist sign in the complete zero-total-Sz bit basis."""
    if (type(length) is not int or type(area) is not int or length < 2
            or length % 2 or area < 1 or area % 2 != 1 or length*area > 16):
        raise ValueError('small even length and odd cross-section required')
    sites = length*area
    count = 0
    for occupied in combinations(range(sites), sites//2):
        mask = sum(1 << i for i in occupied)
        shifted = sum(1 << (((i//area+1) % length)*area+i % area) for i in occupied)
        def exponent(word):
            return sum(x*sum(1-2*((word >> (x*area+j)) & 1) for j in range(area)) for x in range(length))
        difference = exponent(shifted)-exponent(mask)
        if difference % (2*length) != length:
            raise ArithmeticError('translation did not reverse the twist eigenphase')
        count += 1
    return {'length': length, 'odd_cross_section': area, 'zero_Sz_columns': count,
            'twist_translation_factor': '-1'}


def certificate():
    return {'repair': protected_columns(),
            'exact_twist_checks': [twist_translation(L, A) for L, A in ((4, 1), (6, 1), (2, 3), (4, 3))],
            'pi_squared_rational_bound': '(22/7)^2 < 10',
            'spin_gap_bound': 'pi^2*J*A/L; zero if the ground is degenerate',
            'source_gap_bound': 'pi^2*J*A/L+6*tau^3/U^2',
            'uniform_positive_nondegenerate_ground_gap_on_trajectory': False,
            'finite_gap_collapse_proves_massless_particle': False}
