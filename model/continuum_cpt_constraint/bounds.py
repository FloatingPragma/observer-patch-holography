"""Exact moment and symmetry bounds; no numerical continuum-state oracle."""
from fractions import Fraction as F


def rational(value, *, positive=False):
    if type(value) not in (int, F) or value < 0 or (positive and not value):
        raise ValueError('nonnegative exact rational required')
    return F(value)


def gaussian(v, coupling, width):
    v = rational(v, positive=True)
    coupling = rational(coupling)
    width = rational(width, positive=True)
    mean = coupling / 2
    second = 3 * v * v / (2 * width * width) + coupling * coupling / 2
    return {
        'v': v, 'U0': coupling, 'width': width,
        'signed_mean': mean, 'second_moment_upper': second,
        'TV_from_reflection_lower': mean * mean / second,
        'W1_from_reflection_lower': 2 * mean,
        'necessary_charge_clock_subtraction': mean,
        'clipping_radius': second / (2 * mean) if mean else None,
        'clipped_mean_lower': mean / 2,
    }


def measure(pairs):
    result = {}
    for energy, weight in pairs:
        if type(energy) not in (int, F):
            raise ValueError('exact rational energy required')
        weight = rational(weight)
        energy = F(energy)
        if weight:
            result[energy] = result.get(energy, F(0)) + weight
    if sum(result.values()) != 1:
        raise ValueError('probability weights must sum to one')
    return result


def moments(law):
    law = measure(law.items())
    return tuple(sum(weight * energy ** n for energy, weight in law.items())
                 for n in (0, 1, 2))


def reflection(law):
    return measure((-energy, weight) for energy, weight in law.items())


def total_variation(left, right):
    left = measure(left.items()); right = measure(right.items())
    return sum(abs(left.get(x, F(0)) - right.get(x, F(0)))
               for x in set(left) | set(right)) / 2


def wasserstein_one(left, right):
    """Exact one-dimensional CDF integral, independently of the mean bound."""
    left = measure(left.items()); right = measure(right.items())
    nodes = sorted(set(left) | set(right))
    cumulative = F(0); result = F(0)
    for lo, hi in zip(nodes, nodes[1:]):
        cumulative += left.get(lo, F(0)) - right.get(lo, F(0))
        result += abs(cumulative) * (hi - lo)
    return result


def clipped_mean(law, radius):
    radius = rational(radius, positive=True)
    return sum(weight * max(-radius, min(radius, energy))
               for energy, weight in measure(law.items()).items())


def shifted_pair(momentum_energy, scalar):
    """A complete two-band response law; not asserted to be the OPH limit."""
    p = rational(momentum_energy, positive=True)
    scalar = rational(scalar)
    return measure(((scalar + p, F(1, 2)), (scalar - p, F(1, 2))))


def certificate():
    rows = [gaussian(F(1, 8), 1, width) for width in (1, 2, 4, 8)]
    return {
        'component_count': 16,
        'mean_identity': 'integral E dnu = U0/2',
        'second_moment_bound': 'integral E^2 dnu <= v^2 ||gradient phi||^2 + U0^2/2',
        'TV_convention': 'sup_B |nu(B)-rho(B)| = one half total variation norm',
        'TV_bound': 'TV(nu,reflected nu) >= (U0/2)^2/C',
        'W1_bound': 'W1(nu,reflected nu) >= U0',
        'necessary_charge_clock_centering': 'kappa=U0/2',
        'gaussian_rows': [{k: str(v) if v is not None else None for k,v in row.items()}
                          for row in rows],
        'broad_gaussian_TV_lower_limit': '1/2 for U0>0',
        'same_field_same_clock_CPT_vacuum_excluded_when_U0_positive': True,
        'Lorentz_covariance_alone_excluded': False,
        'centered_dynamics_ground_positivity_proved': False,
        'actual_nonzero_coupling_continuum_measure_computed': False,
        'physical_mass_or_species_identified': False,
        'whole_action_or_all_field_dictionaries_excluded': False,
    }
