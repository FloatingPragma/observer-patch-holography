"""Exact rational identities for a projected positive generator."""
from fractions import Fraction as F


def spectral_memory(energies, weights, eta):
    if not energies or len(energies) != len(weights):
        raise ValueError('a complete nonempty spectral measure is required')
    energies, weights, eta = list(map(F, energies)), list(map(F, weights)), F(eta)
    if eta <= 0 or min(energies) < 0 or min(weights) < 0 or sum(weights) != 1:
        raise ValueError('positive generator, normalized measure and eta required')
    mean = sum(e*w for e, w in zip(energies, weights))
    variance = sum(w*(e-mean)**2 for e, w in zip(energies, weights))
    response = sum(w/(e+eta) for e, w in zip(energies, weights))
    memory = mean+eta-1/response
    if not 0 <= memory <= variance/eta:
        raise ArithmeticError('positive complement memory bound failed')
    return {'mean': mean, 'variance': variance, 'resolvent': response, 'memory': memory}


def memory_interval(mean, response, eta):
    lo, hi = map(F, mean)
    gl, gh = map(F, response)
    eta = F(eta)
    if eta <= 0 or lo > hi or not 0 < gl <= gh:
        raise ValueError('ordered mean and strictly positive response intervals required')
    return max(F(), lo+eta-1/gl), max(F(), hi+eta-1/gh)


def readout_memory_error(eta, energy_upper, mean_error, score_error):
    """Exact bound for normalized complete energy scores, without postselection."""
    eta, energy_upper, mean_error, score_error = map(F, (eta, energy_upper, mean_error, score_error))
    if eta <= 0 or min(energy_upper, mean_error, score_error) < 0:
        raise ValueError('positive eta and nonnegative bounds required')
    response_floor = 1/(energy_upper+eta)
    response_error = score_error/eta
    if response_error >= response_floor:
        raise ValueError('score accuracy does not certify a positive resolvent')
    return mean_error+response_error/(response_floor*(response_floor-response_error))


def exact_controls():
    # The one-pole control has zero complement memory. A two-pole response
    # with the same mean has a nonzero memory; neither is a source assumption.
    one = spectral_memory([2], [1], F(1, 4))
    two = spectral_memory([1, 3], [F(1, 2)]*2, F(1, 4))
    assert one['mean'] == two['mean'] == 2
    assert one['memory'] == 0 and two['memory'] == F(4, 9)
    return {'single_pole': one, 'two_pole': two,
            'equal_means_do_not_determine_dynamics': True}
