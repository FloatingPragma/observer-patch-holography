"""Infer the native droplet dispersion from an equal-load pair-binding inbox only."""
import math

import numpy as np

from . import instrument


def clock_radius(apparatus):
    ticks = np.asarray(apparatus['clock_ticks'])
    return math.asin(math.sqrt(2)/apparatus['shots'])*ticks.sum()/(ticks@ticks)


def infer(clock, binding, size, apparatus, numbers):
    rate = instrument.clock_rate(clock,apparatus)
    reading = instrument.read_line(binding,rate,apparatus)
    if reading is None:
        raise ValueError('unresolved binding line')
    clock_error = clock_radius(apparatus)
    center,radius = reading
    if center <= radius or rate <= clock_error:
        raise ValueError('unresolved binding or clock scale')
    low = (center-radius)*rate/(rate+clock_error)
    high = (center+radius)*rate/(rate-clock_error)

    def invert(binding_energy):
        return (binding_energy+2+math.sqrt(binding_energy*(binding_energy+4)))/2

    # A positive trial vector bounds the finite K=0 wrap error by Delta**(3-size).
    # Bootstrap from the universal error <= 1, without feeding the law's coupling to the reader.
    first = invert(max(0,low-1))
    finite_bound = first**(3-size)
    if low <= finite_bound:
        raise ValueError('binding does not resolve an attractive infinite-chain branch')
    deltas = [invert(low-finite_bound),invert(high)]
    predictions = {}
    for n in numbers:
        def inertia(delta):
            previous,current = 0.,1.
            for _ in range(1,n):
                previous,current = current,2*delta*current-previous
            return current
        predictions[str(n)] = [inertia(d) for d in deltas]
    return dict(binding_in_exchange_units=[low,high],finite_ring_bound=finite_bound,
                anisotropy=deltas,inertia=predictions)


def native_transport_bound(chemical, size, n, tick, rate, apparatus, shots):
    """Finite-ring TV enclosure from the chemistry inbox alone (see contrast_memory_spherical_v1/PROOF.md)."""
    if not (2 <= n and 2*n < size and tick >= 0 and shots > 0):
        raise ValueError('flight outside the finite-ring bound domain')
    low,high = chemical['anisotropy']
    if not (1 < low <= high and math.isfinite(high)):
        raise ValueError('invalid attractive branch interval')
    if low <= 2:
        return 1.
    rho = 2/(low+math.sqrt(low*low-1))
    tail = rho**(size-n)
    if (n-1)*tail >= 1:
        return 1.
    residual = math.sqrt(math.comb(size,n)/size)*n*(high+rho)*tail/(1-(n-1)*tail)
    if residual >= low-2:
        return 1.
    hidden_weight = 1/(1+(low-2)**2)
    rlo,rhi = chemical['inertia'][str(n)]
    center = rate/((rlo+rhi)/2)
    dr = clock_radius(apparatus)
    drift = max(center-(rate-dr)/rhi, (rate+dr)/rlo-center)
    phase_error = 2*tick*((rate+dr)*residual+drift)
    return min(1.,math.sqrt(4*hidden_weight+phase_error**2)+(size+1)/(2*shots))
