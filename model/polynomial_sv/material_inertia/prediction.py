"""Infer the native droplet dispersion from an equal-load pair-binding inbox only."""
import math

import numpy as np

from . import instrument


def infer(clock, binding, size, apparatus, numbers):
    rate = instrument.clock_rate(clock,apparatus)
    reading = instrument.read_line(binding,rate,apparatus)
    if reading is None:
        raise ValueError('unresolved binding line')
    ticks = np.asarray(apparatus['clock_ticks'])
    angle_error = math.asin(math.sqrt(2)/apparatus['shots'])
    clock_error = angle_error*ticks.sum()/(ticks@ticks)
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
