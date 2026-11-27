"""Local recorded bath instrument on a native word and an arbitrary spectator."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
import math

import numpy as np

from ... import finite_rods, geometry
from ...archive import source_densities
from ...experiment import _frames
from ...memories import Pair, THETA
from ...native_pair import NativePair, TICKETS

# One representative for each native input word in the eight-ticket register.
ENCODE = (0, 5, 6, 7)
SWEEP = (0, 1, 2)


@lru_cache(maxsize=128)
def native_gate(pair, work):
    return pair.gate(THETA+Pair.from_weighted_work(work).angle_increment)


def _rotate(block, unitary):
    return np.einsum("ai,irjs,bj->arbs", unitary, block, unitary.conj(), optimize=True)


def bath_branch(joint, addresses, outcomes, work, pair=NativePair(0, "right")):
    """Unnormalized material/spectator output for one public displaced-bit record.

    The first factor is the four-word native memory. An unmeasured, independent
    fair bit replaces each measured outgoing bit. All public rows sum to a channel.
    """
    if len(addresses) != len(outcomes) or any(a not in SWEEP for a in addresses):
        raise ValueError("each addressed ticket bit needs one displaced-bit record")
    if any(z not in (0, 1) for z in outcomes):
        raise ValueError("outgoing records must be binary")
    joint = np.asarray(joint, dtype=complex)
    if joint.ndim != 2 or joint.shape[0] != joint.shape[1] or joint.shape[0] % 4 or not joint.size:
        raise ValueError("expected a four-word memory tensored with a finite spectator")
    spectator = joint.shape[0]//4
    q = range(spectator)
    gate = native_gate(pair, work)
    native = _rotate(joint.reshape(4, spectator, 4, spectator), gate.conj().T)
    ticket = np.zeros((8, spectator, 8, spectator), dtype=complex)
    ticket[np.ix_(ENCODE, q, ENCODE, q)] = native
    for a, z in zip(addresses, outcomes, strict=True):
        index = [[i for i in range(8) if ((i >> a) & 1) == bit] for bit in (0, 1)]
        old = ticket[np.ix_(index[z], q, index[z], q)]
        ticket.fill(0)
        for fresh in index:
            ticket[np.ix_(fresh, q, fresh, q)] = old/2
    # A reversible table lookup writes the native word and releases the ticket.
    # The output ticket and its fresh bits are never copied into the public record.
    native.fill(0)
    for t, word in enumerate(TICKETS):
        native[word, :, word, :] += ticket[t, :, t, :]
    return _rotate(native, gate).reshape(joint.shape)


def output_without_conditioning(joint, addresses, work, pair=NativePair(0, "right")):
    return sum(bath_branch(joint, addresses, bits, work, pair)
               for bits in product((0, 1), repeat=len(addresses)))


def ready_records(joint, work, pair=NativePair(0, "right")):
    """The material appointment, including the eight complete public record blocks."""
    return {bits: bath_branch(joint, SWEEP, bits, work, pair) for bits in product((0,1), repeat=3)}


@lru_cache(maxsize=32)
def stock_read(stock, shape, inverse_radius):
    """Use the existing sheet action and finite rod work on the archived source."""
    ordinary, collar = source_densities()[stock]
    ell = Fraction(1, inverse_radius)
    frame = _frames()[shape]
    density = Fraction(ordinary).limit_denominator(33) + Fraction(collar).limit_denominator(33)
    rod = finite_rods.read_ball(finite_rods.angle_gain(ell)*density, ell, frame)
    stress = geometry.project(geometry.energy_tensor(ordinary, collar), frame)
    work = 4*math.pi*float(rod.weighted_over_4pi)*stress
    return dict(stock=stock, shape=shape, inverse_radius=inverse_radius, stress=stress, work=work,
                proper_area=4*math.pi*float(rod.area_over_4pi), dwell_increment=Pair.from_weighted_work(work).angle_increment)


def replacement_horizon(levels):
    """Independent uniform addresses with replacement, cutoff 1/(32*4^(J-1))."""
    if not isinstance(levels, int) or levels < 1:
        raise ValueError("expected a positive integer number of enrolled levels")
    pairs = 7*levels
    # Choose the first s satisfying 3M*(2/3)^s <= ell^8/(M+1)^2, in integers.
    inverse_ell = 32*4**(levels-1)
    target_factor = 3*pairs*(pairs+1)**2*inverse_ell**8
    estimate = math.ceil(math.log(target_factor)/math.log(3/2))
    while 3**estimate < target_factor*2**estimate:
        estimate += 1
    while estimate and 3**(estimate-1) >= target_factor*2**(estimate-1):
        estimate -= 1
    return estimate
