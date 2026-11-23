"""Integer enclosures for finite gateway preparation and unlimited free flight."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from math import ceil, factorial

import numpy as np
import sympy as sp
from scipy import sparse

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import gateway_table


ROUNDS = 1024
ATTEMPTS = 48
SCALE = 1 << 80
SCORE = 10**12


def chart_key(gateways):
    tips = tuple(gateways)
    if (len(tips) != 3 or any(type(v) is not int or not 0 <= v < 12 for v in tips)
            or tuple(sorted(set(tips))) != tips):
        raise ValueError("Expected three distinct, ordered local gateway ports")
    return tips


@lru_cache
def dark_projector(gateways):
    """Largest invariant subspace with zero amplitude at each outward edge."""
    tips = chart_key(gateways)
    adjacency = sp.Matrix(cell().adjacency().toarray().astype(int).tolist())
    eye = sp.eye(12)
    read = eye[list(tips), :]
    if (adjacency-5*eye)*(adjacency+eye)*(adjacency**2-5*eye) != sp.zeros(12):
        raise ArithmeticError("Cell spectrum exceeds the four boundary moments")
    constraints = sp.Matrix.vstack(*(read*adjacency**k for k in range(4)))
    vectors = constraints.nullspace()
    if not vectors:
        raise ArithmeticError("No confined material profile in this chart")
    basis = sp.Matrix.hstack(*vectors)
    proj = basis*(basis.T*basis).inv()*basis.T
    if (proj != proj.T or proj**2 != proj or read*proj != sp.zeros(3, 12)
            or adjacency*proj != proj*adjacency):
        raise ArithmeticError("Invalid confined-space projector")
    return sp.ImmutableMatrix(proj)


def board_charts():
    """Check the native single-load reduction and every home-to-exterior column."""
    local = cell()
    expected = set(local.seams)
    reports = {}
    for level in range(3):
        board = federation(level)
        lap = board.laplacian()
        kinetic = candidate.contrast.generator(board, 0)
        unit = candidate.contrast.generator(board, 1)
        identity = sparse.eye(board.n, format="csr")
        if ((2*kinetic != lap).nnz or (unit-kinetic != 6*identity).nnz
                or (candidate.generator(board) != kinetic+6*candidate.KAPPA*identity).nnz):
            raise ArithmeticError("The registered free generator has a different material action")
        inside = [set() for _ in range(board.cells)]
        tips = gateway_table(board)
        for a, b in board.seams:
            if a//12 == b//12:
                inside[a//12].add((a % 12, b % 12))
        if any(edges != expected for edges in inside):
            raise ArithmeticError("A local exchange or its endpoint is missing")
        counts = Counter(chart_key(row) for row in tips)
        # Every exterior row touches a gateway. P vanishes there, so the
        # embedded range reduces the full Laplacian, including its degree term.
        reports[board.cells] = counts
    return reports


def nearest(numerator, denominator):
    """Round real dyadics to the nearest grid point, with ties rounded up."""
    return (2*numerator + denominator)//(2*denominator)


@lru_cache
def internal_pulse(pulse):
    """Enclose exp(i pulse A/2), omitting only a common phase."""
    if pulse not in (Fraction(1, 2), Fraction(1), Fraction(2)):
        raise ValueError("Pulse outside the three registered preparations")
    adjacency = cell().adjacency().toarray().astype(int).astype(object)
    size = len(adjacency)
    degree = 64
    p, q = pulse.numerator, 2*pulse.denominator
    identity = np.eye(size, dtype=int).astype(object)
    power = identity
    real, imag = identity, 0*identity
    denominator = 1
    for k in range(1, degree+1):
        power = (power @ adjacency)*p
        real, imag = real*(q*k), imag*(q*k)
        denominator *= q*k
        if k % 2:
            imag += (-1)**(k//2)*power
        else:
            real += (-1)**(k//2)*power
    re = nearest(real*SCALE, denominator)
    im = nearest(imag*SCALE, denominator)
    radius = Fraction(5*p, q)
    tail = 3**ceil(radius)*radius**(degree+1)/factorial(degree+1)
    # Each complex entry has rounding error <= sqrt(2)/(2*SCALE).
    units = ceil(tail*SCALE) + size
    return re, im, units


def success_intervals(gateways, pulse):
    """Bound the success probability of each sharp launch using integer matrices."""
    tips = chart_key(gateways)
    re, im, error = internal_pulse(Fraction(pulse))
    re, im = re.copy(), im.copy()
    re[list(tips)] = 0
    im[list(tips)] = 0
    if ROUNDS != 2**(ROUNDS.bit_length()-1):
        raise ValueError("The registered enclosure uses a power-of-two round budget")
    for _ in range(ROUNDS.bit_length()-1):
        real = nearest(re @ re - im @ im, SCALE)
        imag = nearest(re @ im + im @ re, SCALE)
        re, im = real, imag
        # The exact power is a contraction. Squaring accumulates at most
        # 2 delta + delta^2 before rounding, bounded by 3 delta for delta <= 1.
        if error >= SCALE:
            raise ArithmeticError("Power enclosure lost contraction accuracy")
        error = 3*error + len(re)
    delta = Fraction(error, SCALE)
    probability_error = delta*(2+delta)
    answers = []
    for total in (re*re+im*im).sum(axis=0):
        center = Fraction(int(total), SCALE**2)
        low = max(0, (center-probability_error)*SCORE//1)
        high = min(SCORE, ceil((center+probability_error)*SCORE))
        if not 0 <= low <= high <= SCORE:
            raise ArithmeticError("Invalid probability enclosure")
        answers.append((Fraction(low, SCORE), Fraction(high, SCORE)))
    return answers


def certify(gateways, pulse=Fraction(1), *, start=0, direction=1):
    """Joint probability that the source accepts and the terminal home reader fires."""
    tips = chart_key(gateways)
    if type(start) is not int or not 0 <= start < 12 or type(direction) is not int or direction not in (-1, 1):
        raise ValueError("Expected a local starting port and a forward or reverse sweep")
    proj = dark_projector(tips)
    ranges = success_intervals(tips, pulse)
    reach_low = reach_high = Fraction(1)
    protected = Fraction(0)
    for attempt in range(ATTEMPTS):
        j = (start+direction*attempt) % 12
        weight = Fraction(proj[j, j])
        low, high = ranges[j]
        if not 0 <= weight <= high:
            raise ArithmeticError("The invariant weight exceeds the acceptance bound")
        protected += reach_low*weight
        reach_low *= 1-high
        reach_high *= 1-max(weight, low)
    if not 0 <= protected <= 1-reach_low:
        raise ArithmeticError("Invalid unconditional recovery bound")
    return {
        "dimension": int(proj.trace()),
        "joint_home_floor_ppt": int(protected*SCORE),
        "exhaustion_ceiling_ppt": ceil(reach_high*SCORE),
        "one_attempt_success_ppt": [(int(lo*SCORE), int(hi*SCORE)) for lo, hi in ranges],
    }


def least_protected_sweep(weights, upper):
    """Minimize the certified contribution over permutations of one complete sweep."""
    if not weights or len(weights) != len(upper) or any(
            not 0 <= q <= s <= 1 for q, s in zip(weights, upper, strict=True)):
        raise ValueError("Expected confined weights and enclosing acceptance probabilities")
    # Swapping adjacent i,j changes the contribution by reach*(q_i*s_j - q_j*s_i).
    # A port with zero upper acceptance also has zero protected weight.
    order = sorted(range(len(weights)), key=lambda i: weights[i]/upper[i] if upper[i] else 0)
    reach, protected = Fraction(1), Fraction(0)
    for i in order:
        protected += reach*weights[i]
        reach *= 1-upper[i]
    return protected, reach


def arbitrary_sweep_floor(gateways):
    """Unit-pulse floor for four sweeps, each in any fixed permutation of the ports."""
    tips = chart_key(gateways)
    weights = tuple(map(Fraction, dark_projector(tips).diagonal()))
    upper = tuple(hi for _, hi in success_intervals(tips, Fraction(1)))
    one, exhausted = least_protected_sweep(weights, upper)
    sweeps, remainder = divmod(ATTEMPTS, len(weights))
    if remainder:
        raise ArithmeticError("Incomplete source sweep")
    return int(SCORE * one * sum(exhausted**i for i in range(sweeps)))


if __name__ == "__main__":
    boards = board_charts()
    floors = {tips: arbitrary_sweep_floor(tips) for tips in set().union(*boards.values())}
    for size, charts in boards.items():
        print(f"{size} cells, every complete port order: joint home probability >= {min(floors[t] for t in charts)}/{SCORE}")
