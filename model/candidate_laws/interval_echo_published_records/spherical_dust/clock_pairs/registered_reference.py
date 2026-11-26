"""A fixed pair register for physical endpoint signs and conditional blank completions."""
from fractions import Fraction as F
from itertools import product
import math

from .channel import ClockPair, LagGrid, VACUUM_CORRELATION


PAIRS = tuple(product((-1, 1), repeat=2))
ALL = 3


def joint(correlation):
    """Distribution in (--, -+, +-, ++) order."""
    correlation = F(correlation)
    if not -1 <= correlation <= 1:
        raise ValueError("correlation must be in [-1,1]")
    return tuple((1+correlation*a*b)/4 for a, b in PAIRS)


BLANK = joint(VACUUM_CORRELATION)


def censor(pair, mask):
    """Mask bit 0 keeps the A sign; bit 1 keeps B."""
    if type(mask) is not int or not 0 <= mask <= ALL or pair not in PAIRS:
        raise ValueError("two physical signs and a two-bit mask are required")
    return tuple(value if mask & 1 << side else None for side, value in enumerate(pair))


def completion_kernel(raw):
    """The exact blank conditional on the physically returned signs, including an empty mask."""
    if len(raw) != 2 or any(bit not in (-1, 1, None) for bit in raw):
        raise ValueError("expected two signed or missing endpoint records")
    allowed = tuple(all(bit is None or bit == value for bit, value in zip(raw, pair, strict=True))
                    for pair in PAIRS)
    mass = sum(p for keep, p in zip(allowed, BLANK, strict=True) if keep)
    return tuple(p/mass if keep else F(0) for keep, p in zip(allowed, BLANK, strict=True))


def complete(raw, rng):
    """Draw the missing signs with exact integer tickets and retain every physical return."""
    kernel = completion_kernel(raw)
    tickets = math.lcm(*(p.denominator for p in kernel))
    draw = rng.randrange(tickets)
    for pair, p in zip(PAIRS, kernel, strict=True):
        draw -= int(tickets*p)
        if draw < 0:
            return pair
    raise ArithmeticError("incomplete blank conditional")


def restrict_record(pair, observed_mask, target_mask, rng):
    """Restrict a completed record only along an inclusion of actual observed endpoints."""
    if type(observed_mask) is not int or not 0 <= observed_mask <= ALL:
        raise ValueError("invalid original mask")
    raw = censor(pair, target_mask)
    if target_mask & observed_mask != target_mask:
        raise ValueError("requested endpoint was never observed")
    return complete(raw, rng)


def projected_law(law, mask):
    """Push a joint law through the raw aperture and conditional completion."""
    if len(law) != len(PAIRS) or any(p < 0 for p in law) or sum(law) != 1:
        raise ValueError("expected a normalized four-sign distribution")
    rows = (completion_kernel(censor(pair, mask)) for pair in PAIRS)
    return tuple(sum(p*k for p, k in zip(law, col, strict=True)) for col in zip(*rows, strict=True))


def entropy(law):
    return -math.fsum(float(p)*math.log(float(p)) for p in law if p)


def divergence(law, reference):
    if any(p and not q for p, q in zip(law, reference, strict=True)):
        return math.inf
    return math.fsum(float(p)*math.log1p(float((p-q)/q))
                     for p, q in zip(law, reference, strict=True) if p)


def verify():
    """Check the exact kernels on every mask and basis input, then independently check both stocks."""
    checked = 0
    for basis in range(len(PAIRS)):
        point = tuple(F(k == basis) for k in range(len(PAIRS)))
        for mask in range(4):
            if projected_law(BLANK, mask) != BLANK:
                raise ArithmeticError("a blank record depends on its aperture")
            for subset in range(4):
                if subset & mask == subset:
                    if projected_law(projected_law(point, mask), subset) != projected_law(point, subset):
                        raise ArithmeticError("composed clock-register restrictions disagree")
                    checked += 1

    readings = {}
    grid = LagGrid(bits=5, max_lag=24)
    for stock, mass in (("ordinary", F(1, 3)), ("held", F(19, 33))):
        source = ClockPair(grid, F(2, 5)*mass/1024)
        delta = source.shift(F(0), F(1, 32), F(1, 8))
        law = joint(VACUUM_CORRELATION+delta)
        if not delta < 0:
            raise ArithmeticError("the two-clock test source has no material signal")
        for mask in range(4):
            retained = projected_law(law, mask)
            if retained != (law if mask == ALL else BLANK):
                raise ArithmeticError("conditional completion changed a physical sign marginal")
            d = divergence(retained, BLANK)
            s = entropy(retained)-entropy(BLANK)
            slope = -math.atanh(float(VACUUM_CORRELATION))*float(delta) if mask == ALL else 0
            if not math.isclose(s, slope-d, abs_tol=2e-15):
                raise ArithmeticError("the measured entropy differs from the exact pair identity")
        readings[stock] = {"correlation_shift": float(delta),
                           "entropy_change": entropy(law)-entropy(BLANK),
                           "divergence": divergence(law, BLANK)}
    return {"basis_and_inclusion_checks": checked, "stocks": readings}


if __name__ == "__main__":
    print(verify())
