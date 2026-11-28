"""Exact four-site spectrum and outward bounds on visibility and spin retrieval."""

from dataclasses import dataclass
from functools import cache

import sympy as sp

from . import ring


X = sp.Symbol("x")
R = sp.Rational
ERROR = R(1, 10**18)


def matrix(operator):
    return sp.Matrix(operator.toarray())


@dataclass(frozen=True)
class Level:
    factor: sp.Poly
    lo: sp.Rational
    hi: sp.Rational
    multiplicity: int

    @property
    def mid(self):
        return (self.lo + self.hi) / 2


@cache
def sector(up, down):
    h = matrix(ring.hamiltonian(4, up, down))
    characteristic = sp.Poly(h.charpoly(X).as_expr(), X)
    factors = characteristic.factor_list()[1]
    minimal = sp.Poly(1, X)
    levels = []
    for factor, multiple in factors:
        minimal *= factor
        for (lo, hi), degree in factor.intervals(eps=R(1, 10**32)):
            assert degree == 1
            levels.append(Level(factor, lo, hi, multiple))
    levels.sort(key=lambda v: v.lo)
    assert sum(v.multiplicity for v in levels) == h.rows
    assert all(a.hi < b.lo for a, b in zip(levels, levels[1:]))
    return h, levels, minimal


def apply_polynomial(h, polynomial, v):
    result = sp.zeros(h.rows, v.cols)
    for coefficient in polynomial.all_coeffs():
        result = h * result + coefficient * v
    return result


@cache
def projector_polynomial(up, down, index):
    _, levels, minimal = sector(up, down)
    target = levels[index]
    other = minimal.exquo(target.factor)
    projector = (other * sp.invert(other, target.factor)).rem(minimal)
    # Exactly one on this irreducible factor, zero on all the other factors.
    assert projector.rem(target.factor) == sp.Poly(1, X, domain=projector.domain)
    others = [v for v in levels if v.factor == target.factor and v != target]
    denomin = sp.prod(target.mid - v.mid for v in others)
    polynomial = sp.Poly(sp.prod(X - v.mid for v in others) / denomin, X)
    for v in [target, *others]:
        radius = (v.hi - v.lo) / 2
        bounds = [abs(v.mid - z.mid) + radius for z in others]
        derivative = sum(
            sp.prod(bounds[:i] + bounds[i + 1 :]) for i in range(len(bounds))
        ) / abs(denomin)
        assert derivative * radius < ERROR / 20
    return (projector * polynomial).rem(minimal)


def projected(up, down, index, v):
    h = sector(up, down)[0]
    return apply_polynomial(h, projector_polynomial(up, down, index), v)


def interval(lo, hi, digits=10):
    assert lo <= hi
    scale = 10**digits
    return [
        str(sp.Rational(sp.floor(scale * lo), scale)),
        str(sp.Rational(sp.ceiling(scale * hi), scale)),
    ]


def interval_ratio(num, den):
    bounds = [
        a / b for a in (num - ERROR, num + ERROR) for b in (den - ERROR, den + ERROR)
    ]
    assert den > ERROR
    return interval(min(bounds), max(bounds))


def probability(amplitude, prep):
    lo = sp.Max(0, abs(amplitude) - ERROR)
    return interval(
        lo**2 / (4 * (prep + ERROR)),
        (abs(amplitude) + ERROR) ** 2 / (4 * (prep - ERROR)),
        12,
    )


def half_ring_translation(up, down):
    words, indices = ring.sector(4, up, down)
    translation = sp.zeros(len(words))
    for j, word in enumerate(words):
        first, a = ring.shifted(word, 4)
        second, b = ring.shifted(first, 4)
        translation[indices[second], j] = a * b
    assert translation**2 == sp.eye(len(words))
    return translation


def translation_on_factor(h, minimal, level, translation):
    quotient = minimal.exquo(level.factor)
    projector = (quotient * sp.invert(quotient, level.factor)).rem(minimal)
    whole = apply_polynomial(h, projector, sp.eye(h.rows))
    assert whole**2 == whole and whole.T == whole
    turned = translation * whole
    if turned == whole:
        return 1
    assert turned == -whole
    return -1


def site_occupations(counts, x):
    return sp.diag(
        *[((w >> (2 * x)) & 3).bit_count() for w in ring.sector(4, *counts)[0]]
    )


def charge_and_current(sign, counts, columns):
    """Certify what charge the two spin ports exchange and how it reaches site 2."""
    y = 2
    source_number = site_occupations((2, 2), y)
    for spin in (0, 1):
        step = (sign * (spin == 0), sign * (spin == 1))
        target_number = site_occupations((2 + step[0], 2 + step[1]), y)
        port = matrix(ring.field(4, 2, 2, y, spin, sign > 0))
        assert target_number * port - port * source_number == sign * port
        # A physical 2π onsite spin turn multiplies every singly occupied
        # component by -1 and leaves the empty and doubly occupied ones fixed.
        incoming_parity = sp.diag(*[(-1) ** n for n in source_number.diagonal()])
        outgoing_parity = sp.diag(*[(-1) ** n for n in target_number.diagonal()])
        assert outgoing_parity * port * incoming_parity == -port
    number = site_occupations(counts, y)
    h = sector(*counts)[0]
    currents = []
    for x in (1, 3):
        hop = sp.zeros(h.rows)
        for spin in (0, 1):
            mid = (counts[0] - (spin == 0), counts[1] - (spin == 1))
            hop += matrix(ring.field(4, *mid, y, spin, True)) * matrix(
                ring.field(4, *counts, x, spin, False)
            )
        currents.append(hop - hop.T)
    assert sum(currents, sp.zeros(h.rows)) == h * number - number * h
    first, second = columns
    norm = first.dot(first) + second.dot(second)
    # At the quarter beat the upper energy component has relative phase -i.
    quarter = [interval_ratio((first.T * j * second)[0] / 8, norm) for j in currents]
    assert quarter[0] == quarter[1] and all(sign * R(v) > 0 for v in quarter[0])
    return {
        "local_charge_step_per_accepted_port": sign,
        "two_pi_port_phase": -1,
        "quarter_beat_inward_bond_currents_given_branch": quarter,
        "source_bonds": [[1, y], [3, y]],
    }


def half_beat_clock(levels):
    # Machin's formula with alternating-series bounds for each arctangent.
    def atan_bounds(n):
        x = R(1, n)
        base = sum((-1) ** j * x ** (2 * j + 1) / (2 * j + 1) for j in range(14))
        return base, base + x**29 / 29

    a, b = atan_bounds(5), atan_bounds(239)
    pi_lo, pi_hi = 16 * a[0] - 4 * b[1], 16 * a[1] - 4 * b[0]
    low, high = (levels[2].lo - levels[0].hi) / 16, (levels[2].hi - levels[0].lo) / 16
    chosen = R(57316543, 10**6)
    phase_error = max(abs(chosen * low - pi_hi), abs(chosen * high - pi_lo))
    assert phase_error < R(26, 10**9)
    return {
        "half_beat_time": interval(pi_lo / high, pi_hi / low, 11),
        "rounded_wait": str(chosen),
        "relative_phase_error_upper": "26/1000000000",
    }


def resident_bell(up, down, wave, site):
    """Two actual spin columns, one-site occupancy gate, spectators traced out."""
    lowered = matrix(ring.spin_plus(4, up - 1, down + 1)).T * wave
    fibres = {}
    for ref, (counts, column) in enumerate(
        (((up, down), wave), ((up - 1, down + 1), lowered))
    ):
        for i, word in enumerate(ring.sector(4, *counts)[0]):
            local = (word >> (2 * site)) & 3
            if local in (1, 2):
                key = word & ~(3 << (2 * site))
                vector = fibres.setdefault(key, sp.zeros(4, 1))
                vector[2 * ref + (local == 2)] = column[i]
    gram = sum((v * v.T for v in fibres.values()), sp.zeros(4)) / 2
    x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    norm = sp.trace(gram)
    xx = sp.trace(gram * sp.kronecker_product(x, x))
    zz = sp.trace(gram * sp.kronecker_product(z, z))
    assert abs(xx - zz) < ERROR / 10
    assert norm > R(1, 50)
    return gram, norm, xx, zz


def original_start_bell_records(gram, sqrt_two):
    """Retain all four public outcomes and the setting-independent refusal."""
    eye = sp.eye(2)
    x, z = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    raw = gram / 2  # the injection port contributes its own factor of one half
    records = {}
    for name, axis in (("Z", z), ("X", x)):
        for side, pauli in (("Z+X", z + x), ("Z-X", z - x)):
            entries = {}
            exact_sum = 0
            for a in (1, -1):
                for b in (1, -1):
                    left = (eye + a * axis) / 2
                    values = [
                        sp.trace(
                            raw
                            * sp.kronecker_product(left, (eye + b * s * pauli / 2) / 2)
                        )
                        for s in (sqrt_two, sqrt_two + R(1, 10**15))
                    ]
                    exact_sum += values[0]
                    entries[f"{a:+},{b:+}"] = interval(
                        min(values) - ERROR, max(values) + ERROR, 12
                    )
            assert exact_sum == raw.trace()
            records[f"{name},{side}"] = entries
    return {
        "original_start_bins_given_setting": records,
        "all_refusals_per_original": interval(
            1 - raw.trace() - ERROR, 1 - raw.trace() + ERROR, 12
        ),
    }


def run():
    # Dimensions are tiny enough to keep every spin sector, including levels
    # with zero weight when a spinor is applied to the singlet.
    g_level = sector(2, 2)[1][0]
    assert g_level.multiplicity == 1 and g_level.factor not in {
        v.factor for v in sector(3, 1)[1]
    }
    initial = sp.zeros(36, 1)
    neel = sum(1 << (2 * x + x % 2) for x in range(4))
    initial[ring.sector(4, 2, 2)[1][neel]] = 1
    g = projected(2, 2, 0, initial)
    prep = g.dot(g)
    assert prep > R(319, 1000)
    banks = []
    for sign in (1, -1):
        up, down = (3, 2) if sign == 1 else (2, 1)
        h, levels, minimal = sector(up, down)
        h_three, high_spin, _ = sector(up + 1, down - 1)
        raising = matrix(ring.spin_plus(4, up, down))
        assert raising * h == h_three * raising
        assert h * raising.T == raising.T * h_three
        high_counts = {v.factor: v.multiplicity for v in high_spin}
        assert levels[0].multiplicity == levels[2].multiplicity == 2
        assert levels[1].multiplicity == high_counts[levels[1].factor] == 1
        assert all(levels[i].factor not in high_counts for i in (0, 2))
        translation = half_ring_translation(up, down)
        parity = [
            translation_on_factor(h, minimal, levels[i], translation) for i in (0, 2)
        ]
        assert parity[0] == -parity[1]
        probes = [
            matrix(ring.field(4, 2, 2, x, int(sign < 0), sign > 0)) * g
            for x in range(4)
        ]
        cols = [projected(up, down, i, probes[0]) for i in (0, 2)]
        weights = [v.dot(v) for v in cols]
        assert all(v > R(1, 20) for v in weights)
        ranks = [
            {
                "spin_z_half": v.multiplicity,
                "spin_half_copies": v.multiplicity - high_counts.get(v.factor, 0),
                "energy_relative_to_reservoir": interval(
                    (v.lo - g_level.hi) / 16 - sign * R(1, 2),
                    (v.hi - g_level.lo) / 16 - sign * R(1, 2),
                ),
                "polynomial_in_16H": [int(c) for c in v.factor.all_coeffs()],
            }
            for v in levels[:4]
        ]
        words = ring.sector(4, up, down)[0]
        reads = {}
        for label, phase in (("parked", 1), ("half_beat", -1)):
            wave = cols[0] + phase * cols[1]
            norm = wave.dot(wave)
            entries = []
            for y, sink in enumerate(probes):
                # Each local observable has norm at most one on the full sector.
                charge = spin = single = 0
                for i, word in enumerate(words):
                    n = (word >> (2 * y)) & 3
                    value = (n & 1) - (n >> 1)
                    charge += (n.bit_count() - 1) * wave[i] ** 2
                    spin += value * wave[i] ** 2
                    single += value**2 * wave[i] ** 2
                entries.append(
                    {
                        "ground_return_per_start": probability(sink.dot(wave), prep),
                        "signed_excess_given_branch": interval_ratio(charge, norm),
                        "twice_spin_given_branch": interval_ratio(spin, norm),
                        "single_occupation_given_branch": interval_ratio(single, norm),
                        "spin_retention_given_single": interval_ratio(spin, single),
                    }
                )
                if y == 2:
                    gram, bell_norm, xx, zz = resident_bell(up, down, wave, y)
                    # Normalize only after the same spin-independent local gate.
                    # The left/right Pauli axes Z,X and (Z±X)/sqrt(2) give this S.
                    c = sp.Rational(1414213562373095, 10**15)
                    assert c**2 < 2 < (c + R(1, 10**15)) ** 2
                    s_bounds = [
                        a * b / d
                        for a in (c, c + R(1, 10**15))
                        for b in (abs(xx + zz) - 2 * ERROR, abs(xx + zz) + 2 * ERROR)
                        for d in (bell_norm - ERROR, bell_norm + ERROR)
                    ]
                    entries[-1]["local_Bell_absolute_CHSH"] = interval(
                        min(s_bounds), max(s_bounds)
                    )
                    entries[-1]["branch_and_single_per_original"] = interval(
                        (bell_norm - ERROR) / 2, (bell_norm + ERROR) / 2
                    )
                    entries[-1]["bell_outcome_records"] = original_start_bell_records(
                        gram, c
                    )
            reads[label] = entries
        banks.append(
            {
                "sign": sign,
                "receiving_port_and_bond_currents": charge_and_current(
                    sign, (up, down), cols
                ),
                "half_beat_clock": half_beat_clock(levels),
                "lowest_four_charged_energies": ranks,
                "two_site_translation_on_the_two_bright_levels": parity,
                "per_original_branch_acceptance": interval(
                    (sum(weights) - ERROR) / 2, (sum(weights) + ERROR) / 2
                ),
                "local_read_intervals": reads,
            }
        )
    return {
        "ground_level_in_original_units": interval(g_level.lo / 16, g_level.hi / 16),
        "neel_preparation": interval(prep - ERROR, prep + ERROR),
        "projector_operator_error_upper": "1/20000000000000000000",
        "banks": banks,
    }
