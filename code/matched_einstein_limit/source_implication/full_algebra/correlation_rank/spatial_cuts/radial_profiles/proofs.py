"""Finite balance and sign tests for radial momentum-transfer meters."""

import json
from fractions import Fraction
from math import ceil, floor

import numpy as np
from mpmath import iv
from scipy.optimize import linprog

from .. import measure as cuts
from . import pencil


def normalize(rows):
    sizes = np.linalg.norm(cuts.mid(rows), axis=1)
    if np.any(sizes == 0):
        raise ValueError("zero material row")
    exponents = [-int(np.floor(np.log2(n))) for n in sizes]
    return rows * np.array([iv.mpf(2) ** e for e in exponents])[:, None], exponents


def projective(regional, shells, area, volume):
    """Allow arbitrary common signed probe factors on both affine ratio charts."""
    free, _ = normalize(np.vstack((regional, volume)))
    a, ae = normalize(area[None, :])
    # The shell-zero price is redundant on the kernel of the two boundary
    # readings. Keep all five other distances even for arbitrarily large prices.
    tiny = shells[:, 1:]
    sizes = np.linalg.norm(cuts.mid(tiny), axis=2).max(axis=0)
    factors = np.array([iv.mpf(2) ** (-int(np.floor(np.log2(v)))) for v in sizes])
    tiny = tiny * factors[None, :, None]
    result = []
    for first in (0, 1):
        base = np.vstack((free, tiny[first], a))
        slope = np.full(base.shape, iv.mpf(0), dtype=object)
        slope[len(free) : -1] = tiny[1 - first]
        result.append(pencil.audit(base, slope, Fraction(2) ** ae[0]))
    return result


def independent(regional, shells):
    rows = np.vstack((regional, shells[0, 1:], shells[1, 1:]))
    b, _ = normalize(rows)
    # Norm checks concern the rescaled basis. Invertibility carries over by
    # diagonality, even when far-shell responses are far below shot precision.
    right = cuts.fixed(np.linalg.inv(cuts.mid(b)))
    err = cuts.norm(np.eye(18, dtype=int) - b @ right)
    if err.b >= 1:
        raise ValueError("separate radial profile rank is uncertified")
    return {
        "rank": 18,
        "scaled_inverse_defect_ppm": ceil(cuts.survey.dyadic(err._mpi_[1]) * 10**6),
    }


def signed_cone(regional, shells, area, volume, signs):
    free, _ = normalize(np.vstack((regional, volume)))
    # At fixed total flux, each inner cumulative outward transfer equals
    # minus the sum of the remaining outer shell changes, as an exact identity.
    complements = [
        -q * np.flip(np.cumsum(np.flip(s[1:], axis=0), axis=0), axis=0)
        for q, s in zip(signs, shells)
    ]
    obs, scaling = normalize(np.vstack((*complements, area)))
    centers, f = cuts.mid(obs), cuts.mid(free)
    instance = linprog(
        np.r_[np.zeros(18), -1.0],
        A_ub=np.column_stack((-centers, np.ones(len(obs)))),
        b_ub=np.zeros(len(obs)),
        A_eq=np.column_stack((f, np.zeros(len(f)))),
        b_eq=np.zeros(len(f)),
        bounds=[(-1, 1)] * 18 + [(0, None)],
    )
    if not instance.success or instance.x[-1] <= 1e-10:
        return None
    # Project the numerical proposal onto the exact material/volume level.
    direction = cuts.fixed(instance.x[:18] / np.linalg.norm(instance.x[:18]))
    right = cuts.fixed(np.linalg.pinv(f))
    defect = cuts.norm(np.eye(len(f), dtype=int) - free @ right)
    assert defect.b < 1
    correction = cuts.norm(right) * cuts.norm(free @ direction) / (1 - defect)
    order = np.array([cuts.norm(v) for v in obs])
    lower = (obs @ direction - order * correction) / (cuts.norm(direction) + correction)
    if min(x.a for x in lower) <= 0:
        raise ValueError("prescribed sign separation could not be verified")
    entropy = cuts.ri("23.1042361") * lower[-1] / (iv.mpf(2) ** scaling[-1])
    return {
        "vacancy_sign": signs[0],
        "excess_sign": signs[1],
        "minimum_scaled_slope_ppb": min(
            floor(cuts.survey.dyadic(x._mpi_[0]) * 10**9) for x in lower[:-1]
        ),
        "minimum_entropy_loss_piconats": floor(
            cuts.survey.dyadic(entropy._mpi_[0]) * 10**12
        ),
        "projection_correction_ppb": ceil(
            cuts.survey.dyadic(correction._mpi_[1]) * 10**9
        ),
    }


def hidden_slope(rows, area, held=()):
    """Unit exact tangent with signed entropy slope and native held read bounds."""
    b, _ = normalize(rows)
    bm = cuts.mid(b)
    right = np.linalg.pinv(bm)
    price = cuts.mid(area) @ right
    residual = cuts.mid(area) - price @ bm
    picked = cuts.fixed(-residual / np.linalg.norm(residual))
    inverse = cuts.fixed(right)
    defect = cuts.norm(np.eye(len(b), dtype=int) - b @ inverse)
    if defect.b >= 1:
        raise ValueError("hidden tangent constraints lost rank")
    dr = cuts.norm(inverse) * cuts.norm(b @ picked) / (1 - defect)
    error = iv.mpf([-1, 1]) * dr
    # On the exact tangent the price multiple disappears. Subtracting it
    # gives a smaller enclosure for the geometric functional.
    a = area - cuts.fixed(price) @ b
    unit = cuts.norm(picked) + error
    signed = (
        iv.mpf(["23.1042361", "23.1042363"])
        * (-sum(a * picked) + cuts.norm(a) * error)
        / unit
    )
    if signed.a <= 0:
        raise ValueError("entropy gain did not separate the material levels")
    bounds = []
    for h in held:
        # Both forms agree on the exact kernel. Intersecting them can improve
        # the price bound when a very small outer-cell signal is measured.
        res = h - cuts.fixed(cuts.mid(h) @ right) @ b
        intervals = [(sum(z * picked) + cuts.norm(z) * error) / unit for z in (h, res)]
        bounds.append(
            iv.mpf([max(z.a for z in intervals), min(z.b for z in intervals)])
        )
    return signed, bounds


def outer_jump(regional, shells, area, volume):
    # Fix every layer except the contrast of the last two. The whole-domain
    # flux is separately fixed by the independent domain and image readings.
    protected = np.vstack(
        (
            regional,
            volume,
            shells[0, 1:4],
            shells[0, 4] + shells[0, 5],
            shells[1, 1:4],
            shells[1, 4] + shells[1, 5],
        )
    )
    entropy, bounds = hidden_slope(protected, area, shells[:, 5])
    delta = entropy / sum(abs(a).b for a in bounds)
    # Integer output rounded down to three significant digits, outward.
    minimum = floor(cuts.survey.dyadic(delta._mpi_[0]))
    quantum = 10 ** max(0, len(str(minimum)) - 3)
    return {
        "entropy_slope": cuts.survey.pack(entropy),
        "outer_price_jump_lower": minimum // quantum * quantum,
    }


def angular_profile(regional, cells, node_rows):
    # All 12 locations receive one common price for vacancy and excess.
    paired = node_rows[0] + node_rows[1]
    rows = np.vstack((regional, paired))
    # The centre and the smallest outer cell can both be priced at zero.
    # This reduces the raw readout amplification while keeping charge blindness.
    ids = list(range(len(regional))) + [
        len(regional) + i for i, u in enumerate(cells) if u not in (0, 70)
    ]
    b, exponents = normalize(rows[ids])
    left = cuts.fixed(np.linalg.inv(cuts.mid(b)))
    defect = cuts.norm(np.eye(18, dtype=int) - left @ b)
    if defect.b >= 1:
        raise ValueError("common cell weighting is not certified")
    k = cuts.exact_entropy_center()
    scaled_price = cuts.mid(k) @ cuts.mid(left)
    prices = ["0"] * len(rows)
    for i, x, e in zip(ids, scaled_price, exponents):
        prices[i] = f"{np.ldexp(x, e):.16g}"
    frozen = np.array([cuts.ri(x) for x in prices])
    residual = cuts.norm(k - frozen @ rows)
    radius = cuts.ri(
        str(
            json.loads((cuts.survey.HERE.parent / "certificate.json").read_text())[
                "state_HS_error_upper"
            ]
        )
    )
    error = residual + radius / (cuts.ri(".02") - radius)
    assert error.b < cuts.ri("2e-5")
    # For a comparison of angular dependence, a constant may be subtracted
    # from all cell prices and absorbed into the two regional boundary sums.
    gain = sum(abs(x) for x in frozen)
    return dict(
        rank=18,
        charge_blind_cells=list(cells),
        scaled_left_inverse_defect_ppm=ceil(
            cuts.survey.dyadic(defect._mpi_[1]) * 10**6
        ),
        regional_and_cell_prices=prices,
        modular_error_upper=cuts.survey.pack(error)[1],
        raw_price_l1_upper=ceil(cuts.survey.dyadic(gain._mpi_[1])),
    )
