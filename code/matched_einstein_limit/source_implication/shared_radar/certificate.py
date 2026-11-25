"""Which state functions can two charge-resolved four-station surveys share?"""

from fractions import Fraction
from pathlib import Path
import json
import numpy as np
import sympy as sp

from ..null_operator import measure as operator
from ..full_algebra import certificate as source
from ..full_algebra.enclosure import geometry_jet, iv


def left_inverse_defect(numerators, unit, error):
    """Exact infinity norm bound for the defect of a proposed left inverse."""
    rows, cols = numerators.shape
    q = operator.dyadic(np.linalg.pinv(numerators.astype(float) / unit, rcond=1e-16))
    identity = np.eye(cols, dtype=int).astype(object)
    norm_q = Fraction(operator.max_row_sum(q), operator.SCALE)
    defect = (
        Fraction(
            operator.max_row_sum(q @ numerators - identity * (unit * operator.SCALE)),
            unit * operator.SCALE,
        )
        + cols * norm_q * error
    )
    if defect >= 1:
        raise ArithmeticError(f"Uncertified left inverse for shape {rows} by {cols}")
    return float(defect)


def invisible_phase(matrices, traces, unit, error):
    """Construct an exactly volume-matched, marginal-preserving correlation tangent."""
    data = json.loads(
        Path(operator.station.__file__).with_name("results.json").read_text()
    )
    x = np.array(data["x_real"]) + 1j * np.array(data["x_imag"])
    blocks = [source.prior.sym_matrix(b) for b in x]
    blocks = [(b + b.H) / 2 for b in blocks]
    tr = sum(b.trace() for b in blocks)
    blocks = [b / tr for b in blocks]
    re, im = operator.integer_directions()
    basis = [
        [(sp.Matrix(r) + sp.I * sp.Matrix(i)) / 4 for r, i in zip(real, imag)]
        for real, imag in zip(re, im)
    ]
    norms = [sum(source.inner(b, b) for b in row) for row in basis]
    coords = [
        sum(source.inner(b, x) for b, x in zip(row, blocks)) / norm
        for row, norm in zip(basis, norms)
    ]
    pooled = matrices[-1, 6] + matrices[1, 6]
    tracial = traces[-1, 6] + traces[1, 6]
    # Hermitian range effects lie between 0 and 3. The source displacement
    # contributes at most 3 sqrt(8) times its certified HS radius.
    state_radius = 3 * iv.sqrt(8) * source.ri("2e-7")
    ranges = []
    for row, trace in zip(pooled, tracial):
        r = (
            3
            + sp.Rational(int(trace), 4 * unit)
            + sum(sp.Rational(int(z), 2 * unit) * c for z, c in zip(row, coords))
        )
        arithmetic = source.ri(error) * (1 + sum(source.ri(abs(c)) for c in coords))
        ranges.append(source.ri(r) + iv.mpf([-1, 1]) * (state_radius + arithmetic))
    area, vol = geometry_jet(ranges)
    _, _, Y, Z = operator.station.phase.AXES
    indexes = []
    for a, b in ((Y, Z), (Z, Y)):
        v = np.stack([np.zeros((4, 4)), np.kron(a, b) / 2])
        found = np.flatnonzero(
            np.all(
                np.array(re, dtype=float) + 1j * np.array(im, dtype=float) == 4 * v,
                axis=(1, 2, 3),
            )
        )
        assert len(found) == 1
        indexes.append(int(found[0]))
    gradients = []
    for i in indexes:
        dr = [source.ri_band(Fraction(int(z), 2 * unit), error) for z in pooled[:, i]]
        gradients.append(
            (
                sum(a * b for a, b in zip(area.g, dr)),
                sum(a * b for a, b in zip(vol.g, dr)),
            )
        )
    (area1, vol1), (area2, vol2) = gradients
    assert vol2.a > 0
    # Unit tangent proportional to vol1 * D2 - vol2 * D1.
    length = iv.sqrt(vol1**2 + vol2**2)
    area_slope = (vol1 * area2 - vol2 * area1) / length
    entropy_slope = -iv.mpf(["23.1042361", "23.1042363"]) * area_slope
    assert entropy_slope.a > source.ri(".00030")
    # Independent matter bound from the odd head block's logarithm. Both
    # perturbations have zero marginal traces and trace norm two.
    coords, log_error = source.prior.entropy_coordinates(
        blocks[1],
        sp.Rational("2e-7"),
        sp.Rational(".02"),
        [basis[i][1] for i in indexes],
    )
    k1, k2 = [source.ri_band(g, 2 * log_error) for g in coords]
    independent = (vol1 * k2 - vol2 * k1) / length
    assert independent.a > source.ri(".00029")

    return dict(
        components=["odd Y tensor Z / 2", "odd Z tensor Y / 2"],
        adjusting_volume_slope=str(vol2),
        matched_area_slope=str(area_slope),
        matched_collar_slope=str(entropy_slope),
        independent_modular_slope=str(independent),
    )


def run():
    iv.dps = 30
    matrices, traces, unit, error = operator.resolved_range_matrices()
    tracial_geometry = {}
    full_defects = {}
    for charge in (-1, 1):
        values = [
            source.ri_band(3 + Fraction(int(z), 2 * unit), error / 2)
            for z in traces[charge, 6]
        ]
        area, volume = geometry_jet(values)
        # Two entries of their differential suffice to exclude proportionality.
        minor = area.g[0] * volume.g[5] - area.g[5] * volume.g[0]
        assert minor.b < -0.08
        tracial_geometry[str(charge)] = dict(
            area=str(area.v), volume=str(volume.v), area_volume_minor=str(minor)
        )
        family = np.vstack([matrices[charge, c] for c in operator.CAPS])
        full_defects[str(charge)] = left_inverse_defect(family, unit, error)

    stacked = np.vstack([matrices[q, 6] for q in (-1, 1)])
    return dict(
        single_timeout_stacked_rank=24,
        stacked_right_inverse_defect=left_inverse_defect(stacked.T, unit, error),
        tracial_geometry=tracial_geometry,
        three_timeout_single_charge_rank=31,
        invisible_phase=invisible_phase(matrices, traces, unit, error),
        left_inverse_defects=full_defects,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
