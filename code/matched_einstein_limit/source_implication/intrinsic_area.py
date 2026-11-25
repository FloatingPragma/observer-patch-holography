"""Measure the full stored-head entropy and radar area on a common mixed preparation."""

import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, root

from .geometry_gain import pairing, pooled_ranges, ranging, simplex_area_volume


def head_source():
    board = pairing.formation.federation(1)
    joint, _, _ = pairing.prepared_states(board, pairing.formation.region(board), 1)
    joint /= np.trace(joint).real
    product = np.kron(*pairing.formation.marginals(joint))
    return product, joint - product


def full_entropy_slopes(product, delta, odd, strength):
    rho = product + strength * delta
    parity = np.diag([1 - odd, odd])
    log_rho = pairing.log_density(np.kron(parity, rho))
    return np.array(
        [
            -np.trace(np.kron(np.diag([-1, 1]), rho) @ log_rho).real,
            -np.trace(np.kron(parity, delta) @ log_rho).real,
        ]
    )


def geometry_slopes(corners, odd, strength):
    at_parity = corners[:, 0] + strength * (corners[:, 1] - corners[:, 0])
    ranges = (1 - odd) * at_parity[0] + odd * at_parity[1]
    ordinary = at_parity[1] - at_parity[0]
    pairing_change = (1 - odd) * (corners[0, 1] - corners[0, 0]) + odd * (
        corners[1, 1] - corners[1, 0]
    )
    area, volume = simplex_area_volume(ranges)
    jacobian = np.array(
        [
            np.imag(simplex_area_volume(ranges + 1e-16j * dr)) / 1e-16
            for dr in (ordinary, pairing_change)
        ]
    ).T
    area_slopes, volume_slopes = jacobian
    # Holding target volume by one homothety in each arm gives this dimensionless area.
    shape_slopes = (area_slopes - (2 / 3) * area / volume * volume_slopes) / volume ** (
        2 / 3
    )
    # Change the conditional pairing oppositely in the two parity blocks.
    # Their weighted average is fixed and the entropy derivative vanishes at independence.
    adjuster = (corners[0, 1] - corners[0, 0]) - (corners[1, 1] - corners[1, 0])
    adjust_area, adjust_volume = (
        np.imag(simplex_area_volume(ranges + 1e-16j * adjuster)) / 1e-16
    )
    if abs(adjust_volume) < 1e-14:
        raise ValueError("this source adjuster cannot hold survey volume")
    adjust = -volume_slopes / adjust_volume
    constrained = area_slopes + adjust_area * adjust
    return area, volume, area_slopes, volume_slopes, shape_slopes, constrained, adjust


def geometry_entropy(corners, product, delta, odd, even_strength, odd_strength):
    even_ranges = corners[0, 0] + even_strength * (corners[0, 1] - corners[0, 0])
    odd_ranges = corners[1, 0] + odd_strength * (corners[1, 1] - corners[1, 0])
    area, volume = simplex_area_volume((1 - odd) * even_ranges + odd * odd_ranges)
    # The parity blocks have orthogonal support, and the preparation lottery is erased.
    spectrum = np.concatenate(
        (
            (1 - odd) * np.linalg.eigvalsh(product + even_strength * delta),
            odd * np.linalg.eigvalsh(product + odd_strength * delta),
        )
    )
    if np.min(spectrum) <= 0:
        raise ValueError("the two blocks must be faithful")
    return np.array([area, volume, -spectrum @ np.log(spectrum)])


def finite_pairs(corners, product, delta):
    odd, strength, step = 1 / 8, 1 / 4, 1 / 128

    def records(u, s):
        return geometry_entropy(corners, product, delta, u, s, s)

    initial = records(odd, strength)
    pairs = []
    for s in (strength - step, strength + step):
        u = brentq(
            lambda x: records(x, s)[1] - initial[1], odd - 0.01, odd + 0.01, xtol=1e-14
        )
        pairs.append(dict(odd=u, strength=s, records=records(u, s).tolist()))

    # A richer preparation can keep both recorded volume and surface area fixed.
    # Let each pairing arm have its own sharp-read probability. The resulting
    # four corner weights are reproducible iff neither conditional odd rate exceeds 1/2.
    base = np.array([0.25, 0.5, 0.5])
    target = geometry_entropy(corners, product, delta, *base)
    held = []
    for odd_strength in (0.35, 0.50, 0.65):

        def difference(uv):
            return 1e4 * (
                geometry_entropy(corners, product, delta, *uv, odd_strength)[:2]
                - target[:2]
            )

        fit = root(difference, base[:2], tol=1e-10)
        if max(abs(fit.fun)) > 2e-10:
            raise ArithmeticError(
                "failed to match the two independent radar quantities"
            )
        u, s = fit.x
        weights = np.array(
            [
                [(1 - u) * (1 - s), (1 - u) * s],
                [u * (1 - odd_strength), u * odd_strength],
            ]
        )
        conditional = weights[1] / weights.sum(axis=0)
        if np.any(weights <= 0) or max(conditional) >= 0.5:
            raise ArithmeticError("native random dephasing cannot prepare this source")
        held.append(
            dict(
                odd=u,
                even_pairing=s,
                odd_pairing=odd_strength,
                dephase_probabilities=(2 * conditional).tolist(),
                records=geometry_entropy(
                    corners, product, delta, u, s, odd_strength
                ).tolist(),
            )
        )
    # Direct finite ordinary and fixed-energy pairing comparisons on separate copies
    # of the same interior source. Only the independent volume determines the adjuster.
    corrections = {}
    for label, du, ds, step in (("ordinary", 1, 0, 1e-7), ("pairing", 0, 1, 1e-4)):
        arms = []
        for t in (-step, step):
            u = base[0] + du * t
            s = base[1] + ds * t

            def at(c):
                return geometry_entropy(
                    corners, product, delta, u, s + c / (1 - u), s - c / u
                )

            c = brentq(lambda x: at(x)[1] - target[1], -0.01, 0.01, xtol=5e-15)
            rates = 2 * np.array([u * (1 - s) + c, u * s - c]) / np.array([1 - s, s])
            if not np.all((rates > 0) & (rates < 1)):
                raise ValueError("invalid sharp-read probabilities")
            arms.append(
                dict(
                    odd=u,
                    mean_pairing=s,
                    adjuster=c,
                    dephase_probabilities=rates.tolist(),
                    records=at(c).tolist(),
                )
            )
        slope = np.subtract(arms[1]["records"], arms[0]["records"]) / (2 * step)
        corrections[label] = dict(arms=arms, step=step, centered_slopes=slope.tolist())
    return dict(
        equal_volume=pairs, equal_area_and_volume=held, direct_calibration=corrections
    )


def radar_banks():
    archive = json.loads(Path(ranging.__file__).with_name("results.json").read_text())
    return {
        (arm["probe"], arm["left_odd"]): arm["panels"]
        for arm in archive["runs"]
        if arm["contrast"] == 8
    }


def survey_corners(banks, clock, horizon):
    return np.array(
        [
            [
                pooled_ranges(
                    [banks[sign, parity][clock] for sign in (-1, 1)],
                    arm,
                    horizon,
                )
                for arm in ("independent", "formed")
            ]
            for parity in (False, True)
        ]
    )


def run():
    ranging.ordinary_write()
    product, delta = head_source()
    banks = radar_banks()
    panels = []
    for clock in (0, 1):
        for horizon in (2, 4, 6):
            corners = survey_corners(banks, clock, horizon)
            if clock == 1 and horizon == 6:
                finite = finite_pairs(corners, product, delta)
            for odd in (1 / 8, 1 / 4, 3 / 8):
                for strength in (1 / 4, 1 / 2, 3 / 4):
                    entropy = full_entropy_slopes(product, delta, odd, strength)
                    (
                        area,
                        volume,
                        area_slopes,
                        volume_slopes,
                        shape,
                        constrained,
                        adjust,
                    ) = geometry_slopes(corners, odd, strength)
                    # Instead of rescaling any length, change the actual source along dV=0.
                    # Per unit strengthening of pairing, du/ds follows from the measured volume.
                    du = -volume_slopes[1] / volume_slopes[0]
                    along = np.array([du, 1.0])
                    eta = -entropy[0] / shape[0]
                    eta_direct = -entropy[0] / constrained[0]
                    panels.append(
                        dict(
                            tick=banks[1, False][clock]["step"],
                            horizon=horizon,
                            odd=odd,
                            pairing=strength,
                            entropy=entropy.tolist(),
                            area=area,
                            volume=volume,
                            area_slopes=area_slopes.tolist(),
                            volume_slopes=volume_slopes.tolist(),
                            shape=shape.tolist(),
                            ordinary_eta_on_unit_volume_shape=eta,
                            held_entropy_residual=entropy[1] + eta * shape[1],
                            direct_volume_adjusters=adjust.tolist(),
                            direct_volume_area_derivatives=constrained.tolist(),
                            direct_volume_ordinary_eta=eta_direct,
                            direct_volume_pairing_residual=entropy[1]
                            + eta_direct * constrained[1],
                            direct_volume_transfer=-eta_direct
                            * constrained[1]
                            / entropy[1],
                            fixed_volume_ordinary_change=du,
                            fixed_volume_entropy=entropy @ along,
                            fixed_volume_area=area_slopes @ along,
                        )
                    )
    return dict(
        cells=80,
        ports=list(pairing.formation.region(pairing.formation.federation(1))),
        panels=panels,
        finite_pairs=finite,
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
