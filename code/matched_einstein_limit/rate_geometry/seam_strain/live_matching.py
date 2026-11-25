"""Finite native preparation curves at fixed five-station interior volume."""

import argparse
import json
from itertools import combinations

import mpmath as mp
import numpy as np
from scipy.optimize import brentq, root

from . import live_radar, live_source, live_survey
from .source import measure as source


def area_volume(base, change):
    # Tiny source changes can move a macroscopic area by < 10^-12.
    with mp.workdps(42):
        d = mp.matrix(5)
        for i, j in combinations(range(5), 2):
            d[i, j] = d[j, i] = mp.mpf(float(base[i, j])) + mp.mpf(float(change[i, j]))
        g = mp.matrix(4)
        for i in range(4):
            for j in range(4):
                g[i, j] = (d[i, 4] ** 2 + d[j, 4] ** 2 - d[i, j] ** 2) / 2
        area = volume = mp.mpf(0)
        for face in live_survey.FACES:
            i, j, k = face
            product = (d[i, j] ** 2 + d[i, k] ** 2 - d[j, k] ** 2) / 2
            area += mp.sqrt(d[i, j] ** 2 * d[i, k] ** 2 - product**2) / 2
            volume += (
                mp.sqrt(mp.det(mp.matrix([[g[i, j] for j in face] for i in face]))) / 6
            )
        return area, volume


def finite_pairs(panel):
    ports = live_survey.stations()[:5]
    board = source.formation.federation(1)
    rho, product = source.form(board, source.formation.region(board), -1)
    generators = live_survey.rotations(rho)
    e = sum(live_radar.range_effects(1, q, True, ports, 12) for q in (-1, 1)) / 2
    kraus = live_source.head_channel(True)
    # Every input contrast has trace exactly zero. Subtract the common timeout.
    effective = np.einsum(
        "pai,eab,pbj->eij", kraus.conj(), e - 3 * np.eye(8), kraus, optimize=True
    )
    base = np.array(panel["lengths"])[:5, :5]
    zero = np.zeros_like(base)
    a0, v0 = area_volume(base, zero)
    rows = []
    fixed = []
    baseline = panel["centers"][0]
    jac = np.array(
        [
            [baseline["slopes"][side][observable] for side in generators]
            for observable in ("area", "volume")
        ]
    )
    for arm, t in (
        ("ordinary", 1e-6),
        ("decorrelated", 1e-6),
        ("hidden", 1e-8),
        ("hidden", -1e-8),
    ):
        dhead = {
            "ordinary": np.zeros_like(rho),
            "decorrelated": product - rho,
            "hidden": live_source.input_witness(),
        }[arm]
        rail0 = np.diag([2 / 3, 1 / 3])
        drail = np.diag([-1 / 3, 1 / 3]) if arm == "ordinary" else np.zeros((2, 2))
        middle = rho + t * dhead
        if np.linalg.eigvalsh(middle).min() <= 0:
            raise ValueError("inadmissible finite source preparation")

        def changes(angles):
            # Stable evaluation of (U_left U_right - I) at sub-microradian angles.
            dl, dr = [
                -2 * np.sin(z / 4) ** 2 * np.eye(4) - 2j * np.sin(z / 2) * g
                for z, g in zip(angles, generators.values(), strict=True)
            ]
            du = dl + dr + dl @ dr
            rotated = (
                t * dhead
                + du @ middle
                + middle @ du.conj().T
                + du @ middle @ du.conj().T
            )
            input_change = np.kron(rail0 + t * drail, rotated) + np.kron(t * drail, rho)
            ranges = live_survey.read(effective, input_change)
            a, v = area_volume(base, live_survey.square(ranges, 5))
            return np.array([float(a - a0), float(v - v0)])

        for i, side in enumerate(generators):
            axis = np.eye(2)[i]
            proposal = baseline["matched_boundary"][side][arm]["rotation"] * t
            low, high = sorted((0.0, 2 * proposal))
            angle = brentq(lambda z: changes(z * axis)[1], low, high, xtol=1e-18)
            da, dv = changes(angle * axis)
            rows.append(
                dict(
                    arm=arm,
                    side=side,
                    admixture=t,
                    rotation=angle,
                    volume_change=dv,
                    area_change=da,
                    area_slope=da / t,
                )
            )
        slope = baseline["slopes"][arm]
        guess = np.linalg.solve(jac, [-slope["area"], -slope["volume"]])
        sol = root(lambda angles: np.linalg.solve(jac, changes(t * angles) / t), guess)
        residual = changes(sol.x * t)
        if np.max(abs(residual)) > 3e-13 * abs(t) or not np.isfinite(sol.x).all():
            raise ArithmeticError(f"unmatched {arm} at {t}: {residual}, {sol.message}")
        fixed.append(
            dict(
                arm=arm,
                admixture=t,
                angles=(sol.x * t).tolist(),
                area_change=residual[0],
                volume_change=residual[1],
            )
        )
    return rows, fixed


def run():
    panels = json.loads((live_survey.HERE / "live_survey.json").read_text())["panels"]
    primary = next(
        p
        for p in panels
        if p["read"] == "pooled" and p["deadline"] == 12 and p["coupled"]
    )
    rows, fixed = finite_pairs(primary)
    return dict(
        stations=live_survey.stations()[:5],
        finite_pairs=rows,
        matched_area_volume=fixed,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (live_survey.HERE / "live_matching.json").write_text(data)
    else:
        print(data, end="")
