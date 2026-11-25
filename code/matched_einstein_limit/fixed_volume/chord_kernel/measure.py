"""Opposite first exits on a routed cubic register and supplied curved rest charts."""

import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo import candidate

from model.candidate_laws.interval_echo.clock import Registers
from model.experiments.interval_echo_einstein_station.proper_ball_ruler.measure import (
    pulse_calibration,
)
from model.experiments.seam_program.clocked_apparatus import patches
from model.experiments.seam_program.compiler import Compiler


from .collar import held_bridge_menu
from .product_chart import survey as product_survey


HERE = Path(__file__).resolve().parent
AXES = np.eye(3, dtype=int)


def clock_branch(radius, radial, mu, curvature, tick):
    """First appointment reaching the rim, from independent space-form geodesic flow."""

    def s(x):
        k = math.sqrt(abs(curvature))
        if curvature > 0:
            return math.sin(k * x) / k
        if curvature < 0:
            return math.sinh(k * x) / k
        return x

    # q(x)=(1-C_K(x))/K, including its flat limit. This avoids cancellation.
    def q(x):
        return 2 * s(x / 2) ** 2

    initial, rim = q(radial), q(radius)
    cosine = 1 - curvature * initial

    def exited(n):
        t = n * tick
        return initial + cosine * q(t) + s(radial) * s(t) * mu >= rim

    low, high = 0, math.ceil(3 * radius / tick)
    assert not exited(low) and exited(high)
    while high - low > 1:
        mid = (high + low) // 2
        if exited(mid):
            high = mid
        else:
            low = mid
    return high


def curved_panels():
    prior = json.loads((HERE.parent / "measurements.json").read_text())
    out = []
    for panel in prior["panels"]:
        ell, tick = panel["flat_radius"], panel["proper_clock_tick"]
        cases = []
        for source in panel["cases"]:
            k = source["prepared_sectional_curvature"]
            rows = []
            for ratio in (0.25, 0.5, 0.75):
                r = ratio * ell
                for mu in (0.0, 1 / 3, 0.8, 1.0):
                    pair = [
                        clock_branch(ell, r, sign * mu, k, tick) for sign in (1, -1)
                    ]
                    product = math.prod(pair) * tick**2
                    radius_ticks, radial_ticks = round(ell / tick), round(r / tick)
                    measured_r = radial_ticks * tick
                    flat = (radius_ticks * tick) ** 2 - measured_r**2
                    rows.append(
                        dict(
                            radial_fraction=ratio,
                            incidence=mu,
                            exit_ticks=pair,
                            central_radius_ticks=radius_ticks,
                            central_launch_ticks=radial_ticks,
                            relative_weight_error=product / flat - 1,
                            perpendicular_curvature=(
                                3 * (product - flat) / (measured_r**2 * flat)
                                if mu == 0
                                else None
                            ),
                        )
                    )
            area = source["inferred"]["volume_threshold"]["area_einstein_contraction"]
            cases.append(
                dict(
                    sectional_curvature=k,
                    equal_volume_area_contraction=area,
                    rays=rows,
                    maximum_product_relative_error=max(
                        abs(x["relative_weight_error"]) for x in rows
                    ),
                    maximum_sectional_error=max(
                        abs(x["perpendicular_curvature"] - k)
                        for x in rows
                        if x["incidence"] == 0
                    ),
                    maximum_contraction_area_disagreement=max(
                        abs(3 * x["perpendicular_curvature"] - area)
                        for x in rows
                        if x["incidence"] == 0
                    ),
                )
            )
        out.append(dict(radius=ell, tick=tick, cases=cases))
    return out


def flat_quadrature(n):
    """Count every first-exit axial appointment for every interior lattice start."""
    # In mesh units, the first exterior station is +b or -b on each line.
    # Calibrate the physical radius to one for these ell^-4 normalized moments.
    axis = np.arange(-n + 1, n)
    i, j = np.meshgrid(axis, axis, indexing="ij")
    moments = np.zeros(3)
    maximum = 0.0
    for k in axis:
        bound = n * n - j * j - k * k
        b = np.ceil(np.sqrt(np.maximum(bound, 0))).astype(int)
        inside = i * i < bound
        distances = (b - i) * (b + i)
        weight = np.where(inside, distances / (2 * n * n), 0.0)
        exact = np.where(inside, (bound - i * i) / (2 * n * n), 0.0)
        maximum = max(maximum, float(np.max(np.abs(weight - exact))))
        moments += (
            np.array(
                [
                    weight.sum(),
                    (weight * i * i / (n * n)).sum(),
                    (weight * k * k / (n * n)).sum(),
                ]
            )
            / n**3
        )
    return dict(
        points_per_radius=n,
        maximum_weight_error_per_radius=maximum,
        integral_per_radius_fourth=moments.tolist(),
        mass_relative_error=moments[0] / (4 * math.pi / 15) - 1,
    )


class Rail:
    """Positive addressed seam services. One logical coordinate step has fixed hand duration."""

    def __init__(self):
        self.host, self.data, self.physical, front, _ = patches(1)
        # exchange() works between any two data ports, independent of the walk's eight-slot orbital.
        self.compiler = Compiler(self.data, 0, (1, 1, 1, 1))
        self.register = Registers(self.host, 183707)
        self.sites = np.arange(1, 730).reshape(9, 9, 9)
        a = 12 * front
        b = next(b for b in self.host.neighbours[a] if b // 12 == front)
        self.hand = (a, b)
        self.vacancies = tuple(
            next(
                c
                for c in self.host.neighbours[x]
                if c not in self.hand and c // 12 == front
            )
            for x in self.hand
        )
        # A contact path excluding the marker has at most d.n-2 edges; its restored
        # transposition costs at most 2*(d.n-2)-1 microperiods. The enrolled routes
        # below have a tighter verified bound and run in fixed slots of 65.
        self.slot = 65

    def address(self, point):
        return int(self.sites[tuple(x + 4 for x in point)])

    @lru_cache(maxsize=512)
    def schedule(self, a, b):
        instructions = self.compiler.exchange(a, b, 4)
        assert all(q == 4 for _, q in instructions) and len(instructions) <= self.slot
        return [
            tuple(self.physical[i] for i in self.data.seams[e]) for e, _ in instructions
        ]

    def exit(self, point, axis, sign, turn=False, stationary=False, z_scale=1):
        self.register.loads.fill(0)
        self.register.loads[self.hand[0]] = 2
        self.register.loads[self.physical[self.address(point)]] = 2
        trace = []
        switch = 0
        pulses = 0
        while True:
            direction = (axis + 1) % 3 if turn and trace else axis
            dest = tuple(point + sign * AXES[direction])
            events = self.schedule(self.address(point), self.address(dest))

            for period in range(self.slot):
                if period < len(events):
                    a, b = events[period]
                    if not stationary:
                        self.register.swap(a, b, {a // 12, b // 12})
                    pulses += 1
                a, b = self.hand
                self.register.swap(a, b, {a // 12})
                switch = 1 - switch
                terminal = self.hand[switch]
                assert self.register.flag(terminal, self.vacancies[switch], {a // 12})
            terminal = self.physical[self.address(dest)]
            # Fresh 0-load reference chosen away from both candidate tag addresses.
            aux = next(
                c
                for c in self.host.neighbours[terminal]
                if c != self.physical[self.address(point)]
                and c // 12 != self.hand[0] // 12
            )
            found = self.register.flag(terminal, aux, {terminal // 12, aux // 12})
            trace.append([list(map(int, dest)), found])
            if (
                not found
                or sum(int(v * z) ** 2 for v, z in zip(dest, (1, 1, z_scale))) >= 3.5**2
            ):
                return dict(
                    hand_reads=len(trace) * self.slot, data_pulses=pulses, trace=trace
                )
            point = np.array(dest)


def native_rails():
    rail = Rail()
    calibration = pulse_calibration()
    rows = []
    for site in ((0, 0, 0), (1, 1, 0), (1, -1, 2), (2, 0, 1)):
        for axis in range(3):
            runs = [rail.exit(np.array(site), axis, sign) for sign in (1, -1)]
            legs = [arm["hand_reads"] // rail.slot for arm in runs]
            assert all(record[1] for arm in runs for record in arm["trace"])
            pairs = sum(u < legs[0] and v < legs[1] for u in range(7) for v in range(7))
            assert pairs == math.prod(legs)
            rows.append(
                dict(
                    site=list(site),
                    axis=axis,
                    legs=legs,
                    product_clock_weight=pairs / 7,
                    flat_weight=(3.5**2 - sum(x * x for x in site)) / 7,
                    trials=runs,
                )
            )
    altered = [rail.exit(np.zeros(3, dtype=int), 0, s, turn=True) for s in (1, -1)]
    anisotropic = [rail.exit(np.zeros(3, dtype=int), 2, s, z_scale=2) for s in (1, -1)]
    # Same data schedule; the contrast-only control preserves each blank endpoint word.
    law = candidate.native_law
    _, agreement = law.local_term_parts()
    changed = expm(-1j * np.pi * candidate.KAPPA * agreement)
    for src in ((2, 0), (0, 2), (0, 0)):
        index = law.SRC.index(src)
        np.testing.assert_allclose(
            np.abs(changed[:, index]), np.eye(law.DIM)[:, index], atol=1e-13
        )
    no_transfer = rail.exit(np.zeros(3, dtype=int), 0, 1, stationary=True)
    assert not no_transfer["trace"][-1][1]
    return dict(
        chart_width=9,
        ball_radius=3.5,
        hand_periods_per_mesh_edge=rail.slot,
        native_hand_single_period=calibration["decoded_single_step"],
        straight=rows,
        turn_after_first=altered,
        anisotropic_aperture=anisotropic,
        contrast_only_data=no_transfer,
    )


def run():
    return dict(
        native_chart=native_rails(),
        flat_mesh=[flat_quadrature(n) for n in (8, 16, 32, 64)],
        comparison_charts=curved_panels(),
        product_charts=[
            product_survey(ell, 2**-18 * ell**5, k)
            for ell in (0.5, 0.25, 0.125, 0.0625)
            for k in (-1.0, 1.0)
        ],
        formed_collar=held_bridge_menu(),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    data = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
