"""Stationarity of the complete collar entropy under a correlated preparation lottery."""

import json

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import brentq, root

from .intrinsic_area import head_source, radar_banks, survey_corners
from .geometry_gain import simplex_area_volume


class Station:
    """Three independent corner probabilities with the lottery erased on every trial."""

    def __init__(self, corners, product, difference):
        self.base = corners[0, 0]
        self.dr = np.array(
            [
                corners[1, 0] - corners[0, 0],
                corners[0, 1] - corners[0, 0],
                corners[1, 1] - corners[1, 0] - corners[0, 1] + corners[0, 0],
            ]
        )
        empty = np.zeros_like(product)
        self.initial = np.stack([product, empty])
        # x=(odd probability, formed probability, odd-and-formed probability).
        self.density_directions = np.array(
            [[-product, product], [difference, empty], [-difference, difference]]
        )

    def geometry(self, x):
        ranges = self.base + np.einsum("i,ijk->jk", x, self.dr)
        values = np.array(simplex_area_volume(ranges))
        gradients = np.array(
            [
                np.imag(simplex_area_volume(ranges + 1e-16j * direction)) / 1e-16
                for direction in self.dr
            ]
        ).T
        return values, gradients

    def entropy(self, x):
        blocks = self.initial + np.einsum("i,ibjk->bjk", x, self.density_directions)
        entropy, first, second = 0.0, np.zeros(3), np.zeros((3, 3))
        for i, block in enumerate(blocks):
            lam, q = np.linalg.eigh(block)
            if min(lam) <= 0:
                raise ValueError("entropy reference must be faithful on both blocks")
            energy = np.log(lam)
            d = q.conj().T @ self.density_directions[:, i] @ q
            entropy -= lam @ energy
            first -= np.einsum("j,ijj->i", energy, d).real
            gap = lam[:, None] - lam[None, :]
            divided = np.divide(
                energy[:, None] - energy[None, :],
                gap,
                out=np.broadcast_to(1 / lam, gap.shape).copy(),
                where=(gap != 0),
            )
            second -= np.einsum("mij,nij,ij->mn", d.conj(), d, divided).real
        return entropy, first, second

    def read(self, x):
        _, (da, dv) = self.geometry(x)
        _, ds, _ = self.entropy(x)
        # Both comparisons use the same dx_3 adjustment read from the radar volume.
        a = da[:2] - da[2] * dv[:2] / dv[2]
        s = ds[:2] - ds[2] * dv[:2] / dv[2]
        eta = -s[0] / a[0]
        return dict(
            area=a.tolist(),
            entropy=s.tolist(),
            eta=eta,
            residual=s[1] + eta * a[1],
            forced_fraction=-eta * a[1] / s[1],
        )

    def stability(self, x, eta, step=1e-4):
        _, (da, dv) = self.geometry(x)
        _, ds, d2s = self.entropy(x)
        price = (ds + eta * da) @ dv / (dv @ dv)

        def derivative(pos):
            _, (a, v) = self.geometry(pos)
            return eta * a - price * v

        geometric_hessian = np.array(
            [
                (derivative(x + axis * step) - derivative(x - axis * step)) / (2 * step)
                for axis in np.eye(3)
            ]
        )
        basis = null_space(dv.reshape(1, 3))
        h = basis.T @ (d2s + geometric_hessian) @ basis
        return dict(
            gradient_norm=float(np.linalg.norm(basis.T @ (ds + eta * da))),
            hessian_eigenvalues=np.linalg.eigvalsh((h + h.T) / 2).tolist(),
            volume_multiplier=price,
        )

    def finite_maximum(self, x, eta):
        base, (da, dv) = self.geometry(x)
        start = self.entropy(x)[0]
        basis = null_space(dv.reshape(1, 3))
        results = []
        for first, second in ((1, 0), (0, 1), (2**-0.5, 2**-0.5)):
            for sign in (-1, 1):
                moved = x + sign * 0.002 * (basis @ np.array([first, second]))
                u, s, v = moved

                def records(odd):
                    pos = np.array([odd, s, v])
                    return self.geometry(pos)[0]

                odd = brentq(
                    lambda odd: records(odd)[1] - base[1],
                    u - 0.001,
                    u + 0.001,
                    xtol=3e-14,
                )
                pos = np.array([odd, s, v])
                current = records(odd)
                rates = dephase_rates(pos)
                if min(rates) <= 0 or max(rates) >= 1:
                    raise ValueError("unavailable source lottery")
                results.append(
                    dict(
                        tangent=[first, second],
                        sign=sign,
                        volume_error=current[1] - base[1],
                        change_generalized_entropy=self.entropy(pos)[0]
                        - start
                        + eta * (current[0] - base[0]),
                    )
                )
        return results


def dephase_rates(x):
    u, s, v = x
    return [2 * (u - v) / (1 - s), 2 * v / s]


def pre_calibrated_selection(station):
    # Sharp-sign probabilities 31/64 and 33/64 on fresh and formed lotteries,
    # respectively. Only the ordinary tangent on THIS source fixes eta.
    calibration = np.array([0.25, 0.5, 0.125 + 1 / 256])
    ordinary = station.read(calibration)
    eta = ordinary["eta"]
    vacuum_comparison = np.array([0.25, 0.5, 0.125])
    target_volume = station.geometry(vacuum_comparison)[0][1]

    def target(x):
        measured = station.read(x)
        return np.r_[
            (station.geometry(x)[0][1] - target_volume) / 0.05,
            np.asarray(measured["entropy"]) + eta * np.asarray(measured["area"]),
        ]

    solution = root(target, vacuum_comparison, tol=1e-10)
    if not solution.success or max(abs(target(solution.x))) > 2e-8:
        raise ValueError(
            "the specified source plate did not yield an interior stationary point"
        )
    x = solution.x
    trials = []
    for i, t in enumerate((1e-7, 1e-5)):
        samples = []
        for scale in (-1, 1):
            u, s = x[:2] + scale * t * np.eye(2)[i]

            def at(v):
                pos = np.array([u, s, v])
                return station.geometry(pos)[0], station.entropy(pos)[0]

            v = brentq(
                lambda v: at(v)[0][1] - target_volume,
                x[2] - 0.005,
                x[2] + 0.005,
                xtol=1e-14,
            )
            geometry, entropy = at(v)
            rates = dephase_rates([u, s, v])
            if min(rates) <= 0 or max(rates) >= 1:
                raise ValueError("infeasible lottery")
            samples.append(
                dict(
                    probabilities=[u, s, v],
                    dephase=rates,
                    area=geometry[0],
                    volume=geometry[1],
                    entropy=entropy,
                )
            )
        slope = np.array(
            [
                (samples[1][key] - samples[0][key]) / (2 * t)
                for key in ("area", "entropy")
            ]
        )
        trials.append(
            dict(
                tangent=("ordinary", "pairing")[i],
                step=t,
                arms=samples,
                measured_slope=slope.tolist(),
                fraction=-eta * slope[0] / slope[1],
            )
        )
    return dict(
        calibration_coordinates=calibration.tolist(),
        reference_rates=dephase_rates(calibration),
        ordinary_area=ordinary["area"][0],
        ordinary_entropy=ordinary["entropy"][0],
        eta=eta,
        target_volume=target_volume,
        selected_coordinates=x.tolist(),
        selected_rates=dephase_rates(x),
        stability=station.stability(x, eta),
        held_finite_differences=trials,
    )


def run():
    p, d = head_source()
    banks = radar_banks()
    u, s = 0.25, 0.5

    def pos(c):
        return np.array([u, s, u * s - c])

    locus = []
    for menu in (0, 1):
        for horizon in (2, 4, 6):
            instrument = Station(survey_corners(banks, menu, horizon), p, d)
            bracket = (-0.02, 0.0)
            lower, upper = (instrument.read(pos(c))["residual"] for c in bracket)
            root = brentq(
                lambda c: instrument.read(pos(c))["residual"], *bracket, xtol=1e-15
            )
            x = pos(root)
            result = instrument.read(x)
            locus.append(
                dict(
                    tick=banks[1, False][menu]["step"],
                    horizon=horizon,
                    correlation=root,
                    dephasing=dephase_rates(x),
                    bracket_residuals=[lower, upper],
                    exchange_energy=16 + u,
                    **result,
                    **instrument.stability(x, result["eta"]),
                    hessian_smaller_step=instrument.stability(
                        x, result["eta"], step=2e-5
                    )["hessian_eigenvalues"],
                )
            )
    target = locus[-1]
    x = pos(target["correlation"])
    eta = target["eta"]
    held = []
    for menu, horizon in ((0, 6), (1, 2), (1, 4), (1, 6)):
        station = Station(survey_corners(banks, menu, horizon), p, d)
        read = station.read(x)
        held.append(
            dict(
                tick=banks[1, False][menu]["step"],
                horizon=horizon,
                ordinary_fraction=-eta * read["area"][0] / read["entropy"][0],
                pairing_fraction=-eta * read["area"][1] / read["entropy"][1],
            )
        )
    selected = Station(survey_corners(banks, 1, 6), p, d)
    control = selected.read(pos(0))
    probe = selected.finite_maximum(x, eta)
    return dict(
        odd=u,
        mean_pairing=s,
        stations=locus,
        held_reference_and_eta=held,
        independent_lottery_at_selected_eta=dict(
            ordinary_residual=control["entropy"][0] + eta * control["area"][0],
            pairing_residual=control["entropy"][1] + eta * control["area"][1],
        ),
        selected_finite_neighbours=probe,
        independent_ordinary_calibration=pre_calibrated_selection(selected),
    )


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, allow_nan=False))
