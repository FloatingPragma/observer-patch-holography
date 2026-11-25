"""Compare the four-hinge coefficients with an independent flight and energy evolution."""

import json

import numpy as np
from scipy.integrate import quad

from .. import aperture, flight, hinge
from ... import live_radar, live_source, live_survey
from ...source import measure as source
from . import native


def individual(alpha):
    # Position-space wedge quadrature. At distance u from the axis integrate
    # the volume missing from the unit ball over forbidden azimuths.
    beta = alpha / 2

    def moment(u):
        h = np.sqrt(1 - u * u)

        def slices(phi):
            c = -u * np.cos(phi)
            psi = np.arctan2(h, c)
            return (
                2 * h**3 / 3 - c * ((h * h + c * c) * psi - c * h),
                2 * (h - c * psi),
            )

        v = 4 * np.pi / 3 - 2 * quad(lambda phi: slices(phi)[0], beta, np.pi)[0]
        a = 4 * np.pi - 2 * quad(lambda phi: slices(phi)[1], beta, np.pi)[0]
        return u * (4 * np.pi / 3) ** (4 / 3) * a * a / (3 * v ** (7 / 3))

    return alpha * (2 * np.pi - quad(moment, 0, 1)[0])


def test_native_coefficients_and_implicit_controls():
    saved = json.loads((native.HERE / "native.json").read_text())
    rebuilt = native.run()
    radar = flight.pooled(json.loads((flight.HERE / "flight.json").read_text()), 16)
    lengths = radar["reference"][:5, :5]
    alphas = 2 * np.pi - hinge.bend(lengths)[1]
    coeff = np.array([individual(a) for a in alphas])
    step = 1e-4
    derivative = np.array(
        [(individual(a + step) - individual(a - step)) / (2 * step) for a in alphas]
    )
    for record in (saved, rebuilt):
        for row, a, v, dv in zip(
            record["radial_edges"], alphas, coeff, derivative, strict=True
        ):
            assert row["angle"][0] < a < row["angle"][1]
            lo, hi = row["individually_fixed_volume_area_coefficient"]
            assert lo < v < hi
            lo, hi = row["individually_fixed_volume_angle_derivative"]
            assert lo < dv < hi
            assert row["mean_fixed_volume_area_deficit_coefficient"][1] < v
        for tag in native.ARMS:
            d = radar[tag][:5, :5]
            da = -hinge.bend(lengths + 1e-19j * d)[1].imag / 1e-19
            slope = (d[4, :4] @ coeff + (lengths[4, :4] * derivative) @ da) / 3
            lo, hi = record["individual_area_slopes"][tag]
            assert lo < slope < hi
        for name in (
            "equal_mean_area_coefficient_and_regional_number",
            "equal_individual_area_coefficient_and_regional_number",
        ):
            trial = record[name]
            assert trial["control_determinant"][0] > 0
            assert trial["matched_slopes"]["energy"][1] < -0.2
            assert trial["matched_slopes"]["regional_entropy"][0] > 0.01
        for protocol in record["ordinary_calibration"].values():
            for prediction in protocol.values():
                assert prediction["measured_minus_predicted"][0] > 0.001
        direct = record["individual_area_slopes"]["decorrelated"]
        entropy_prediction = record["ordinary_calibration"]["individual"]["regional_entropy"]
        assert direct[0] / entropy_prediction["predicted_repairing"][1] > 59
        assert record["radial_edges"][0]["axis_area_slopes_per_radius_squared"][
            "decorrelated"
        ][1] < -0.00376
        trial = record["equal_regional_energy_and_number"]
        for name in ("mean_area_coefficient", "individual_area_coefficient"):
            assert trial["matched_slopes"][name][0] > 0.001
        assert trial["matched_slopes"]["regional_entropy"][0] > 0.19


def test_regional_energy_and_number_from_independent_evolution():
    record = native.regional_energy()
    board = source.formation.federation(1)
    route = source.formation.region(board)
    inside = live_radar.station_rows(board, aperture.enroll(board, 1)[0]).ravel()
    w = (
        live_source.head_channel()
        .reshape(-1, 2, 4, 8)
        .transpose(1, 0, 2, 3)
        .reshape(-1, 8)
    )
    h = source.system(board, route, -1)[0]
    hw = h @ w / 40
    ph = w[inside].conj().T @ hw[inside]
    effects = {
        "number": w[inside].conj().T @ w[inside],
        "energy": (ph + ph.conj().T) / 2,
    }
    ref, product = source.form(board, route, -1)
    gens = live_survey.rotations(ref)
    inputs = {
        "reference": ([2 / 3, 1 / 3], ref),
        "ordinary": ([-1 / 3, 1 / 3], ref),
        "decorrelated": ([2 / 3, 1 / 3], product - ref),
    }
    for name, g in gens.items():
        inputs[name] = ([2 / 3, 1 / 3], -1j * (g @ ref - ref @ g))
    for name, observable in effects.items():
        for kind, (weights, rho) in inputs.items():
            x = np.trace(observable @ np.kron(np.diag(weights), rho)).real
            lo, hi = record[name][kind]
            # Floating regeneration of the pump, formation and trace are separate
            # from the integer enclosures. Its own machine precision is finite.
            assert lo - 2e-11 < x < hi + 2e-11
