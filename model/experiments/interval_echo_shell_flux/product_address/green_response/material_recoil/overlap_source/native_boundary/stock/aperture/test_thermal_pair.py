import json

import numpy as np
from scipy import sparse
from mpmath import mp

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.aperture import (
    thermal_energy as energy,
    thermal_pair as experiment,
    thermal_phase_survey as combined,
)


def test_thermal_enclosures_against_direct_noncommuting_matrices():
    kernel = experiment.kernel
    h = (
        np.diag([2.0, 1.0, 3.0, 2.0, 1.0])
        - np.diag([0.5] * 4, 1)
        - np.diag([0.5] * 4, -1)
    )
    sparse_h = sparse.csr_matrix(h)
    sites = [0, 3]
    heat = kernel.mixed_bounds(sparse_h, sites, 3)

    def direct(a, b):
        k = h.copy()
        k[0, 0] -= a
        k[3, 3] -= b
        with mp.workdps(55):
            mat = mp.matrix(k.tolist())
            w = mp.expm(-mp.mpf(kernel.BETA) * mat)
            return tuple(
                np.array([float(a[i, i]) for i in range(len(h))]) for a in (w, mat * w)
            )

    z, y = direct(0, 0)
    both, _ = direct(32, 33)
    da, _ = direct(32, 0)
    db, _ = direct(0, 33)
    r = both - da - db + z
    assert (
        min(r) >= -2e-14
        and heat["lower_at_first_site"] <= r[0] <= heat["upper_each_diagonal"]
    )
    for e, target in enumerate((z, y)):
        approx, err = kernel.baseline(sparse_h, energy=bool(e))
        assert max(abs(approx - target)) <= err
        increments, cert = kernel.single_updates(sparse_h, sites, energy=bool(e))
        for i, row in enumerate(increments):
            for u, v in row.items():
                full = direct(u, 0) if i == 0 else direct(0, u)
                np.testing.assert_allclose(v, full[e] - target, atol=3e-13, rtol=2e-13)
    trace = energy.trace_data(sparse_h, sites, kernel.BETA, heat)
    read = energy.trace_data(sparse_h, sites, kernel.BETA, heat, energy=True)
    assert (
        float(energy.log_mix(trace).a)
        <= sum(
            np.log(
                direct(a, b)[0].sum()
                * z.sum()
                / (direct(a, 0)[0].sum() * direct(0, b)[0].sum())
            )
            / 16
            for a in (32, 33)
            for b in (32, 33)
        )
        <= float(energy.log_mix(trace).b)
    )

    def U(a, b):
        w, ew = direct(a, b)
        return sum(ew) / sum(w)

    exact = sum(
        (U(a, b) - U(a, 0) - U(0, b) + U(0, 0)) / 16 for a in (32, 33) for b in (32, 33)
    )
    lo, hi = experiment.interval_output(energy.energy_mix(trace, read))
    assert lo <= exact <= hi


def test_phase_flags_use_the_same_parity_as_the_material_reading():
    # Make a sharp address reveal both parities. Averaging over parity before
    # interrogating the heads preserves the marginals but loses this relation.
    values = experiment.STRENGTHS
    sectors = {}
    for i, a in enumerate(values):
        for j, b in enumerate(values):
            sectors[f"{-a},{-b}"] = np.eye(9)[3 * i + j]
    joint = combined.thermal_flags(sectors).reshape(4, 4, 9)
    assert not np.any(joint[3, 3, [0, 1, 2, 3, 6]])
    assert np.all(joint[3, 3, [4, 5, 7, 8]] > 0)
    assert not np.any(joint[0, 1:])
    np.testing.assert_allclose(joint.sum(axis=(1, 2)), 1)


def test_source_energy_and_two_calibrated_receivers_have_different_implications():
    def archive(name):
        return json.loads((experiment.HERE / (name + ".json")).read_bytes())

    e = archive("thermal_energy")
    assert min(e["ordinary_internal_energy"]) > max(e["isothermal_free_energy"]) > 0
    assert min(e["energy_detector"]["held"]) > 0
    lab = archive("thermal_pair")
    assert max(lab["certified_work"]["enrolled"]["2_6"]) < 0
    assert lab["single_filter_gain_upper"] < 1e-26
    ordinary = archive("thermal_ordinary")
    assert not ordinary["common_ratio"]
    assert (
        ordinary["tests"]["work_0_1"]["ratio"][1]
        < 0
        < ordinary["tests"]["release_2"]["ratio"][0]
    )
