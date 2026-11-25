import json
from fractions import Fraction
import numpy as np

from . import measure as run


def test_energy_resolution_reconstructs_original_unfiltered_channel_and_reference():
    data = json.loads((run.HERE / "results.json").read_text())
    arrays = np.load(run.HERE / "lines.npz")
    original = json.loads((run.trimer.HERE / "results.json").read_text())
    for sign, label in [(1, "matched"), (-1, "matched_odd")]:
        e = arrays[f"energies_{sign}"]
        effects = arrays[f"source_weights_{sign}"]
        mats = arrays[f"amplitudes_{sign}"]
        np.testing.assert_allclose(effects.sum(axis=0), np.eye(2), atol=2e-12)
        np.testing.assert_allclose(mats.sum(axis=0), 0, atol=1e-12)
        # With one exact energy every coin-resolved matrix has the same
        # source bra; summing the positive separable Choi matrices remains PPT.
        dephased = np.zeros((4, 4))
        for effect, mat in zip(effects, mats, strict=True):
            assert abs(np.linalg.det(effect)) < 1e-13
            null = np.linalg.eigh(effect)[1][:, 0]
            assert np.linalg.norm(mat @ null) < 1e-11
            bell = mat.transpose(0, 2, 1).reshape(8, 4) / np.sqrt(2)
            dephased += bell.T @ bell
        partial = dephased.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
        assert np.linalg.eigvalsh(partial)[0] > -1e-12
        all_lines = np.einsum("k,kcij->cij", np.exp(-1j * run.FLIGHT * e), mats)
        mean, pt = run.process(all_lines)
        old = original["settings"][label]["reads"][str(run.FLIGHT)]
        np.testing.assert_allclose(
            [mean, pt],
            [
                old["mixed_code_acceptance"],
                old["choi_partial_transpose_eigenvalues"][0],
            ],
            atol=2e-8,
        )
    assert data["cross_parity_mod_p_gcd_ascending"] == [1]


def test_certificates_current_and_independent_spectral_beat():
    fresh, _ = run.run()
    rec = fresh["cases"]["1"]
    for sign, det in [("1", 59318), ("-1", 19952)]:
        r = fresh["cases"][sign]
        assert r["cyclicity"]["det_mod_prime"] == det
        assert (
            r["cyclicity"]["second_rail_det_mod_prime"]
            == {"1": 36521, "-1": 25337}[sign]
        )
        assert r["cyclicity"]["dimension"] == 168
        assert float(Fraction(r["ordered_certificate"]["ordered_error_bound"])) < 3e-10
    lo, hi = rec["main_gap_inclusion"]
    assert 0 < lo < hi and hi - lo < 1e-9
    raw = rec["count_spectroscopy"]
    assert np.all(np.sum(raw["first_and_second_tip_else"], axis=1) == run.SHOTS)
    for f in (raw["first_half_frequency"], raw["held_half_frequency"]):
        assert abs(f - (lo + hi) / 2) < 1.1e-9
    assert (
        raw["max_observed_residual"]
        < raw["other_tone_absolute_bound"] + raw["simultaneous_99_shot_bound"]
    )
    current = rec["integrated_boundary_current"]
    assert current[0]["observed_net_charge"] == 0
    assert current[-1]["observed_net_charge"] < -0.88
    assert all(sum(r["signed_load_bins"]) == run.SHOTS for r in current)
    assert abs(fresh["dimer_clock"]["angular_frequency"] - 1) < 1e-6
