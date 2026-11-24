import json
import numpy as np

from . import dispersion, probe


def test_signed_contact_replay_and_complete_enrollments():
    packet = json.loads((probe.HERE / "probe.json").read_text())
    for study in packet["studies"]:
        curves = probe.instrument(study["dwell"], study["strength_at_measured_contact"])
        for sign in range(2):
            np.testing.assert_allclose(
                curves[sign], curves[1 - sign, ::-1, ::-1], atol=3e-14
            )
        cal = study["calibration_counts"]
        unit = (
            packet["studies"][0]["calibration_counts"]
            if study["strength_at_measured_contact"] == 0
            else cal
        )
        for panel in study["panels"]:
            for settings in panel["records"].values():
                for row in settings["probes"]:
                    assert (
                        row["filter_rejected"]
                        + row["no_arrival"]
                        + np.sum(row["terminal_counts"])
                        == probe.N
                    )
                    assert all(c % 5 == 2 for c in row["availability_codes"])
            assert probe.decode(panel["records"], cal, unit) == panel["read"]
    for row in packet["studies"][0]["panels"][:3]:
        read = row["read"]
        assert abs(read["whole_odd_even"][0]) + 5 * read["sampling_se"][0] < 0.1
        assert read["whole_odd_even"][1] > 0.19
        assert read["classes"]["5"]["odd_even_in_probe_units"][0] > 0.8
        assert read["classes"]["6"]["odd_even_in_probe_units"][0] < -0.8


def test_graph_wavelength_and_timing_forecast():
    packet = json.loads((dispersion.HERE / "dispersion.json").read_text())
    assert packet["rest_cut_inertia"]["upper_lines_below"] == 2
    for row in packet["sources"]:
        for shell in row["shells"]:
            assert (
                shell["rest_preparation_failed"] + shell["passed_rest"]
                == dispersion.SHOTS
            )
            assert (
                shell["grating_or_final_energy_filter_failed"]
                + sum(shell["elastic_counts"])
                + sum(shell["inelastic_counts"])
                == shell["passed_rest"]
            )
            assert dispersion.read(shell["inelastic_counts"]) == shell["read"]
        assert dispersion.fit(row["shells"]) == row["fit"]
        assert not any(r["passes"] for r in row["fit"]["held"])
        assert row["shells"][1]["any_three_bins_upper_99_with_numerical_error"] < 0.5
    assert all(r["passes"] for r in packet["scalar_reader_control"]["fit"]["held"])
