import numpy as np

from model.candidate_laws.interval_echo import candidate
from model.executable_model import cell, federation
from .experiment import dark_bands, gateway_table, geometry, port_reader, prepare


def test_four_port_source_and_local_occupation_reader():
    local, axes = geometry()
    car = federation(0)
    tips = gateway_table(car)[0]
    available = [pair for pair in axes if set(tips).isdisjoint(pair)]
    target = np.zeros(12)
    target[list(available[0])] = 0.5
    target[list(available[1])] = -0.5
    source, _ = prepare(local, tips, target)
    for port in range(12):
        _, effect = port_reader(local, port)
        np.testing.assert_allclose(
            effect @ abs(source) ** 2, target[port] ** 2, atol=1e-12
        )
    action = candidate.generator(car)[:, :12] @ source
    np.testing.assert_allclose(action[:12], 51 * source, atol=1e-12)
    np.testing.assert_allclose(action[12:], 0, atol=1e-12)


def test_internal_beat_requires_two_distinct_dark_energies():
    local_h = candidate.generator(cell()).toarray()
    rich = dark_bands(local_h, [0, 1, 2])
    ordinary = dark_bands(local_h, [1, 3, 10])
    assert len(ordinary) == 1
    assert len(rich) == 3
    np.testing.assert_allclose(rich[-1][0] - rich[0][0], np.sqrt(5), atol=1e-12)


def test_anchored_vacancy_additivity_and_opposite_sign_control():
    from .anchored_vacancy import run

    result = run()
    assert result["source"]["minimum_port_distance"] == 12
    for reference, expected in (("0", (35, 67, 70, 102)), ("1", (51, 51, 102, 102))):
        panel = result["records_by_degree_reference"][reference]
        for key, energy in zip(
            ("one_vacancy", "one_excess", "two_vacancies", "neutral_pair"),
            expected,
            strict=True,
        ):
            assert abs(panel[key]["ramsey_estimate"] - energy) < 1e-5
            assert panel[key]["eigen_residual"] < 1e-12
        assert panel["adjacent_sharp_two_vacancies_minus_singles"] == -17


def test_same_load_material_clock_reads_transfer_and_resolves_changed_pair_interaction():
    from . import material_reference as material

    result = material.panel(federation(0), 0, np.random.default_rng(178031))
    for degree, expected in (("5", (43, 43, 86, 86)), ("6", (35, 67, 70, 102))):
        for name, freq in zip(
            ("vacancy", "excess", "two_vacancies", "opposite_signs"),
            expected,
            strict=True,
        ):
            record = result["panels"][degree][name]
            counts = np.array(record["detector_counts"])
            assert (counts.sum(axis=-1) == material.SHOTS).all()
            signs = np.array([1, -1, 0])
            weights = signs if len(counts) == 1 else np.kron(signs, signs)
            data = counts[:, :, :-1] @ weights / material.SHOTS
            signal = data[0] if len(counts) == 1 else data[0] - data[1]
            recovered = material.fit(signal, record["timer_swaps"])[
                "frequency_per_clock_gap"
            ]
            assert abs(recovered - freq) < 0.08
            assert record["max_probability_change_under_sector_shift"] < 2e-13
    changed = result["separated_interaction_control"]
    assert (
        changed["frequency_per_clock_gap"] - 70
        > 6 * changed["frequency_standard_error"]
    )
    assert (
        abs(result["local_phase_flip"]["-1"]["mean_energy_change_from_counts"] + 2 / 3)
        < 0.004
    )
    assert result["all_single_excess_formation_energy_bounds"] == [48.0, 69.0]


def test_released_native_vacancy_has_no_anchored_recoil_line():
    from . import localized_release

    record = localized_release.run()
    native = record["studies"]["registered"]
    for board in record["boards"]:
        lower, upper = board["degree_block_enclosures"]
        assert lower[1] < 35 * np.sqrt(2) < upper[0]
        assert lower[1] < board["vacancy_floor_lower_bound"] * np.sqrt(2)
        assert board["vacancy_floor_trial_upper_bound"] * np.sqrt(2) < upper[0]
    assert native["number_of_eigenvalues_in_one_percent_window"] == 0
    assert native["rest_energy_component_outside_home"] < 1e-22
    for setting in native["flight"]:
        assert sum(setting["counts"]) == record["shots"]
        if not setting["pulse"]:
            assert setting["counts"][1] == 0
    assert native["flight"][-1]["counts"][1] > 0
    assert record["studies"]["half_contrast"]["kicked_probability_in_one_percent_window"] > 0
