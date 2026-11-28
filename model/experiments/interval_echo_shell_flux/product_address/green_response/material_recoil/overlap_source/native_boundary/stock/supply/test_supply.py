"""Independent neighboring-donor contraction, sector bounds and released local support."""

from itertools import product

import numpy as np

from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock import (
    measure as native,
)
from model.experiments.interval_echo_shell_flux.product_address.green_response.material_recoil.overlap_source.native_boundary.stock.supply import (
    measure,
    module,
)


def test_neighboring_fueled_donors_and_coherent_address_selection():
    # Two six-port modules linked by an extra real blank, read with all
    # eleven native terms. This checks the cross-donor term, which an
    # isolated one-source transfer cannot test.
    edges = module.EDGES + [(a + 6, b + 6) for a, b in module.EDGES] + [(3, 8)]
    angles = [1, 0.5j, 0.3 - 0.2j, -0.8]
    psi = np.array(angles) / np.linalg.norm(angles)
    words = []
    for hot in [(2,), (3,), (8,), (9,)]:
        word = [0] * (12 + len(edges))
        for v in hot:
            word[v] = 1
        words.append(tuple(word))
    two = {}
    for left, right in product((2, 3), (8, 9)):
        word = [0] * (12 + len(edges))
        word[left] = word[right] = 1
        two[tuple(word)] = 0.5
    _, _, h4, _ = native.kernels()
    h = h4 / 4
    e, u = np.linalg.eigh(h)
    swap = (u * np.exp(-1j * np.pi * e)) @ u.T
    echo = native.LAW.echo()
    donor = np.array([[0, 1], [2, 3]])
    other = np.array([[1, 2]])
    expected = measure.survey(
        psi, donor, other, np.array([0, 1]), np.random.default_rng(917)
    )
    expected_many = measure.many([True, True], donor, other)

    def apply(state, gate, i):
        return module.apply(state, gate, i, edges=edges, ports=12)

    for state, theoretical in [
        (dict(zip(words, psi, strict=True)), expected),
        (two, expected_many),
    ]:
        for batch, value in zip(
            [
                [],
                [(echo, 0), (echo, 5)],
                [(swap, i) for i in [1, 2, 6, 7]],
                [(swap, i) for i in [3, 4, 8, 9]],
                [(echo, i) for i in [1, 2, 6, 7]],
            ],
            theoretical["stage_energies"],
            strict=True,
        ):
            for gate, i in batch:
                state = apply(state, gate, i)
            energy = sum(
                module.inner(state, apply(state, h, i)) for i in range(len(edges))
            )
            np.testing.assert_allclose(energy, value, atol=1e-12)
        for word in state:
            assert sum(word[:12]) == sum(native.MAG[word[12 + i]] for i in [0, 5])
            assert sum(word[i] for i in (2, 3, 8, 9)) == 0
        # Once released each term's component that exits the sector of
        # parked ordinary charges and empty core vanishes on these states.
        for i in range(len(edges)):
            evolved = apply(state, h, i)
            assert all(not any(word[j] for j in (2, 3, 8, 9)) for word in evolved)
    assert np.isclose(expected_many["stage_energies"][0], 37.25)


def test_three_physical_profiles_and_bounds():
    r = measure.run()
    assert r["transferred_even_odd_eigenvalues"] == [8, 9]
    assert r["ordinary_parity_readout_checked"]
    assert (
        r["single_load_unmodified_blank_energy_interval"][1]
        < r["one_separated_empty_mark_plus_original_rest_energy_floor"]
    )
    assert [b["reachable"] for b in r["boards"]] == [64, 8, 8]
    for board in r["boards"]:
        assert board["full_native_degree_max"] == 6
        assert board["private_reference_analyzers"] == 2 * board["owners"]
        assert board["unaddressed_no_capture_unit_marks"] == 0
        assert board["many"]["whole_board"]["unit_marks"] == 2 * board["owners"]
        for t in board["traces"]:
            assert 8 <= t["native"]["record_mean"] <= 9
            assert (
                sum(np.asarray(t["native"]["address_energy_counts"]).ravel())
                == measure.COUNT
            )
            assert t["zero_exchange"]["radius_first"] == 0
            assert t["native"]["net_work"] < 0


def test_charge_from_cold_and_hot_ordinary_dimer():
    m = module.native_microscope()
    assert m["six_port_stage_energies"]["vacuum"] == [0] * 5
    assert m["six_port_stage_energies"]["even"] == [16.5, 8.5, 25, 16.5, 16.5]
    assert m["six_port_stage_energies"]["odd"] == [17.5, 8.5, 26, 17.5, 17.5]
