"""Verify the receiver over all reachable words and independently propagate a release."""

import json
import numpy as np
from scipy.linalg import expm
from . import outlet, outlet_clock


def receipt(name):
    return json.loads((outlet.HERE / (name + ".json")).read_text())


def test_site_hamiltonian_and_charge_instrument():
    board, dimers, path, cell, internal = outlet.outlet_geometry()
    assert path == (75, 80, 192) and set(p // 12 for p in cell) == {16}
    assert set(cell).isdisjoint(set(sum(dimers, ())))
    sectors = outlet.receiving_cell(board, cell, internal)
    assert outlet.inventory_reader(board, cell, sectors)[
        "inventory_zero_one_two_to_pointer"
    ] == [0, 1, 2]
    for closed in (False, True):
        words, h = outlet.material(board, dimers, path, cell, internal, closed)
        sites = tuple(sorted({*dimers[-2], *dimers[-1], *path, *cell}))
        # Direct full physical loads on all active ports and the occupied pin.
        configurations = [{dimers[-2][1 - d], *word} for d in (0, 1) for word in words]
        lookup = {frozenset(q): i for i, q in enumerate(configurations)}
        blank = [*internal, *dimers[-2:], (80, 192)]
        if not closed:
            blank.append((75, 80))
        direct = np.zeros_like(h)
        for i, word in enumerate(configurations):
            direct[i, i] = sum(8 * (6 - board.degree[p]) for p in sites if p in word)
            for a, b in blank:
                if (a in word) != (b in word):
                    direct[i, i] += 8.5
                    direct[lookup[frozenset(word ^ {a, b})], i] -= 0.5
            active = word | {79}
            for a, b, cost in ((74, 75, 33), (74, 76, 32), (78, 79, 32)):
                direct[i, i] += cost * (1 - int(a in active and b in active))
        np.testing.assert_allclose(
            direct - h, np.eye(len(h)) * (direct - h)[0, 0], atol=3e-14, rtol=0
        )


def test_direct_flights_and_clock_support():
    first = receipt("outlet")
    repeat = receipt("outlet_clock")
    board, dimers, path, cell, internal = outlet.outlet_geometry()
    sectors = outlet.receiving_cell(board, cell, internal)
    stocks, _ = outlet.source_stocks()
    _, one, v = sectors[1]
    gaps = one - one[0]
    assert np.count_nonzero((gaps >= 0.25) & (gaps < first["cutoff"])) == 0
    low = (gaps > 1e-8) & (gaps < first["cutoff"])
    assert low.sum() == 2 and np.ptp(gaps[low]) < 1e-12
    assert min(gaps[gaps >= first["cutoff"]]) > 6.8
    assert (v[0, low] ** 2).sum() > 0.5
    assert np.max(np.abs(np.imag(v))) < 1e-13
    # Exact energy projector on the primary source propagated by a direct exponential;
    # compare to the separately archived low-line weight, before the finite clock.
    words, h = outlet.material(board, dimers, path, cell, internal, False)
    projector = v[:, low] @ v[:, low].T
    for arm in ("pulse", "idle"):
        state = expm(-32j * h) @ outlet.inject(
            stocks["native", arm], words, dimers, cell, v[:, :1]
        )
        groups = {}
        for d in (0, 1):
            for i, word in enumerate(words):
                inside = [p for p in word if p in cell]
                if len(inside) == 1:
                    outside = next(p for p in word if p not in cell)
                    s = groups.setdefault((d, outside), np.zeros((12, 8), complex))
                    s[cell.index(inside[0])] = state[d * len(words) + i]
        true = sum(np.vdot(block, projector @ block).real for block in groups.values())
        panel = repeat["panels"]["prompt"]["settings"]["native"][arm][4]
        old = first["panels"]["native"][arm][4]
        np.testing.assert_allclose(
            true,
            np.asarray(panel["core_by_sector_energy"]["1"])[:, 1].sum(),
            atol=2e-12,
        )
        np.testing.assert_allclose(
            old["core_by_sector_energy"]["1"],
            panel["core_by_sector_energy"]["1"],
            atol=2e-12,
        )
    for label, host in repeat["panels"].items():
        assert host["ready_load_in_cell"] > 0.9999
        for name, cohorts in host["settings"].items():
            for p, i in zip(cohorts["pulse"], cohorts["idle"], strict=True):
                excess = (np.array(p["core_by_outlet"]) - i["core_by_outlet"]).sum(0)
                if name != "native":
                    np.testing.assert_allclose(excess, 0, atol=2e-12)
                else:
                    assert abs(excess[2]) < 0.01


def test_finite_clock_and_ready_ensemble():
    wave = receipt("outlet_clock")
    board, _, path, cell, edges = outlet.outlet_geometry()
    sec = outlet.receiving_cell(board, cell, edges)
    ports, hosts = outlet_clock.ready_inlet(board, path, cell, edges, sec)
    _, h = outlet.binary_generator(board, ports, [*edges, (80, 192)], 1)
    # Sum the uniform exposures by dyadic doubling of the density operator.
    rho = hosts["prompt"] @ hosts["prompt"].T
    u = expm(-2j * np.pi * h / outlet.PERIOD)
    for _ in range(int(np.log2(outlet_clock.LABELS))):
        rho = (rho + u @ rho @ u.conj().T) / 2
        u = u @ u
    np.testing.assert_allclose(
        rho, hosts["pre_exposed"] @ hosts["pre_exposed"].conj().T, atol=3e-12
    )
    assert wave["low_line_acceptance_floor"] > 0.997
    assert wave["false_light_ceiling"] < 0.00008
    assert (
        min(wave["sharp_address_clock_other_light_heavy"][cell.index(path[-1])][1:])
        > 0.02
    )


def test_steady_interface_stimulus_and_all_time_envelope():
    from . import transmission

    record = receipt("transmission")
    board, dimers, path, cell, edges = outlet.outlet_geometry()
    cells = outlet.receiving_cell(board, cell, edges)
    words, h = outlet.material(board, dimers, path, cell, edges, False)
    n = len(words)
    phase = np.diag(np.r_[np.ones(n), -np.ones(n)])
    _, v = np.linalg.eigh(outlet.hamiltonians(4)[1])
    guide = np.zeros((4, 2, 2, 2))
    guide[0, 0] = v[:, 0].reshape(2, 2)
    initial = outlet.inject(guide, words, dimers, cell, cells[1][2][:, :1])[:, 0]
    rho = np.outer(initial, initial.conj())
    u = expm(-2j * np.pi * h / outlet.PERIOD)
    for _ in range(int(np.log2(outlet_clock.LABELS))):
        rho = (rho + u @ rho @ u.conj().T) / 2
        u = u @ u
    delta = phase @ rho @ phase - rho
    native = record["settings"]["native"]
    assert abs(np.trace(h @ delta).real - native["source_energy"]) < 2e-11
    one = cells[1][1]
    response = outlet.clock(
        one,
        one[0],
        receipt("outlet_clock")["parent_gap"] / 10,
        labels=outlet_clock.LABELS,
        low_edge=outlet.PERIOD / outlet_clock.LABELS,
    )
    effects = transmission.outlet_effects(
        outlet.read_groups(words, cells, cell), cells, response
    )
    for t in (0, 32, 64):
        u = expm(-1j * t * h)
        evolved = u @ delta @ u.conj().T
        measured = np.einsum("aij,ji->a", effects, evolved).real
        idx = record["waits"].index(t)
        target = (
            np.array(native["panels"]["pulse"][idx]) - native["panels"]["idle"][idx]
        )
        np.testing.assert_allclose(measured, target, atol=2e-11, rtol=0)
    # Independent diagonalizer and density preparation reproduce the bound.
    values, v = np.linalg.eigh(np.round(2 * h) / 2)
    assert min(np.diff(values)) > 3e-9
    assert np.linalg.norm(h @ v - v * values, 2) < 1e-12
    coeff = (v.T @ effects[2] @ v) * (v.T @ delta @ v).T
    bound = np.trace(coeff).real + abs(coeff - np.diag(coeff.diagonal())).sum()
    assert abs(bound - native["all_time_upper"]) < 1e-8 and bound < 0.00445
