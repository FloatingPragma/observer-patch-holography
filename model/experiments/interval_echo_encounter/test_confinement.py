"""Independent global-sector tests of local closure, sign read and controls."""
from copy import deepcopy
from itertools import product
import json

import numpy as np
import pytest
from scipy.linalg import expm

from model.candidate_laws.contrast_memory_spherical_v1 import law
from model.candidate_laws.contrast_memory_spherical_v1.sectors import PairSector
from model.executable_model.carrier import Carrier, federation
from . import check_transport, confinement, conjugation, sign_read, transport


def test_exact_spaces_and_cell_counts():
    actual = confinement.certificate()
    assert actual == json.loads(confinement.RECEIPT.read_text())
    compact = [sum(profile["cells_by_level"].get(str(level), 0) *
                   profile["windows"]["single"]["0"]["dimension"]
                   for profile in actual["profiles"]) for level in range(3)]
    assert compact == [40, 232, 952]
    assert actual["distinct_cell_product_dimensions"] == {"0": 760, "1": 26500, "2": 451420}
    assert actual["quiet_quadratic_coefficients"] == {"home": "1/4", "seams": "5/4"}
    for p in actual["profiles"]:
        for window in ("home", "seams"):
            assert [case["dimension"] for case in p["windows"][window].values()] == [0]*3
    # An isolated cell has a wholly invariant home space. A solver always
    # returning zero would pass every connected two-load negative case.
    _, matrix, outward = confinement.blocks((), "home", 8)
    _, basis = confinement.maximal_space(matrix, outward.values())
    np.testing.assert_array_equal(basis, np.eye(66))


def test_global_sector_includes_every_exit_and_diagonal():
    car = federation(0)
    pair = PairSector(car)
    h = 2 * pair.hamiltonian(8)
    cell = 13
    gates = confinement.gate_profiles(car)[cell]
    outside = [next(n for n in car.neighbours[12*cell+g] if n//12 != cell) for g in gates]
    address = list(range(12*cell, 12*cell+12)) + outside
    for window in ["home", "seams"]:
        states, inside, exits = confinement.blocks(gates, window, 8)
        index = [pair.index(address[a], address[b]) for a,b in states]
        columns = h[:, index]
        np.testing.assert_array_equal(columns[index].toarray(), inside)
        for (a,b), row in exits.items():
            np.testing.assert_array_equal(columns[pair.index(address[a],address[b])].toarray()[0], row)
        assert columns.count_nonzero() == np.count_nonzero(inside) + np.count_nonzero(list(exits.values()))


def test_maximality_by_independent_modular_observability_gramian():
    # A nonzero minor modulo a prime also survives over the rationals.
    # The rational invariant witnesses supply the opposite dimension bound.
    prime = 1009
    for profile in json.loads(confinement.RECEIPT.read_text())["profiles"]:
        for window, cases in profile["windows"].items():
            for strength, case in cases.items():
                _, a, exits = confinement.blocks(profile["gates"], window, int(strength))
                rows = np.array(list(exits.values()))
                gram = np.zeros_like(a)
                for _ in range(len(a)):
                    gram = (gram + rows.T @ rows) % prime
                    rows = (rows @ a) % prime
                rank = 0
                for col in range(len(a)):
                    pivots = np.flatnonzero(gram[rank:,col])
                    if not pivots.size:
                        continue
                    pivot = rank + pivots[0]
                    gram[[rank,pivot]] = gram[[pivot,rank]]
                    factor = pow(int(gram[rank,col]), -1, prime)
                    gram[rank] = (gram[rank]*factor) % prime
                    gram[rank+1:] = (gram[rank+1:] -
                                     gram[rank+1:,col,None]*gram[rank]) % prime
                    rank += 1
                assert rank == len(a) - case["dimension"]


def test_distinct_neighbor_cell_product_in_global_hamiltonian():
    car = federation(0)
    pair = PairSector(car)
    i,j = (p//12 for p in car.external[0])
    profiles = confinement.gate_profiles(car)
    modes = []
    for cell in (i,j):
        _, a, out = confinement.blocks(profiles[cell], "single", 0)
        _, v = confinement.maximal_space(a, out.values())
        modes.append(v.astype(float))
    factors = np.kron(*modes)
    indices = [pair.index(12*i+a,12*j+b) for a in range(12) for b in range(12)]
    columns = np.zeros((pair.P, factors.shape[1]))
    columns[indices] = factors
    for strength in (-8,0,8):
        np.testing.assert_array_equal(pair.hamiltonian(strength) @ columns, 6*columns)


def test_local_pulses_on_full_native_collar_and_registered_sign_read():
    k, v = law.local_term_parts()
    binary = [law.SRC.index((a,b)) for a in range(2) for b in range(2)]
    for strength in (8,0,-8):
        native = expm(-0.5j*np.pi*(k+strength*v))
        reduced = expm(-0.5j*np.pi*law.seam_term(strength))
        np.testing.assert_allclose(native[:, binary], np.eye(law.DIM,9) @ reduced[:,binary], atol=2e-13)
    result = sign_read.campaign()
    expected = json.loads(sign_read.RECEIPT.read_text())
    sign_read.verify(result, expected)
    for strength, arm in result["arms"].items():
        assert arm["slopes"] == [int(strength)/2, -1/8]
    def readings(strength):
        return np.array(result["arms"][strength]["readings"]["0.0625"])[:,0]
    coherent, dephased = readings("8") - readings("-8")
    assert coherent > 0.10
    assert abs(dephased) < 0.001
    corrupted = deepcopy(expected)
    corrupted["arms"]["8"]["readings"]["0.0625"][0][-1] += 0.01
    with pytest.raises(ValueError):
        sign_read.verify(result, corrupted)


def test_published_release_and_complementary_outcomes():
    data = json.loads(transport.RECEIPT.read_text())
    # Sums of 30 rounded probabilities in the M10.3 deposited receipt.
    published = {"8": [0.990982505,0.989653268,0.980271946],
                 "0": [0.450001679,0.275768631,0.202905801]}
    for strength, known in published.items():
        hits = [sum(data["arms"][strength][t][:30]) for t in ("1","7","37")]
        np.testing.assert_allclose(hits, known, rtol=0, atol=2e-8)
    for readings in data["arms"].values():
        for bins in readings.values():
            assert len(bins) == 34 and min(bins) >= 0 and abs(sum(bins)-1) < 2e-11
    corrupt = deepcopy(data)
    corrupt["arms"]["-8"]["148"].pop()
    with pytest.raises(ValueError):
        transport.verify(data, corrupt)


def test_load_conjugation_on_all_qutrit_states_and_terminal_records():
    flip = np.eye(9)[::-1]
    seam = law.seam_term(8)
    np.testing.assert_array_equal(flip @ seam @ flip, seam)
    for ref, (a,b) in product((0,1), law.SRC):
        # One missing contact at a, none at b. The identity extends by summation.
        difference = 8*((2-a-ref)**2-(a-ref)**2)
        assert difference == 32*(1-ref)*(1-a)
    result = conjugation.campaign()
    conjugation.verify(result, json.loads(conjugation.RECEIPT.read_text()))
    for board in result.values():
        port = board["ports"][0] % 12
        native, exchange, centered, singleton = board["arms"]
        assert native["excess"]["slope"] - native["vacancy"]["slope"] == 8
        assert (native["excess"]["readings"]["0.0625"][port] -
                native["vacancy"]["readings"]["0.0625"][port]) > 0.10
        for law_control in (exchange, centered):
            assert law_control["excess"] == law_control["vacancy"]
        assert singleton["excess"] == singleton["vacancy"] == native["vacancy"]
    corrupted = deepcopy(result)
    corrupted["2"]["arms"][-1]["vacancy"]["readings"]["0.125"][-1] += 0.001
    with pytest.raises(ValueError):
        conjugation.verify(result, corrupted)


def test_singleton_reference_preserves_binary_and_changes_doubled_sector():
    profile = conjugation.singleton_profile()
    assert profile == (0,1,0) and profile == profile[::-1]
    car = federation(0)
    missing = 6 - car.degree
    pair = PairSector(car)
    # Compare every matrix entry on the published release board. One assembly
    # walks each physical seam; the native PairSector walks configuration neighbors.
    native_pair = pair.hamiltonian(8)
    assert (conjugation.binary_generator(pair, 8, "singleton") - native_pair).count_nonzero() == 0
    centered = conjugation.binary_generator(pair, 8, 1)
    difference = centered - native_pair
    np.testing.assert_array_equal(difference.diagonal(), -16*(missing[pair.aa] + missing[pair.bb]))
    assert difference.count_nonzero() == np.count_nonzero(difference.diagonal())
    # Every two-single configuration, including cross-cell and separated pairs.
    delta_binary = (profile[1]-profile[0]-1)*(missing[pair.aa]+missing[pair.bb])
    assert np.all(delta_binary == 0)
    # Every sharp doubled-load state. Only its missing-seam terms change.
    delta_double = (profile[2]-profile[0]-4)*8*missing
    assert set(delta_double.tolist()) == {-32, 0}
    from scipy.sparse import diags, eye
    from model.candidate_laws.contrast_memory_spherical_v1.sectors import doubled_hamiltonian
    native = conjugation.defect_generator(car, 2, 8, 0, background=0)
    modified = conjugation.defect_generator(car, 2, 8, "singleton", background=0)
    assert (native - doubled_hamiltonian(car,8) - 96*eye(car.n)).count_nonzero() == 0
    assert (modified - native - diags(delta_double, dtype=float)).count_nonzero() == 0


def test_full_qutrit_path_including_sector_complements():
    car = Carrier("three-port path", 3, [(0, 1), (1, 2)])
    states = list(product(range(3), repeat=3))
    seam = law.seam_term(8)
    physical = np.kron(seam, np.eye(3)) + np.kron(np.eye(3), seam)
    binary = [i for i, state in enumerate(states) if max(state) <= 1]
    native = None
    for reference, profile in ((0, (0, 1, 4)), (1, (1, 0, 1)), ("singleton", (0, 1, 0))):
        correction = [sum((6-d) * profile[n] for d, n in zip(car.degree, state)) for state in states]
        full = physical + 8*np.diag(correction)
        if reference == 0:
            native = full
            assert not np.array_equal(full[::-1, ::-1], full)
        else:
            np.testing.assert_array_equal(full[::-1, ::-1], full)
        if reference == "singleton":
            np.testing.assert_array_equal(full[:, binary], native[:, binary])
        for background, change in ((1, 1), (1, -1), (0, 2)):
            basis = [states.index(tuple(background + change*(i == j) for j in range(3)))
                     for i in range(3)]
            embedding = np.eye(len(states))[:, basis]
            medium = states.index((background,)*3)
            expected = full @ embedding - full[medium, medium]*embedding
            reduced = conjugation.defect_generator(car, change, 8, reference, background)
            np.testing.assert_array_equal(embedding @ reduced.toarray(), expected)


def test_probability_and_registration_verifiers_reject_invalid_inputs(tmp_path):
    for vector in ([True, 0.0], [-1e-12, 1.0], [0.0, float("nan")], [0.6, 0.6], [0.0]):
        with pytest.raises(ValueError):
            transport.probability_vector(vector, 2)
    path = tmp_path / transport.PROTOCOL.name
    original = transport.PROTOCOL.read_bytes().replace(b"\r\n", b"\n")
    path.write_bytes(original.replace(b"\n", b"\r\n"))
    assert transport.read_protocol(path) == transport.read_protocol(transport.PROTOCOL)
    path.write_bytes(original.replace(b"148", b"147"))
    with pytest.raises(ValueError, match="Registered protocol"):
        transport.read_protocol(path)


def test_terminal_decoder_resolves_each_address_and_each_quiet_alternative():
    sector = PairSector(federation(0))
    internal, _ = transport.release(sector, 13)
    labels = check_transport.terminal_labels(sector, 13, internal)
    assert set(labels) == set(range(34))
    for label in range(34):
        state = np.zeros(sector.P)
        state[np.flatnonzero(labels == label)[0]] = 1
        expected = np.zeros(34)
        expected[label] = 1
        np.testing.assert_array_equal(transport.census(sector, 13, internal, state), expected)


def test_centered_release_against_complete_archived_first_read():
    from scipy.sparse.linalg import expm_multiply
    pair = PairSector(federation(0))
    internal, source = transport.release(pair, 13)
    start = np.zeros(pair.P, complex)
    start[pair.index(*source)] = 1
    final = expm_multiply(-1j * conjugation.binary_generator(pair, 8, 1), start)
    bins = transport.census(pair, 13, internal, final)
    archived = json.loads((transport.HERE / "centered_pair.json").read_text())
    np.testing.assert_allclose(bins, transport.probability_vector(archived["1"], 34), rtol=0, atol=1e-10)
