"""Check the detector against an independent three-load operator and orbital projector."""

from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from model.candidate_laws.contrast_memory_spherical_v1.sectors import TripleSector
from model.candidate_laws.interval_echo.openwave.quark_winding.bound_face import ground_singlet, measure, routings
from model.executable_model import cell

HERE = Path(__file__).parent


def compare(actual, archived):
    if isinstance(archived, dict):
        assert actual.keys() == archived.keys()
        for k in archived:
            compare(actual[k], archived[k])
    elif isinstance(archived, list):
        for a, b in zip(actual, archived, strict=True):
            compare(a, b)
    elif isinstance(archived, float):
        assert actual == pytest.approx(archived, abs=1e-9, rel=0)
    else:
        assert actual == archived


@pytest.fixture(scope='module')
def data():
    result = json.loads((HERE / 'results.json').read_text())
    compare(json.loads(json.dumps(measure.run())), result)
    return result


def test_ground_and_flight(data):
    triple = TripleSector(cell())
    for case in data['cases']:
        start = np.zeros(triple.dim)
        start[triple.rank(*case['birth_face'])] = 1
        for run in case['runs']:
            h = triple.hamiltonian(run['kappa'])
            if run['source'] == 'ground_filter':
                _, vectors = np.linalg.eigh(h.toarray())
                psi = vectors[:, 0]
                assert 1 - run['source_failure'] == pytest.approx(abs(psi @ start) ** 2)
            else:
                psi = start if run['source'] == 'sharp' else expm_multiply(-1j * h, start)
            assert run['any_face_given_source'] == pytest.approx(sum(abs(psi[triple.faces]) ** 2))
            for port in run['target_ports_given_source']:
                p = sum(abs(psi[(triple.tri == port['port']).any(axis=1)]) ** 2)
                assert port['probes']['0']['signed_echo']['1'] == pytest.approx(p)
            for port in run['target_ports_given_birth_face']:
                assert port['probes']['0']['signed_echo']['1'] == pytest.approx(1)
            assert run['source_failure'] + sum(run['occupancy_per_launch'].values()) == pytest.approx(1)
            colour = run['colour_per_launch']
            assert colour['accepted'] == pytest.approx(run['occupancy_per_launch']['1,1,1'])
            assert sum(colour['triplet_count_weights']) == pytest.approx(colour['accepted'])
            dark = colour['ports']['centre']['dark']
            ideal = .75 * sum(colour['triplet_count_weights'][1:3])
            assert abs(dark - ideal) <= 1.5 * data['colour_calibration']['compiled_centre_error'] + 1e-9


def test_discriminating_controls(data):
    for n in range(3):
        for r in range(3):
            record = data['signed_calibration']['8.0'][f'{n},{r}']
            assert record['signed_echo'][str(n-r)] == 1
            # Exact blank-contact diagonal amplitude includes the exchange path.
            swap = 1 if n == r else (1 + np.exp(-1j * np.pi / 32)) / 2
            overlap = np.exp(-1j * 8 * (n-r)**2 * np.pi / 32) * swap
            assert record['bright_y'] == pytest.approx((1 + overlap.imag) / 2)
    controls = data['receiver_controls']
    assert all(ports['dark'] < 1e-8 for ports in controls['antisymmetric']['ports'].values())
    assert controls['three_identical']['ports']['half_turn']['dark'] > .99999
    for case in data['cases']:
        sharp = case['runs'][0]['colour_per_launch']
        assert sharp['accepted'] == 1
        assert sharp['ports']['centre']['dark'] == pytest.approx(27 / 64, abs=2e-6)


def test_all_routings():
    full = json.loads((HERE / 'routings.json').read_text())
    assert json.loads(json.dumps(routings.run())) == full
    assert full['partitions'] == 5775
    assert sum(full['face_count_histogram'].values()) == 5775
    faces = cell().faces()
    eigenvalues, vectors = np.linalg.eigh(cell().adjacency().toarray())
    triplet = vectors[:, np.isclose(eigenvalues, np.sqrt(5))]
    for f in faces:
        assert np.linalg.det(triplet[list(f)]) ** 2 / 6 == pytest.approx(
            full['common_frame_singlet_probability_per_face'])
    P = triplet @ triplet.T
    centre = np.eye(12) + np.expm1(2j * np.pi / 3) * P
    least = 1.
    maximum_singlet_return = 0.
    for groups in routings.partitions():
        slot = {p: i for i, group in enumerate(groups) for p in group}
        eligible = [sorted(f, key=slot.__getitem__) for f in faces if len({slot[p] for p in f}) == 3]
        cols = np.array(eligible).T
        rotation = np.prod([centre[np.ix_(c, c)] for c in cols], axis=0)
        effect = (2 * np.eye(len(eligible)) - rotation - rotation.conj().T) / 4
        least = min(least, np.linalg.eigvalsh(effect)[0])
        total = sum(np.linalg.det(triplet[f]) ** 2 / 6 for f in eligible)
        maximum_singlet_return = max(maximum_singlet_return, total)
    assert least >= full['worst']['bound_in_256ths'] / 256
    assert full['worst']['bound_in_256ths'] == 63
    assert maximum_singlet_return == pytest.approx(full['max_face_probability_common_frames'])
    assert maximum_singlet_return < Fraction(full['face_probability_bound_any_sharp_destination_frames'])


def test_ground_mirrors_and_control():
    full = json.loads((HERE / 'ground.json').read_text())
    compare(json.loads(json.dumps(ground_singlet.run())), full)
    car = cell()
    val, basis = np.linalg.eigh(car.adjacency().toarray())
    v = basis[:, np.isclose(val, np.sqrt(5))]
    for (a, b), i in zip(combinations(range(12), 2),
                         full['witness_index_for_each_lexicographic_port_pair'], strict=True):
        perm = full['reflections_of_triplet_trace_one'][i]
        assert (perm[a], perm[b]) == (a, b)
        assert np.linalg.det(v.T @ v[perm]) == pytest.approx(-1)
    assert all(full['common_ground_amplitudes_below_1e_minus_10'].values())
    weights = full['max_per_preparation_singlet_weight']
    assert weights['common']['ground_k8'] == 0
    assert weights['common']['flight_0'] > .0019
    assert weights['rotate_first']['ground_k8'] > .0008
