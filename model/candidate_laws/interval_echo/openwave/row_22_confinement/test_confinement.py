"""Independent term-by-term checks on bit words, trace subtraction, and raw records."""

import json
from itertools import product

import numpy as np

from model.candidate_laws.contrast_memory_spherical_v1.sectors import BinarySector
from model.candidate_laws.interval_echo import currents
from model.executable_model import carrier

from . import measure, spectrum


def records():
    return json.loads((measure.HERE / 'records.json').read_text())


def test_full_seams_on_a_coherent_three_load_wave():
    board = carrier.federation(0)
    graph, ports = measure.skeleton(board)
    data = records()['levels']['0']['colour']
    b, c = data['antitriplet_owners']
    a = data['rows'][0]['owner']
    owners = (a, b, c)
    band = spectrum.physical_frame(currents._service(carrier.cell(), False, 1.)[4][0])
    frames = (measure.aligned(band, ports[a, b] % 12),
              measure.aligned(band, ports[b, a] % 12), band)
    rng = np.random.default_rng(771)
    z = rng.normal(size=27) + 1j*rng.normal(size=27)
    z /= np.linalg.norm(z)
    wave = np.einsum('pi,qj,rk,ijk->pqr', *frames, z.reshape((3,)*3)).ravel()
    words = [tuple(12*x+n for x,n in zip(owners, w, strict=True))
             for w in product(range(12), repeat=3)]
    amplitudes = {frozenset(w): v for w, v in zip(words, wave, strict=True)}
    # Apply native swaps to full bit words, including moves to unoccupied cells.
    for k in (8, 0, -8):
        diagonal = off = 0.
        for word, v in amplitudes.items():
            for p in word:
                for q in board.neighbours[p]:
                    if q in word:
                        diagonal -= k*abs(v)**2
                    else:
                        diagonal += .5*abs(v)**2
                        moved = (word - {p}) | {q}
                        off -= .5*(v.conjugate()*amplitudes.get(moved, 0)).real
        base, contact = spectrum.fixed_source_matrices(graph, ports, owners, frames, k)
        np.testing.assert_allclose(diagonal + off, np.vdot(z, (base+contact) @ z), atol=1e-10)


def test_relaxed_source_spectra_and_raw_counts():
    r = records()
    for level in r['levels'].values():
        rows = level['colour']['rows']
        assert rows[0]['fixed_source_spectrum']['all']['operator_norm'] == 17/32
        for row in rows[1:]:
            assert all(v['ground_shift'] == v['operator_norm'] == 0
                       for v in row['fixed_source_spectrum'].values())
        for row in rows:
            for label, prep in row['preparations'].items():
                assert prep['reader']['BC_antisymmetric'] == r['colour_shots']
                assert prep['reader']['joint_determinant'] == (label == 'singlet')*r['colour_shots']
                terms = np.array(prep['odd_both_by_seam_and_arm'])
                linear = terms[:,:,0]-16*terms[:,:,1]
                measured = (linear @ [1,-1,-1,1]).sum() / r['colour_shots']
                np.testing.assert_allclose(prep['connected_energy']['8']['measured'], measured, atol=1e-11)


def test_full_guide_words_and_physical_spectral_bound():
    for lev, experiment in records()['levels'].items():
        board = carrier.federation(int(lev))
        med = experiment['medium']
        for mode, rows in med['preparations'].items():
            base = set(med['aligned_ports'])
            if mode != 'removed_spectators':
                base.update(med['spectators'])
            def energy(word):
                degree = sum(board.degree[p] for p in word)
                adjacent = sum(i in word and j in word for i,j in board.seams)
                return degree/2 - 17*adjacent + 48*len(word)
            reference = energy(base)
            for row in rows:
                m = row['length']
                moved = {0,m-1} if mode == 'outer_records_only' else set(range(m))
                word = (base - {med['aligned_ports'][i+1] for i in moved}) | {
                    med['forward_ports'][i] for i in moved}
                assert energy(word)-reference == row['energy']['8']['exact']
    cell = BinarySector(carrier.cell(), 3)
    eigenvalues = np.linalg.eigvalsh(cell.hamiltonian(8)) + 144
    assert min(eigenvalues) >= 96-1e-10 and max(eigenvalues) <= 162+1e-10


def test_controlled_native_pulses():
    rec = records()
    for k, arms in rec['pulse_controls'].items():
        for name, energy in (('odd', int(k)+1), ('even', int(k)), ('both', 0)):
            phase = np.exp(-1j*np.pi*energy/32)
            expected = [(1+phase.real)/2, (1+phase.imag)/2, 1]
            np.testing.assert_allclose(arms[name]['probabilities'], expected, atol=1e-10)
            counts = np.array(arms[name]['x_y_blank']) / rec['phase_shots']
            assert np.max(abs(counts-expected)) < .003


def test_entangled_dressing_against_native_four_load_swaps():
    from . import dressing

    board = carrier.federation(0)
    graph, ports = measure.skeleton(board)
    layouts = records()['levels']['0']['colour']
    b, c = layouts['antitriplet_owners']
    a = layouts['rows'][0]['owner']
    owners = a, b, c
    band = spectrum.physical_frame(currents._service(carrier.cell(), False, 1.)[4][0])
    frames = measure.aligned(band, ports[a, b] % 12), measure.aligned(band, ports[b, a] % 12), band
    mobile, live, bridges = dressing.geometry(board, graph, ports, owners)
    h = dressing.full_source_h(graph, ports, owners, frames, mobile, bridges, 8, 'ABC')
    source = sum(spectrum.fixed_source_matrices(graph, ports, owners, frames, 8))
    embed = dressing.source_embedding('ABC')
    offset = np.linalg.eigvalsh(embed.T @ source @ embed)[0]
    port = int(live[bridges[0][2]])
    other = next(q for q in board.neighbours[port] if q//12 == port//12)
    rng = np.random.default_rng(118)
    coeff = np.zeros((9, len(live)))
    coeff[:, np.isin(live, (port, other))] = rng.normal(size=(9, 2))
    coeff /= np.linalg.norm(coeff)
    amplitudes = {}
    for m in (port, other):
        color = (embed @ coeff[:, np.where(live == m)[0][0]]).reshape(3, 3, 3)
        amplitudes_for_source = np.einsum('pi,qj,rk,ijk->pqr', *frames, color)
        for idx in product(range(12), repeat=3):
            word = frozenset([m]+[12*o+j for o, j in zip(owners, idx, strict=True)])
            amplitudes[word] = amplitudes_for_source[idx]
    native = 0.
    for word, amp in amplitudes.items():
        for p in word:
            for q in board.neighbours[p]:
                if q in word:
                    native -= 8*amp**2
                else:
                    native += .5*amp*(amp-amplitudes.get((word-{p})|{q}, 0.))
    np.testing.assert_allclose(native, coeff.ravel() @ (h @ coeff.ravel()) + offset, atol=1e-10)


def test_dressing_receivers_and_binding_budget():
    data = json.loads((measure.HERE / 'dressing.json').read_text())
    for rows in data['levels'].values():
        for row in rows:
            for trial in row['trials'].values():
                for label, arm in trial['arms'].items():
                    x, y = 2*np.array(arm['phase_read']['x_y_counts'])/data['shots']-1
                    np.testing.assert_allclose(-np.arctan2(y,x), arm['phase_read']['energy'], atol=1e-10)
                    assert arm['residual_upper_bound'] <= 1e-8, label
                for key in ('connected', 'beyond_best'):
                    error = trial['observer'][key+'_se']
                    assert abs(trial['exact'][key] - trial['observer'][key]) < 6*error
    rows = data['levels']['1']
    assert -.110 < rows[1]['trials']['singlet_8']['exact']['beyond_best'] < -.109
    assert abs(rows[-1]['trials']['singlet_8']['exact']['beyond_best']) < 1e-8
    assert rows[-1]['trials']['singlet_8']['exact']['connected'] > .07
