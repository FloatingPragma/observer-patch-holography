import itertools
import json
import math

import numpy as np
import pytest
from scipy import sparse as sp
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from ... import contrast_replay
from ...contrast_memory_spherical_v1 import collar, law
from .. import compiler
from ..cage import construction as cage, instruments
from ..instruments import prepare_source
from . import __main__ as replay, pair
from .law import collar_states, interaction, memory_write, ramp, ramsey


def test_native_factors_and_visible_internal_levels():
    x, square = law.local_term_parts()
    reverse = contrast_replay.seam_orientation_chart(True)
    np.testing.assert_allclose(square[np.ix_(reverse, reverse)], square, atol=2e-14)
    np.testing.assert_allclose(square[:9, :9], np.diag([(a-b)**2 for a, b in law.SRC]), atol=2e-14)
    np.testing.assert_allclose(square[9:, :9], 0, atol=2e-14)
    for sign in (1, -1):
        arms = collar_states(sign)
        np.testing.assert_allclose(square @ arms, arms @ np.diag([0, 1]), atol=2e-14)
        np.testing.assert_allclose(x @ arms, arms @ np.diag([0, (1-sign)/2]), atol=2e-14)
    car = cage.cell()
    Q = car['Q']
    L = 10 * np.eye(12) - Q @ Q.T
    expected = np.repeat([0., 5-math.sqrt(5), 6., 5+math.sqrt(5)], [1, 3, 5, 3])
    np.testing.assert_allclose(np.linalg.eigvalsh(L), expected, atol=1e-13)
    _, _, complete = collar.build(1, 8)
    native = complete[::collar.M, ::collar.M]
    np.testing.assert_allclose(np.linalg.eigvalsh(native), 48 + expected/2, atol=2e-13)
    # In scalar vacuum the contact vanishes even when the load contrast is nonzero.
    g = compiler.build()
    np.testing.assert_allclose(interaction(g, (0, 2), 1).diagonal().sum(), 1)


def test_rank_one_contact_and_all_ring_placements():
    car = cage.cell()
    edges = car['edges']
    rings = cage.rings(car)
    incidence = np.array([cage.ring_profile(car, r) != 0 for r in rings])
    assert incidence.shape == (30, 30)
    np.testing.assert_array_equal(incidence.sum(axis=0), np.full(30, 4))
    constraints = []
    for mask in incidence:
        support = np.flatnonzero(mask)
        constraints.extend(np.eye(len(edges))[support[1:]] - np.eye(len(edges))[support[0]])
    assert cage.rank(np.asarray(constraints, int).tolist()) == 29
    D = cage.law(car)
    for ring, mask in zip(rings, incidence):
        modes, _, _ = cage.modes(car, cage.ring_profile(car, ring))
        # Columns of this matrix are the off-pair coupling for each of the 30 placements.
        received = np.eye(len(D))[:, :len(edges)] * modes[:len(edges), 0]
        spill = received - modes @ (modes.T @ received)
        np.testing.assert_allclose((spill**2).sum(axis=0), mask * (3/16), atol=1e-14)
    marked = car['index'][(0, 2)]
    for g, multiplicity in ((0, 18), (1, 17), (-1, 17)):
        H = D.copy()
        H[marked, marked] += g
        values = np.linalg.eigvalsh(H)
        for level in (-math.sqrt(6), math.sqrt(6)):
            assert np.count_nonzero(abs(values-level) < 1e-11) == multiplicity


def test_ramsey_against_the_full_native_memory_and_finite_write():
    for g in (0, 1, -1):
        actual = ramsey(8, g)
        for n, (j, sign) in itertools.product((0, 1), enumerate((1, -1))):
            np.testing.assert_allclose(actual[n, j], collar.reader(0, 8+n*g, sign), atol=1e-12)
    rabi = math.hypot(64, 4)
    transfer = (64/rabi * math.sin(math.pi*rabi/128))**2
    np.testing.assert_allclose(memory_write(8), [1-transfer, transfer], atol=1e-14)


def test_dense_conditional_ramp_and_capture_against_forward_evolution():
    spec = replay.protocol()
    car = cage.cell()
    D = cage.law(car)
    modes, _, _ = cage.modes(car, cage.ring_profile(car, cage.RING))
    e = modes[:, 0]
    frame, _ = prepare_source(spec['ramp_midpoints_per_star'])
    np.testing.assert_allclose(ramp(car, D*0, 16), frame @ frame.conj().T, atol=2e-11)
    term = np.zeros_like(D)
    term[car['index'][(0, 2)], car['index'][(0, 2)]] = 1
    h = D + term
    state = ramp(car, term, 16)
    calculated, _ = instruments.persistence_law(D, e, modes, state, term=term)
    _, potential, duration = instruments.PREPARATIONS[1]
    _, setting, time = instruments.ANALYZERS[2]
    project = np.outer(e, e)
    for source, action in ((h, cage.STEP), (h + potential*project, duration),
                           (h, cage.STEP*(1+instruments.WAITS[3])), (h + setting*project, time), (h, cage.STEP)):
        U = expm(-1j * action * source)
        state = U @ state @ U.conj().T
    total = np.zeros((151, 151), complex)
    total[:-1, :-1] = state
    H = np.zeros_like(total)
    H[:-1, :-1] = h
    H[-1, :-1] = H[:-1, -1] = cage.STRENGTH * e
    U = expm(-1j * math.pi/(2*cage.STRENGTH) * H)
    assert abs((U @ total @ U.conj().T)[-1, -1] - calculated[1, 3, 2]) < 2e-11


@pytest.mark.parametrize('level', [0, 1, 2])
def test_same_local_spill_on_a_degree_five_ring_after_refinement(level):
    car = cage.federation(level)
    ports = np.flatnonzero(car['degree'] == 5)[:12]
    ring = next(r for r in cage.rings(car, ports) if all(car['degree'][p] == 5 for p in r))
    profile = cage.ring_profile(car, ring)
    e, _, _, native_error = cage.pair(car, profile)
    assert native_error < 1e-13
    for s in np.flatnonzero(profile):
        local = np.zeros_like(e)
        local[s] = e[s]
        spill = local - e * np.dot(e, local)
        assert abs(np.dot(spill, spill) - 3/16) < 1e-14


def test_mobile_pointers_by_evolving_the_two_bits_directly():
    # Independent basis order: scalar mode, load/fuel, then two explicit pointer bits.
    car = cage.cell()
    D = sp.csr_matrix(cage.law(car))
    spec = json.loads((replay.HERE / 'pair_protocol.json').read_text())
    vals, states, _ = contrast_replay.cell(8)
    start_pair = states.index(tuple(spec['emitter']))
    s = car['index'][tuple(spec['scalar_seed_and_reader'])]
    n = len(states) + 1
    H = sp.block_diag((sp.csr_matrix(vals), np.array([[vals[start_pair,start_pair]]])), format='csr')
    source = sp.kronsum(H, D, format='csr')
    shifts = np.zeros((D.shape[0], n))
    for k, (a,b) in enumerate(car['edges']):
        shifts[k, :-1] = [(loads.count(a)-loads.count(b))**2 for loads in states]
    source += sp.diags(shifts.ravel())
    T = sp.lil_matrix((n,n))
    T[start_pair,-1] = T[-1,start_pair] = 64
    psi = np.zeros(source.shape[0])
    psi[s*n+n-1] = 1
    pulse = math.pi/128
    psi = expm_multiply(-1j*pulse*(source+sp.kron(sp.eye(D.shape[0]), T)), psi)
    psi = expm_multiply(-1j*source, psi)
    X = np.array([[0, 1], [1, 0]])
    paired = np.tile(np.arange(n) == start_pair, D.shape[0]).astype(float)
    occupied = np.repeat(np.arange(D.shape[0]) == s, n).astype(float)
    read = (sp.kron(source, sp.eye(4), format='csr')
            + 64*sp.kron(sp.diags(paired), np.kron(X, np.eye(2)), format='csr')
            + 64*sp.kron(sp.diags(occupied), np.kron(np.eye(2), X), format='csr'))
    after = expm_multiply(-1j*pulse*read, np.kron(psi, [1,0,0,0]))
    probabilities = (abs(after.reshape(-1, 4))**2).sum(axis=0)
    saved = json.loads((replay.HERE / 'pair_records.json').read_text())
    np.testing.assert_allclose(probabilities, np.ravel(saved['pointer_laws'][1][3]), atol=2e-11)


def test_binding_bound_inputs_on_the_complete_cell():
    H, states, adjacent = contrast_replay.cell(8)
    binary = np.array([a != b for a, b in states])
    separate = binary & ~adjacent
    largest_adjacent = np.linalg.eigvalsh(H[np.ix_(adjacent, adjacent)])[-1]
    lowest_separate = np.linalg.eigvalsh(H[np.ix_(separate, separate)])[0]
    np.testing.assert_allclose(lowest_separate - largest_adjacent, 15-math.sqrt(6), atol=2e-13, rtol=0)
    np.testing.assert_allclose(np.linalg.norm(H[np.ix_(separate, adjacent)], 2), 2, atol=2e-14)
    assert not np.any(H[np.ix_(~binary, binary)])
    emitter = states.index((0, 2))
    assert np.linalg.norm(np.delete(H[:, emitter], emitter))**2 == pytest.approx(2)
    levels = np.linalg.eigvalsh(cage.law(cage.cell()))
    np.testing.assert_allclose(levels[-1] - levels[0], 2*math.sqrt(6), atol=1e-13)


def test_mobile_receipt_rejects_false_summaries_and_forbidden_weights():
    spec = replay.load_json(replay.HERE / 'pair_protocol.json')
    result = replay.load_json(replay.HERE / 'pair_records.json')
    independent = np.full((3, 6, 2, 2), .25).tolist()
    with pytest.raises(ValueError, match='differs from the joint'):
        pair.assess(spec, dict(result, pointer_laws=independent))
    with pytest.raises(ValueError, match='lost its joint pointer signal'):
        pair.assess(spec, dict(result, pointer_laws=independent, total_variation_from_product=[[0]*6]*3))
    with pytest.raises(ValueError, match='wrong all-time bounds'):
        pair.assess(spec, dict(result, all_time_bounds=[1, 1, 1]))
    for value in (False, -1e-12, 1e-12, float('nan')):
        result['doubled'][0][0] = value
        with pytest.raises(ValueError, match='invalid probability|sector bound'):
            pair.assess(spec, result)


def test_incomplete_freeze_and_redigested_draw_are_rejected(tmp_path):
    manifest = tmp_path / 'sources.json'
    manifest.write_text('{}')
    with pytest.raises(ValueError, match='source list'):
        replay.check_sources(manifest, replay.SOURCE_FILES)
    manifest.write_text('{"a": "first", "a": "second"}')
    with pytest.raises(ValueError, match='duplicate JSON key'):
        replay.check_sources(manifest, {'a'})
    spec = replay.protocol()
    laws = {'ramsey': np.array([[.5, .5]])}
    altered = replay.draw_counts(spec, laws)
    altered['ramsey'][0] += [1, -1]
    receipt = tmp_path / 'counts.npz'
    np.savez_compressed(receipt, **altered)
    receipt.with_suffix('.sha256').write_text(replay.digest(receipt))
    with pytest.raises(ValueError, match='differs from the seeded acquisition'):
        replay.read_counts(receipt, spec, laws)


def test_optimized_replays_and_false_records(monkeypatch):
    monkeypatch.setattr('sys.argv', ['contact'])
    replay.main()
    pair.main()
    spec = replay.protocol()
    with np.load(replay.HERE / 'counts.npz', allow_pickle=False) as z:
        counts = {key: z[key] for key in z.files}
    # Even if purported probabilities exactly match the substituted counts,
    # identical laws at all couplings cannot support the observer claims.
    for records in counts.values():
        records[1:] = records[0]
    fake = {key: c / spec['shots_per_setting'] for key, c in counts.items()}
    with pytest.raises(ValueError, match='collar missed'):
        replay.assess(spec, fake, counts)
    bad = counts['collision'].copy()
    bad[0, 0, 0, 0, 0] -= 1
    bad[0, 0, 0, 0, 1] += 1
    with pytest.raises(ValueError, match='forbidden'):
        replay.assess(spec, fake, dict(counts, collision=bad))
    # Shift every failed write to the written branch, preserving trial totals.
    bad = counts['collision'].copy()
    bad[..., :4] += bad[..., 4:]
    bad[..., 4:] = 0
    with pytest.raises(ValueError, match='variance bound|sampling band'):
        replay.assess(spec, fake, dict(counts, collision=bad))
