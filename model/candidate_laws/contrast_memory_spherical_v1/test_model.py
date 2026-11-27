"""Tests of the A-contrast-memory-spherical-v1 model: local identities, carriers, the theorem's inequalities, the headline checks,
the controls, mutations that each headline check must notice, and the default verification."""
from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm

from model.executable_model import carrier as cr
from . import law, stability, theorem, verify
from .sectors import PairSector, emit, free_run, record_bound, theorem_bound

KAPPA = 8.0


@pytest.fixture(scope='module')
def result():
    """One default verification: recomputes the headline subset and compares it with RESULTS.json."""
    return verify.main([])


# ------------------------------------------------------------------------------------------- local law, carrier

def test_local_identities_exact():
    out = law.exact_checks()
    assert all(out['checks'].values())
    assert out['local_dimension'] == 99


def test_local_reduction_numeric():
    out = law.numeric_checks()
    assert out['U_unitary'] and out['reduction_residual'] < 1e-12


def test_written_memory_transfer_and_free_invariants_on_all_99_inputs():
    h0, h1 = law.local_term_parts()
    tags = [None, *law.LABELS]
    inputs = [(tag, x, y) for tag in tags for x, y in law.SRC]
    for kappa in (0.0, 1.5, 8.0):
        h = h0 + kappa * h1
        for values in ([tag is None for tag, x, y in inputs],
                       [0 if tag is None else tag[1]**2 for tag, x, y in inputs],
                       [x + y for _, x, y in inputs],
                       [x*x + y*y for _, x, y in inputs]):
            charges = np.asarray(values, dtype=float)
            np.testing.assert_allclose(h * (charges[:, None] - charges), 0, atol=2e-14)
        written_zero = [i for i, (tag, x, y) in enumerate(inputs) if tag and tag[1] == 0]
        np.testing.assert_allclose(h[:, written_zero], 0, atol=2e-14)
        initial = inputs.index(((0, 1), 0, 1))
        final = inputs.index(((1, 1), 1, 0))
        assert h[final, initial] == pytest.approx(-kappa/2 - 0.25)
        assert h[inputs.index((None, 2, 0)), inputs.index((None, 0, 2))] == pytest.approx(-0.5)


def test_strict_classical_write_decreases_the_endpoint_cost_on_both_coins():
    for a, b in law.SRC:
        for coin in (0, 1):
            x, y = law.native_event(a, b, coin)
            assert a + b == x + y
            decrease = a*a + b*b - x*x - y*y
            assert decrease == ((a-b)**2 - ((a+b) % 2))/2
            assert (decrease > 0) == (abs(a-b) == 2)


def test_changed_write_breaks_the_reduction(monkeypatch):
    # an output that keeps the order of the inputs: the write stays an isometry, the reduction is lost
    monkeypatch.setattr(law, 'native_event', lambda a, b, c: (a, b))
    assert law.numeric_checks()['reduction_residual'] > 0.1


def test_carriers():
    res = cr.run(levels=(0, 1))
    assert cr.check(res) == []
    assert (res['level0']['ports'], res['level0']['seams']) == (240, 630)
    assert (res['level1']['ports'], res['level1']['seams']) == (960, 2520)


def test_broken_carrier_is_noticed():
    car = cr.federation(0)
    a = car.external[0][0]
    c = next(v for v in range(car.n) if v // 12 != a // 12 and v not in car.neighbours[a])
    broken = cr.Carrier('broken', car.n, car.seams + [(min(a, c), max(a, c))])
    entry = cr.federation_checks(broken)
    assert not entry['external_seams_form_a_matching']


# ---------------------------------------------------------------------------------------------------- theorem

def test_operator_inequalities_20_cells():
    car = cr.federation(0)
    for kappa in (4.0, 8.0):
        entry = theorem.operator_inequalities(car, kappa)
        assert entry['holds'], entry
        assert entry['PHP_largest'] <= 10 * kappa + 6
        assert entry['QHQ_smallest'] >= 12 * kappa
        assert entry['QHP_norm'] <= math.sqrt(10)


def test_single_cell_constants():
    s = theorem.single_cell()
    assert s['amax_is_5'] and s['beta_is_2'] and s['separate_minimum_is_4_minus_sqrt6']
    assert (s['adjacent_pairs'], s['distance_two_pairs'], s['antipodal_pairs']) == (30, 30, 6)


def test_emitter_rationals():
    e = theorem.emitter_bounds()
    assert e['no_pair_doubled_or_fuel'] == '64009/25600000'
    assert e['outside_band_or_fuel'] == '35721/25600000'


@pytest.mark.parametrize('key', ['cell', 0])
def test_probability_bound_against_a_run(key):
    car = cr.carrier(key)
    ps = PairSector(car)
    c = theorem.constants(car)
    seam = cr.EMITTER_SEAM[key]
    state, Ef = emit(ps.hamiltonian(KAPPA), ps.index(*seam))
    lines, _ = free_run(ps, ps.hamiltonian(KAPPA), state[1:], stability.TIMES, shift=Ef)
    worst = max(r['no_pair'] for r in lines)
    da, db = (int(car.degree[v]) for v in seam)
    assert worst <= record_bound(KAPPA, c['adjacent_block_largest_level'], c['coupling_largest_singular_value'], da, db)
    assert worst <= theorem_bound(KAPPA, 64.0)


def test_complete_repeated_readings_can_exceed_the_free_emitter_bound():
    sector = PairSector(cr.cell())
    h = sector.hamiltonian(KAPPA)
    state, _ = emit(h, sector.index(0, 1), 64.0)
    fuel = abs(state[0]) ** 2  # The fuel decouples after the emitter pulse.
    rho = np.outer(state[1:], state[1:].conj())
    flight = expm(-1j * (np.pi / 16) * h.toarray())
    # Distinct seam outcomes are orthogonal. All separated pairs share one quiet outcome.
    quiet = ~sector.adjacent
    same_outcome = np.eye(sector.P, dtype=bool) | np.outer(quiet, quiet)
    readings = {}
    for n in range(1, 201):
        rho = (flight @ rho @ flight.conj().T) * same_outcome
        if n in (1, 100, 200):
            readings[n] = fuel + rho.diagonal()[quiet].real.sum()
    assert [readings[n] for n in (1, 100, 200)] == pytest.approx(
        [0.012802793824, 0.448793890194, 0.506936158729], abs=1e-10, rel=0)

    # Omit exchange during each finite window. The diagonal evolution preserves each outcome
    # of the dephased material state. On one cell 0 <= K <= 10, also with the idle fuel level,
    # so centering K at 5 costs at most 5*tau in trace distance per window (global phase).
    tau = 1 / 16384
    assert readings[200] - 200 * 5 * tau > theorem_bound(KAPPA, 64.0) + 12 * tau


# ------------------------------------------------------------------------------------- headline and controls

def test_default_verification_passes(result):
    assert result['mode'] == 'default'
    assert all(result['particle_stability']['verdict'].values())


def test_controls(result):
    ch = result['particle_stability']['changed_laws']
    for carrier_rows in ch.values():
        assert carrier_rows['exchange alone']['max_no_pair'] > 0.5
        assert carrier_rows['exchange alone']['lowest_level'] > -1e-9
        assert carrier_rows['agreement energy reversed']['lowest_level'] > 0
        assert carrier_rows['law']['lowest_level'] < -1.0 and carrier_rows['law']['pair_in_lowest_level'] > 0.5


def test_face_object(result):
    face = result['particle_stability']['face_object']['one cell']
    assert face['law']['max_no_face'] < 0.05 < face['exchange alone']['max_no_face']


def test_lowest_level_against_dense_diagonalization():
    """The sparse path of `lowest_level` returns the bottom of the spectrum for the law and both changed laws."""
    for key in ('cell', 0):
        ps = PairSector(cr.carrier(key))
        for kind in ('law', 'exchange alone', 'agreement energy reversed'):
            H = stability.law(ps, KAPPA, kind)
            level, weight, deg = stability.lowest_level(H, ps.adjacent, dense_limit=0)
            if key == 'cell':
                w = np.linalg.eigvalsh(H.toarray())
                assert abs(level - w[0]) < 1e-9 and deg == int((w - w[0] < 1e-8).sum())
                assert abs(weight - stability.lowest_level(H, ps.adjacent)[1]) < 1e-9
            if kind == 'exchange alone':
                # the uniform state is annihilated by every exchange term: level 0, pair weight = seams / pairs
                assert abs(level) < 1e-9 and abs(weight - len(ps.car.seams) / ps.P) < 1e-9
            if kind == 'agreement energy reversed':
                assert level > 0


# ------------------------------------------------------------------------------------------------- mutations

def four_checks(H, H_control, p, ps, bound):
    """The four headline checks on one cell for a given law, start state and control."""
    state, Ef = emit(H, p)
    lines, _ = free_run(ps, H, state[1:], stability.TIMES, first=p if ps.adjacent[p] else None, shift=Ef)
    lost = max(r['no_pair'] for r in lines)
    moved = max(r['pair_on_another_seam'] if ps.adjacent[p] else 1 - r['no_pair'] for r in lines)
    level, pair_in_it, _ = stability.lowest_level(H, ps.adjacent)
    s0, E0 = emit(H_control, p)
    lost0 = max(r['no_pair'] for r in free_run(ps, H_control, s0[1:], stability.TIMES, shift=E0)[0])
    return {'kept within the bound': lost <= bound and lost < 0.05, 'read on other seams': moved > 0.05,
            'level below two free loads': level < -1.0 and pair_in_it > 0.5, 'control loses the pair': lost0 > 0.5}


@pytest.fixture(scope='module')
def cell_setup():
    ps = PairSector(cr.cell())
    p = ps.index(0, 1)
    bound = record_bound(KAPPA, 5.0, 2.0, 5, 5)
    return ps, p, bound


def test_emitted_flight_observes_each_requested_exposure_on_its_own_copy():
    sector = PairSector(cr.cell())
    h = sector.hamiltonian(1)
    launch = sector.index(0, 1)
    state, energy = emit(h, launch)
    rows, fuel, last, _ = stability.emitted_run(sector, h, (1, 0), (1, 0, 1, 0.5))
    for row, t in zip(rows, (1, 0, 1, 0.5), strict=True):
        evolved = expm(-1j*t*(h.toarray()-energy*np.eye(sector.P))) @ state[1:]
        expected = fuel + np.linalg.norm(evolved[~sector.adjacent])**2
        assert row["t"] == t and row["no_pair"] == pytest.approx(expected, abs=1e-14)
    np.testing.assert_allclose(last, evolved, rtol=0, atol=1e-14)
    for invalid in ((0, 0), (-1, 0), (0, sector.n), (0, 2)):
        with pytest.raises(ValueError, match="seam"):
            stability.emitted_run(sector, h, invalid, (0,))
    for invalid in (-1, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="exposures"):
            stability.emitted_run(sector, h, (0, 1), (invalid,))


def test_pair_absence_retains_small_fuel_and_separated_branches(monkeypatch):
    sector = PairSector(cr.cell())
    state = np.zeros(1+sector.P, complex)
    state[0] = 1e-12
    state[1+sector.index(0, 1)] = 1
    state[1+sector.index(0, 2)] = 2e-12
    state[1+sector.index(0, 5)] = 3e-12
    monkeypatch.setattr(stability, "emit", lambda h, p: (state.copy(), 0.))
    rows, _, _, _ = stability.emitted_run(sector, sector.hamiltonian(1), (0, 1), (0,))
    assert rows[0]["no_pair"] == pytest.approx(5e-24, rel=1e-14, abs=0)
    assert rows[0]["pair_on_another_seam"] == pytest.approx(9e-24, rel=1e-14, abs=0)


def test_the_law_passes_the_four_checks(cell_setup):
    ps, p, bound = cell_setup
    assert all(four_checks(ps.hamiltonian(KAPPA), ps.hamiltonian(0.0), p, ps, bound).values())


def test_every_check_fails_under_some_mutation(cell_setup):
    ps, p, bound = cell_setup
    H, H0 = ps.hamiltonian(KAPPA), ps.hamiltonian(0.0)
    half = ps.adjacent & ((ps.aa + ps.bb) % 2 == 0)
    apart = next(i for i in range(ps.P) if ps.aa[i] == 0 and not ps.adjacent[i])
    mutations = {
        'agreement energy removed': (H0, H0, p),
        'agreement energy reversed': (ps.hamiltonian(-KAPPA), H0, p),
        'exchange removed': (sparse.diags(H.diagonal()).tocsr(), H0, p),
        'agreement energy on half of the seams': ((H + sparse.diags(2 * KAPPA * (ps.adjacent & ~half))).tocsr(), H0, p),
        'start on two ports that share no seam': (H, H0, apart),
        'control equal to the law': (H, H, p),
    }
    noticed = {}
    for name, (Hm, Hc, pm) in mutations.items():
        for check, ok in four_checks(Hm, Hc, pm, ps, bound).items():
            if not ok:
                noticed.setdefault(check, []).append(name)
    assert set(noticed) == {'kept within the bound', 'read on other seams', 'level below two free loads',
                            'control loses the pair'}, noticed


def test_changed_reader_is_noticed(cell_setup):
    ps, p, bound = cell_setup
    H = ps.hamiltonian(KAPPA)
    state, Ef = emit(H, p)
    lines, psi = free_run(ps, H, state[1:], stability.TIMES, shift=Ef)
    seam_reader = np.abs(psi[ps.seam_states]) ** 2
    # a reader that reports "both ends loaded" on pairs of ports that share no seam
    wrong = np.flatnonzero(~ps.adjacent)[:len(ps.seam_states)]
    assert 1 - seam_reader.sum() < 0.05
    assert 1 - (np.abs(psi[wrong]) ** 2).sum() > 0.5


def test_comparison_names_the_failing_key(result):
    import copy
    stored = copy.deepcopy(result)
    stored['theorem']['level0']['adjacent_block_largest_level'] += 1e-3
    failures = []
    verify.compare(stored, result, '', failures, False)
    assert failures and failures[0].startswith('theorem.level0.adjacent_block_largest_level')
