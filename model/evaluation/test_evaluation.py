"""The evaluation rejects laws off the fixed carrier, nonlocal or asymmetric laws and hardcoded outcomes."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from model.executable_model import CELL_SEAMS, carrier, laws

from . import evaluate as ev


def toy(H, seams=CELL_SEAMS, parameters=None, response=None):
    """A one-port-per-register law on the cell with the generator H."""
    cell = carrier.cell()

    def represent(rotation, car=None):
        P = np.zeros((12, 12))
        for p, q in enumerate(rotation):
            P[q, p] = 1.0
        return P

    return SimpleNamespace(NAME='toy', LABELLING=seams, basis=lambda car=None: [('port', p) for p in range(12)],
                           generator=lambda car=None: H if car is None or car.cells == 1 else
                           (_ for _ in ()).throw(NotImplementedError('one cell')),
                           response=lambda car=None: H if response is None else response,
                           represent=represent, PARAMETERS=parameters or {}, cell=cell)


def laplacian():
    return carrier.cell().laplacian().toarray()


def test_the_cell_laplacian_passes_every_carrier_check():
    checks = ev.carrier_checks(toy(laplacian() / 2))
    assert {k: c['status'] for k, c in checks.items()} == {
        'native_cell': 'pass', 'hermitian': 'pass', 'locality': 'pass', 'rotations': 'pass',
        'rank_three_response': 'pass', 'federation': 'partial'}


def test_a_supplied_lattice_is_not_the_native_cell():
    square = [(i, j) for i in range(12) for j in (i + 1, i + 4) if j < 12 and (j != i + 1 or (i + 1) % 4)][:30]
    assert ev.carrier_checks(toy(np.zeros((12, 12)), seams=square))['native_cell']['status'] == 'fail'


def test_a_coupling_between_distant_ports_breaks_locality():
    H = laplacian()
    a, b = next((a, b) for a in range(12) for b in range(12) if a < b and H[a, b] == 0)
    H[a, b] = H[b, a] = -1.0
    assert ev.carrier_checks(toy(H))['locality']['status'] == 'fail'


def test_a_preferred_port_breaks_the_rotations():
    H = laplacian()
    H[0, 0] += 0.5
    assert ev.carrier_checks(toy(H))['rotations']['status'] == 'fail'


def test_the_slowest_response_must_have_rank_three():
    flat = np.diag([0.0] + [1.0] * 11)
    assert ev.carrier_checks(toy(laplacian(), response=flat))['rank_three_response']['status'] == 'fail'


@pytest.mark.parametrize('quantity', ['electron particle mass', 'spatial dimension', 'spacetime metric'])
def test_parameters_may_not_hardcode_an_observed_outcome(quantity):
    cand = toy(laplacian(), parameters={'x': {'value': 3, 'provenance': 'designed', 'quantity': quantity}})
    assert ev.parameter_check(cand)['status'] == 'fail'


def test_a_calibrated_parameter_names_its_observation():
    cand = toy(laplacian(), parameters={'x': {'value': 1, 'provenance': 'calibrated', 'quantity': 'service rate'}})
    assert ev.parameter_check(cand)['status'] == 'fail'
    cand.PARAMETERS['x']['observation'] = 'a declared clock comparison'
    assert ev.parameter_check(cand)['status'] == 'pass'


def test_every_candidate_declares_known_premises_and_runs_on_the_cell():
    for name in laws.names():
        result = ev.evaluate(name)
        assert result['checks']['native_cell']['status'] == 'pass'
        assert result['checks']['parameters']['status'] == 'pass'
        assert len(result['premises']) == len(ev.json.loads(ev.REGISTER.read_text())['rows'])


def test_written_evaluations_are_current():
    order = sorted(laws.names(), key=lambda n: (laws.load(n).ROLE != 'main', n))
    stale = [p for p, text in ev.outputs([ev.evaluate(n) for n in order]).items()
             if not p.exists() or p.read_text(encoding='utf-8') != text]
    assert not stale, f'run python -m model evaluate --write: {stale}'
