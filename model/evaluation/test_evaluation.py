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


def test_rotations_move_the_native_registers():
    cand = toy(laplacian())
    cand.represent = lambda rotation, car: np.eye(12)
    assert ev.carrier_checks(cand)['rotations']['status'] == 'fail'


def test_complex_generators_use_the_adjoint_and_retain_imaginary_couplings():
    H = laplacian().astype(complex)
    a, b = next((a, b) for a in range(12) for b in range(a + 1, 12) if H[a, b] == 0)
    H[a, b], H[b, a] = 1j, -1j
    checks = ev.carrier_checks(toy(H, response=laplacian()))
    assert checks['hermitian']['status'] == 'pass'
    assert checks['locality']['status'] == 'fail'
    H[b, a] = 1j
    assert ev.carrier_checks(toy(H, response=laplacian()))['hermitian']['status'] == 'fail'


def test_federation_requires_the_native_basis_in_each_cell():
    cand = toy(laplacian())
    cand.generator = lambda car: laplacian()
    assert ev.carrier_checks(cand)['federation']['status'] == 'fail'


def test_federation_rejects_a_link_between_unjoined_cells():
    cand = toy(laplacian())
    cand.basis = lambda car: [('port', p) for p in range(car.n)]

    def generator(car):
        H = car.laplacian().tolil()
        if car.cells > 1:
            a, b = next((0, b) for b in range(12, car.n) if b not in car.neighbours[0])
            H[a, b] = H[b, a] = 1.0
        return H.tocsr()

    cand.generator = generator
    check = ev.carrier_checks(cand)['federation']
    assert check['status'] == 'fail'
    assert 'nonlocal couplings 2' in check['detail']


def test_the_slowest_response_must_have_rank_three():
    flat = np.diag([0.0] + [1.0] * 11)
    assert ev.carrier_checks(toy(laplacian(), response=flat))['rank_three_response']['status'] == 'fail'


def test_a_nonhermitian_response_cannot_pass_via_its_lower_triangle():
    response = laplacian().astype(complex)
    response[0, 1] += 1j
    assert ev.carrier_checks(toy(laplacian(), response=response))['rank_three_response']['status'] == 'fail'


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


def test_every_existing_law_has_an_explicit_assessment_for_each_required_id():
    from . import assessment

    ids = set()
    for filename in assessment.FILES:
        for group in assessment.read_json(assessment.HERE / filename)['groups']:
            for row in group['criteria']:
                assert set(assessment.SHORT) <= set(row['results']), row['id']
                assert row['id'] not in ids
                ids.add(row['id'])
    assert {'A1', 'A2', 'A3', 'FT-BAND', 'FT-CONFLUENCE', 'GR-NULL'} <= ids
    assert {f'G{i}' for i in range(1, 7)} <= ids
    assert {key for key in ids if key.startswith('OBS-')} == {f'OBS-{i:02}' for i in range(1, 8)}
    ow = assessment.read_json(assessment.HERE / 'openwave_criteria.json')['groups'][0]['criteria']
    assert [(row['id'], row['upstream_number']) for row in ow] == [(f'OW-{i:02}', i) for i in range(1, 32)]

    for name in laws.names():
        particles = assessment.species_group(laws.load(name))['criteria']
        assert len(particles) == len({r['id'] for r in particles}) == 17


def test_readout_matches_the_original_measurement_and_scope_is_preserved():
    result = ev.evaluate('interval_echo')
    rows = {r['id']: r for group in result['requirements'] for r in group['criteria']}
    m = result['measurement']
    original = ev.json.loads((ev.ROOT / m['source']).read_text())
    measured = original['cohorts']['panels']['alternating_2_6']['three_bits']
    assert (m['value'], m['standard_error']) == (measured['ratio'], measured['standard_error'])
    assert rows['PATCH-02']['status'] == 'demonstrated'
    assert rows['OBS-01']['status'] == 'partial'
    assert rows['SM-e']['status'] == 'open'
    assert ev.criterion(result, 'JOINT-relativistic-matter')['status'] == 'fails'
    assert ev.criterion(result, 'pr-43')['status'] == 'fulfilled'
    assert ev.criterion(result, 'unregistered') is None


def test_every_scorecard_requirement_links_to_its_own_evidence_heading():
    import re

    scorecard = ev.SCORECARD.read_text()
    evidence = (ev.HERE / 'EVIDENCE.md').read_text()
    headings = {title.casefold() for title in re.findall(r'^### (.+)$', evidence, re.M)}
    targets = re.findall(r'\]\(EVIDENCE\.md#([^)]*)\)', scorecard)
    assert 'ow-31' in targets and 'sm-h' in targets and 'joint-one-law' in targets
    assert set(targets) - {'full-premise-register'} <= headings


def test_openwave_rows_distinguish_free_stability_from_electron_identity_and_driven_waves_from_free_waves():
    result = ev.evaluate('interval_echo')
    scorecard = ev.report.scorecard([result])
    ledger = ev.report.evidence_ledger([result])
    for key in ('OW-04', 'OW-30'):
        row = ev.criterion(result, key)
        measurement = row['demonstrated_measurement']
        assert measurement['source'] in row['evidence']
        link = f"[{measurement['label']}]({ev.report.local_link(measurement['source'])})"
        grade = f"**Demonstrated** · {link}"
        if key == 'OW-30':
            assert row['status'] == 'partial'
            grade += ". Full criterion · Partial"
        else:
            assert row['status'] == 'demonstrated'
            assert ev.criterion(result, 'SM-e')['status'] == 'open'
            filing = row['filing']
            assert filing['decision'] in row['sources'] and filing['review'] in row['evidence']
            grade += (f". OpenWave · [⚠️ Partial]({filing['decision']})"
                      f". Resubmission · [HOLD]({ev.report.local_link(filing['review'])})")
        scorecard_row = next(line for line in scorecard.splitlines() if line.startswith(f'| [{key}]'))
        assert grade in scorecard_row
        assert grade in ledger.split(f'### {key}\n', 1)[1].split('\n### ', 1)[0]
