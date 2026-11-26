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


@pytest.mark.parametrize('phase_on_all_registers', [True, False])
def test_rotation_products_allow_common_phases_and_reject_relative_ones(phase_on_all_registers):
    from scipy.linalg import block_diag
    from model.executable_model import rotations

    labels = [('port', p) for p in range(12)] + [('seam', edge) for edge in CELL_SEAMS]
    positions = {label: i for i, label in enumerate(labels)}
    involution = next(g for g in rotations() if g != tuple(range(12)) and all(g[g[p]] == p for p in range(12)))
    candidate = toy(block_diag(laplacian(), 24 * np.eye(30)))
    candidate.basis = lambda car: labels

    def represent(rotation, car):
        targets = [positions[ev.register_at(label, rotation.__getitem__)] for label in labels]
        matrix = np.eye(len(labels), dtype=complex)[:, targets]
        if tuple(rotation) == involution:
            matrix[:, 0 if phase_on_all_registers else 12:] *= 1j
        return matrix

    candidate.represent = represent
    result = ev.carrier_checks(candidate)['rotations']
    assert result['status'] == ('pass' if phase_on_all_registers else 'fail')
    if not phase_on_all_registers:
        assert 'relative register phases' in result['detail']


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


def test_an_isolated_owner_fails_the_live_carrier_command_and_candidate_evaluation(monkeypatch, capsys):
    from model.__main__ import main
    from model.executable_model.validation import require_three_gateways

    fed = carrier.federation(0)
    isolated = [(fed.n + a, fed.n + b) for a, b in CELL_SEAMS]
    board = carrier.Carrier('one isolated owner', fed.n + 12, fed.seams + isolated)
    monkeypatch.setattr(carrier, 'federation', lambda level: board)
    assert main(['carrier']) == 1
    assert 'owner 20 has 0 external ports' in capsys.readouterr().out
    check = ev.carrier_checks(toy(laplacian()))['federation']
    assert check['status'] == 'fail'
    assert 'owner 20 has 0 external ports' in check['detail']
    with pytest.raises(ValueError, match='owner 0 has 0 external ports'):
        require_three_gateways(carrier.cell())


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
    assert ev.carrier_checks(toy(flat))['rank_three_response']['status'] == 'fail'


@pytest.mark.parametrize(('unit', 'origin'), [(1e-10, 0), (1, 1e12)])
def test_energy_units_and_origin_preserve_the_threefold_band(unit, origin):
    response = unit * laplacian() + origin * np.eye(12)
    assert ev.carrier_checks(toy(response))['rank_three_response']['status'] == 'pass'


def test_small_energy_units_preserve_nonlocal_couplings_and_response_mismatch():
    h = laplacian()
    a, b = next((0, b) for b in range(1, 12) if h[0, b] == 0)
    h[a, b] = h[b, a] = -1
    checks = ev.carrier_checks(toy(1e-15 * h, response=1e-15 * laplacian()))
    assert checks['locality']['status'] == 'fail'
    assert checks['rank_three_response']['status'] == 'fail'


def test_a_zero_generator_cannot_borrow_the_cell_laplacian_response():
    check = ev.carrier_checks(toy(np.zeros((12, 12)), response=laplacian()))['rank_three_response']
    assert check['status'] == 'fail'
    assert 'differs from the registered' in check['detail']


def test_scalar_response_is_reconstructed_from_the_registered_seam_wedge_couplings(monkeypatch):
    scalar = laws.load('scalar_channel')
    checks = ev.carrier_checks(scalar)
    assert checks['rank_three_response']['status'] == 'pass'
    assert checks['federation']['status'] == 'partial'
    generate = scalar.generator
    monkeypatch.setattr(scalar, 'generator', lambda car: 2 * generate(car))
    check = ev.carrier_checks(scalar)['rank_three_response']
    assert check['status'] == 'fail'
    assert 'seam-wedge Gram' in check['detail']


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
        assert len(result['premises']) == len(ev.premise_register())


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
    assert {'A1', 'A2', 'A3', 'FT-BAND', 'FT-CONFLUENCE'} <= ids
    gravity = next(group for group in assessment.definitions() if group['id'] == 'einstein')
    assert [row['id'] for row in gravity['criteria']] == [f'G{i}' for i in range(1, 7)]
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
            assert 'unread connector memories' in measurement['label']
            assert any(path.endswith('row_30_cautions/README.md') for path in row['evidence'])
            grade = f"**Partial**. Demonstrated measurement · {link}"
            assert measurement['source'].endswith('row_30_collective/README.md')
            filing = row['filing']
            assert filing['decision'] in row['sources'] and filing['review'] in row['evidence']
            grade += (f". OpenWave · [🚧 Not yet tested]({filing['decision']})"
                      f". Full-credit submission · [HOLD]({ev.report.local_link(filing['review'])})")
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


def test_written_field_comparison_is_scoped_to_the_separate_generator():
    result = ev.evaluate('interval_echo')
    rows = next(group['criteria'] for group in result['requirements'] if group['id'] == 'openwave')
    electron = ev.criterion(result, 'SM-e')
    assert electron['status'] == 'open'
    study = electron['comparison']
    assert study['law'] == 'Written-field Laplacian extension'
    assert 'Strong candidate for the electron’s charge and field response' in ev.report.comparison_link(study)
    assert study['source'] not in electron['evidence']
    assert study['source'].endswith('interval_echo/openwave/electron_field/README.md')
    for row in rows:
        if 'comparison' in row:
            assert row['comparison']['source'] not in row['evidence']
            assert ev.report.criterion_grade(row) == ev.report.criterion_grade(
                {key: value for key, value in row.items() if key != 'comparison'}
            )
    assert all('comparison' not in row for row in ev.assessment.species_group(laws.load('scalar_channel'))['criteria'])
    ledger = ev.report.evidence_ledger([result])
    electron_entry = ledger.split('### SM-e\n', 1)[1].split('\n### ', 1)[0]
    assert study['finding'] in electron_entry
    assert ev.report.local_link(study['source']) in electron_entry
    assert ev.report.comparison_link(study) in ev.report.scorecard([result])


def test_pinned_dark_contact_control_uses_all_native_load_columns():
    from model.candidate_laws.interval_echo.openwave.electron_field import dark_contact

    dark_contact.run()


def test_photon_comparisons_keep_each_construction_and_evidence_visible():
    result = ev.evaluate('interval_echo')
    photon = ev.criterion(result, 'SM-photon')
    studies = ev.assessment.comparison_studies(photon)
    assert photon['status'] == 'open'
    assert len(studies) == len({s['law'] for s in studies}) == 2
    scorecard = ev.report.scorecard([result])
    evidence = ev.report.evidence_ledger([result])
    for study in studies:
        assert ev.report.comparison_link(study) in scorecard
        assert study['finding'] in evidence
        assert study['source'] not in photon['evidence']


def test_measured_premises_preserve_the_candidate_and_reject_duplicate_sources(tmp_path, monkeypatch):
    import json

    directory = tmp_path / 'model' / 'candidate_laws' / 'test_source'
    directory.mkdir(parents=True)
    (directory / 'receipt.json').write_text('{}')
    path = directory / 'premise_results.json'
    path.write_text(json.dumps({'PR-1': {'status': 'fulfilled', 'evidence': ['receipt.json'],
                                        'note': 'Finite run with a complete public record.'}}))
    candidate = SimpleNamespace(NAME='test_source', PREMISES={})
    premise = {'id': 'PR-1', 'name': 'record', 'statement': 'A record is emitted.',
               'premise_kind': 'derivable'}
    roles = {'types': {}, 'overrides': {}, 'defaults': {'derivable': 'law'}}
    monkeypatch.setattr(ev, 'ROOT', tmp_path)
    assert ev.premise_rows(candidate, [premise], roles)[0]['status'] == 'fulfilled'
    assert candidate.PREMISES == {}
    candidate.PREMISES['PR-1'] = {'status': 'open'}
    with pytest.raises(SystemExit, match='duplicate premise assessments'):
        ev.premise_rows(candidate, [premise], roles)

    # A new observation may supersede exactly the historical assessment it names.
    from hashlib import sha256
    old = {'status': 'partial', 'evidence': ['first.json'], 'note': 'First source only.'}
    (directory / 'first.json').write_text('{}')
    candidate.PREMISES['PR-1'] = old
    measured = json.loads(path.read_text())
    measured['PR-1']['supersedes_sha256'] = sha256(
        json.dumps(old, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    path.write_text(json.dumps(measured))
    updated = ev.premise_rows(candidate, [premise], roles)[0]
    assert updated['status'] == 'fulfilled' and updated['evidence'] == ['first.json', 'receipt.json']
    assert updated['note'] == 'First source only. Finite run with a complete public record.'
    assert candidate.PREMISES['PR-1'] is old
    candidate.PREMISES['PR-1'] = dict(old, status='fails')
    with pytest.raises(SystemExit, match='prior premise assessment has changed'):
        ev.premise_rows(candidate, [premise], roles)
    del candidate.PREMISES['PR-1']
    with pytest.raises(SystemExit, match='no prior premise assessment'):
        ev.premise_rows(candidate, [premise], roles)
