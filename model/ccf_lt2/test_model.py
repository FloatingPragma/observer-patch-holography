"""Integrity tests for finite carriers, exact dressing and strict reproduction."""
import json
import numpy as np
import pytest
import sympy as sy
from . import core, source, spatial, verify


def test_hard_wall_retains_integer_flux_and_does_not_wrap():
    cap = 3
    raising = spatial.shift_matrix(cap).toarray()
    electric = np.diag(np.arange(-cap, cap+1))
    assert np.array_equal(electric@raising-raising@electric, raising)
    assert np.count_nonzero(raising[:, -1]) == 0
    assert np.array_equal(spatial.shift_matrix(cap, -3).toarray(),
                          np.linalg.matrix_power(raising.T, 3))
    assert not np.array_equal(raising@raising.T, np.eye(2*cap+1))


def test_wrong_dressing_is_detected_without_changing_hermiticity():
    hop = core.mono(((129, True), (33, False)), sy.Rational(1, 8))
    dressed = core.dress(hop, (3,))
    good = core.add(dressed, core.gadjoint(dressed))
    assert core.check_gauged(good, 2, ((0, 1),))['all_local_gauss_identities']
    mutant = core.dress(hop, (2,))
    mutant = core.add(mutant, core.gadjoint(mutant))
    assert mutant == core.gadjoint(mutant)
    with pytest.raises(AssertionError):
        core.check_gauged(mutant, 2, ((0, 1),))


def test_complete_native_column_keeps_high_occupation_and_finite_wall():
    graph = spatial.edge()
    hop = core.dress(core.mono(((129, True), (33, False)), -1), (3,))
    law = {'graph': graph, 'matter': core.add(hop, core.gadjoint(hop)),
           'field': {'electric': sy.Integer(1), 'magnetic': sy.Integer(1),
                     'B': graph.incidence()[0], 'C': graph.incidence()[1]}}
    # Highest native coordinate is a retained occupied spectator.
    mask = (1 << 191) | (1 << 33)
    column = source.apply_source_column(law, mask, (0,), cap=3)
    destination = (1 << 191) | (1 << 129)
    assert column == {(destination, (3,)): -1}
    wall = source.apply_source_column(law, mask, (1,), cap=3)
    assert wall == {(mask, (1,)): sy.Rational(1, 2)}
    shards = core.add(*(source.apply_source_column(law, mask, (0,), cap=3, vertex=v)
                        for v in range(2)))
    assert shards == column


def test_source_commitment_rejects_modified_snapshot(tmp_path, monkeypatch):
    from . import snapshot
    payload = json.loads((snapshot.HERE/'source_snapshot.json').read_text())
    payload['terms'][0]['weight'] = '2'
    (tmp_path/'source_snapshot.json').write_text(json.dumps(payload))
    snapshot.snapshot.cache_clear()
    with monkeypatch.context() as context:
        context.setattr(snapshot, 'HERE', tmp_path)
        with pytest.raises(ValueError, match='snapshot has changed'):
            snapshot.snapshot()
    snapshot.snapshot.cache_clear()
    assert len(snapshot.snapshot()['terms']) == 114


@pytest.mark.parametrize('actual', [
    {'passed': False, 'p': .25}, {'passed': True},
    {'passed': True, 'p': float('nan')}, {'passed': True, 'p': .2501},
])
def test_result_comparison_rejects_changed_gate_or_probability(actual):
    with pytest.raises(AssertionError):
        verify.compare({'passed': True, 'p': .25}, actual)


def test_radial_carrier_has_closed_face_boundaries():
    graph = spatial.radial_shell()
    boundary, faces = graph.incidence()
    assert graph.nvertices == 8
    assert not np.any(boundary@faces.T)
    directions, lengths = graph.edge_geometry()
    assert len(directions) == len(graph.edges)
    assert all(length.is_positive for length in lengths)


@pytest.mark.parametrize('flux,cap,vertex', [((.5,), 3, None), ((0,), 3.5, None), ((0,), 3, 2)])
def test_column_rejects_labels_outside_integer_carrier(flux, cap, vertex):
    graph = spatial.edge()
    law = {'graph': graph, 'matter': {},
           'field': {'electric': sy.Integer(1), 'magnetic': sy.Integer(1),
                     'B': graph.incidence()[0], 'C': graph.incidence()[1]}}
    with pytest.raises(ValueError):
        source.apply_source_column(law, 0, flux, cap=cap, vertex=vertex)


def test_phase_response_includes_magnetic_loop_and_its_contact():
    graph = spatial.triangle()
    boundary, curl = graph.incidence()
    law = {'graph': graph, 'matter': {},
           'field': {'electric': sy.Integer(1), 'magnetic': sy.Rational(3, 7),
                     'B': boundary, 'C': curl}}
    angle = sy.symbols('angle', real=True)
    shift = tuple(map(int, curl[0]))
    loop_coefficient = -sy.Rational(3, 14)*sy.exp(sy.I*shift[0]*angle)
    current = source.phase_derivative(law, 0)
    contact = source.phase_derivative(law, 0, 0)
    assert current[((), shift)] == sy.diff(loop_coefficient, angle).subs(angle, 0)
    assert contact[((), shift)] == sy.diff(loop_coefficient, angle, 2).subs(angle, 0)
    assert not source.phase_derivative(law, 0, matter_only=True)
    assert source.phase_response_audit(law)['exact_native_charge_continuity']


@pytest.mark.parametrize('payload', ['{"ok":true,"ok":false}', '{"p":NaN}', '{"p":1e999}'])
def test_receipt_loader_rejects_ambiguous_or_nonfinite_json(tmp_path, payload):
    path = tmp_path/'result.json'
    path.write_text(payload)
    with pytest.raises(ValueError):
        verify.read_result(path)


def test_optimized_python_cannot_report_success():
    import subprocess, sys
    process = subprocess.run([sys.executable, '-B', '-O', '-m', 'model.ccf_lt2.verify'],
                             capture_output=True, text=True)
    assert process.returncode != 0
    assert 'requires enabled assertions' in process.stderr
    assert 'all_passed' not in process.stdout


def test_nested_phase_derivative_matches_complete_hamiltonian_column():
    graph = spatial.edge()
    hop = core.dress(core.mono(((129, True), (33, False)), -1), (3,))
    law = {'graph': graph, 'matter': core.add(hop, core.gadjoint(hop)),
           'field': {'electric': sy.Integer(1), 'magnetic': sy.Integer(1),
                     'B': graph.incidence()[0], 'C': graph.incidence()[1]}}
    key = ((0,0), (1, ((1<<191)|(1<<33), (0,))))
    full = source.apply_total_column(law, key, cap=3)
    current = source.apply_total_phase_column(law, key, 0, cap=3)
    contact = source.apply_total_phase_column(law, key, 0, 0, cap=3)
    expected_current = core.clean({out: sy.I*out[1][1][1][0]*value for out,value in full.items()})
    expected_contact = core.clean({out: -out[1][1][1][0]**2*value for out,value in full.items()})
    assert current == expected_current and current
    assert contact == expected_contact and contact


def test_local_capacity_factors_exactly_on_disconnected_components():
    from .local_capacity_checks import disconnected_factorization
    assert disconnected_factorization()['exact_component_sum']


def test_common_capacity_reference_and_complete_dimension():
    from .local_capacity_checks import reference_convention
    from .local_capacity import carrier_dimension
    assert reference_convention()['all_passed']
    assert carrier_dimension(2,1)==2**192*13*3**2*16**2


def test_zero_geometry_edge_weight_is_admitted():
    from .local_capacity import geometry_shard_column
    graph={'nvertices':2,'edges':((0,1,0),)}
    empty={'nvertices':2,'edges':()}
    for index in range(9):
        for vertex in range(2):
            assert geometry_shard_column(index,vertex,graph)==geometry_shard_column(index,vertex,empty)


@pytest.mark.parametrize('face', [(0,1,0,2),(0,1),(0,1,2,0)])
def test_repeated_or_backtracking_face_is_rejected(face):
    from dataclasses import replace
    graph=replace(spatial.triangle(),faces=(face,))
    with pytest.raises(ValueError,match='distinct vertices'):
        graph.incidence()
