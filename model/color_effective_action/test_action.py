"""Independent source solves, complete complements and falsification controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy import sparse
from scipy.linalg import expm, sqrtm
from scipy.sparse.linalg import spsolve

from model.interacting_extension.higgs_join.interacting_carrier import build_join, NativeCarrier
from . import action, link, packet


@pytest.fixture(scope='module')
def receipt():
    return packet.build()


def numeric(matrix):
    return np.asarray([[float(x) for x in row] for row in matrix])


def full_source(control=action.Controls()):
    return build_join(electric=float(control.electric), magnetic=float(control.magnetic),
        hopping=float(control.hopping), higgs_coupling=float(control.higgs),
        scalar_energy=float(control.scalar))['sectors'][0]['H'].real.tocsr()


def direct_schur(mu, control=action.Controls()):
    h = full_source(control)
    a, c, d = h[:6, :6].toarray(), h[6:, :6].toarray(), h[6:, 6:]
    inverse_c = spsolve(mu*sparse.eye(432, format='csr')+d, c)
    return mu*np.eye(6)+a-c.T@inverse_c, (a, c, d, inverse_c)


def test_committed_receipt_and_full_source_replay(receipt):
    assert packet.canonical(receipt) == packet.canonical(json.loads((packet.HERE/'results.json').read_text()))
    assert receipt['literal_source_replay']['operator_columns_replayed'] == 20778
    assert receipt['arbitrary_control_component_identities']['occupied_scalar_Higgs_columns'] == 432
    assert receipt['scope']['physical_alpha_s'] is None
    assert receipt['scope']['physical_particle_masses'] is None
    assert not receipt['scope']['spatial_refining_source_family_derived']
    assert not receipt['scope']['physical_scale_or_QCD_scheme_matching_derived']
    assert not receipt['scope']['filter_instrument_native_pulse_execution_claimed']


@pytest.mark.parametrize('control', (action.Controls(),
    action.Controls(F(3, 2), F(2, 3), F(5, 4), F(7, 5), F(3, 2)),
    action.Controls(F(2), F(-1, 2), F(-3, 4), F(-1, 3), F(5, 4))))
def test_source_coefficients_match_full_independent_generator(control):
    model = action.Elimination(control)
    h = full_source(control)
    a, c, d = h[:6, :6].toarray(), h[6:, :6].toarray(), h[6:, 6:]
    assert np.max(np.abs(a-numeric(model.a))) < 1e-14
    for k, expected in enumerate(model.moments(2)):
        computed = c.copy()
        for _ in range(k):
            computed = d@computed
        assert np.max(np.abs(c.T@computed-numeric(expected))) < 3e-12
    action.exact_moment_certificate(control)


@pytest.mark.parametrize('mu,order', ((64, 6), (64, 10), (128, 8), (256, 6)))
def test_rational_operator_error_bounds_full_438_state_schur_solve(mu, order):
    model = action.Elimination()
    approx, error = model.enclosure(F(mu), order)
    actual, _ = direct_schur(mu)
    # The coarse Neumann errors exceed binary64 roundoff in these tests.
    assert float(error) > 1e-11
    assert np.linalg.norm(actual-numeric(approx), 2) < float(error)
    h = full_source()
    projected_inverse = spsolve(mu*sparse.eye(438, format='csr')+h, np.eye(438, 6))[:6, :]
    assert np.linalg.norm(projected_inverse-np.linalg.inv(actual), 2) < 2e-16


def test_same_source_minimization_monotonicity_and_reference_shift():
    mu = 64
    effective, (a, c, d, inverse_c) = direct_schur(mu)
    rng = np.random.default_rng(7021)
    phi = rng.normal(size=6)+1j*rng.normal(size=6)
    xi = -inverse_c@phi
    h = full_source().toarray()+mu*np.eye(438)
    minimizer = np.concatenate([phi, xi])
    minimum = np.vdot(phi, effective@phi)
    assert abs(np.vdot(minimizer, h@minimizer)-minimum) < 2e-12
    for _ in range(4):
        delta = rng.normal(size=432)+1j*rng.normal(size=432)
        perturbed = minimizer+np.concatenate([np.zeros(6), delta])
        difference = np.vdot(perturbed, h@perturbed)-minimum
        expected = np.vdot(delta, (mu*sparse.eye(432)+d)@delta)
        assert abs(difference-expected) < 1e-9
        assert difference.real > 0
    derivative = np.eye(6)+inverse_c.T@inverse_c
    curvature = -2*inverse_c.T@spsolve(mu*sparse.eye(432, format='csr')+d, inverse_c)
    assert np.linalg.eigvalsh(derivative).min() >= 1
    assert np.linalg.eigvalsh(curvature).max() < 0


def test_clock_rescaling_preserves_dimensionless_full_action():
    original = action.Elimination()
    changed = action.Elimination(action.Controls(3, 3, 3, 3, 3))
    a, _ = original.enclosure(F(64), 12)
    b, _ = changed.enclosure(F(192), 12)
    assert b == action.matrix_scale(a, 3)
    original_moments = original.moments(4)
    for k, moment in enumerate(changed.moments(4)):
        assert moment == action.matrix_scale(original_moments[k], 3**(k+2))
    base = action.uniform_family_remainders(F(64), original.d_bound, F(1), F(1))
    scaled = action.uniform_family_remainders(F(192), 3*original.d_bound, F(3), F(3))
    assert base == scaled
    assert F(base['normalized_remainder_after_first_two_moments']) > 0


def test_effective_action_has_source_generated_terms_outside_bare_ansatz(receipt):
    for row in receipt['effective_actions']:
        for name in ('induced_vacuum_to_weak_flux', 'induced_vacuum_to_charged_singlet_flux'):
            lo, hi = map(F, row[name])
            assert lo < hi < 0
        assert row['scalar_complement_states_eliminated_exactly_in_the_formula'] == 432
        assert not row['single_Yang_Mills_coupling_ansatz_used']
    # A term omitted by the original electric/magnetic ansatz is present in
    # an independently computed full-generator inverse.
    actual, _ = direct_schur(64)
    assert actual[0, 2] < -4e-5 and actual[0, 5] < -3e-5
    for matrix in (action.identity(6), action.ELECTRIC, action.MAGNETIC):
        assert matrix[0][2] == matrix[0][5] == 0
    ratio_intervals = [tuple(map(F, row['projected_electric_over_positive_magnetic']))
                       for row in receipt['effective_actions']]
    assert ratio_intervals[0][0] > ratio_intervals[1][1] > ratio_intervals[2][1] > 1


def test_higgs_and_hopping_controls_remove_their_own_generated_term():
    for control in (action.Controls(higgs=0), action.Controls(hopping=0)):
        model = action.Elimination(control)
        effective, error = model.enclosure(F(64), 12)
        assert effective[0][2] == effective[0][5] == 0
        if control.higgs == 0:
            assert error == 0
            assert effective == action.matrix_add(model.a, action.identity(6), 64)
    base = action.Elimination()
    assert base.moments(2)[2][0][2] == 12


def test_component_certificate_rejects_false_hopping_gram(monkeypatch):
    damaged = deepcopy(action.INDUCED)
    damaged[0][2] = damaged[2][0] = 0
    monkeypatch.setattr(action, 'INDUCED', damaged)
    with pytest.raises(ArithmeticError, match='hopping Gram'):
        action.component_certificate()


def test_induced_positive_operator_has_exact_rank_and_characteristic_polynomial():
    import sympy as s
    matrix = s.Matrix(action.INDUCED)
    z = s.Symbol('z')
    assert matrix.charpoly(z).as_expr().expand() == (z**4*(z*z-54*z+432)).expand()
    assert matrix.rank() == 2
    assert matrix.is_positive_semidefinite


def test_same_hamiltonian_transfer_has_positive_reflected_kernel_and_complete_filter():
    model = action.Elimination()
    h = full_source().toarray()+float(model.source.norm_bound)*np.eye(438)
    rng = np.random.default_rng(77)
    # Complex vectors also test the sesquilinear signs in the reflection.
    v = rng.normal(size=(438, 3))+1j*rng.normal(size=(438, 3))
    v /= np.linalg.norm(v, axis=0)
    times = (.005, .011, .017)
    states = [expm(-t*h)@v[:, i] for i, t in enumerate(times)]
    reflected = np.array([[np.vdot(v[:, i], expm(-(t+s)*h)@v[:, j])
        for j, s in enumerate(times)] for i, t in enumerate(times)])
    direct = np.array([[np.vdot(a, b) for b in states] for a in states])
    assert np.max(np.abs(reflected-direct)) < 5e-15
    assert np.linalg.eigvalsh(reflected).min() > 0
    ok = expm(-.01*h)
    fail = sqrtm(np.eye(438)-ok@ok)
    assert np.max(np.abs(ok.T@ok+fail.conj().T@fail-np.eye(438))) < 4e-14
    entangled = v[:, :2]/np.sqrt(2)
    assert abs(np.linalg.norm(ok@entangled)**2+np.linalg.norm(fail@entangled)**2-1) < 2e-14


def test_native_link_defect_on_full_register_and_spectator():
    receipt = link.certificate()
    carrier = NativeCarrier()
    for r, block in enumerate(carrier.data['blocks']):
        d = len(block['indices'])
        u = np.zeros((d*60, d*60))
        for a, ia in enumerate(block['indices']):
            for b, ib in enumerate(block['indices']):
                u[a*60+carrier.reverse[r, ia, ib], b*60] = 1/np.sqrt(d)
        assert np.max(abs(u@u)) == 0
        defect = u.T@u-np.eye(d*60)
        assert abs(np.linalg.norm(defect, 2)-1) < 3e-15
        assert abs(np.linalg.eigvalsh((u+u.T)/2).max()-.5) < 1e-15
        assert np.linalg.svd(np.eye(d*60)-u, compute_uv=False).min() > .6
        assert receipt['blocks'][r]['normalized_partial_isometry_unitarity_defect_norm'] == '1'


def test_exact_control_and_resolution_domain_failures():
    model = action.Elimination()
    with pytest.raises(ValueError, match='exceed'):
        model.enclosure(model.source.norm_bound)
    with pytest.raises(TypeError, match='exact rational'):
        model.enclosure(64.)
    with pytest.raises(ValueError, match='integer'):
        model.moments(-1)
    with pytest.raises(ValueError, match='moment list'):
        model.enclosure(F(64), 2, model.moments(1))


def test_exact_outward_interval_and_fit_normalization():
    lo, hi = map(F, action.interval_packet(F(1, 3), F(1, 7), bits=8))
    assert lo <= F(1, 3)-F(1, 7) <= F(1, 3)+F(1, 7) <= hi
    basis = [action.identity(6), action.ELECTRIC, action.matrix_scale(action.MAGNETIC, -1)]
    for i, weights in enumerate(action.linear_fit_weights().values()):
        assert [sum(weights[a][b]*matrix[a][b] for a in range(6) for b in range(6))
                for matrix in basis] == [int(i == j) for j in range(3)]


def test_complete_energy_scores_reconstruct_the_resolvent_and_action():
    mu = 64
    h = full_source().toarray()
    energies, vectors = np.linalg.eigh(h)
    score = 1/(mu+energies)
    direct = (vectors[:6]*score)@vectors[:6].T
    measured = np.zeros((6, 6))
    for i in range(6):
        measured[i, i] = np.sum(vectors[i]**2*score)
        for j in range(i):
            plus = (vectors[i]+vectors[j])/np.sqrt(2)
            minus = (vectors[i]-vectors[j])/np.sqrt(2)
            assert abs(np.sum(plus**2)-1) < 2e-15
            assert abs(np.sum(minus**2)-1) < 2e-15
            measured[i, j] = measured[j, i] = (np.sum(plus**2*score)-np.sum(minus**2*score))/2
    assert np.max(np.abs(measured-direct)) < 1e-17
    effective, _ = direct_schur(mu)
    assert np.linalg.norm(np.linalg.inv(measured)-effective, 2) < 2e-12
    bound = action.Elimination().source.norm_bound
    delta = F(1, 10**7)
    error = action.inverse_readout_error(mu, bound, delta)
    changed = direct+float(delta)*np.ones((6, 6))
    assert np.linalg.norm(np.linalg.inv(changed)-np.linalg.inv(direct), 2) < float(error)
    with pytest.raises(ValueError, match='invertibility'):
        action.inverse_readout_error(mu, bound, F(1))


def test_requested_receipt_path_is_checked_and_corruption_is_rejected(tmp_path, receipt, monkeypatch):
    target = tmp_path/'requested.json'
    target.write_text(packet.canonical(receipt))
    monkeypatch.setattr(packet, 'build', lambda:receipt)
    packet.verify(target)
    altered = deepcopy(receipt)
    altered['scope']['physical_alpha_s'] = '1'
    target.write_text(packet.canonical(altered))
    with pytest.raises(ValueError, match='scope differs'):
        packet.verify(target)
    with pytest.raises(FileNotFoundError):
        packet.verify(tmp_path/'missing.json')


def test_source_bindings_cover_full_producer_runtime_imports():
    script = '''
from pathlib import Path
import sys
from model.color_effective_action import packet
receipt=packet.build()
root=packet.ROOT.resolve()
loaded=set()
for module in tuple(sys.modules.values()):
    value=getattr(module,"__file__",None)
    if not value: continue
    path=Path(value).resolve()
    if path.suffix==".py" and path.is_relative_to(root): loaded.add(path.relative_to(root).as_posix())
assert loaded<=set(receipt["source_bindings"]),sorted(loaded-set(receipt["source_bindings"]))
assert "model/color_loop_readout/source.py" in loaded
assert "model/color_effective_action/action.py" in loaded
print("complete runtime source closure verified")
'''
    run = subprocess.run([sys.executable, '-c', script], cwd=packet.ROOT,
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout+run.stderr
    assert 'complete runtime source closure verified' in run.stdout
