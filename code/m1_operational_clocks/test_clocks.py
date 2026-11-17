"""Independent mathematics, executed countercontrols and hostile verifier inputs."""

import copy
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from m1_operational_clocks import check, model, verify


@pytest.fixture(scope='module')
def expected():
    return check.reconstruct()


@pytest.fixture(scope='module')
def candidate():
    return model.candidate()


def test_independent_complete_replay(candidate, expected):
    check.same(candidate, expected)


def test_committed_receipt(expected):
    packet = verify.load(verify.HERE/'receipt.json')
    check.same(packet['evidence'], expected)
    check.same(packet['sources'], verify.pins())
    check.same(packet['claims'], verify.claim_pins())
    assert packet['schema'] == 'oph-m1-operational-clocks-v1'


def test_symbolic_full_clifford_and_mass():
    j, root = sp.I, sp.sqrt(3)
    c = sp.Matrix([[0, 1, 1, 1], [1, 0, j, -j], [1, -j, 0, j], [1, j, -j, 0]])/root
    assert c*c == sp.eye(4)
    points = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    xs = []
    for i in range(3):
        b = sp.diag(*(p[i] for p in points))
        xs.append(root*(b+c*b*c)/2)
    for i, x in enumerate(xs):
        assert sp.simplify(x*c-c*x) == sp.zeros(4)
        for k, y in enumerate(xs):
            assert sp.simplify(x*y+y*x) == 2*(i == k)*sp.eye(4)
    beta = sp.kronecker_product(sp.Matrix([[0, 1], [1, 0]]), sp.eye(4))
    alpha = [sp.kronecker_product(sp.diag(1, -1), x) for x in xs]
    assert beta*beta == sp.eye(8)
    for x in alpha:
        assert sp.simplify(beta*x+x*beta) == sp.zeros(8)


@pytest.mark.parametrize('p', [(0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1), (3, 2, 1)])
def test_exact_massive_characteristic_polynomial(p):
    j, root = sp.I, sp.sqrt(3)
    c = sp.Matrix([[0, 1, 1, 1], [1, 0, j, -j], [1, -j, 0, j], [1, j, -j, 0]])/root
    signs = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    w = sp.diag(*[(-j)**sum(a*b for a, b in zip(s, p)) for s in signs])*c
    # Exactly rational cosine/sine: no rounded eigenvalue can prove this test.
    mass = sp.kronecker_product(sp.Matrix([[sp.Rational(3, 5), -4*j/5], [-4*j/5, sp.Rational(3, 5)]]), sp.eye(4))
    u = mass*sp.diag(w, w.conjugate().T)
    z = sp.Symbol('lambda')
    cos_eps_sq = 1-sp.Rational(sum(a % 2 for a in p), 3)
    expected = sp.expand(((z*z+1)**2-4*sp.Rational(9, 25)*cos_eps_sq*z*z)**2)
    assert sp.simplify(u.charpoly(z).as_expr()-expected) == 0
    assert sp.simplify(u.conjugate().T*u) == sp.eye(8)


@pytest.mark.parametrize('a', [.02, .005])
@pytest.mark.parametrize('k', [(0., 0., 0.), (.9, -.4, .7)])
def test_actual_gate_sequence_and_inverse(a, k):
    u = model.walk(k, a, .7)
    f = np.diag(np.exp(-1j*a*(model.SIGNS @ np.array(k))))
    pre, post = np.eye(8, dtype=complex), np.eye(8, dtype=complex)
    pre[:4, :4], post[4:, 4:] = model.C, model.C
    flight = np.zeros((8, 8), complex)
    flight[:4, :4], flight[4:, 4:] = f, f.conj()
    mu = np.sqrt(3)*a*.7/3
    mass = np.cos(mu)*np.eye(8)-1j*np.sin(mu)*model.BETA
    assert np.linalg.norm(mass@post@flight@pre-u) < 2e-15
    # Execute each reversed instruction, not just an abstract unitarity check.
    assert np.linalg.norm(pre.conj().T@flight.conj().T@post.conj().T@mass.conj().T@u-np.eye(8)) < 3e-15


def test_wrong_inverse_flight_is_detectable():
    u = model.walk([.7, -.2, .4], .5, 1.)
    w = model.walk([.7, -.2, .4], .5, 0.)[:4, :4]
    bad = np.zeros((8, 8), complex)
    bad[:4, :4], bad[4:, 4:] = w, w.conj()  # Entrywise conjugation is not the inverse.
    mu = np.sqrt(3)*.5/3
    bad = (np.cos(mu)*np.eye(8)-1j*np.sin(mu)*model.BETA)@bad
    assert np.linalg.norm(bad-u) > .5
    assert np.max(abs(np.sort(abs(np.angle(np.linalg.eigvals(bad))))-
                      np.sort(abs(np.angle(np.linalg.eigvals(u)))))) > .01


def test_odd_carrier_and_retiming_are_not_optional():
    m, a, n = 1., .01, 1
    tau = np.sqrt(3)*a/3
    correct = model.R@(np.cos(tau*m)*np.eye(8)-1j*np.sin(tau*m)*model.BETA)
    no_carrier = np.cos(tau*m)*np.eye(8)-1j*np.sin(tau*m)*model.BETA
    assert np.linalg.norm(model.walk([0, 0, 0], a, m)-correct) < 1e-14
    assert np.linalg.norm(correct-no_carrier) > 3
    n = 200
    actual = np.linalg.matrix_power(model.walk([0, 0, 0], a, m), n)
    wrong_time = np.cos(2*n*tau*m)*np.eye(8)-1j*np.sin(2*n*tau*m)*model.BETA
    assert np.linalg.norm(actual-wrong_time) > 1


def test_finite_clock_witness_is_informative(expected):
    witness = expected['analytic_clock']
    assert witness['probability_bound'] < .04
    assert witness['visibility'] > .74
    assert witness['cycle']/witness['gamma'] == pytest.approx(2*np.pi/.1)
    assert witness['sampling_error'] > 0 and witness['tail_norm'] > 0
    assert witness['interval_certificate']['observable_swing_lower'] == '0.6753'
    assert witness['visibility']-2*witness['probability_bound'] > .6753


def test_common_dirac_velocity_is_not_assumed_exact_on_the_lattice():
    # Different masses at common reference velocity have unequal native
    # velocities at finite spacing; both approach the reference under refinement.
    def speed(a, mass):
        tau = a/np.sqrt(3)
        k = 1.25*mass*.6
        s = np.sin(a*k)**2/3
        omega = np.arccos(np.cos(mass*tau)*np.sqrt(1-s))
        return (a/tau)*np.cos(mass*tau)*np.sin(2*a*k)/(6*np.sqrt(1-s)*np.sin(omega))
    coarse = [speed(.1, m) for m in (4., 4.2)]
    fine = [speed(.001, m) for m in (4., 4.2)]
    assert abs(coarse[0]-coarse[1]) > 1e-4
    assert max(abs(x-.6) for x in coarse) > .001
    assert max(abs(x-.6) for x in fine) < 1e-5


def test_finite_packets_keep_coherence_and_refine(expected):
    rows = expected['packets']['rows']
    coarse = max(abs(p-r) for row in rows[:3] for p, r in zip(row['probabilities'], row['rigid']))
    fine = max(abs(p-r) for row in rows[-3:] for p, r in zip(row['probabilities'], row['rigid']))
    assert coarse > .02 and fine < .002 and fine < coarse/10
    for row in rows:
        assert abs(row['probabilities'][0]-row['probabilities'][1]) > .1
        assert row['total_modes'] == 8*row['sites']
        assert max(row['reverse_errors']) < 2e-10
        assert row['final_energies'] == pytest.approx(row['initial_energies'], abs=2e-10)
        assert all(4 < e < 8 for e in row['initial_energies'])


@pytest.mark.parametrize('overlap', [0., .01, .3, .8, 1.])
def test_state_uniform_car_control(overlap):
    # All four Fock sectors of two modes, including fermionic signs.
    annihilators = []
    for mode in range(2):
        a = np.zeros((4, 4), complex)
        for state in range(4):
            if state & (1 << mode):
                a[state ^ (1 << mode), state] = (-1)**((state & ((1 << mode)-1)).bit_count())
        annihilators.append(a)
    f = annihilators[0]
    g = overlap*annihilators[0]+np.sqrt(1-overlap**2)*annihilators[1]
    nf, ng = f.conj().T@f, g.conj().T@g
    comm = np.linalg.norm(nf@ng-ng@nf, 2)
    assert comm <= 2*overlap+1e-14
    encoding = np.eye(4)-2*ng  # exp(i pi N_g)
    contrast = np.linalg.norm(encoding.conj().T@nf@encoding-nf, 2)
    assert contrast == pytest.approx(2*overlap*np.sqrt(1-overlap**2), abs=1e-14)


def test_finite_site_record_energy_cost():
    # A uniform Fourier density, and arbitrary internal polarization.
    # The lower bound uses the complete spectrum, not just an occupied band.
    for q in (8, 16):
        a, mass, c = 1/q, .5, 1.
        tau = np.sqrt(3)*a/c
        low = []
        for k in __import__('itertools').product(range(q), repeat=3):
            s = sum(np.sin(2*np.pi*j/q)**2 for j in k)/3
            low.append(np.arccos(np.cos(tau*mass)*np.sqrt(max(0, 1-s)))/tau)
        assert np.mean(low) >= np.arcsin(1/np.sqrt(6))/(2*tau)


def test_bounded_clock_visibility_is_not_a_global_phase():
    # Dephasing the two mass labels removes the measured quadrature difference.
    psi = np.array([1., np.exp(-.7j)])/np.sqrt(2)
    coherent = np.outer(psi, psi.conj())
    dephased = np.diag(np.diag(coherent))
    effects = [np.outer(np.array([1., np.exp(1j*t)])/np.sqrt(2),
                        np.array([1., np.exp(1j*t)]).conj()/np.sqrt(2)) for t in (0., np.pi/2)]
    signals = [np.trace(f@coherent).real for f in effects]
    no_signal = [np.trace(f@dephased).real for f in effects]
    assert abs(signals[0]-signals[1]) > .5
    assert no_signal == pytest.approx([.5, .5])


@pytest.mark.parametrize('mu', [.1, .5, 1.2])
def test_all_band_speed_limit_is_sharp(mu):
    # Analytic ray approaching a retained middle crossing, no low-band filter.
    values = []
    for delta in (.1, .01, .001):
        s = np.cos(delta)**2
        speed = np.cos(mu)*np.sqrt(s/(1-np.cos(mu)**2*(1-s)))
        values.append(speed)
        assert speed < np.cos(mu)
    assert np.cos(mu)-values[-1] < 1e-6
    assert values == sorted(values)


def test_compact_smear_signals_vanish_with_resolution():
    # 1D invariant waveguide sector, all lattice modes retained. The analytic
    # theorem handles compact support in 3D; this is a finite execution control.
    values = []
    for count in (256, 512, 1024, 2048, 4096):
        a = 64/count
        x = (np.arange(count)-count//2)*a
        f = np.maximum(0, 1-x*x)**3
        g = np.maximum(0, 1-(x-6.25)**2)**3
        f, g = f/np.linalg.norm(f), g/np.linalg.norm(g)
        spin = np.eye(8)[:, 0]
        initial = f[:, None]*spin
        kh = np.zeros((count, 3))
        kh[:, 0] = 2*np.pi*np.fft.fftfreq(count, d=a)
        n = round(4/(a/np.sqrt(3)))
        assert 6.25-2 > n*a/np.sqrt(3)  # Outside the field cone, inside native cone.
        final = check.binary_apply(check.gates(kh, a, 1., 3.), np.fft.fft(initial, axis=0, norm='ortho'), n)
        # Retain the known local carrier: an orthogonal fixed spin probe
        # would give a vacuous zero at every odd tick.
        probe = check.frames()[3]@spin if n % 2 else spin
        z = abs(np.vdot(g[:, None]*probe, np.fft.ifft(final, axis=0, norm='ortho')))
        values.append(z)
    assert values[0] > 1e-3 and values[-1] < 1e-8
    assert values[-1] < values[0]/10


MUTATIONS = ['missing_band', 'missing_grid', 'wrong_mode_count', 'wrong_mass', 'wrong_time',
             'wrong_clock', 'hidden_pi', 'erased_phase', 'fake_inverse', 'superluminal',
             'empty_catalog', 'bool_count', 'string_probability', 'extra_key', 'nan', 'infinity',
             'zero_spacing', 'zero_tail', 'zero_error_bound', 'fake_contrast']


def corrupt(evidence, name):
    if name == 'missing_band': evidence['spectra'][0]['energies'][0].pop()
    elif name == 'missing_grid': evidence['spectra'].pop()
    elif name == 'wrong_mode_count': evidence['spectra'][0]['modes'] //= 2
    elif name == 'wrong_mass': evidence['spectra'][0]['mass'] = 0.
    elif name == 'wrong_time': evidence['packets']['rows'][0]['time'] *= 2
    elif name == 'wrong_clock': evidence['analytic_clock']['gamma'] = .8
    elif name == 'hidden_pi': evidence['spectra'][0]['energies'][0][-1] = 0.
    elif name == 'erased_phase': evidence['packets']['rows'][0]['probabilities'] = [.5, .5]
    elif name == 'fake_inverse': evidence['packets']['rows'][0]['reverse_errors'][0] = .1
    elif name == 'superluminal': evidence['velocities'][0]['speeds'][0] = 3.
    elif name == 'empty_catalog': evidence['dynamics'] = []
    elif name == 'bool_count': evidence['dynamics'][0]['ticks'] = False
    elif name == 'string_probability': evidence['fast_records'][0]['empty_probability'] = '0.0'
    elif name == 'extra_key': evidence['passed'] = True
    elif name == 'nan': evidence['analytic_clock']['visibility'] = float('nan')
    elif name == 'infinity': evidence['analytic_clock']['visibility'] = float('inf')
    elif name == 'zero_spacing': evidence['analytic_clock']['a'] = 0.
    elif name == 'zero_tail': evidence['analytic_clock']['tail_norm'] = 0.
    elif name == 'zero_error_bound': evidence['analytic_clock']['probability_bound'] = 0.
    elif name == 'fake_contrast': evidence['analytic_clock']['interval_certificate']['observable_swing_lower'] = '0.7457'
    else: raise ValueError(name)


@pytest.mark.parametrize('name', MUTATIONS)
def test_hostile_catalogs_rejected(expected, name):
    bad = copy.deepcopy(expected)
    corrupt(bad, name)
    with pytest.raises(ValueError):
        check.same(bad, expected)


@pytest.mark.parametrize('name', ['missing_band', 'wrong_time', 'wrong_clock', 'erased_phase', 'bool_count', 'nan', 'zero_spacing', 'zero_tail'])
def test_optimized_cli_rejects_forgery(tmp_path, name):
    bad = verify.load(verify.HERE/'receipt.json')
    corrupt(bad['evidence'], name)
    path = tmp_path/'bad.json'
    path.write_text(json.dumps(bad), encoding='utf-8')
    run = subprocess.run([sys.executable, '-O', str(verify.HERE/'verify.py'), str(path)], capture_output=True, text=True)
    assert run.returncode != 0
    assert 'ValueError' in run.stderr


@pytest.mark.parametrize('text', ['{}', '[]', 'null', '{"schema":1,"schema":2}', '{"x":NaN}', '{"x":Infinity}'])
def test_junk_envelopes_fail_closed(tmp_path, text):
    path = tmp_path/'bad.json'
    path.write_text(text, encoding='utf-8')
    with pytest.raises((ValueError, TypeError)):
        verify.verify(path)


def test_source_tamper_rejected(tmp_path):
    packet = verify.load(verify.HERE/'receipt.json')
    packet['sources']['code/m1_operational_clocks/model.py'] = '0'*64
    path = tmp_path/'bad.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(ValueError, match='sources'):
        verify.verify(path)


def test_checker_does_not_import_producer():
    import ast
    tree = ast.parse((verify.HERE/'check.py').read_text(encoding='utf-8'))
    imports = [x.module for x in ast.walk(tree) if isinstance(x, ast.ImportFrom)]
    names = [n.name for x in ast.walk(tree) if isinstance(x, ast.Import) for n in x.names]
    assert not any(x and 'model' in x for x in imports+names)
    assert 'np.linalg.eig' not in (verify.HERE/'check.py').read_text(encoding='utf-8')


