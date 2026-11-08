"""Independent complete-source and error-bound controls for dynamic readouts."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import types

import numpy as np
import pytest
from scipy.linalg import expm
import sympy as sp

from . import finite, packet


@pytest.fixture(scope='module')
def data():
    return finite.source_frames()


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def numeric_detector(data):
    blocks = []
    for columns, norms, frame in data['frames']:
        base = np.column_stack([np.array(c, dtype=float)/np.sqrt(float(n))
                                for c, n in zip(columns, norms, strict=True)])
        rotation = np.array([[float((v.lo+v.hi)/2) for v in row] for row in frame])
        blocks.append(base@rotation)
    return np.column_stack([sum((-1)**((k & x).bit_count())*blocks[k] for k in range(8))/np.sqrt(8)
                            for x in range(8)])


def test_exact_source_parity_and_all_character_ranks():
    h, scale, bases = finite.finite_basis()
    complement = list(range(255, -1, -1))
    assert np.array_equal(h[np.ix_(complement, complement)], h)
    assert tuple(len(columns) for columns, _ in bases) == (16,)*8
    assert scale > 0
    # Direct trace calculation gives16 independently of the chosen bases.
    actions = finite.algebra.reflections()
    for k in range(8):
        rank = F(sum((-1)**((g & k).bit_count())
                     *(sum(action[m] == m for m in range(256))
                       -sum((action[m] ^ 255) == m for m in range(256)))
                     for g, action in enumerate(actions)), 16)
        assert rank == 16


def test_complete_detector_covariance_preserves_source_first_axes(data):
    eta = numeric_detector(data)
    odd = (np.eye(256)-np.eye(256)[list(range(255, -1, -1))])/2
    assert np.linalg.norm(eta.T@eta-np.eye(128)) < 1e-12
    assert np.linalg.norm(eta@eta.T-odd) < 1e-12
    for shift, action in enumerate(finite.algebra.reflections()):
        order = [16*(x ^ shift)+j for x in range(8) for j in range(16)]
        assert np.linalg.norm(eta[list(action)]-eta[:, order]) < 1e-12
    probe = np.array(data['probe'], dtype=float)/np.sqrt(float(data['norm']))
    decoded = (eta.T@probe).reshape(8, 16)
    assert np.linalg.norm(decoded[:, 1:]) < 1e-12
    assert abs(np.sum(abs(decoded[:, 0])**2)-1) < 1e-12


def test_actual_clock_keeps_every_response_in_detection_sector(data, result):
    h = np.array(data['h'], dtype=float)/data['scale']
    eta = numeric_detector(data)
    odd = eta@eta.T
    assert np.linalg.norm(h@odd-odd@h) < 1e-12
    probe = np.array(data['probe'], dtype=float)/np.sqrt(float(data['norm']))
    retained = result['finite_source_readout']['same_clock_position_readouts']
    for t, certificate in zip((0., .125), retained, strict=True):
        evolved = expm(-1j*t*h)@probe
        assert np.linalg.norm(evolved-odd@evolved) < 1e-12
        probabilities = np.sum(abs((eta.T@evolved).reshape(8, 16))**2, axis=1)
        for probability, interval in zip(probabilities, certificate['position_probability_intervals'], strict=True):
            assert float(interval[0]) <= probability <= float(interval[1])
        assert certificate['dark_even_parity_probability_exact'] == '0'
    # The source time changes actual complete readout probabilities.
    lo0 = F(retained[0]['position_probability_intervals'][0][0])
    hi1 = F(retained[1]['position_probability_intervals'][0][1])
    assert lo0-hi1 > F(1, 250)


def test_original_first_mode_code_leaks_but_complete_response_does_not(data):
    eta = numeric_detector(data)
    first = eta[:, ::16]
    probe_projection = first@first.T
    h = np.array(data['h'], dtype=float)/data['scale']
    probe = np.array(data['probe'], dtype=float)/np.sqrt(float(data['norm']))
    assert np.linalg.norm((np.eye(256)-probe_projection)@h@probe) > .05
    assert data['rank8_code_squared_energy_leakage_lower'] > F(1, 100)


def test_source_position_effect_cannot_be_supported_on_proper_subset(data, result):
    eta = numeric_detector(data)[:, :16]
    entry = float((eta@eta.T)[1, 254])
    control = result['finite_source_readout']['single_site_locality_countercontrol']
    assert control['source_basis_masks'] == [1, 254]
    assert control['differing_sites'] == 8
    assert float(control['true_effect_matrix_entry_interval'][0]) <= entry <= float(control['true_effect_matrix_entry_interval'][1])
    assert abs(entry) > .3


@pytest.mark.parametrize('left,right', [((F(-3, 2), F(2, 3)), (F(-7, 4), F(5, 3))),
                                       ((F(1, 3), F(1, 2)), (F(2, 7), F(3, 4))),
                                       ((F(-2), F(-1)), (F(-1, 7), F(-1, 9)))])
def test_directed_interval_products_and_sums_enclose_exact_corners(left, right):
    a, b = finite.Interval(*left), finite.Interval(*right)
    product, total = a*b, a+b
    for x in left:
        for y in right:
            assert product.lo <= x*y <= product.hi
            assert total.lo <= x+y <= total.hi
    assert all((a-b).lo <= x-y <= (a-b).hi for x in left for y in right)


def test_interval_inverse_and_square_handle_zero_correctly():
    a = finite.Interval(F(-2), F(3))
    assert a.square().lo == 0 and a.square().hi == 9
    with pytest.raises(ValueError):
        a.reciprocal()
    b = finite.Interval(F(-3), F(-2)).reciprocal()
    assert b.lo <= F(-1, 2) <= F(-1, 3) <= b.hi
    with pytest.raises(ValueError):
        finite.Interval(F(2), F(1))


@pytest.mark.parametrize('value', [F(0), F(2), F(7, 11), F(10**30)])
def test_rational_root_intervals_bound_squared_value(value):
    interval = finite.root(value)
    assert interval.lo**2 <= value <= interval.hi**2


def test_taylor_remainder_bounds_actual_source_state(data):
    state, error = finite.taylor_state(data, F(1, 8), 6)
    numeric = np.array([complex(float((a.lo+a.hi)/2), float((b.lo+b.hi)/2)) for a, b in state])
    h = np.array(data['h'], dtype=float)/data['scale']
    probe = np.array(data['probe'], dtype=float)/np.sqrt(float(data['norm']))
    truth = expm(-.125j*h)@probe
    assert np.linalg.norm(truth-numeric) <= float(error)
    assert error > 0


def test_covariant_povm_variance_includes_multiplicity_noise():
    k, lam = sp.symbols('k lam', real=True)
    f = sp.exp(-k*k/2)
    vector = sp.Matrix([sp.cos(lam*k)*f, sp.sin(lam*k)*f])
    derivative = vector.diff(k)
    noise = sp.simplify((derivative.T*derivative)[0]-sp.diff(f, k)**2)
    assert noise == lam*lam*f*f


def test_receipt_does_not_promote_response_to_particle_mass(result):
    assert result['scope']['sharp_particle_mass_shell_or_observed_species_identified'] is False
    assert result['scope']['native_position_instrument_tape_executed'] is False
    changed = deepcopy(result)
    changed['scope']['sharp_particle_mass_shell_or_observed_species_identified'] = True
    with pytest.raises(ValueError, match='does not reproduce'):
        packet.validate(changed, result)
    changed = deepcopy(result)
    changed['scope']['dark_outcome_commutes_with_actual_time_generator'] = 1
    with pytest.raises(ValueError, match='does not reproduce'):
        packet.validate(changed, result)


def test_explicit_check_path_and_duplicate_keys_fail_before_replay(tmp_path):
    path = tmp_path/'supplied.json'
    path.write_text('{"schema":"oph.spin_dynamical_localization.v1","schema":"forged"}')
    bad = subprocess.run([sys.executable, '-m', 'model.spin_dynamical_localization', '--check', str(path)],
                         cwd=packet.ROOT, capture_output=True, text=True)
    assert bad.returncode != 0
    assert 'duplicate JSON key' in bad.stderr
    path.unlink()
    missing = subprocess.run([sys.executable, '-m', 'model.spin_dynamical_localization', '--check', str(path)],
                             cwd=packet.ROOT, capture_output=True, text=True)
    assert missing.returncode != 0
    assert str(path) in missing.stderr


def test_source_bindings_cover_actual_local_imports(monkeypatch):
    # Other owners in a combined pytest invocation may load unrelated modules
    # from a shared package. Their imports are not this producer's inputs.
    unrelated = types.ModuleType('oph_test_unrelated_spectrum_comparison')
    unrelated.__file__ = str(packet.ROOT/'model/centered_particle_spectrum/comparison.py')
    monkeypatch.setitem(sys.modules, unrelated.__name__, unrelated)
    code = '''
import sys
from pathlib import Path
from model.spin_dynamical_localization import packet

# Execute the complete owner before inspecting runtime imports: lazy imports
# used by the scientific computation must be included, not just packet imports.
pins = packet.produce()['source_bindings']
root = packet.ROOT.resolve()
runtime = set()
for module in tuple(sys.modules.values()):
    filename = getattr(module, '__file__', None)
    if filename:
        path = Path(filename).resolve()
        if path.is_relative_to(root) and path.suffix == '.py':
            runtime.add(path.relative_to(root).as_posix())
assert not runtime.difference(pins), sorted(runtime.difference(pins))
dependency = 'model/centered_particle_spectrum/source.py'
assert dependency in runtime
missing_dependency = dict(pins)
del missing_dependency[dependency]
assert runtime.difference(missing_dependency) == {dependency}
assert not any(name.startswith(('docs/', 'paper/', 'flagship/')) for name in pins)
'''
    completed = subprocess.run([sys.executable, '-c', code], cwd=packet.ROOT,
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout+completed.stderr
