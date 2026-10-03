"""Scientific and hostile controls; a digest match is not a scientific test."""
from copy import deepcopy
from fractions import Fraction as F
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import sympy as sp

from . import source, source_check, observables, observable_check, verify
from .format import load, rational, equal


@pytest.fixture(scope='module')
def packet():
    return verify.verify()


@pytest.fixture(scope='module')
def rows():
    return source_check.model()


def test_complete_independent_replay_and_rebuild(packet):
    actual = source.build()
    equal(actual, packet['evidence']['source'], 'deterministic full source reproduction')
    equal(observables.build(verify.HERE/'measurements.json', actual),
          packet['evidence']['observables'], 'deterministic complete observable reproduction')
    assert len(actual) == 16
    assert all(x['events'] == 2112 and x['reads'] == 11264 for x in actual)


def test_all_normalized_positive_actions_include_non_einstein_members():
    # Exact independent examples off the four displayed branches, including k!=1.
    for gamma in (F(-7, 3), F(-1, 5), F(0), F(7, 2)):
        for kappa in (F(1, 7), F(3, 2), F(9)):
            transformation = sp.Matrix([[1, 0], [-gamma, 1]])
            matrix = transformation.T*sp.diag(1, kappa)*transformation
            inverse = matrix.inv()
            assert matrix[0, 0] > 0 and matrix.det() > 0
            assert inverse[0, 0] == 1 and inverse[1, 0] == gamma
            assert -matrix[1, 0]/matrix[1, 1] == gamma
            assert F(1+gamma, 2) != 1


def test_nonlinear_source_and_schedule_controls_are_real(packet):
    items = packet['evidence']['source']
    assert items[0]['u'] != items[1]['u']  # nonlinear matter changes the field
    assert items[0]['projection'] != items[4]['projection']  # actual schedule
    assert items[0]['u'] == items[2]['u']  # w perturbation cannot change u
    assert set(items[0]['w']) == {'0'} and set(items[2]['w']) != {'0'}
    for item in items:
        assert max(map(F, item['u'])) > 0
        assert len(set(item['chains'].values())) == 4


def test_coercivity_does_not_imply_a_normalizable_gibbs_law():
    # A smooth confining energy tending to infinity still has a divergent
    # partition function. This defeats the tempting but false generalization.
    x, limit = sp.symbols('x limit', real=True, positive=True)
    energy = sp.log(1+x*x)/2
    assert sp.limit(energy, x, sp.oo) == sp.oo
    partition = 2*sp.integrate(sp.exp(-energy), (x, 0, limit))
    assert sp.simplify(partition-2*sp.asinh(limit)) == 0
    assert sp.limit(partition, limit, sp.oo) == sp.oo
    # The actual one-coordinate resting-source term instead has a Gaussian
    # lower bound after adding j, for every positive stiffness and source.
    stiffness, j = sp.symbols('stiffness j', positive=True)
    resting = stiffness*x*x/2+j*(sp.exp(-x)-1)
    assert sp.simplify(resting+j-stiffness*x*x/2) == j*sp.exp(-x)


def test_compact_density_has_only_the_claimed_regularity():
    r = sp.symbols('r', real=True)
    density = (1-r*r)**2  # zero continuation beyond the unit source radius
    assert density.subs(r, 1) == sp.diff(density, r).subs(r, 1) == 0
    assert sp.diff(density, r, 2).subs(r, 1) != 0


def test_finite_optical_separation_cannot_be_hidden_by_consistent_errors(packet):
    evidence = deepcopy(packet['evidence'])
    state = evidence['source'][0]
    reads = evidence['observables']['finite_reads'][0]
    state['error_intervals'][1] = ['0', '1/1000']
    eps_u, eps_w = [F(a[1]) for a in state['error_intervals']]
    for read in reads['optical']:
        read['log_stationary_error'] = str(abs(1+read['gamma'])*eps_u+eps_w)
    # Keep every propagated error internally consistent; rejection must be
    # the scientific failure to resolve the pair, not custody or bookkeeping.
    with pytest.raises(ValueError, match='stationary optical intervals must be disjoint'):
        observable_check.verify(verify.data(), evidence['observables'], evidence['source'])


@pytest.mark.parametrize('attack', ['clock', 'redshift_sign', 'metric', 'normalization', 'target', 'tensor',
                                  'drop_branch', 'sample_clock', 'source_writer', 'error',
                                  'nonlinear_ray', 'systematic', 'physical', 'boolean'])
def test_semantic_observable_attacks(packet, attack):
    evidence = deepcopy(packet['evidence'])
    p = evidence['observables']
    if attack == 'clock':
        p['finite_reads'][0]['clock_ratio'] = ['1', '1']  # conserved coordinate frequency, not local clock
    elif attack == 'redshift_sign':
        p['branches'][0]['clock_exponent'] = '1/20000'
    elif attack == 'metric':
        p['branches'][3]['optical_exponent'] = '1/10000'  # temporal potential alone
    elif attack == 'normalization':
        p['branches'][0]['stiffness'][0][0] = '3'
    elif attack == 'target':
        p['branches'][0]['comparisons'][0]['prediction'] = '24997/25000'
    elif attack == 'tensor':
        p['branches'][0]['exterior_tensor_coefficient'] = '0'
    elif attack == 'drop_branch':
        p['branches'].pop(0)
    elif attack == 'sample_clock':
        p['finite_reads'][0]['clock_exponent'] = str(F(2112, 11264))
    elif attack == 'source_writer':
        p['finite_reads'][0]['emit'][1] = 32
    elif attack == 'error':
        p['finite_reads'][0]['log_stationary_clock_error'] = '0'
    elif attack == 'nonlinear_ray':
        p['nonlinear_rays'][-1]['numerical_angle'] = p['nonlinear_rays'][-1]['linear_angle']
    elif attack == 'systematic':
        p['redshift']['profiled_statistical_residual'] = '5/38'
    elif attack == 'physical':
        p['physical_establishment'] = True
    elif attack == 'boolean':
        p['branches'][2]['gamma'] = True
    with pytest.raises(ValueError):
        observable_check.verify(verify.data(), p, evidence['source'])


@pytest.mark.parametrize('attack', ['strength', 'sampling', 'writer', 'fields', 'transverse',
                                  'residual', 'law', 'exponent', 'extra'])
def test_source_semantics_without_custody_shortcut(packet, rows, attack):
    item = deepcopy(packet['evidence']['source'][0])
    if attack == 'strength':
        item['strength'] = '1/16'
    elif attack == 'sampling':
        item['events'] *= 2
    elif attack == 'writer':
        item['chains']['-1'] = item['chains']['1']
    elif attack == 'fields':
        item['u'][21] = '0'
    elif attack == 'transverse':
        item['w'][21] = '1/10000'
    elif attack == 'residual':
        item['error_intervals'][0] = ['0', '1/10000000000000000']
    elif attack == 'law':
        item['nonlinear'] = True
    elif attack == 'exponent':
        item['denominator'] = 10**90  # rejected before allocation or iteration
    elif attack == 'extra':
        item['successful_outcomes_only'] = True
    with pytest.raises(ValueError):
        source_check.verify_case(item, F(1, 32), False, False, False, rows)


def test_self_consistent_changed_source_action_is_rejected(rows, monkeypatch):
    original_rows, mass = source.data()
    changed = [{j: 2*a for j, a in row.items()} for row in original_rows]
    monkeypatch.setattr(source, 'data', lambda: (changed, mass))
    # The producer recomputes all writes, reads, hashes and residuals coherently.
    # The checker must still reconstruct the declared primitive action.
    item = source.execute(F(1, 32), False, False, False)
    with pytest.raises(ValueError, match='source record'):
        source_check.verify_case(item, F(1, 32), False, False, False, rows)


@pytest.mark.parametrize('bad', ['', '1.0', '1e999999999', '01', '-0', '2/2', '1/0',
                                '1/01', True, 1, '9'*161])
def test_bounded_exact_number_parser(bad):
    with pytest.raises(ValueError):
        rational(bad)


@pytest.mark.parametrize('bad', ['{"x":1,"x":2}', '{"x":NaN}', '{"x":1.0}',
                                '{"x":'+('9'*101)+'}', ' '*1_000_001],
                         ids=['duplicate', 'nonfinite', 'float', 'large_integer', 'oversize'])
def test_strict_json(bad, tmp_path):
    p = tmp_path/'bad.json'
    p.write_text(bad, encoding='ascii')
    with pytest.raises(ValueError):
        load(p)


def test_data_target_substitution_rejected(tmp_path, monkeypatch):
    p = verify.data()
    p['deflection']['rows'][0]['gamma_minus_one'] = '-2'
    (tmp_path/'measurements.json').write_text(json.dumps(p), encoding='ascii')
    monkeypatch.setattr(verify, 'HERE', tmp_path)
    with pytest.raises(ValueError, match='primary measurement'):
        verify.data()


def test_custody_rejects_new_claim_or_source(packet, tmp_path):
    for key in ('claim', 'sources'):
        row = deepcopy(packet)
        row[key] = '0'*64 if key == 'claim' else {}
        p = tmp_path/(key+'.json')
        p.write_text(json.dumps(row), encoding='ascii')
        with pytest.raises(ValueError, match='custody'):
            verify.verify(p)


def test_optimized_verifier_without_producer_imports():
    script = '''
import sys, importlib.abc, importlib.machinery
from pathlib import Path
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in {'paired_gravity.source','paired_gravity.observables','paired_gravity.build',
                        'source_scalar_execution.scalar_execution_algebra'}:
            raise RuntimeError('producer imported: '+fullname)
sys.meta_path.insert(0, Block())
original_exec = importlib.machinery.SourceFileLoader.exec_module
blocked_files = {Path(p).resolve() for p in [
    'code/paired_gravity/source.py','code/paired_gravity/observables.py',
    'code/paired_gravity/build.py','code/source_scalar_execution/scalar_execution_algebra.py']}
def checked_exec(loader, module):
    if Path(loader.path).resolve() in blocked_files:
        raise RuntimeError('producer file loaded: '+loader.path)
    return original_exec(loader, module)
importlib.machinery.SourceFileLoader.exec_module = checked_exec
from paired_gravity.verify import verify, data
from paired_gravity.observable_check import verify as semantic
from paired_gravity.format import rational
p=verify()
p['evidence']['observables']['branches'][0]['clock_exponent']='1/20000'
for f in [lambda: semantic(data(),p['evidence']['observables'],p['evidence']['source']),
          lambda: rational('1e999999999')]:
    try: f()
    except ValueError: pass
    else: raise RuntimeError('hostile input accepted under -O')
print('producer-free optimized replay and rejections passed')
'''
    result = subprocess.run([sys.executable, '-O', '-c', script], cwd=verify.ROOT,
                            env={**os.environ, 'PYTHONPATH': str(verify.ROOT/'code')},
                            capture_output=True, text=True, timeout=300)
    assert result.returncode == 0, result.stdout+result.stderr


def test_legacy_module_import_order_is_preserved():
    script = '''
import sys
from pathlib import Path
sys.path.insert(0, str(Path('code/source_scalar_execution').resolve()))
import source_scalar_execution as legacy
from paired_gravity import source, source_check
if sys.modules['source_scalar_execution'] is not legacy:
    raise RuntimeError('legacy producer was replaced')
rows, mass = source.data()
if len(rows) != 64 or len(source_check.model()) != 64:
    raise RuntimeError('primitive support not reconstructed')
print('legacy producer and independent package imports coexist')
'''
    result = subprocess.run([sys.executable, '-c', script], cwd=verify.ROOT,
                            env={**os.environ, 'PYTHONPATH': str(verify.ROOT/'code')},
                            capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout+result.stderr
