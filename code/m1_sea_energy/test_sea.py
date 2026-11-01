"""Semantic, hostile and representation-equivalence controls for the sea result."""

from copy import deepcopy
import fnmatch
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy.linalg import expm
from . import build, verify, spectrum, spectrum_check, preparation, preparation_check
from . import clock_check, bounds_check
from .format import exact, unpack, pack


@pytest.fixture(scope='module')
def evidence():
    return verify.load(verify.HERE/'receipt.json')['evidence']


def test_source_pinned_receipt_and_all_parents():
    verify.verify()


def test_independently_replay_fresh_producer():
    verify.verify_evidence(build.candidate())


def test_all_checkers_work_without_producers_and_under_optimization():
    script = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname in {
            'm1_sea_energy.spectrum', 'm1_sea_energy.preparation',
            'm1_sea_energy.clock', 'm1_sea_energy.bounds', 'm1_sea_energy.build',
            'm1_operational_clocks.model', 'm1_fermionic_source.model',
            'm1_source_realization.model', 'm1_fermionic_source.walk'}:
            raise RuntimeError('forbidden producer import: '+fullname)
sys.meta_path.insert(0, Block())
from m1_sea_energy import verify
verify.verify()
bad = verify.load(verify.HERE/'receipt.json')['evidence']['bounds']
bad['storage_time'][0] = 0
try:
    verify.bounds_check.verify(bad)
except ValueError:
    pass
else:
    raise RuntimeError('optimized interpreter accepted false lifetime')
'''
    env = dict(os.environ, PYTHONPATH=str(verify.ROOT/'code'), OPENBLAS_NUM_THREADS='1')
    subprocess.run([sys.executable, '-O', '-c', script], check=True, env=env, timeout=300)


@pytest.mark.parametrize('group,check', [('spectrum', spectrum_check.verify),
    ('preparation', preparation_check.verify), ('clock', clock_check.verify), ('bounds', bounds_check.verify)])
@pytest.mark.parametrize('trash', [None, False, 0, '', [], {}])
def test_rejects_empty_or_wrong_type_roots(group, check, trash):
    with pytest.raises(ValueError):
        check(trash)


def changed(evidence, path, value):
    row = deepcopy(evidence[path[0]])
    target = row
    for key in path[1:-1]:
        target = target[key]
    target[path[-1]] = value
    return row


MUTATIONS = [
    (('spectrum', 'blocks'), []),
    (('spectrum', 'blocks', 0, 'a'), 0.),
    (('spectrum', 'blocks', 0, 'mass'), True),
    (('spectrum', 'blocks', 0, 'energies'), [0.]*8),
    (('spectrum', 'blocks', 0, 'p', 0, 0, 0), 0.),
    (('spectrum', 'blocks', 0, 'h', 0, 1, 0), 7.),
    (('spectrum', 'blocks', 0, 'logz', 0), 0.),
    (('spectrum', 'blocks', 0, 'phase_energy'), 0.),
    (('spectrum', 'blocks', 0, 'excitations'), 0.),
    (('spectrum', 'blocks', 0, 'unitary', 0, 0, 0), float('nan')),
    (('spectrum', 'boundary'), []),
    (('spectrum', 'boundary', 1, 'gap'), 0.),
    (('spectrum', 'boundary', 1, 'central_columns', 0, 0, 0), .25),
    (('spectrum', 'boundary', 1, 'modes'), 1),
    (('spectrum', 'projection'), []),
    (('spectrum', 'projection', 0, 'error'), 0.),
    (('spectrum', 'projection', 0, 'bound'), 10.),
    (('spectrum', 'projection', 0, 'ordinary_dirac_error'), 0.),
    (('preparation', 0, 'branches'), []),
    (('preparation', 0, 'branches', 0, 'probability'), 1.),
    (('preparation', 0, 'branches', 1, 'outcomes'), [0]),
    (('preparation', 0, 'branches', 1, 'density', 0, 0, 0), .5),
    (('preparation', 0, 'gates'), []),
    (('preparation', 0, 'gates', 0, 'rotations', 0, 'theta'), .73),
    (('preparation', 0, 'gates', 0, 'scalar'), [0., 0.]),
    (('preparation', 0, 'read'), []),
    (('preparation', 0, 'qnd'), []),
    (('preparation', 0, 'kraus'), []),
    (('preparation', 0, 'kraus', 1, 0, 0, 0), .71),
    (('preparation', 0, 'probabilities'), [0., 1.]),
    (('preparation', 0, 'mean_read_energy'), 0.),
    (('preparation', 0, 'pulse_rotations'), 0),
    (('clock', 'packets', 'rows'), []),
    (('clock', 'packets', 'rows', 0, 'background'), 0.),
    (('clock', 'packets', 'rows', 0, 'blocked'), 0.),
    (('clock', 'packets', 'rows', 0, 'energy'), -1.),
    (('clock', 'packets', 'rows', 0, 'probabilities'), [1., 0.]),
    (('clock', 'packets', 'rows', 1, 'lab_probabilities'), [1., 0.]),
    (('clock', 'packets', 'rows', 1, 'lab_rigid'), [1., 0.]),
    (('clock', 'packets', 'rows', 0, 'norm'), .5),
    (('clock', 'witness', 'a'), 0.),
    (('clock', 'witness', 'boundary_error'), 0.),
    (('clock', 'witness', 'boundary_margin'), 1),
    (('clock', 'witness', 'interval', 'probability_upper'), '1'),
    (('clock', 'witness', 'interval', 'observable_swing_lower'), '1'),
    (('bounds', 'storage_time'), [1, 1, 1]),
    (('bounds', 'joint_error'), [-64, 0, 0]),
    (('bounds', 'energy_error'), [-55, 0, 0]),
    (('bounds', 'schedules'), []),
    (('bounds', 'schedules', 0, 'levels'), 1),
    (('bounds', 'schedules', 0, 'word_bits'), 0),
    (('bounds', 'routing', 0, 'phases'), False),
    (('bounds', 'routing', 4, 'swaps_upper'), 0),
]


@pytest.mark.parametrize('path,value', MUTATIONS)
def test_semantic_corruptions_fail(evidence, path, value):
    check = dict(spectrum=spectrum_check.verify, preparation=preparation_check.verify,
                 clock=clock_check.verify, bounds=bounds_check.verify)[path[0]]
    with pytest.raises(ValueError):
        check(changed(evidence, path, value))


@pytest.mark.parametrize('rank', [2, 3])
def test_occupied_orbital_gauge_is_accepted(rank):
    # Different orbital representatives and tapes, same complete Slater state.
    generator = np.zeros((4, 4), complex)
    generator[:rank, :rank] = .23*np.ones((rank, rank))
    generator[rank:, rank:] = -.41*np.eye(4-rank)
    row = preparation.case(rank, expm(1j*generator))
    preparation_check.verify_case(row, rank)


def test_isometry_global_phase_is_accepted(evidence):
    rows = deepcopy(evidence['preparation'])
    for row in rows:
        row['isometry'] = pack(np.exp(.73j)*unpack(row['isometry'], (16, 8)))
    preparation_check.verify(rows)


def test_same_spectrum_wrong_eigenvectors_are_rejected(evidence):
    row = deepcopy(evidence['spectrum'])
    h = unpack(row['blocks'][0]['h'], (8, 8))
    # Spectral scalars alone would miss this corruption of the generator.
    row['blocks'][0]['h'] = pack(np.diag(np.linalg.eigvalsh(h)))
    with pytest.raises(ValueError):
        spectrum_check.verify(row)


def test_high_band_cannot_be_replaced_by_ordinary_dirac_sea(evidence):
    row = deepcopy(evidence['spectrum'])
    first = row['blocks'][0]
    h0 = __import__('m1_operational_clocks.model', fromlist=['hamiltonian']).hamiltonian(first['k'], first['mass'])
    first['p'] = pack((np.eye(8)-h0/first['mass'])/2)
    with pytest.raises(ValueError):
        spectrum_check.verify(row)


@pytest.mark.parametrize('raw', ['{"schema":1,"schema":2}', '{"value":NaN}', '{"value":Infinity}', '{"value":-Infinity}'])
def test_json_duplicates_and_nonfinite_tokens_rejected(tmp_path, raw):
    path = tmp_path/'bad.json'
    path.write_text(raw, encoding='utf-8')
    with pytest.raises(ValueError):
        verify.load(path)


@pytest.mark.parametrize('field', ['sources', 'claims'])
def test_custody_is_checked_before_semantics(tmp_path, field):
    packet = verify.load(verify.HERE/'receipt.json')
    packet[field].pop(next(iter(packet[field])))
    path = tmp_path/'bad.json'
    path.write_text(json.dumps(packet), encoding='utf-8')
    with pytest.raises(ValueError, match='custody'):
        verify.verify(path)


@pytest.mark.parametrize('x', [False, float('nan'), float('inf'), 0., 1e-30])
def test_positive_small_inputs_use_relative_checks(x):
    with pytest.raises(ValueError):
        exact(x, 1e-10)


@pytest.mark.parametrize('levels', [64, 2**20])
def test_protection_exponent_is_bounded_before_exponentiation(evidence, levels):
    row = deepcopy(evidence['bounds'])
    row['schedules'][0]['levels'] = levels
    # Keep even a regressed implementation's test allocation small, while
    # requiring rejection at the guard before the ceiling calculation.
    with pytest.raises(ValueError, match='bounded quantum exponent'):
        bounds_check.verify(row)


def test_boundary_weighted_resolvent_estimate_on_actual_walk():
    u = spectrum_check.closed_walk(3, .13)
    k = (u.conj().T-u)/(2j)
    delta = np.sin(.7*.13/np.sqrt(3))
    s = np.log1p(delta/2)
    distances = np.repeat([max(abs(x-1), abs(y-1), abs(z-1)) for x in range(3) for y in range(3) for z in range(3)], 16)
    weight = np.exp(s*distances)
    moved = weight[:, None]*k/weight[None, :]
    assert np.linalg.norm(moved-k, 2) <= np.expm1(s)+1e-12
    for t in (0., .1, 1.):
        inverse = np.linalg.inv(moved-1j*t*np.eye(len(k)))
        assert np.linalg.norm(inverse, 2) <= 2/np.sqrt(delta*delta+t*t)+1e-10


def test_fock_density_matches_phase_energy_formula():
    from m1_operational_clocks.check import gates
    u = gates([.2, .1, -.3], .17, .8, 3.)
    h, p = spectrum_check.sea(u, .17/np.sqrt(3))
    vals, vecs = np.linalg.eigh(h)
    bits = sum(1 << j for j in range(8) if vals[j] < 0)
    vacuum = spectrum_check.exterior(vecs)[:, bits]
    f = np.exp(.31j*np.arange(8))/np.sqrt(8)
    v = np.eye(8)+(np.exp(.61j)-1)*np.outer(f, f.conj())
    excited = spectrum_check.exterior(v)@vacuum
    generator = spectrum_check.occupations(h)-np.trace(h@p)*np.eye(256)
    energy = np.vdot(excited, generator@excited).real
    assert abs(energy-np.trace(h@(v@p@v.conj().T-p)).real) < 1e-8
    assert energy > 0


def test_clock_witness_is_an_interval_certificate():
    assert clock_check.interval_certificate() == dict(probability_upper='0.038', observable_swing_lower='0.923')


def test_old_envelope_only_probe_fails_the_corrected_read_contract(evidence):
    old = deepcopy(evidence['clock'])
    for row in old['packets']['rows']:
        row['probabilities'] = row['lab_probabilities'][:]
    with pytest.raises(ValueError, match='complete sea-read probabilities'):
        clock_check.verify(old)


def test_finite_clock_control_distinguishes_both_phase_laws(evidence):
    rows = evidence['clock']['packets']['rows']
    final = rows[-1]
    assert final['sites'] == 1024
    assert abs(final['probabilities'][0]-final['rigid'][0]) < .006
    assert abs(final['lab_probabilities'][0]-final['rigid'][0]) > .6
    assert abs(final['lab_probabilities'][0]-final['lab_rigid'][0]) < .025


def test_uv_branch_shift_preserves_smooth_energy_and_its_limit():
    from m1_operational_clocks.check import gates, frames
    _, _, _, carrier, _ = frames()
    f = (np.eye(8)+carrier)[:, 0]
    f /= np.linalg.norm(f)
    for a in (.2, .1, .05, .025):
        tau = a/np.sqrt(3)
        h, p = spectrum_check.sea(gates([.3, -.4, .7], a, .7, 3.), tau)
        ev, z = np.linalg.eigh(h)
        high = (z*(abs(ev) > np.pi/(2*tau)))@z.conj().T
        uv = h+2*np.pi/tau*(np.eye(8)-2*p)@high
        excess = float((2*np.pi/tau)*np.vdot(high@f, high@f).real)
        assert 0 < excess <= (8*tau/np.pi)*np.linalg.norm(h@f)**2
        assert np.linalg.norm(uv@f) <= 5*np.linalg.norm(h@f)
        assert np.linalg.norm(h@f) < 2.


def test_ci_covers_every_pinned_input_on_both_platforms():
    import yaml
    workflow = yaml.load((verify.ROOT/'.github/workflows/m1-sea-energy.yml').read_text(encoding='utf-8'), Loader=yaml.BaseLoader)
    for event in ('push', 'pull_request'):
        patterns = workflow['on'][event]['paths']
        assert all(any(fnmatch.fnmatchcase(path, pattern) for pattern in patterns) for path in verify.SOURCES)
    job = workflow['jobs']['controls']
    assert job['strategy']['matrix']['os'] == ['ubuntu-latest']
    assert job['steps'][-1]['run'] == 'python -m pytest -q code/m1_sea_energy'
    assert not job.get('continue-on-error')
    assert not any(step.get('continue-on-error') for step in job['steps'])
