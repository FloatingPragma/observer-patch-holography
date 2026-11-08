"""Full-source, complete-instrument and interval-inversion controls."""
from copy import deepcopy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from scipy import sparse
from scipy.sparse.linalg import expm_multiply

from model.interacting_extension.higgs_join.interacting_carrier import build_join
from . import algebra, packet, readout, source


@pytest.fixture(scope='module')
def receipt():
    return packet.build()


def normalized_matrix(sector):
    rows, cols, values = [], [], []
    for (i, j), value in sector.numerators.items():
        rows.append(i)
        cols.append(j)
        values.append(float(value)/np.sqrt(float(sector.norms[i]*sector.norms[j])))
    return sparse.csr_matrix((values, (rows, cols)), shape=(sector.dimension,)*2)


def test_committed_receipt_replays_every_source_and_complete_record(receipt):
    saved = json.loads((packet.HERE/'results.json').read_text())
    assert packet.canonical(saved) == packet.canonical(receipt)
    replay = receipt['literal_source_replay']
    assert replay['basis_vectors_rebuilt'] == 2034
    assert replay['operator_columns_replayed'] == 20778
    assert replay['seed_pair_is_actual_reference_neutral_addition']
    for experiment, dimension in zip(receipt['finite_time_complete_experiments'], (438, 1596), strict=True):
        assert experiment['dimension'] == experiment['complete_record_rank_sum'] == dimension
        assert len(experiment['all_record_probabilities']) == 6
        for table in experiment['all_record_probabilities']:
            assert len(table['records']) == experiment['complete_record_count']
            bounds = [tuple(map(F, row['probability'])) for row in table['records']]
            assert sum(x[0] for x in bounds) <= 1 <= sum(x[1] for x in bounds)
        bound = list(map(F, experiment['bare_coefficient_readout']['electric_over_abs_magnetic']))
        assert bound[0] < 1 < bound[1] and bound[1]-bound[0] < F(1, 10**10)
    assert receipt['scope']['physical_alpha_s'] is None
    assert receipt['scope']['physical_particle_masses'] is None
    assert not receipt['scope']['flux_instrument_compiled_from_bare_repair_records']
    assert not receipt['scope']['colored_chiral_matter_and_refining_3plus1_field_constructed']


@pytest.mark.parametrize('q', (0, 1))
@pytest.mark.parametrize('control', (source.Controls(), source.Controls(F(3, 2), F(2, 3), F(5, 4), F(7, 5), F(3, 2))))
def test_full_exact_generator_matches_independent_retained_sparse_source(q, control):
    sector = source.Sector(q, control)
    actual = normalized_matrix(sector)
    expected = build_join(electric=float(control.electric), magnetic=float(control.magnetic),
        hopping=float(control.hopping), higgs_coupling=float(control.higgs),
        scalar_energy=float(control.scalar))['sectors'][q]['H']
    difference = actual-expected
    assert np.max(np.abs(difference.data), initial=0) < 4e-15
    assert float(np.max(np.asarray(abs(actual).sum(axis=1)))) <= float(sector.norm_bound)


@pytest.mark.parametrize('q', (0, 1))
def test_every_probability_encloses_full_scipy_evolution(q):
    sector = source.Sector(q)
    h = build_join()['sectors'][q]['H']
    seed = np.zeros(sector.dimension, complex)
    seed[sector.seed] = 1
    evolved = expm_multiply(-.125j*h, seed)
    intervals, radius = readout.complete_probabilities(sector, F(1, 8), order=20)
    assert radius > 0
    for record, indices in sector.records.items():
        p = float(np.vdot(evolved[indices], evolved[indices]).real)
        bound = intervals[record]
        assert float(bound.lo)-1e-14 <= p <= float(bound.hi)+1e-14
    assert abs(np.vdot(evolved, evolved)-1) < 2e-15
    # This entire probability lies outside the four-dimensional flux plane.
    assert sum(float(intervals[key].lo) for key in intervals if key[0] == 1) > .02


@pytest.mark.parametrize('q', (0, 1))
def test_complete_instrument_error_includes_an_entangled_reference(q):
    sector = source.Sector(q)
    h = build_join()['sectors'][q]['H']
    rng = np.random.default_rng(23+q)
    x = rng.standard_normal((sector.dimension, 2))+1j*rng.standard_normal((sector.dimension, 2))
    x /= np.linalg.norm(x)
    exact = expm_multiply(-.125j*h, x)
    polynomial, term = x.copy(), x.copy()
    for k in range(1, 13):
        term = (-.125j/k)*(h@term)
        polynomial += term
    trace_error = 0.
    for indices in sector.records.values():
        a, b = exact[indices].ravel(), polynomial[indices].ravel()
        aa, bb, ab = np.vdot(a, a).real, np.vdot(b, b).real, np.vdot(a, b)
        trace_error += np.sqrt(max(0., (aa+bb)**2-4*abs(ab)**2))
    r = float((sector.norm_bound/8)**13/factorial(13))
    assert trace_error <= 2*r+r*r+2e-7
    assert abs(np.linalg.norm(exact)-1) < 1e-14


def test_higgs_outcomes_and_sixth_order_interaction_are_not_deleted(receipt):
    for row, lost, sixth in zip(receipt['full_space_moment_identities'], (108, 90), (F(-81, 5), F(-66, 5)), strict=True):
        assert row['full_Higgs_loop_commutator_columns'] == row['dimension']
        assert F(row['quartic_probability_lost_if_Higgs_outcomes_removed']['0']) == lost
        assert F(row['Higgs_induced_sixth_order_color_change']['0']) == sixth
        assert row['loop_sign_involution_checked_on_every_generator_entry']


@pytest.mark.parametrize('q', (0, 1))
def test_clock_rescaling_preserves_exact_bare_ratio(q):
    original = source.Sector(q)
    rescaled = source.Sector(q, source.Controls(3, 3, 3, 3, 3))
    for block, _, _ in algebra.BLOCKS:
        a, b = original.moments(block), rescaled.moments(block)
        assert all(b[k] == 3**k*a[k] for k in a)
    def exact_ratio(sector):
        a, b = sector.moments(0), sector.moments(3)
        return -3*(a[4]/9-b[4]/9)/(7*(a[2]/9)**2)
    assert exact_ratio(original) == exact_ratio(rescaled) == 1


def test_exact_even_interpolation_annihilates_every_unwanted_power():
    for count in range(2, 9):
        for target, row in enumerate(readout.interpolation_weights(count), start=1):
            for power in range(1, count+1):
                assert sum(w*j**(2*power) for j, w in enumerate(row, start=1)) == int(power == target)


def test_nonunit_controls_are_inferred_from_complete_time_data():
    control = source.Controls(F(3, 2), F(2, 3), F(5, 4), F(7, 5), F(3, 2))
    result = readout.witness(1, controls=control)
    bounds = result['bare_coefficient_readout']
    for name, exact in [('magnetic_squared', F(4, 9)), ('electric_squared', F(9, 4)),
                        ('electric_over_abs_magnetic', F(9, 4))]:
        lo, hi = map(F, bounds[name])
        assert lo < exact < hi and hi-lo < F(1, 10**7)


@pytest.mark.parametrize('damage', ('labels', 'electric', 'norms', 'dimension'))
def test_raw_source_replay_rejects_comutated_metadata(monkeypatch, damage):
    data = deepcopy(source.inputs())
    if damage == 'labels':
        data[0]['0']['labels'][0]['link0_block'] = 0
    elif damage == 'electric':
        data[1]['q0_diagonal']['electric_native_trace'][0] = '1'
    elif damage == 'norms':
        data[0]['0']['raw_norms'][0] = '2'
    else:
        data[0]['0']['dimension'] += 1
    monkeypatch.setattr(source, 'inputs', lambda:data)
    with pytest.raises(ArithmeticError, match='source'):
        source.source_replay()


def test_zero_magnetic_or_nonexact_controls_do_not_produce_a_ratio():
    with pytest.raises(ValueError, match='nonzero'):
        readout.witness(0, controls=source.Controls(magnetic=0))
    with pytest.raises(TypeError, match='exact'):
        source.Controls(electric=1.)
    with pytest.raises(ValueError, match='positive denominator'):
        readout.Interval(1, 2).divide_positive(readout.Interval(-1, 1))


def test_outward_serialization_and_sampling_resource_bound():
    interval = readout.Interval(F(-1, 3), F(7, 11))
    lo, hi = map(F, interval.packet(bits=12))
    assert lo <= interval.lo <= interval.hi <= hi
    assert interval.lo-lo < F(1, 4096) and hi-interval.hi < F(1, 4096)
    tolerance = F(1, 2**80)
    n = readout.sample_budget(3, 12, tolerance, F(1, 10**6))
    assert 2*n*tolerance*tolerance >= 27


def test_changed_complete_record_or_physical_scope_is_rejected(receipt, tmp_path, monkeypatch):
    monkeypatch.setattr(packet, 'build', lambda:receipt)
    for mutation in ('drop', 'promote'):
        bad = deepcopy(receipt)
        if mutation == 'drop':
            bad['finite_time_complete_experiments'][0]['all_record_probabilities'][0]['records'].pop()
        else:
            bad['scope']['renormalized_QCD_coupling_identified'] = True
        path = tmp_path/(mutation+'.json')
        path.write_text(packet.canonical(bad))
        with pytest.raises(ValueError, match='complete probabilities'):
            packet.verify(path)


def test_bindings_cover_all_actual_local_runtime_imports():
    code = '''
from pathlib import Path
import sys
from model.color_loop_readout import packet
r=packet.build()
root=packet.ROOT.resolve()
actual=set()
for module in tuple(sys.modules.values()):
    name=getattr(module,'__file__',None)
    if name:
        p=Path(name).resolve()
        if p.is_relative_to(root) and p.suffix=='.py': actual.add(p.relative_to(root).as_posix())
assert 'model/interacting_extension/higgs_join/interacting_carrier.py' in actual
assert not actual-set(r['source_bindings']), sorted(actual-set(r['source_bindings']))
'''
    result = subprocess.run([sys.executable, '-c', code], cwd=packet.ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout+result.stderr
