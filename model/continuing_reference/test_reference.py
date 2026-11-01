from fractions import Fraction as F
from itertools import product
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import mpmath as mp
import numpy as np
import pytest

from ..certified_reference import source
from ..charged_defect import dictionary
from ..field_symmetry.symmetry import operators
from ..full_fock_field.field import FieldGeometry
from ..reference_recovery.recovery import execute_complete
from . import algebra, packet, preparation as p


def test_independent_high_precision_complete_source_response():
    with mp.workdps(85):
        control = tuple(F(x) for x in p.one_mode_certificate()['literal_controls'])
        spec = source.specification()
        rational = lambda x: mp.mpf(x.numerator)/x.denominator
        matrix = mp.matrix(6)
        for i, j in product(range(6), repeat=2):
            for c, g in zip(control, spec, strict=True):
                re, im = g[i][j]
                matrix[i, j] += rational(c)*(rational(re[0])+mp.sqrt(5)*rational(re[1])
                                              +1j*(rational(im[0])+mp.sqrt(5)*rational(im[1])))
        unitary = mp.expm(matrix)
        a, b = abs(unitary[0, 0])**2, abs(unitary[1, 0])**2
        bounds = p.one_mode_certificate()
        assert rational(bounds['a_lower']) < a < rational(bounds['a_upper'])
        assert rational(bounds['b_lower']) < b < rational(bounds['b_upper'])
        assert sum(abs(unitary[j, 0])**2 for j in range(2, 6)) < rational(bounds['one_mode_abort_probability_upper'])
        assert rational(bounds['gamma_lower']) < mp.log(b/a) < rational(bounds['gamma_upper'])
        lo, hi = p.logarithm_three()
        assert rational(lo) < mp.log(3) < rational(hi)


@pytest.mark.parametrize('count', [1, 2, 3, 4])
def test_all_product_records_and_abort_mass(count):
    weights = np.array([np.vdot(v, v).real for v in p.numerical_exports()])
    mass, success = 0., 0.
    for records in product(range(6), repeat=count):
        address, amplitude = p.numerical_branch(records)
        probability = abs(amplitude)**2
        assert probability == pytest.approx(np.prod(weights[list(records)]), rel=2e-14, abs=1e-300)
        assert address == sum(min(r, 2)*3**j for j, r in enumerate(records))
        mass += probability
        if max(records) < 2:
            success += probability
    assert mass == pytest.approx(weights.sum()**count, abs=2e-14)
    assert success == pytest.approx((weights[0]+weights[1])**count, abs=2e-14)


def test_local_bound_has_no_total_volume_factor():
    a, b = p.resource_certificate(4, 4), p.resource_certificate(10000, 4)
    for key in ('initial_local_modular_operator_error', 'initial_local_proper_block_minimum_weight',
                'initial_local_complete_output_trace_error'):
        assert a[key] == b[key]
    assert a['abort_mass_upper'] < b['abort_mass_upper']
    assert a['source_pulses'] == 4 and b['source_pulses'] == 10000


@pytest.mark.parametrize('value', [0, -1, True, 1.5, 1000001])
def test_invalid_mode_counts(value):
    with pytest.raises(ValueError):
        p.resource_certificate(value)


@pytest.mark.parametrize('records', [[], [True], [6], [-1], [0.5]])
def test_invalid_record_labels(records):
    with pytest.raises(ValueError):
        p.numerical_branch(records)


def test_other_public_domain_boundaries():
    with pytest.raises(ValueError): p.resource_certificate(4, 5)
    with pytest.raises(ValueError): p.resource_certificate(200, 129)
    with pytest.raises(ValueError): p.modular_word(4097)
    with pytest.raises(ValueError): p.modular_word(4, True)
    with pytest.raises(ValueError): p.modular_word(4, float('nan'))


def test_exact_stationarity_premises_and_full_space_diagnostics():
    cert = algebra.exact_number_premises()
    assert cert['exact_H_and_K_number_commutators_zero']
    assert all(row['all_30_seams_unital'] for row in cert['complete_native_laws'])
    assert max(row['Frobenius_residual'] for row in algebra.numerical_stationarity()) < 1e-15


def test_nonuniform_preparation_is_not_automatically_stationary():
    h = operators(FieldGeometry((2, 1, 1))).H.matrix().toarray()
    local = [np.diag([.25, .75]), np.diag([.5, .5]),
             np.diag([.25, .75]), np.diag([.25, .75])]
    rho = np.array([[1.]])
    for value in reversed(local): rho = np.kron(rho, value)
    assert np.linalg.norm(h@rho-rho@h) > .01


def test_outside_outcome_is_required_for_complete_reference():
    rho = np.diag([3**j.bit_count()/4**8 for j in range(256)])
    ops = dictionary.full_extension(2)
    incomplete = sum(op@rho@op.conj().T for label, op in ops if label != ('outside',))
    assert np.trace(incomplete).real == pytest.approx(36/65536)
    full = incomplete+ops[-1][1]@rho@ops[-1][1].conj().T
    assert np.linalg.norm(full-rho) < 1e-15


def test_conditioned_reference_is_not_unconditioned_stationarity():
    rho = np.diag([3**j.bit_count()/4**8 for j in range(256)])
    op = dict(dictionary.full_extension(2))[('native', 0, 0)]
    output = op@rho@op.conj().T
    probability = np.trace(output).real
    distance = np.sum(abs(np.linalg.eigvalsh(output/probability-rho)))/2
    assert probability == pytest.approx(9/32768)
    assert distance == pytest.approx(16375/16384)


def test_complete_feedback_counterexample():
    cert = algebra.feedback_counterexample()
    assert cert['complete_retained_record_count'] == 33
    assert cert['reference_trace_distance'] == F(1, 12)
    assert cert['unitality_numerator_over4'] != cert['trace_preserving_Gram_numerator_over4']


def test_fixed_point_dirichlet_identity_independent_matrix_example():
    rng = np.random.default_rng(128)
    x = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4)); x += x.conj().T
    a, _ = np.linalg.qr(rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4)))
    b, _ = np.linalg.qr(rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4)))
    kraus = [a/math.sqrt(3), b*math.sqrt(2/3)]
    phi = sum(k@x@k.conj().T for k in kraus)
    left = sum(np.linalg.norm(k@x-x@k)**2 for k in kraus)
    right = 2*(np.linalg.norm(x)**2-np.vdot(x, phi).real)
    assert left == pytest.approx(right, rel=1e-13)


def test_exact_active_projector_and_energy_obstruction():
    result = algebra.fixed_reference_classification()
    assert result['active_K_trace'] == -135
    assert result['active_K_squared_trace'] == F(1815, 8)
    assert result['active_K_variance'] == F(13, 48) > F(1, 4)
    assert result['commutant_complex_dimension'] == 198
    assert sum(result['dark_energy_block_dimensions']) == 45


@pytest.mark.parametrize('parameter', [-1., 0., 1.])
def test_modular_clock_complete_all_input_source_execution(parameter):
    word = p.modular_word(2, parameter)
    result = execute_complete(word, 2, np.eye(4))
    target = np.diag(np.exp(1j*parameter*math.log(3)*np.array([0, 1, 1, 2])))
    assert np.linalg.norm(result.success-target, 2) < 2e-14
    assert len(result.failure_weights) == 4
    assert np.linalg.norm(result.success.conj().T@result.success+result.failure_gram-np.eye(4), 2) < 2e-14
    cert = p.modular_certificate(2, parameter)
    assert cert['complete_all_input_isometry_error'] < F(1, 10**12)
    assert all(len(result.failure(j)) == 4 for j in range(4))


def test_charge_clock_is_trivial_at_fixed_number():
    phases = [np.exp(1j*.37*j.bit_count()) for j in range(1 << 8) if j.bit_count() == 4]
    assert np.max(abs(np.array(phases)-phases[0])) == 0
    assert p.one_mode_certificate()['gamma_lower'] > 1
    assert algebra.aligned_spin_reference()['exact_target_fixed_by_field_and_each_aligned_complete_repair']


def test_receipt_custody_and_replay():
    retained = json.loads(packet.OUTPUT.read_text())
    fresh = packet.run()
    from ..observer.scalar_response.experiment import compare_receipt
    compare_receipt(retained, fresh)
    for name, digest in retained['source_sha256'].items():
        assert hashlib.sha256((packet.ROOT/name).read_bytes()).hexdigest() == digest


def test_cli_rejects_modified_receipt(tmp_path):
    retained = json.loads(packet.OUTPUT.read_text())
    retained['successful_field_substate_silently_normalized'] = True
    bad = tmp_path/'bad.json'; bad.write_text(json.dumps(retained))
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    result = subprocess.run([sys.executable, '-m', 'model.continuing_reference', '--check', str(bad)],
                            cwd=packet.ROOT, env=env, capture_output=True, text=True)
    assert result.returncode != 0
