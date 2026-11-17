"""Replay the prepared source/FIFO/modular-service composition.

Long preparation is executed as a sparse quantum channel with an entangled
reference. Exact channel and mixing certificates, rather than numerical
convergence, justify the uniform finite error bound. The exported archive is
specified by the complete instrument, never replaced by independent records.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

from tools import strict_json
from . import local_bias as local, mixing, continuing, modular_clock

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CONTRACT = HERE/'observer_campaign.json'
OUTPUT = HERE/'observer_results.json'
CAMPAIGN = 'OPH-LOCAL-REPAIR-CONTINUING-CLOCK-2026-10-28-v1'
# The finite service class changes only through a versioned contract revision.
CONTRACT_SHA256 = '85d58c58b3b8981f97c3101033ca4aeb83cd9cf5ce9807d0596f6b73bcfb79d7'


def contract(path=CONTRACT):
    value = strict_json.load(Path(path))
    if hashlib.sha256(Path(path).read_bytes()).hexdigest() != CONTRACT_SHA256:
        raise ValueError('the frozen continuing-observer contract changed')
    local.contract()
    return value


def source_channel(parameters):
    """Same complete branches as the actual FIFO, on all source matrix units."""
    catalogue = continuing.branches(parameters)
    population = [[F() for _ in range(66)] for _ in range(66)]
    rows, columns, weights = [], [], []
    for branch in catalogue:
        for i, j in branch.pairs:
            population[i][j] += branch.weight
            if branch.weight:
                for k, l in branch.pairs:
                    rows.append(66*j+l)
                    columns.append(66*i+k)
                    weights.append(float(branch.weight))
    if tuple(map(tuple, population)) != mixing.probability_kernel(parameters):
        raise ArithmeticError('FIFO instrument and preparation channel disagree')
    return coo_matrix((weights, (rows, columns)), shape=(66**2, 66**2)).tocsr()


def preparation_execution():
    """An actual 1600+1 attempt marginal, retaining an initial quantum R.

    Export tracing is explicit. Only nonzero blocks need storage; the complete
    zero/raw alphabet remains in continuing.branches. This finite arithmetic
    diagnostic is separate from the exact arbitrary-input error certificate.
    """
    p = mixing.canonical_parameters()
    states = local.basis()
    i = states.index(tuple(int(k in (0, 2)) for k in range(12)))
    j = states.index(tuple(int(k in (3, 4)) for k in range(12)))
    density = np.zeros((66, 66, 2, 2), complex)
    density[i, i, 0, 0] = density[j, j, 1, 1] = .5
    density[i, j, 0, 1], density[j, i, 1, 0] = -.5j, .5j
    channel = source_channel(p)
    for _ in range(1600):
        density = (channel@density.reshape(66**2, 4)).reshape(66, 66, 2, 2)
    target = np.einsum('ij,ab->ijab', np.diag(list(map(float, local.reference(p)))), np.eye(2)/2)
    difference = density-target
    cq_distance = 0.
    mass = 0.
    for branch in continuing.branches(p):
        if not branch.weight or not branch.pairs:
            continue
        inputs = [a for a, _ in branch.pairs]
        delta = float(branch.weight)*difference[inputs][:, inputs]
        # Output indices form an injection; relabeling its supported corner
        # does not change trace norm. R remains quantum within every block.
        matrix = delta.transpose(0, 2, 1, 3).reshape(2*len(inputs), 2*len(inputs))
        cq_distance += np.linalg.svd(matrix, compute_uv=False).sum()/2
        mass += float(branch.weight)*sum(np.trace(density[a, a]).real for a in inputs)
    if (not np.isfinite(cq_distance) or cq_distance > 1e-10
            or not abs(mass-1) < 1e-10
            or not np.allclose(np.einsum('iiab->ab', density), np.eye(2)/2, atol=1e-10, rtol=0)):
        raise ArithmeticError('executed quantum preparation or untouched reference failed')
    return {'source_input': '(|02,0>+i|34,1>)/sqrt(2)',
            'repair_attempts': 1600, 'FIFO_filling_attempts': 1,
            'complete_raw_slots_per_attempt': len(continuing.branches(p)),
            'source_and_FIFO_half_trace_diagnostic_upper': '1/10000000000',
            'normalization_and_reference_residual_upper': '1/10000000000',
            'arithmetic_diagnostic_is_the_uniform_proof': False,
            'exports_traced_only_in_stationary_comparison': True}


def composition_budget(clock_budget):
    if type(clock_budget) is not dict:
        raise ValueError('an original clock budget mapping is required')
    values = []
    for key in ('complete_instrument_halfdiamond_upper', 'channel_norm_upper'):
        value = clock_budget.get(key)
        if type(value) is not str:
            raise ValueError('canonical rational clock bounds required')
        number = F(value)
        if str(number) != value:
            raise ValueError('canonical rational clock bounds required')
        values.append(number)
    compilation, gain = values
    if not 0 <= compilation <= 1 or gain < 1:
        raise ValueError('nonnegative clock error and norm gain at least one required')
    delta = F(mixing.canonical_certificate()['half_diamond_upper'])
    minimum = min(local.reference(mixing.canonical_parameters()))
    epsilon = delta/minimum
    if not 0 <= epsilon < 1:
        raise ArithmeticError('preparation does not certify a faithful reference')
    modular = epsilon/(1-epsilon)
    return {'preparation_half_trace_upper': str(delta),
            'reference_minimum_eigenvalue': str(minimum),
            'relative_operator_sandwich_epsilon': str(epsilon),
            'actual_source_and_FIFO_reference_faithful_on_target_support': True,
            'reference_comparison_traces_preexisting_R_and_exports': True,
            'complete_protocol_target': '(Lambda tensor id_FIFO,R)(Omega_1 tensor rho_R); FIFO correlations are retained',
            'clock_source_readout_is_marginal_of_correlated_FIFO_extension': True,
            'normalized_rare_record_division_required': False,
            'modular_unitary_error_per_abs_parameter_upper': str(modular),
            'complete_nominal_clock_readout_half_trace_upper': str(gain*delta+compilation),
            'complete_actual_source_reference_clock_comparison_at_s1_upper': str(gain*(delta+modular)+compilation),
            'FIFO_modular_stability_domain': 'Only the unchanged supported direct-sum algebra; a probe leaving a record corner is not an implemented flow on that corner.',
            'physical_service_error_included': False,
            'physical_allowance_rule': 'Add declared joint apparatus error; pulse bound alone does not include load/export, storage, selection or record hardware.'}


def report():
    frozen = contract()
    mix, fifo, clock = mixing.report(), continuing.report(), modular_clock.report()
    # These are semantic joins, not matches between similarly named receipts.
    p = mixing.canonical_parameters()
    if p.port_weights != tuple(frozen['canonical']['port_weights']):
        raise ArithmeticError('prepared and declared reference differ')
    for c in (F(0), F(1, 2), F(1)):
        member = replace(p, c_weight=c)
        source_channel(member)
        continuing.certificate(member)
    paths = {'model/local_repair_family/observer.py', 'model/local_repair_family/test_observer.py',
             'model/local_repair_family/observer_campaign.json', 'model/local_repair_family/README.md',
             'tools/strict_json.py'}
    bindings = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in sorted(paths)}
    for child in (mix, fifo, clock):
        for name, digest in child['source_sha256'].items():
            if name in bindings and bindings[name] != digest:
                raise ArithmeticError('component source custody differs')
            bindings[name] = digest
    return {'schema': 'oph.local_repair_family.observer.v1', 'campaign_id': CAMPAIGN,
            'contract_sha256': CONTRACT_SHA256, 'source_sha256': dict(sorted(bindings.items())),
            'mixing': mix, 'continuing': fifo, 'modular_clock': clock,
            'composition_budget': composition_budget(clock['source_budget']),
            'preparation_execution': preparation_execution(),
            'same_channel_endpoint_and_tagged_mixture_join': True,
            'record_resources_through_preparation': continuing.resource_counts(1601),
            'record_resources_through_two_repairs_and_detector': continuing.resource_counts(1604),
            'scope': {'integrated_gate_accepted': False, 'physical_comparison_count_added': 0,
                      'physical_duration_or_energy_normalization_derived': False,
                      'coherent_services_derived_from_passive_repair': False,
                      'weights_source_selected': False}, 'remaining': frozen['remaining']}


def _compare(path, expected):
    actual = strict_json.load(Path(path))
    # Exact JSON comparison also distinguishes false from0 and rejects NaN.
    if json.dumps(actual, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
        raise ValueError('continuing observer receipt differs from complete replay')


def verify(path=OUTPUT):
    _compare(path, report())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--write', action='store_true')
    modes.add_argument('--check', nargs='?', const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.write:
        OUTPUT.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    if args.check:
        _compare(args.check, value)
    print('Prepared continuing observer: exact mixing, complete FIFO and source modular clock verified.')


if __name__ == '__main__':
    main()
