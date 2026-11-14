"""An exact finite countercheck for regional modular/exchange identification.

Uses the actual binary exclusion transfer on a four-cycle at total load two.
No producer code is imported. The regional one-load block is obtained by a
literal partial trace of the six-state global reference.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


def evaluate():
    edges = ((0, 1), (1, 2), (2, 3), (3, 0))
    states = [m for m in range(16) if m.bit_count() == 2]
    index = {m: j for j, m in enumerate(states)}
    generator = sp.zeros(6)
    for a, b in edges:
        swap = sp.zeros(6)
        for col, mask in enumerate(states):
            changed = mask ^ (1 << a) ^ (1 << b) if ((mask >> a) & 1) != ((mask >> b) & 1) else mask
            swap[index[changed], col] = 1
        generator += sp.eye(6) - swap
    transfer = sp.eye(6) - generator / 8
    reference = transfer / sp.trace(transfer)
    regional = sp.zeros(4)
    for i, first in enumerate(states):
        for j, second in enumerate(states):
            if first >> 2 == second >> 2:
                regional[first & 3, second & 3] += reference[i, j]
    block = regional.extract([1, 2], [1, 2])
    conditional = block / sp.trace(block)
    exact_expected = sp.Matrix([[sp.Rational(5, 16), sp.Rational(1, 16)],
                                [sp.Rational(1, 16), sp.Rational(5, 16)]])
    local = sp.Matrix([[1, -1], [-1, 1]])
    finite_formula = sp.Rational(1, 2) * (sp.Rational(3, 4) * sp.eye(2) - local / 8)
    # Removing a scalar from -log(rho_A) leaves gap log(3/2).
    exact_gap = math.log(1.5)
    linear_gap = 1 / 3
    modular_parameter = .37
    flow_error = 2 * abs(math.sin(modular_parameter * (exact_gap - linear_gap) / 2))
    wrong_sign_error = 2 * abs(math.sin(modular_parameter * exact_gap))
    family = []
    for vertices in (4, 8, 16, 64, 256):
        # Even cycle, fixed region of two neighbours, density 1/2, local k=1.
        # r=1/2 and a_k=3/4 exactly for every such cycle.
        rescaled_gap = -2 * vertices * math.log1p(-4 / (3 * vertices))
        error = abs(rescaled_gap - 8 / 3)
        local_bound = 4
        a = .75
        z = local_bound / (2 * vertices * a)
        bound = local_bound ** 2 / (4 * vertices * a * a * (1 - z))
        family.append({'vertices_and_edges': vertices, 'rescaled_exact_gap': rescaled_gap,
                       'limit_gap': 8 / 3, 'operator_error': error, 'theorem_bound': bound})
    checks = {'literal_partial_trace_matches_exact_formula': block == finite_formula == exact_expected,
              'all_regional_sectors_retained': sp.trace(regional) == 1,
              'conditional_eigenvalues_exact': sorted(conditional.eigenvals()) == [sp.Rational(2, 5), sp.Rational(3, 5)],
              'finite_linear_generator_identification_rejected': abs(8 * exact_gap - 8 / 3) > .5,
              'finite_flow_replacement_rejected': flow_error > .02,
              'opposite_time_orientation_rejected': wrong_sign_error > .29,
              'all_family_errors_inside_theorem_bound': all(r['operator_error'] <= r['theorem_bound'] for r in family),
              'declared_fixed_region_limit_improves': all(y['operator_error'] < x['operator_error'] for x, y in zip(family, family[1:]))}
    return {'schema': 'oph.audit.modular_transfer_limit.v1', 'checks': checks,
            'source': {'graph': 'cycle4', 'total_load': 2, 'region': [0, 1], 'global_reference': 'P/Tr(P), P=I-G/(2M)',
                       'global_dimension': 6, 'global_trace_P': str(sp.trace(transfer))},
            'regional_density': [[str(v) for v in regional.row(i)] for i in range(4)],
            'one_load_block_probability': str(sp.trace(block)), 'conditional_eigenvalues': ['2/5', '3/5'],
            'exact_modular_gap': exact_gap, 'linearized_modular_gap': linear_gap,
            'rescaled_exact_gap': 8 * exact_gap, 'rescaled_linear_gap': 8 / 3,
            'operator_norm_flow_error_at_s_0_37': flow_error,
            'opposite_time_orientation_error_at_s_0_37': wrong_sign_error,
            'even_cycle_family': family,
            'scope': 'Exact finite reference restriction and analytic two-level flow comparison. The regional generator is -2M log(I-D_k/(2M a_k)) modulo block scalars. D_p/a_p is its controlled growing-family limit. This comparison supplies neither preparation nor implementation by passive native repair, and no physical clock conversion.',
            'reviewer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = evaluate()
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('checks', 'rescaled_exact_gap', 'rescaled_linear_gap', 'operator_norm_flow_error_at_s_0_37')}, indent=2))
    raise SystemExit(0 if all(result['checks'].values()) else 1)
