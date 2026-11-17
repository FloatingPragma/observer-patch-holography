"""Reproducible intervention design; no empirical outcomes are consumed."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from ..observer.source_scalar import compiler
from . import response, source


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    'model/identification_design/response.py',
    'model/identification_design/source.py',
    'model/identification_design/experiment.py',
    'model/observer/source_scalar/compiler.py',
    'model/observer/scalar_response/experiment.py',
    'code/source_selection_model/response.json',
    'code/m1_source_realization/receipt.json',
)


def run():
    coherent = [source.execute(extra) for extra in (0., 1.)]
    committed = [source.execute(extra, commit_before_read=True) for extra in (0., 1.)]
    reference = [source.independent_reference(extra) for extra in (0., 1.)]
    residuals = [float(np.max(np.abs(np.array(row['complete_read_probabilities'])-law)))
                 for row, law in zip(coherent, reference, strict=True)]
    tv = lambda rows: float(np.sum(np.abs(np.array(rows[0]['complete_read_probabilities'])
                                         -np.array(rows[1]['complete_read_probabilities'])))/2)
    design = response.complete_design()
    background = source.background_slope_check()
    checks = {
        'exact_full_nongauge_rank_after_selected_probes': design['after']['rank'] == 5,
        'initial_calibration_has_two_invisible_action_directions': design['before']['rank'] == 3,
        'native_reference_agreement': all(
            np.isfinite(error) and 0 <= error < 2e-12 for error in residuals),
        'complete_native_evolution_matches_independent_H': all(
            np.isfinite(row['full_phase_evolution_target_residual'])
            and 0 <= row['full_phase_evolution_target_residual'] < 2e-12
            for row in coherent+committed),
        'complete_probability_laws': all(abs(row['complete_probability_mass']-1) < 2e-12
                                         for row in coherent+committed),
        'coherent_read_separates_candidates': tv(coherent) > .7,
        'commitment_control_erases_separation': tv(committed) < 2e-12,
        'background_slope_agrees': all(abs(row['twice_derivative_difference']
                                           -row['exact_binomial_response']) < 2e-12
                                      for row in background['records']),
    }
    if not all(checks.values()):
        raise ArithmeticError('Intervention design check failed')
    return {
        'schema': 'oph.identification_design.v1',
        'scope': 'component intervention design for supplied finite action candidates; no empirical observations or source-law selection',
        'supplied': ['eight-level code and coherent preparation/read services',
                     'known action background and binomial correction family through order five',
                     'calibrated model-coordinate duration',
                     'candidate pair c4=0 versus c4=1 for the locked finite-time witness',
                     'independent fresh trials for the sampling theorem'],
        'derived': ['exact response rank and nullspace, retaining global-energy gauge',
                    'minimum occupation additions within the declared nested family',
                    'complete native and independent reference detector probabilities',
                    'conditional finite-shot discrimination bound'],
        'input_sha256': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                         for path in SOURCE_PATHS},
        'source_controls': compiler.source_summary(),
        'exact_zero_compressions_at_cutoff_four': response.exact_blind_compression(),
        'response_design': design,
        'nonconstant_gap_inverse_exact': response.finite_difference_inverse(),
        'native_coherent_trials': coherent,
        'native_committed_controls': committed,
        'independent_reference_laws': [row.tolist() for row in reference],
        'maximum_source_reference_probability_errors': residuals,
        'coherent_total_variation': tv(coherent),
        'committed_total_variation': tv(committed),
        'sampling_design': response.bernoulli_design(source.TIME, 1.),
        'sampling_with_assumed_one_percent_probability_error': response.bernoulli_design(
            source.TIME, 1., trials=20, probability_error=.01),
        'short_time_derivative_sampling_control': {
            'time_parameter': 'supplied model coordinate', 'model_time': 1/128,
            'fresh_independent_trials': 200, 'failure_probability': .01,
            'twice_derivative_sampling_radius': response.derivative_noise_radius(1/128, 200, .01),
            'finite_time_bias_included': False,
            'purpose': 'rank does not itself establish finite-shot identification',
        },
        'noncommuting_background_derivative_check': background,
        'spatial_metric_selection': 'unchanged; the planar cube counterexample remains',
        'checks': checks, 'all_pass': all(checks.values()),
    }
