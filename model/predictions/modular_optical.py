"""Conditional assignment audit of native three-mode gaps to optical clocks.

The source generator is H=-log(T3). One shared positive conversion factor
cancels in frequency ratios, including a common Gibbs beta. All nine ordered
single-gap assignments are checked. This is a retrospective test of a stated
minimal attachment, not an atomic model or an exclusion of OPH.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from ..transfer_modular_probe import construction as source
from ..repair_identification import core as exact
from tools import strict_json

ROOT = source.ROOT
DATA = 'model/predictions/data/modular_optical_ratios.json'
SOURCE_PATHS = tuple(sorted(set((*source.SOURCE_PATHS, *exact.SOURCE_PATHS,
    DATA, 'model/predictions/__init__.py', 'model/predictions/modular_optical.py',
    'model/repair_identification/__init__.py',
    # Transitive initialization modules executed by the source-spectrum CLI.
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/claims.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/constants/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/constants/oph_pixel.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/array_geometry.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/echosahedral_dynamics.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/echosahedral_federation.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/graph.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/icosahedral.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/patch_state.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/patchnet.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/pixel_scale.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/screen_microphysics.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/core/screen_ports.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/anneal.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/consensus_certificate.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/kernel_dispatcher.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/modular.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/positive_geometry.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/repair.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/bundle.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/cross_repo_artifacts.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/hashes.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/particle_input_policy.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/production_envelope.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/evidence/verifier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/finite_groups.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/a5_sm_certificate.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/covariant_overlap.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/higgs_carrier.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/mar_sieve.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/physical_a5_sm_requirements.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/repair_projection.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/source_to_eft.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/gauge/yang_mills_gap.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/__init__.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/base.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/clock.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/registry.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/s3.py',
    'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/groups/z2.py',
    'model/__init__.py',
    'model/anharmonic_transfer/__init__.py',
    'model/common_source/__init__.py',
    'model/common_source/composition.py',
    'model/common_source/native_join.py',
    'model/common_source/operational.py',
    'model/common_source/records.py',
    'model/common_source/thermal_clock.py',
    'model/native_clock/__init__.py',
    'model/observer/__init__.py',
    'model/observer/p_response/__init__.py',
    'model/observer/p_response/experiment.py',
    'model/observer/scalar_response/__init__.py',
    'model/observer/source_discrimination/__init__.py',
    'model/observer/source_scalar/__init__.py',
    'model/observer/source_scalar/experiment.py',
    'model/predictions/protocol.py',
    'model/primitive_observer/__init__.py',
    'model/record_laws/__init__.py',
    'model/relativistic_scalar/__init__.py',
    'model/relativistic_scalar/lattice.py',
    'model/scalar_observations/__init__.py',
    'model/scalar_observations/experimental_readout/__init__.py',
    'model/scalar_observations/lattice.py',
    'model/scalar_observations/source.py',
))))
LABELS = ('E_plus', 'E_minus', 'E_minus-E_plus')
ERROR_BUDGET = F(1, 10**12)


def _fraction(value, *, positive=False):
    value = exact._fraction(value)
    if value < 0 or (positive and value == 0):
        raise ValueError('a canonical nonnegative exact fraction is required')
    return value


def _interval(value):
    if type(value) is not tuple or len(value) != 2:
        raise ValueError('an ordered positive exact interval is required')
    lo, hi = (_fraction(x, positive=True) for x in value)
    if hi < lo:
        raise ValueError('an ordered positive exact interval is required')
    return lo, hi


def outward(interval, digits=12):
    if type(digits) is not int or not 1 <= digits <= 30:
        raise ValueError('an original decimal precision from one through thirty is required')
    lo, hi = _interval(interval)
    scale = 10**digits
    lo, hi = lo*scale, hi*scale
    return (F(lo.numerator//lo.denominator, scale),
            F(-((-hi).numerator//(-hi).denominator), scale))


def gaps():
    """Outward enclosures of every positive gap, derived before data access."""
    energies = tuple(tuple(exact._fraction(x) for x in interval) for interval in source.energies())
    if energies[0] != (F(0), F(0)):
        raise ArithmeticError('the source transfer ground energy changed')
    plus, minus = energies[1:]
    difference = minus[0]-plus[1], minus[1]-plus[0]
    return tuple(outward(x) for x in (plus, minus, difference))


def assignment_intervals():
    """Nine assignments, with exact same-gap cancellation and no extra gains."""
    energy_gaps = gaps()
    return tuple((a, b, (F(1), F(1)) if i == j else outward((x[0]/y[1], x[1]/y[0])))
                 for i, (a, x) in enumerate(zip(LABELS, energy_gaps, strict=True))
                 for j, (b, y) in enumerate(zip(LABELS, energy_gaps, strict=True)))


def compare_ratio(measured, *, error=ERROR_BUDGET):
    """Deterministic interval compatibility under a declared total ratio error.

    The error is a hypothesis combining observation and attachment errors; a
    quoted standard uncertainty does not by itself imply this hard enclosure.
    This function reports no tail probability or statistical significance.
    """
    measured, error = _fraction(measured, positive=True), _fraction(error)
    if error >= measured:
        raise ValueError('the declared ratio error must be smaller than the central value')
    rows = []
    for numerator, denominator, (lo, hi) in assignment_intervals():
        separation = max(F(0), lo-(measured+error), measured-error-hi)
        rows.append({'numerator_gap': numerator, 'denominator_gap': denominator,
            'predicted_ratio_enclosure': [str(lo), str(hi)],
            'minimum_absolute_ratio_separation': str(separation),
            'intervals_disjoint': separation > 0})
    return {'measured_ratio': str(measured), 'declared_total_absolute_ratio_error': str(error),
        'all_nine_assignments_incompatible': all(x['intervals_disjoint'] for x in rows),
        'minimum_absolute_ratio_separation': str(min(F(x['minimum_absolute_ratio_separation']) for x in rows)),
        'assignments': rows}


def run():
    # Source gaps are evaluated without reading any observation or choosing a scale.
    energy_gaps = gaps()
    data = strict_json.load(ROOT/DATA)
    expected = {40: ('1.05776938758748094', '0.00000000000000011'),
                36: ('1.05776646273518748', '0.00000000000000013')}
    if (data['schema'] != 'oph.optical_absolute_ratios.v1'
            or [(x['isotope'], x['ratio'], x['standard_uncertainty']) for x in data['rows']]
            != [(k, *v) for k, v in expected.items()]):
        raise ValueError('the two declared primary-table observations changed')
    results = []
    for row in data['rows']:
        result = compare_ratio(F(row['ratio']))
        if not result['all_nine_assignments_incompatible'] or F(result['minimum_absolute_ratio_separation']) <= F(57, 1000):
            raise ArithmeticError('the minimal single-generator attachment comparison changed')
        results.append({'isotope': row['isotope'], 'quoted_standard_uncertainty': row['standard_uncertainty'], **result})
    return {'schema': 'oph.modular_optical_assignment.v1',
        'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCE_PATHS},
        'source_generator': 'H=-log(T3), spectrum {0,E_plus,E_minus}',
        'positive_energy_gap_enclosures': {key: list(map(str, value)) for key, value in zip(LABELS, energy_gaps, strict=True)},
        'attachment': {'frequency_relation': 'nu_A/nu_B=gap_A/gap_B plus total absolute ratio error',
            'one_common_positive_frequency_scale': True,
            'common_inverse_temperature_cancels': True,
            'same_gap_assignments_included': True,
            'independent_species_gains_offsets_harmonics_or_Floquet_log_branches': False,
            'atomic_spectrum_or_transition_dictionary_derived': False,
            'source_preparation_identified_with_measured_ions': False},
        'observation': {key: value for key, value in data.items() if key != 'rows'},
        'comparisons': results,
        'interpretation': {'status': 'retrospective conditional assignment audit',
            'each_ratio_compared_separately_no_independence_assumed': True,
            'quoted_standard_uncertainties_treated_as_hard_bounds': False,
            'p_value_or_sigma_significance_claimed': False,
            'frozen_prospective_prediction': False,
            'minimal_three_mode_direct_transition_attachment_incompatible': True,
            'OPH_or_all_modular_clocks_excluded': False}, 'all_pass': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    result = run()
    if args.check:
        if strict_json.load(args.check) != result:
            raise ValueError('the conditional modular optical assignment receipt does not replay')
        print('All direct native-gap optical assignments and source enclosures replayed.')
    else:
        print(json.dumps(result, sort_keys=True, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
