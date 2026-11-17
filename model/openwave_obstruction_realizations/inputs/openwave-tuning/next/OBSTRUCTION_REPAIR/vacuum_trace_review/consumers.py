"""Additive stronger cut and same-source cross-preparation fidelity bounds."""
from pathlib import Path
from fractions import Fraction as F
import json
import math
import mpmath as mp

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'species_leptons/chiral_ring4_fidelity'
mp.iv.dps = 45
iv = mp.iv
DELTA = F(1, 6)  # Verification cut, not a Hamiltonian coefficient.


def q(x):
    return F(str(x))


def interval(x):
    return iv.mpf(x.numerator)/x.denominator


def up(x):
    return math.nextafter(float(x.b), math.inf)


def main():
    own = json.loads((HERE/'RESULT.json').read_text())
    prepared = json.loads((SOURCE/'inputs/prepared_results.json').read_text())
    coulomb = json.loads((SOURCE/'inputs/coulomb_results.json').read_text())
    old = json.loads((SOURCE/'RESULT.json').read_text())
    # The interval trace has already been calculated for offset1/8. An
    # increased cut only multiplies that rigorous bound by exp(beta dc).
    multiplier = iv.exp(iv.mpf(5)/4*interval(DELTA-F(1, 8)))
    witnesses = {}
    for r in own['vacuum_counts']:
        count = up(interval(F.from_float(r['exact_boundary_derivative_count_upper']))*multiplier)
        assert count < 2
        case = next(x for x in prepared if x['case'] == r['case'])
        shift = 4*abs(q(case['mass']))
        cut = q(case['variational_energy_upper'])+DELTA
        for row in coulomb['every_number_sector']:
            if row['N'] != 8:
                assert q(row['complete_source_lower'])-(0 if row['N'] in (0, 16) else shift) > cut
        witnesses[r['case']] = count
    outputs = []
    for case in prepared:
        family = 'zero' if case['mass'] == 0 else 'register_electron'
        sharp = next(x for x in prepared if x['case'] == family+'_K4')
        sigma4 = q(sharp['full_preparation_spectral_residual'])+F(1, 10**8)
        low = q(sharp['variational_energy'])-F(1, 10**6)-sigma4*sigma4/DELTA
        high = q(sharp['variational_energy_upper'])
        cut = high+DELTA
        mu_upper = q(case['variational_energy_upper'])
        sigma = q(case['full_preparation_spectral_residual'])+F(1, 10**8)
        # The residual bound is relative to THIS shared certified cut.
        residual_weight = min(F(1), sigma*sigma/(cut-mu_upper)**2)
        # E0 is in[low,high]. (mu-E0)/(cut-E0) decreases as E0 grows,
        # so the lower endpoint supplies the conservative weight bound.
        energy_weight = min(F(1), (mu_upper-low)/(cut-low))
        weight = min(residual_weight, energy_weight)
        shift = 4*abs(q(case['mass']))
        charged = min(q(r['complete_source_lower'])-shift for r in coulomb['every_number_sector'] if r['N'] in (7, 9))
        gap = charged-mu_upper
        assert gap > 0 and 0 <= weight <= 1
        prep = next(r for r in coulomb['prepared_response_consumers'] if r['case'] == case['case'])
        state_error = 2*iv.sqrt(interval(weight))/interval(gap)
        energy_error = interval(mu_upper-low)/interval(gap*gap)
        total = iv.mpf(str(prep['prepared_Green_error_upper']))+state_error+energy_error
        previous = next(r for r in old['vacuum_consumers'] if r['case'] == case['case'])
        assert up(total) < previous['actual_vacuum_Green_error_upper']
        outputs.append({'case': case['case'], 'larger_cut_eigen_count_upper': witnesses[case['case']],
            'shared_K4_cut_exact': str(cut), 'neutral_gap_lower_exact': str(DELTA),
            'same_source_ground_energy_interval_exact': [str(low), str(high)],
            'residual_weight_upper_exact': str(residual_weight), 'energy_weight_upper_exact': str(energy_weight),
            'ground_fidelity_lower_exact': str(1-weight), 'ground_fidelity_lower': float(1-weight),
            'actual_vacuum_Green_error_upper': up(total),
            'previous_frozen_Green_error_upper': previous['actual_vacuum_Green_error_upper']})
    result = {'schema': 'oph.root.vacuum_trace_consumers.v1',
        'larger_cut_offset_exact': str(DELTA), 'source_parameters_changed': False,
        'complete_source_shared_between_K2_and_K4': True,
        'all_four_counts_pass': True, 'consumers': outputs,
        'physical_lepton_pole_or_ratio_claim': False}
    (HERE/'CONSUMERS.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
