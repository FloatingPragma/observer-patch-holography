"""Complete energy-bin laws on the original interacting eight-cell source."""
from fractions import Fraction as F
from functools import lru_cache

from model.native_charge_resolution.bounds import data as native_data, frac, root, round_up
from model.coherent_defect_wave.chiral_joint.interacting_ground_preparation.model import certificate as preparation

from .local import CHANNELS
from .response import positive_response


def down(value):
    return -round_up(-F(value))


def encoded(pair):
    return [str(down(pair[0])), str(round_up(pair[1]))]


@lru_cache(maxsize=None)
def windows(family, multiplier=128):
    source = native_data()
    local = positive_response(family)
    prep = preparation()['finite_energy_agreement_filter']
    eta = round_up(source['eta'])
    v = source['v']
    elo, ehi = source['e']
    delta = round_up(F(prep['conditional_one_source_trace_distance_upper']))
    qlo = down(F(prep['success_probability_lower']))
    radius = multiplier*v
    energies = [(down(frac(row['energy_interval'][0])+7*elo+4),
                 round_up(frac(row['energy_interval'][1])+7*ehi+4))
                for row in local['records']]
    weights = [(down(frac(row['weight_interval'][0])), round_up(frac(row['weight_interval'][1])))
               for row in local['records']]
    bins = []
    lowers = []
    uppers = []
    for i, ((lo, hi), (wl, wu)) in enumerate(zip(energies, weights, strict=True)):
        center = (lo+hi)/2
        window = (center-radius, center+radius)
        inner_distance = radius-(hi-lo)/2
        assert inner_distance > v
        selected_lower = root(wl*(1-(v/inner_distance)**2))
        leakage = F(0)
        for j, ((other_lo, other_hi), (_, other_w)) in enumerate(zip(energies, weights, strict=True)):
            if j == i:
                continue
            distance = max(window[0]-other_hi, other_lo-window[1])
            assert distance > 0
            leakage += round_up(v*root(other_w, True)/distance)
        lower = max(F(0), selected_lower-leakage-eta)**2
        upper = min(F(1, 2), (root(wu, True)+leakage+eta)**2)
        lower, upper = down(lower), round_up(upper)
        lowers.append(lower)
        uppers.append(upper)
        prepared = (max(F(0), lower-delta), min(F(1), upper+delta))
        excitation = (window[0]-4-8*ehi,
                      window[1]-4-8*elo+source['neutral_shift'])
        bins.append(dict(index=i, absolute_active_source_plus_reference_bin=encoded(window),
                         source_excitation_relative_actual_ground=encoded(excitation),
                         exact_ground_input_probability=encoded((lower, upper)),
                         prepared_conditional_probability=encoded(prepared),
                         original_attempt_probability_lower=str(down(qlo*prepared[0])),
                         other_visible_atom_leakage_upper=str(leakage),
                         diagnostic_source_interval=list(map(float, excitation)),
                         diagnostic_original_attempt_lower=float(qlo*prepared[0])))
    assert all(F(a['absolute_active_source_plus_reference_bin'][1])
               < F(b['absolute_active_source_plus_reference_bin'][0])
               for a, b in zip(bins, bins[1:]))
    outside = (max(F(0), F(1, 2)-sum(uppers)), min(F(1, 2), F(1, 2)-sum(lowers)))
    assert 0 <= outside[0] <= outside[1]
    return dict(family=family, radius_in_hopping_norm_units=multiplier,
                full_hopping_norm=str(v), exact_ground_trial_vector_error=str(eta),
                actual_preparation_trace_distance_upper=str(delta),
                preparation_success_lower=str(qlo),
                exact_ground_pump_success='1/2', exact_ground_pump_failure='1/2',
                prepared_pump_failure_interval=encoded((F(1, 2)-delta, F(1, 2)+delta)),
                records=bins, exact_ground_outside_all_bins=encoded(outside),
                prepared_outside_all_bins=encoded((max(F(0), outside[0]-delta),
                                                  min(F(1), outside[1]+delta))),
                total_energy_pointer_dimension=len(bins)+1,
                success_reference_energy='4', failure_reference_energy='0',
                full_outside_bin_spectrum_identified=False,
                exact_source_spectral_service_declared=True,
                efficient_native_gate_synthesis_claimed=False,
                complete_law=['1-q', 'q tr(K_fail rho K_fail^dagger)',
                              'q tr(P_i B rho B^dagger), every listed bin i',
                              'q tr((I-sum_i P_i) B rho B^dagger)'])


if __name__ == '__main__':
    for family in CHANNELS:
        result = windows(family)
        print(family, len(result['records']), result['records'][0]['diagnostic_source_interval'],
              result['records'][0]['diagnostic_original_attempt_lower'], flush=True)
