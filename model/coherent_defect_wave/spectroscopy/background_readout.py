"""Complete source-relative phase records from a stationary background.

Every face projector and both clock-filter outcomes are retained. The active
filter annihilates the stationary component without selecting an energy log.
"""
from fractions import Fraction as F
from functools import lru_cache

from model.coherent_defect_wave.model import JointObserver
from model.coherent_defect_wave.background.model import background
from .readout import powers, overlap, _real_norm, filter_clock_branches


@lru_cache(maxsize=1)
def background_experiment():
    observer = JointObserver()
    omega = background(observer, internal=0)
    if observer.step(omega) != omega or _real_norm(omega) != 1:
        raise ArithmeticError('stationary source reference changed')
    omega_one = background(observer, internal=1)
    if (observer.step(omega_one) != omega_one or _real_norm(omega_one) != 1
            or overlap(omega, omega_one) != 0):
        raise ArithmeticError('orthogonal stationary energy-branch witness changed')
    records, summaries = [], []
    expected_second = {0:F(107,120), 2:F(29,30), 5:F(107,120),
                       16:F(29,30), 17:F(641,720)}
    for face_label, state in observer.region_face_measurement(omega):
        weight = _real_norm(state)
        if weight != F(1,5):
            raise ArithmeticError('complete face preparation weight changed')
        evolved = powers(observer, state)
        moments = tuple(overlap(state, v)/weight for v in evolved)
        face = face_label['face_index']
        if moments[2] != expected_second[face]:
            raise ArithmeticError('background-relative source response changed')
        rows = filter_clock_branches(evolved)
        active = F(0)
        for label, output in rows:
            if not output.gauss_positive():
                raise ArithmeticError('background readout left the Gauss sector')
            probability = output.norm_squared()
            nonflat = label['filter'] == 'nonflat_filter'
            orthogonal = (overlap(omega, output.real) == 0
                          and overlap(omega, output.imag) == 0)
            if nonflat:
                if not orthogonal:
                    raise ArithmeticError('active readout contains the stationary reference')
                active += probability
            records.append({**face_label, **label, 'probability':str(probability),
                            'positive_constraint_preserved':True,
                            'orthogonal_to_stationary_reference':orthogonal})
        if sum(output.norm_squared() for _,output in rows) != weight:
            raise ArithmeticError('a prepared face lost a complete clock outcome')
        if active != weight*(1-moments[2])/2:
            raise ArithmeticError('filter output disagrees with source moment')
        summaries.append({'face_index':face, 'face':face_label['face'],
                          'preparation_probability':str(weight),
                          'conditional_moments':list(map(str,moments)),
                          'nonflat_probability_per_background_trial':str(active)})
    if len(records) != 40 or sum(F(r['probability']) for r in records) != 1:
        raise ArithmeticError('background-relative experiment is incomplete')
    active_total = sum(F(r['probability']) for r in records
                       if r['filter'] == 'nonflat_filter')
    if active_total != F(283,7200):
        raise ArithmeticError('total background-relative active weight changed')
    return {'initial_reference':'Omega_0: uniform links, ports and addresses; internal0; parity0',
            'conditioned_on_successful_background_preparation':True,
            'background_preparation_evidence':'model/coherent_defect_wave/background',
            'source_tick':'the unchanged joint U',
            'maximum_actual_U_power':5,
            'complete_original_background_trial_records':records,
            'face_preparations':summaries,
            'total_nonflat_probability_per_background_trial':str(active_total),
            'all_active_outputs_orthogonal_to_stationary_reference':True,
            'same_source_positive_unique_ground_branch_control':{
                'base_generator':'H0: U spectral phases in [0,2pi), with U=exp(-i H0)',
                'generators':'Hk = H0 + 2pi*k*(I-|Omega0><Omega0|), k=1,2',
                'same_actual_joint_U':True,
                'same_unique_ground_vector':'Omega0',
                'orthogonal_gap_eigenvector':'Omega1',
                'exact_gaps_in_pi_units':['2','4'],
                'nonlocal_countercontrols_not_selected_physical_generators':True},
            'scope':{'stationary_source_reference_used':True,
                     'all_five_face_outcomes_retained':True,
                     'physical_ground_state_or_energy_gap_identified':False,
                     'particle_species_identified':False}}
