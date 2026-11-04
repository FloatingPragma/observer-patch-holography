"""Full-Fock weak-contact consequences of an exhausted inversion-even kernel."""
from fractions import Fraction as F

from ..native_mass_channels.algebra import rational


def bounds(case, coupling):
    coupling = rational(coupling)
    if coupling <= 0:
        raise ValueError('positive exact contact coupling required')
    result = {'coupling': str(coupling), 'input_one_body_kernel': case['exact_kernel_dimension'],
              'bound_proves_ground_and_charged_pair': False,
              'bound_failure_is_a_physical_exclusion': False}
    if case['degree'] % 2 == 0 or case['exact_kernel_dimension'] != 2:
        return dict(result, reason='requires the exactly two-dimensional inversion-even kernel')
    gap = F(case['nonzero_gap_lower'])
    volume = (2*case['radius']+1)**3
    w = coupling*volume/4
    ipr = F(1, volume)
    if case['mass'] > 2:
        alpha = 1/(1+F(3, 2*(case['mass']-2)**2))
        ipr = max(ipr, alpha*alpha)
    result.update(sites=volume, full_Fock_dimension=2**(2*volume),
                  contact_norm=str(w), one_body_gap_lower=str(gap),
                  inverse_participation_lower=str(ipr),
                  full_Fock_matrix_executed=False)
    if 4*w >= gap:
        return dict(result, reason='complete-Fock contact norm is outside the sufficient gap domain')
    epsilon, angle = w*w/(gap-2*w), w/(gap-2*w)
    lower = coupling*ipr/2-2*epsilon
    result.update(ordered_low_energy_error=str(epsilon), projection_angle_upper=str(angle),
                  candidate_gap_lower_bound=str(lower))
    if lower <= 0 or 3*angle >= F(1, 2):
        return dict(result, reason='the sufficient ground-ordering or residue inequality is unresolved')
    return dict(result, bound_proves_ground_and_charged_pair=True,
                ground_particle_number=volume, exact_ground_multiplicity=2,
                charged_particle_numbers=[volume-1, volume+1],
                each_charged_state_multiplicity=1,
                charged_gap_interval=[str(lower), str(coupling/2+2*epsilon)],
                other_charged_threshold_lower=str(gap-2*w),
                each_pole_matrix_rank=2,
                each_pole_matrix_error=str(3*angle),
                pole_gap_is_relativistic_defect_rest_mass=False,
                volume_uniform_interacting_gap=False)
