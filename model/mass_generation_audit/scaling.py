"""Quantitative charge-threshold scaling, without an inserted Dirac mass."""
from fractions import Fraction as F

from ..common_source.mass_response import gap_certificate as old


I = old.RationalInterval


def gap_envelope(hopping, repulsion):
    if type(hopping) not in (int, F) or type(repulsion) not in (int, F) or min(hopping, repulsion) <= 0:
        raise ValueError('strictly positive exact source couplings required')
    t, u = I.enclose(hopping), I.enclose(repulsion)
    pi = old.pi_bounds(); a = 2*pi*t/u
    exponential = old.exp_negative_bounds(a)
    leading = 8/pi*old.sqrt_bounds(t*u)*exponential
    ratio = (1+3/(8*a))/(1-exponential**2)
    return I(leading.lower, (leading*ratio).upper), leading, ratio


def certificate():
    speed, repulsion = F(1, 8), F(1)
    rows = []
    for power in range(9):
        spacing = F(1, 2**power)
        hopping = speed/(2*spacing)
        gap, _, _ = gap_envelope(hopping, repulsion)
        rows.append({'spacing': str(spacing), 'hopping': str(hopping),
                     'fixed_repulsion': str(repulsion), 'charge_gap_envelope': gap.as_json()})
    running = []
    for a in (4, 8, 16, 32, 64):
        x = I.enclose(a); exponential = old.exp_negative_bounds(x)
        relative = (1+3/(8*x))/(1-exponential**2)
        running.append({'A': a, 'half_gap_over_target_interval': I(F(1), relative.upper).as_json()})
    return {'Bessel_bound': 'base(z)<=K1(z)<=base(z)(1+3/(8z)),base=sqrt(pi/(2z))*exp(-z)',
            'threshold_envelope': 'Delta0<=Delta<=Delta0*(1+3/(8A))/(1-exp(-2A)),A=2pi*t/U,Delta0=(8/pi)sqrt(Ut)exp(-A)',
            'fixed_U_refinement': 't=v/(2a),U fixed>0 implies Delta->0 as a->0',
            'fixed_U_decay': 'Delta0=(8/pi)sqrt(Uv/(2a))*exp[-pi*v/(aU)]',
            'explicit_fixed_U_enclosures': rows,
            'running_condition': 'For chosen M>0,v>0: a(A)=4v*exp(-A)/(M*sqrt(2pi*A)), t=v/(2a), U=pi*v/(a*A)',
            'running_half_gap_limit': '1<=Delta/(2M)<=(1+3/(8A))/(1-exp(-2A))->1',
            'running_relative_bounds': running,
            'target_M_determined_by_fixed_source': False,
            'running_coupling_is_additional_family_input': True,
            'quantitative_charge_threshold_only': True,
            'fixed_lattice_Lorentz_dispersion_proved': False,
            'interacting_three_dimensional_mass_generation_proved': False,
            'reference': 'https://arxiv.org/html/cond-mat/9911281'}
