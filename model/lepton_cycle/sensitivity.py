"""Architecture sensitivity and retrospective residuals, without parameter fitting."""
from __future__ import annotations
import mpmath as mp

# The retained comparison's rounded CODATA central ratio is comparison-only.
MUON_REFERENCE = '206.768283'


def spectrum(theta, balance, normalization=1):
    return [normalization * (1 + balance * mp.cos(theta + 2 * mp.pi * k / 3)) for k in range(3)]


def ratios(theta, balance=None, normalization=1):
    balance = mp.sqrt(2) if balance is None else balance
    values = spectrum(theta, balance, normalization)
    if values[1] == 0:
        raise ValueError('electron-character eigenvalue is zero')
    return values[2] ** 2 / values[1] ** 2, values[0] ** 2 / values[1] ** 2


def logarithmic_jacobian(theta, balance):
    angles = [theta + 2 * mp.pi * k / 3 for k in range(3)]
    base = [1 + balance * mp.cos(a) for a in angles]
    if any(v == 0 for v in base):
        raise ValueError('logarithmic sensitivity is singular at a zero eigenvalue')
    return mp.matrix([[2 * balance * (mp.sin(angles[1]) / base[1] - mp.sin(angles[k]) / base[k]),
                       2 * (mp.cos(angles[k]) / base[k] - mp.cos(angles[1]) / base[1])]
                      for k in (2, 0)])


def finite_phase(n, gap=1):
    if type(n) is not int or n < 1 or type(gap) is not int or gap < 1:
        raise ValueError('positive integer round count and current gap required')
    p = mp.mpf(2) / 9
    return -n * mp.arg(1 - p + p * mp.exp(-1j * gap / n)) / gap


def report():
    with mp.workdps(60):
        theta, balance = mp.mpf(2) / 9, mp.sqrt(2)
        muon, tau = ratios(theta, balance)
        jac = logarithmic_jacobian(theta, balance)
        residual = muon / mp.mpf(MUON_REFERENCE) - 1
        fmt = lambda x: mp.nstr(x, 25)
        finite = []
        for n in (1, 10, 100, 1000, 10000):
            phase = finite_phase(n)
            rm, rt = ratios(phase, balance)
            finite.append({'rounds': n, 'gap_one_phase_per_unit': fmt(phase),
                           'gap_two_phase_per_unit': fmt(finite_phase(n, 2)),
                           'muon_ratio_if_one_assigns_gap_one_phase': fmt(rm),
                           'tau_ratio_if_one_assigns_gap_one_phase': fmt(rt),
                           'muon_ratio_fractional_bias_from_infinite_round_value': fmt(rm / muon - 1)})
        return {'theta': '2/9', 'balance_2b_over_a': 'sqrt(2)',
                'muon_over_electron': fmt(muon), 'tau_over_electron': fmt(tau),
                'log_ratio_jacobian_columns_theta_balance': [[fmt(jac[i, j]) for j in range(2)] for i in range(2)],
                'jacobian_determinant': fmt(mp.det(jac)),
                'common_normalization_log_sensitivity': ['0', '0'],
                'exclusive_count_input_derivatives': {'dtheta_d_forward_probability_at_independence':'2/3',
                    'dtheta_d_written_probability_at_independence':'5/6',
                    'dtheta_d_joint_probability_at_fixed_marginals':'-2'},
                'marginals_without_independence': {'forward_probability':'1/12','written_probability':'1/6',
                    'allowed_joint_probability':['0','1/12'],'allowed_exclusive_phase':['1/12','1/4'],
                    'independent_joint_probability':'1/72','independent_exclusive_phase':'2/9'},
                'positive_chamber_phase_margin': fmt(mp.pi / 12 - theta),
                'koide_positive_chamber': '1/3 + balance^2/6',
                'comparison_only': {'reference_muon_ratio': MUON_REFERENCE,
                    'reference_convention': 'Rounded CODATA central value used by the retained paper table; no uncertainty supplied here.',
                    'relative_residual_ppm': fmt(residual * 10**6),
                    'linear_phase_displacement_for_residual_diagnostic': fmt(-mp.log(1 + residual) / jac[0, 0]),
                    'phase_displacement_applied': False, 'reference_used_to_select_parameters': False},
                'finite_round_diagnostics': finite,
                'interpretation': 'Squared-eigenvalue architecture ratios; no physical pole-mass or absolute-unit identification.'}
